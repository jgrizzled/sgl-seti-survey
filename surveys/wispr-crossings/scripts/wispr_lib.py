"""PSP/WISPR-I measurement chain (hypotheses freeze D3).

Per frame: fetch (NRL L3 synoptic WISPR-I; size-checked) -> load ->
frame gate (InnerFFV synoptic full field, 960 x 1024, NSUMEXP >= 5) ->
rescale L3 -> L2 units by (DSUN_OBS / 0.2 AU)^-2.3 (the header value the
L3 pipeline used) -> header celestial WCS ('A' RA/DEC-ZPN) refined by a
Hipparcos match (translation; >= 30 stars, rms <= 1.5 px) -> 15-px
median high-pass -> 2-D colour-corrected star ZP (quadratic in x, y;
Hipparcos V <= 8) -> top-hat aperture photometry (r 3 px, annulus 5-9)
at ICRS-fixed patches.

Flux units: MSB (L2) normalised to the frame's ZP at the patch position,
F_norm = F * 10^(-0.4 (ZP(x,y) - ZP_REF)) with ZP_REF = 8.65 = the recon
B-V = 0.65 zero point (-21.35) + 30 mag, i.e. one flux unit is 1e-12 of
the MSB-based scale and V_eq = 8.65 - 2.5 log10 F_norm (a V = 8.65
solar-colour star reads 1.0).
"""

from __future__ import annotations

import hashlib
import time
import urllib.error
import urllib.request
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from scipy import ndimage
from scipy.spatial import cKDTree

REPO = Path(__file__).resolve().parents[3]
FRAMES = REPO / "runs" / "wispr-crossings" / "frames"
STARCAT = REPO / "runs" / "heliospheric-crossings" / "starcat" / "hip_v95.npz"
TYCHO = REPO / "runs" / "heliospheric-crossings" / "starcat" / "tyc2_v11.npz"

NRL = "https://wispr.nrl.navy.mil/data/rel/fits"
FRAME_BYTES = 3_954_240            # 2x2-binned 960x1024 float32 (every orbit but 02)
FRAME_BYTES_UNBINNED = 15_750_720  # orbit 02 (E2): unbinned 1920x2048, binned on load
AU_M = 1.495978707e11

NX, NY = 960, 1024
PIX_ARCSEC = 152.3
AP_R, ANN_IN, ANN_OUT = 3.0, 5.0, 9.0
HP_SIZE = 15
ZP_REF = 8.65                   # V-equivalent of one flux unit; unit = 1e-12 of the MSB-based scale (frame ZPs are ~ -21.35)
COLOUR_COEFF = 0.27            # mag/(B-V), recon mean of three frames; re-measured at dev
BV_REF = 0.65
PSF_FWHM_PX = 1.6
MIN_STARS, MAX_RMS_PX = 30, 1.5   # rms after the affine refinement
CALIB_VMAX, MIN_CALIB, MAX_ZP_MAD, CALIB_SNR_MIN = 8.0, 30, 9.9, 15.0   # v1.1: S/N >= 15 calibrators, >= 30 (linear ZP surface below 60); v1.5: the MAD gate is replaced by MAX_ZP_UNC
MAX_ZP_UNC = 0.10               # v1.5: ZP uncertainty gate, MAD / sqrt(n_calib) <= 0.10 mag
QUAD_MIN = 60

_cat = {k: np.array(v) for k, v in np.load(STARCAT).items()}   # eager (fork-safe)


# ----------------------------------------------------------------- fetch

def _sha(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def fetch_frame(orbit_day: str, fname: str, dest_dir: Path = FRAMES) -> tuple[Path | None, str]:
    """Fetch an L3 frame `fname` from `L3/<orbit>/<day>/`. Returns
    (path, source). Size-checked (3,954,240 bytes) and FITS-signature
    checked on every fetch. `requests` with connect/read timeouts and a
    streamed body with a wall-clock cap (a stalled TLS handshake or read
    raises instead of hanging — the first dev run hung in do_handshake)."""
    import requests

    orbit, day = orbit_day.split("/")
    d = dest_dir / day
    d.mkdir(parents=True, exist_ok=True)
    p = d / fname
    if p.exists() and p.stat().st_size in (FRAME_BYTES, FRAME_BYTES_UNBINNED):
        return p, "cached"
    url = f"{NRL}/L3/{orbit}/{day}/{fname}"
    for attempt in range(3):
        try:
            with requests.get(url, stream=True, timeout=(15, 45)) as r:
                if r.status_code == 404:
                    return None, "404"
                r.raise_for_status()
                chunks, t0 = [], time.time()
                for ch in r.iter_content(262144):
                    chunks.append(ch)
                    if time.time() - t0 > 90:
                        raise IOError("slow body")
                raw = b"".join(chunks)
            if len(raw) not in (FRAME_BYTES, FRAME_BYTES_UNBINNED) or not raw.startswith(b"SIMPLE"):
                raise IOError(f"bad payload {len(raw)} bytes")
            p.write_bytes(raw)
            return p, "nrl"
        except Exception:
            time.sleep(2 + 3 * attempt)
    return None, "failed"


# ----------------------------------------------------------------- frame

@dataclass
class Frame:
    path: str
    mjd_avg: float
    xposure: float
    nsumexp: int
    gain: str
    code: str
    r_au: float                     # header DSUN (for the rescale only)
    scale: float                    # (r/0.2)^2.3 divided out
    wcs: WCS
    raw: np.ndarray                 # rescaled L3 (L2 units)
    img: np.ndarray                 # high-passed
    translation: np.ndarray = field(default_factory=lambda: np.zeros(2))
    affine: np.ndarray | None = None      # 2x3: [x', y'] = A @ [x, y, 1] (refinement on header pixels)
    n_match: int = 0
    astrom_rms: float = np.inf
    zp_coef: np.ndarray | None = None
    zp_scatter: float = np.inf
    n_calib: int = 0
    calibrators: list = field(default_factory=list)

    def zp_at(self, x, y) -> float:
        if self.zp_coef is None:
            return float("nan")
        u, v = (x - NX / 2) / (NX / 2), (y - NY / 2) / (NY / 2)
        return float(np.array([1, u, v, u * u, u * v, v * v]) @ self.zp_coef)


def frame_gate(h: fits.Header) -> str | None:
    if (h.get("NAXIS1"), h.get("NAXIS2")) not in ((NX, NY), (2 * NX, 2 * NY)):
        return f"shape_{h.get('NAXIS1')}x{h.get('NAXIS2')}"
    if str(h.get("OBJECT", "")).strip() != "InnerFFV":
        return f"object_{h.get('OBJECT')}"
    if int(h.get("NSUMEXP", 0)) < 5:
        return f"nsumexp_{h.get('NSUMEXP')}"
    if "CRVAL1A" not in h:
        return "no_celestial_wcs"
    if str(h.get("LEVEL", "")).strip() != "L3":
        return f"level_{h.get('LEVEL')}"
    return None


def load_frame(path: Path) -> Frame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        hdu = fits.open(path)[0]
        h = hdu.header
        wcs = WCS(h, key="A")
    r_au = float(h["DSUN_OBS"]) / AU_M
    scale = (r_au / 0.2) ** 2.3
    raw = np.asarray(hdu.data, np.float64) / scale
    if raw.shape == (2 * NY, 2 * NX):            # orbit 02: unbinned -> 2x2 mean (MSB is per-pixel brightness)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            raw = np.nanmean(raw.reshape(NY, 2, NX, 2), axis=(1, 3))
            wcs = wcs.slice((slice(0, None, 2), slice(0, None, 2)))
    fill = np.nan_to_num(raw, nan=float(np.nanmedian(raw)))
    img = raw - ndimage.median_filter(fill, size=HP_SIZE)
    name = Path(path).name
    return Frame(path=str(path), mjd_avg=Time(h["DATE-AVG"], scale="utc").mjd,
                 xposure=float(h["XPOSURE"]), nsumexp=int(h.get("NSUMEXP", 0)),
                 gain=str(h.get("GAINMODE", "")), code=name.split("_")[-1][:4],
                 r_au=r_au, scale=scale, wcs=wcs, raw=raw, img=img)


# ----------------------------------------------------------------- photometry

def aper_flux(img: np.ndarray, x: float, y: float) -> tuple[float, float, float]:
    ny, nx = img.shape
    ix, iy = int(round(x)), int(round(y))
    r = int(np.ceil(ANN_OUT))
    if not (r < ix < nx - r and r < iy < ny - r):
        return np.nan, np.nan, 0.0
    sub = img[iy - r: iy + r + 1, ix - r: ix + r + 1]
    yy, xx = np.mgrid[-r:r + 1, -r:r + 1]
    rr = np.hypot(xx - (x - ix), yy - (y - iy))
    ann = sub[(rr >= ANN_IN) & (rr <= ANN_OUT)]
    ann = ann[np.isfinite(ann)]
    if len(ann) < 30:
        return np.nan, np.nan, 0.0
    bg = float(np.median(ann))
    noise = float(1.4826 * np.median(np.abs(ann - bg))) + 1e-30
    a = rr <= AP_R
    vals = sub[a]
    fin = np.isfinite(vals)
    valid = float(np.mean(fin))
    if fin.sum() == 0:
        return np.nan, np.nan, 0.0
    return float(np.sum(vals[fin] - bg)), float(noise * np.sqrt(a.sum())), valid


def sky_to_pix(fr: Frame, ra, dec) -> tuple[np.ndarray, np.ndarray]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        px, py = fr.wcs.all_world2pix(np.atleast_1d(ra), np.atleast_1d(dec), 0)
    if fr.affine is not None:
        A = fr.affine
        return A[0, 0] * px + A[0, 1] * py + A[0, 2], A[1, 0] * px + A[1, 1] * py + A[1, 2]
    return px + fr.translation[0], py + fr.translation[1]


def detect(fr: Frame, snr_min: float = 6.0):
    fin = np.isfinite(fr.img)
    hp = np.where(fin, fr.img, 0.0)
    sig = 1.4826 * np.median(np.abs(hp[fin])) + 1e-30
    pk = (ndimage.maximum_filter(hp, size=5) == hp) & (hp > snr_min * sig) & fin
    pk[:6, :] = pk[-6:, :] = pk[:, :6] = pk[:, -6:] = False
    ys, xs = np.nonzero(pk)
    return np.c_[xs, ys].astype(float)


def catalog_at(mjd: float):
    dyr = (mjd - 48348.5625) / 365.25
    ra = _cat["ra"] + _cat["pmra_masyr"] * dyr / 3.6e6 / np.cos(np.radians(_cat["dec"]))
    dec = _cat["dec"] + _cat["pmde_masyr"] * dyr / 3.6e6
    return ra, dec


def fit_frame(fr: Frame) -> bool:
    """Verify + refine the header WCS against Hipparcos; fit the 2-D
    colour-corrected ZP. Returns astrometric validity."""
    ra, dec = catalog_at(fr.mjd_avg)
    # ZPN inverse can wrap far-off-axis stars into the frame: keep the
    # 35-deg cone about the frame centre
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c_ra, c_dec = fr.wcs.all_pix2world([[NX / 2, NY / 2]], 0)[0]
    cv = np.array([np.cos(np.radians(c_dec)) * np.cos(np.radians(c_ra)),
                   np.cos(np.radians(c_dec)) * np.sin(np.radians(c_ra)), np.sin(np.radians(c_dec))])
    sv = np.c_[np.cos(np.radians(dec)) * np.cos(np.radians(ra)),
               np.cos(np.radians(dec)) * np.sin(np.radians(ra)), np.sin(np.radians(dec))]
    near = sv @ cv > np.cos(np.radians(35))
    fr.translation = np.zeros(2)
    fr.affine = None
    px, py = sky_to_pix(fr, ra[near], dec[near])
    inside = np.isfinite(px) & (px > 10) & (px < NX - 10) & (py > 10) & (py < NY - 10)
    if inside.sum() < MIN_STARS:
        return False
    px, py = px[inside], py[inside]
    vmag, bv = _cat["vmag"][near][inside], _cat["bv"][near][inside]
    det = detect(fr)
    if len(det) < MIN_STARS:
        return False
    tree = cKDTree(det)
    bright = vmag <= 7.5
    d, i = tree.query(np.c_[px[bright], py[bright]], distance_upper_bound=3.0)
    ok = np.isfinite(d)
    if ok.sum() < MIN_STARS:
        return False
    # translation, then a 2x3 affine on the matched pairs (ZPN residual
    # scale/rotation), then a tighter re-match
    off = np.median(det[i[ok]] - np.c_[px[bright][ok], py[bright][ok]], axis=0)
    fr.translation = off
    d, i = tree.query(np.c_[px[bright] + off[0], py[bright] + off[1]], distance_upper_bound=3.0)
    ok = np.isfinite(d)
    if ok.sum() >= MIN_STARS:
        P = np.c_[px[bright][ok], py[bright][ok], np.ones(ok.sum())]
        Q = det[i[ok]]
        A, *_ = np.linalg.lstsq(P, Q, rcond=None)
        fr.affine = A.T
        px2, py2 = sky_to_pix(fr, ra[near][inside][bright], dec[near][inside][bright])
        d, i = tree.query(np.c_[px2, py2], distance_upper_bound=2.5)
        ok = np.isfinite(d)
    fr.n_match = int(ok.sum())
    fr.astrom_rms = float(np.sqrt(np.mean(d[ok] ** 2))) if ok.any() else np.inf
    if fr.n_match < MIN_STARS or fr.astrom_rms > MAX_RMS_PX:
        return False
    rows = []
    pxr, pyr = sky_to_pix(fr, ra[near][inside], dec[near][inside])
    for k in np.nonzero(vmag <= CALIB_VMAX)[0]:
        if not (-0.3 <= bv[k] <= 1.6) or vmag[k] < 2.5:
            continue
        x, y = float(pxr[k]), float(pyr[k])
        f, e, val = aper_flux(fr.img, x, y)
        if np.isfinite(f) and f > 0 and e > 0 and f / e >= CALIB_SNR_MIN and val >= 0.99:
            zp = vmag[k] - COLOUR_COEFF * (bv[k] - BV_REF) + 2.5 * np.log10(f)
            rows.append((x, y, zp, vmag[k], bv[k], f))
    if len(rows) < MIN_CALIB:
        fr.zp_coef = None
        return True
    R = np.array(rows)
    u, v = (R[:, 0] - NX / 2) / (NX / 2), (R[:, 1] - NY / 2) / (NY / 2)
    A6 = np.c_[np.ones(len(R)), u, v, u * u, u * v, v * v]
    ncol = 6 if len(R) >= QUAD_MIN else 3          # linear surface on sparse (perihelion) frames
    A = A6[:, :ncol]
    coef, *_ = np.linalg.lstsq(A, R[:, 2], rcond=None)
    resid = R[:, 2] - A @ coef
    mad = 1.4826 * np.median(np.abs(resid))
    keep = np.abs(resid) < 3 * max(mad, 0.05)
    if keep.sum() >= MIN_CALIB:
        coef, *_ = np.linalg.lstsq(A[keep], R[keep, 2], rcond=None)
        resid = R[:, 2] - A @ coef
        mad = 1.4826 * np.median(np.abs(resid[keep]))
    coef = np.r_[coef, np.zeros(6 - ncol)]
    fr.zp_coef, fr.zp_scatter, fr.n_calib = coef, float(mad), int(keep.sum())
    fr.calibrators = [(float(a), float(b), float(c), float(dd), float(ee), float(ff))
                      for (a, b, c, dd, ee, ff), kp in zip(R, keep) if kp]
    return True


def measure_patch(fr: Frame, ra: float, dec: float) -> dict:
    x, y = sky_to_pix(fr, ra, dec)
    x, y = float(x[0]), float(y[0])
    f, e, val = aper_flux(fr.img, x, y)
    zp = fr.zp_at(x, y)
    scale = 10 ** (-0.4 * (zp - ZP_REF)) if np.isfinite(zp) else np.nan
    return {"mjd": fr.mjd_avg, "x": round(x, 2), "y": round(y, 2),
            "flux": f * scale if np.isfinite(f) else np.nan,
            "err": e * scale if np.isfinite(e) else np.nan,
            "flux_msb": f, "valid": round(val, 3), "zp": round(zp, 4) if np.isfinite(zp) else None}


def in_footprint(x: float, y: float, margin_px: float = 12.0) -> bool:
    return margin_px <= x <= NX - 1 - margin_px and margin_px <= y <= NY - 1 - margin_px


# ----------------------------------------------------------------- template

def encircled_fraction(sep_px: np.ndarray, r: float = AP_R, fwhm: float = PSF_FWHM_PX) -> np.ndarray:
    """Fraction of a Gaussian PSF (FWHM fwhm) centred sep_px from the
    aperture centre that falls inside radius r (numeric, 0.1-px grid)."""
    sig = fwhm / 2.355
    g = np.arange(-8, 8.01, 0.1)
    yy, xx = np.meshgrid(g, g)
    out = np.empty(len(sep_px))
    for k, s in enumerate(np.atleast_1d(sep_px)):
        psf = np.exp(-0.5 * ((xx - s) ** 2 + yy ** 2) / sig ** 2)
        out[k] = psf[np.hypot(xx, yy) <= r].sum() / psf.sum()
    return out


def template_stars(ra0: float, dec0: float, mjd: float, tyc, tyc_tree, hip_tree,
                   max_sep_px: float = AP_R + 2.0, extra: list | None = None,
                   exclude_sep_arcsec: float = 0.0) -> list:
    """Catalogue stars near a patch centre: [(V, B-V, sep_px)] from the
    union of Tycho-2 (VT <= 11; V ~ VT - 0.06, B-V from a Hipparcos match
    within 5 arcsec else 0.65) and Hipparcos (V <= 9.5, all stars — Tycho-2
    lacks the brightest and the high-proper-motion stars), de-duplicated
    within 10 arcsec. `extra` = explicit (V, B-V, sep_px) rows (the S2
    target star from the config table); catalogue stars within
    `exclude_sep_arcsec` of the centre are dropped when `extra` is given."""
    v0 = np.array([np.cos(np.radians(dec0)) * np.cos(np.radians(ra0)),
                   np.cos(np.radians(dec0)) * np.sin(np.radians(ra0)), np.sin(np.radians(dec0))])
    lim = 2 * np.sin(np.radians(max_sep_px * PIX_ARCSEC / 3600) / 2)
    dup = 2 * np.sin(np.radians(10 / 3600) / 2)
    out, vecs = [], []
    hra, hdec = catalog_at(mjd)
    for j in hip_tree.query_ball_point(v0, lim):
        sv = np.array([np.cos(np.radians(hdec[j])) * np.cos(np.radians(hra[j])),
                       np.cos(np.radians(hdec[j])) * np.sin(np.radians(hra[j])), np.sin(np.radians(hdec[j]))])
        sep_as = 2 * np.degrees(np.arcsin(np.linalg.norm(sv - v0) / 2)) * 3600
        if extra and sep_as < exclude_sep_arcsec:
            continue
        bv = float(_cat["bv"][j]) if np.isfinite(_cat["bv"][j]) and _cat["bv"][j] < 90 else BV_REF
        out.append((float(_cat["vmag"][j]), bv, float(sep_as / PIX_ARCSEC)))
        vecs.append(sv)
    for j in tyc_tree.query_ball_point(v0, lim):
        ra, dec = float(tyc["ra"][j]), float(tyc["dec"][j])
        sv = np.array([np.cos(np.radians(dec)) * np.cos(np.radians(ra)),
                       np.cos(np.radians(dec)) * np.sin(np.radians(ra)), np.sin(np.radians(dec))])
        if any(np.linalg.norm(sv - w) < dup for w in vecs):
            continue
        sep_as = 2 * np.degrees(np.arcsin(np.linalg.norm(sv - v0) / 2)) * 3600
        if extra and sep_as < exclude_sep_arcsec:
            continue
        d, i = hip_tree.query(sv, distance_upper_bound=2 * np.sin(np.radians(5 / 3600) / 2))
        bv = float(_cat["bv"][i]) if np.isfinite(d) and np.isfinite(_cat["bv"][i]) and _cat["bv"][i] < 90 else BV_REF
        out.append((float(tyc["vt"][j]) - 0.06, bv, float(sep_as / PIX_ARCSEC)))
        vecs.append(sv)
    if extra:
        out.extend(tuple(map(float, e)) for e in extra)
    return out


def template_flux(stars: list, colour_coeff: float = COLOUR_COEFF) -> float:
    """ZP_REF-normalised aperture flux predicted for the catalogue stars."""
    if not stars:
        return 0.0
    s = np.array(stars, float)
    v_inst = s[:, 0] - colour_coeff * (s[:, 1] - BV_REF)
    f = 10 ** (-0.4 * (v_inst - ZP_REF))
    return float(np.sum(f * encircled_fraction(s[:, 2])))
