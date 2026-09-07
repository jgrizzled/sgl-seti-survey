"""Solar Orbiter / SoloHI inner-tile measurement chain (hypotheses freeze D3).

Per frame: fetch (NRL L2 tree or the SOAR data service; size-checked)
-> load -> frame gate (inner tile 1ft/2ft, 1024 x 960, NBIN 4, NSUMEXP
5, L2, celestial WCS present) -> header celestial WCS ('A' RA/DEC-ZPN)
refined by a Hipparcos match (translation + 2x3 affine; >= 30 stars,
rms <= 1.5 px) -> 25-px median high-pass -> 2-D colour-corrected star
ZP (quadratic in x, y; Hipparcos V <= 9) -> top-hat aperture photometry
(r 3 px, annulus 5-9) at ICRS-fixed patches.

Flux units: MSB normalised to the frame's ZP at the patch position,
F_norm = F * 10^(-0.4 (ZP(x,y) - ZP_REF)) with ZP_REF = 10.05 = the recon
B-V = 0.65 zero point (-19.95) + 30 mag, i.e. one flux unit is 1e-12 of
the MSB-based scale and V_eq = 10.05 - 2.5 log10 F_norm.
"""

from __future__ import annotations

import hashlib
import time
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
FRAMES = REPO / "runs" / "solohi-crossings" / "frames"
STARCAT = REPO / "runs" / "heliospheric-crossings" / "starcat" / "hip_v95.npz"

NRL = "https://solohi.nrl.navy.mil/so_data/L2"
SOAR = "https://soar.esac.esa.int/soar-sl-tap/data"
FRAME_BYTES = 3_960_000            # 2x2-binned 1024x960 float32
AU_M = 1.495978707e11

NX, NY = 1024, 960
PIX_ARCSEC = 74.2
AP_R, ANN_IN, ANN_OUT = 3.0, 5.0, 9.0
HP_SIZE = 25
ZP_REF = 10.05                  # V-equivalent of one flux unit; unit = 1e-12 of the MSB scale (frame ZPs ~ -19.95)
COLOUR_COEFF = 0.38             # mag/(B-V), recon mean of 14 frames; re-measured at dev
BV_REF = 0.65
PSF_FWHM_PX = 1.8
MIN_STARS, MAX_RMS_PX = 30, 1.5
CALIB_VMAX, MIN_CALIB, CALIB_SNR_MIN = 9.0, 30, 15.0
MAX_ZP_UNC = 0.10               # ZP uncertainty gate, MAD / sqrt(n_calib) <= 0.10 mag (WISPR v1.5)
SAT_FRAC = 0.85                 # v1.1: pixels >= SAT_FRAC * DSATVAL are masked (the F-corona plateau at the
                                # sunward edge in the >= 45 s exposure regime sits at 0.90-0.92 DSATVAL)
QUAD_MIN = 60

_cat = {k: np.array(v) for k, v in np.load(STARCAT).items()}   # eager (fork-safe)


# ----------------------------------------------------------------- fetch

def _sha(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def fetch_frame(day: str, fname: str, dest_dir: Path = FRAMES, source: str = "nrl") -> tuple[Path | None, str]:
    """Fetch an L2 frame `fname` (day YYYYMMDD) from the NRL tree or the
    SOAR data service (data_item_id = stem without _V0N.fits). Returns
    (path, source). Size-checked (3,960,000 bytes) and FITS-signature
    checked on every fetch; `requests` with timeouts and a streamed body
    with a wall-clock cap."""
    import requests

    d = dest_dir / day
    d.mkdir(parents=True, exist_ok=True)
    p = d / fname
    if p.exists() and p.stat().st_size == FRAME_BYTES:
        return p, "cached"
    for attempt in range(3):
        src = source if attempt == 0 else ("soar" if source == "nrl" else "nrl")
        try:
            if src == "nrl":
                req = requests.get(f"{NRL}/{day}/{fname}", stream=True, timeout=(15, 60))
            else:
                req = requests.get(SOAR, params={"retrieval_type": "LAST_PRODUCT",
                                                 "data_item_id": fname[:-9], "product_type": "SCIENCE"},
                                   stream=True, timeout=(15, 60))
            with req as r:
                if r.status_code == 404:
                    continue
                r.raise_for_status()
                chunks, t0 = [], time.time()
                for ch in r.iter_content(262144):
                    chunks.append(ch)
                    if time.time() - t0 > 120:
                        raise IOError("slow body")
                raw = b"".join(chunks)
            if len(raw) != FRAME_BYTES or not raw.startswith(b"SIMPLE"):
                raise IOError(f"bad payload {len(raw)} bytes")
            p.write_bytes(raw)
            return p, src
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
    tile: str
    mode: str
    roll: float
    r_au: float
    wcs: WCS
    raw: np.ndarray
    img: np.ndarray
    dbox: tuple = (1, NX, 1, NY)    # DSTART1, DSTOP1, DSTART2, DSTOP2 (1-based valid region)
    nx: int = NX                    # tile 1 is 1024 x 960, tile 2 is 960 x 1024 (RECTROTA 5 / 6)
    ny: int = NY
    sat_frac_masked: float = 0.0
    translation: np.ndarray = field(default_factory=lambda: np.zeros(2))
    affine: np.ndarray | None = None
    n_match: int = 0
    astrom_rms: float = np.inf
    zp_coef: np.ndarray | None = None
    zp_scatter: float = np.inf
    n_calib: int = 0
    calibrators: list = field(default_factory=list)

    def zp_at(self, x, y) -> float:
        if self.zp_coef is None:
            return float("nan")
        u, v = (x - self.nx / 2) / (self.nx / 2), (y - self.ny / 2) / (self.ny / 2)
        return float(np.array([1, u, v, u * u, u * v, v * v]) @ self.zp_coef)


def frame_gate(h: fits.Header) -> str | None:
    if (h.get("NAXIS1"), h.get("NAXIS2")) not in ((NX, NY), (NY, NX)):
        return f"shape_{h.get('NAXIS1')}x{h.get('NAXIS2')}"
    if str(h.get("DETECTOR", "")).strip() not in ("1", "2"):
        return f"detector_{h.get('DETECTOR')}"
    if int(h.get("NSUMEXP", 0)) < 5:
        return f"nsumexp_{h.get('NSUMEXP')}"
    if int(h.get("NBIN", 0)) != 4:
        return f"nbin_{h.get('NBIN')}"
    if "CRVAL1A" not in h:
        return "no_celestial_wcs"
    if str(h.get("LEVEL", "")).strip() != "L2":
        return f"level_{h.get('LEVEL')}"
    return None


def load_frame(path: Path) -> Frame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        hdu = fits.open(path)[0]
        h = hdu.header
        wcs = WCS(h, key="A")
    raw = np.asarray(hdu.data, np.float64)
    satval = float(h.get("DSATVAL", np.inf))
    if np.isfinite(satval) and satval > 0:
        raw = np.where(raw >= SAT_FRAC * satval, np.nan, raw)          # v1.1 saturation mask
    fill = np.nan_to_num(raw, nan=float(np.nanmedian(raw)))
    img = raw - ndimage.median_filter(fill, size=HP_SIZE)
    return Frame(path=str(path), mjd_avg=Time(h["DATE-AVG"], scale="utc").mjd,
                 xposure=float(h["XPOSURE"]), nsumexp=int(h.get("NSUMEXP", 0)),
                 gain=str(h.get("GAINMODE", "")), tile=str(h["DETECTOR"]).strip(),
                 mode=str(h.get("OBS_MODE", "")), roll=float(h.get("SC_ROLL", np.nan)),
                 r_au=float(h["DSUN_OBS"]) / AU_M, wcs=wcs, raw=raw, img=img,
                 dbox=(int(h.get("DSTART1", 1)), int(h.get("DSTOP1", h["NAXIS1"])),
                       int(h.get("DSTART2", 1)), int(h.get("DSTOP2", h["NAXIS2"]))),
                 nx=int(h["NAXIS1"]), ny=int(h["NAXIS2"]),
                 sat_frac_masked=float(np.mean(~np.isfinite(raw))))


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
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c_ra, c_dec = fr.wcs.all_pix2world([[fr.nx / 2, fr.ny / 2]], 0)[0]
    cv = np.array([np.cos(np.radians(c_dec)) * np.cos(np.radians(c_ra)),
                   np.cos(np.radians(c_dec)) * np.sin(np.radians(c_ra)), np.sin(np.radians(c_dec))])
    sv = np.c_[np.cos(np.radians(dec)) * np.cos(np.radians(ra)),
               np.cos(np.radians(dec)) * np.sin(np.radians(ra)), np.sin(np.radians(dec))]
    near = sv @ cv > np.cos(np.radians(25))
    fr.translation = np.zeros(2)
    fr.affine = None
    px, py = sky_to_pix(fr, ra[near], dec[near])
    inside = np.isfinite(px) & (px > 10) & (px < fr.nx - 10) & (py > 10) & (py < fr.ny - 10)
    if inside.sum() < MIN_STARS:
        return False
    px, py = px[inside], py[inside]
    vmag, bv = _cat["vmag"][near][inside], _cat["bv"][near][inside]
    det = detect(fr)
    if len(det) < MIN_STARS:
        return False
    tree = cKDTree(det)
    bright = vmag <= 8.0
    d, i = tree.query(np.c_[px[bright], py[bright]], distance_upper_bound=3.5)
    ok = np.isfinite(d)
    if ok.sum() < MIN_STARS:
        return False
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
        if not (-0.3 <= bv[k] <= 1.6) or vmag[k] < 3.0:
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
    u, v = (R[:, 0] - fr.nx / 2) / (fr.nx / 2), (R[:, 1] - fr.ny / 2) / (fr.ny / 2)
    A6 = np.c_[np.ones(len(R)), u, v, u * u, u * v, v * v]
    ncol = 6 if len(R) >= QUAD_MIN else 3
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


def in_footprint(fr: Frame, x: float, y: float, margin_px: float = 12.0) -> bool:
    x0, x1, y0, y1 = fr.dbox
    return (x0 - 1 + margin_px) <= x <= (x1 - 1 - margin_px) and (y0 - 1 + margin_px) <= y <= (y1 - 1 - margin_px)


# ----------------------------------------------------------------- template (Gaia DR3)

def encircled_fraction(sep_px: np.ndarray, r: float = AP_R, fwhm: float = PSF_FWHM_PX) -> np.ndarray:
    sig = fwhm / 2.355
    g = np.arange(-8, 8.01, 0.1)
    yy, xx = np.meshgrid(g, g)
    out = np.empty(len(np.atleast_1d(sep_px)))
    for k, s in enumerate(np.atleast_1d(sep_px)):
        psf = np.exp(-0.5 * ((xx - s) ** 2 + yy ** 2) / sig ** 2)
        out[k] = psf[np.hypot(xx, yy) <= r].sum() / psf.sum()
    return out


def gaia_cone(ra0: float, dec0: float, radius_deg: float, gmax: float = 15.0) -> list:
    """Gaia DR3 sources in a cone: [(ra, dec, pmra, pmdec, G, bp_rp)];
    the archive TAP (sync), three attempts."""
    import requests
    q = ("SELECT ra, dec, pmra, pmdec, phot_g_mean_mag, bp_rp FROM gaiadr3.gaia_source WHERE "
         f"1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra0:.6f}, {dec0:.6f}, {radius_deg:.4f})) "
         f"AND phot_g_mean_mag < {gmax}")
    for attempt in range(3):
        try:
            r = requests.get("https://gea.esac.esa.int/tap-server/tap/sync",
                             params={"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": q},
                             timeout=180)
            r.raise_for_status()
            out = []
            for line in r.text.splitlines()[1:]:
                p = line.split(",")
                if len(p) < 6:
                    continue
                out.append(tuple(float(x) if x not in ("", "null") else np.nan for x in p[:6]))
            return out
        except Exception:
            time.sleep(3 + 5 * attempt)
    raise RuntimeError("gaia query failed")


def gaia_to_v_bv(g: float, bp_rp: float) -> tuple[float, float]:
    """Gaia EDR3 -> Johnson: G - V = -0.02704 + 0.01424 c - 0.2156 c^2 + 0.01426 c^3
    (Riello+2021); B-V ~ 0.75 c - 0.10 (rough, solar 0.82 -> 0.52..0.65)."""
    c = bp_rp if np.isfinite(bp_rp) else 0.82
    v = g - (-0.02704 + 0.01424 * c - 0.2156 * c ** 2 + 0.01426 * c ** 3)
    return float(v), float(np.clip(0.75 * c - 0.10, -0.3, 1.8))


def template_stars_gaia(ra0: float, dec0: float, mjd: float, gaia: list,
                        max_sep_px: float = AP_R + 2.0, extra: list | None = None,
                        exclude_sep_arcsec: float = 0.0) -> list:
    """[(V, B-V, sep_px)] of the Gaia stars (proper motion propagated from
    J2016.0 to `mjd`) within max_sep_px of the patch centre; `extra` =
    explicit (V, B-V, sep_px) rows (the S2 target star), with catalogue
    stars within `exclude_sep_arcsec` of the centre dropped."""
    v0 = np.array([np.cos(np.radians(dec0)) * np.cos(np.radians(ra0)),
                   np.cos(np.radians(dec0)) * np.sin(np.radians(ra0)), np.sin(np.radians(dec0))])
    dyr = (mjd - 57388.5) / 365.25          # J2016.0 = MJD 57388.5 (approx.)
    out = []
    for ra, dec, pmra, pmdec, g, c in gaia:
        if np.isfinite(pmra):
            ra = ra + pmra * dyr / 3.6e6 / np.cos(np.radians(dec))
            dec = dec + pmdec * dyr / 3.6e6
        sv = np.array([np.cos(np.radians(dec)) * np.cos(np.radians(ra)),
                       np.cos(np.radians(dec)) * np.sin(np.radians(ra)), np.sin(np.radians(dec))])
        sep_as = 2 * np.degrees(np.arcsin(np.linalg.norm(sv - v0) / 2)) * 3600
        if sep_as / PIX_ARCSEC > max_sep_px:
            continue
        if extra and sep_as < exclude_sep_arcsec:
            continue
        v, bv = gaia_to_v_bv(g, c)
        out.append((v, bv, float(sep_as / PIX_ARCSEC)))
    if extra:
        out.extend(tuple(map(float, e)) for e in extra)
    return out


def template_flux(stars: list, colour_coeff: float = COLOUR_COEFF) -> float:
    if not stars:
        return 0.0
    s = np.array(stars, float)
    v_inst = s[:, 0] - colour_coeff * (s[:, 1] - BV_REF)
    f = 10 ** (-0.4 * (v_inst - ZP_REF))
    return float(np.sum(f * encircled_fraction(s[:, 2])))
