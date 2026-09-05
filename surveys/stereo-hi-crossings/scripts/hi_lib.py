"""STEREO-A HI-1 level-1 measurement chain (hypotheses freeze D3).

Per frame: fetch (RAL CGI, SSC L2 fallback; size-checked) -> load ->
frame gate (synoptic 30-sum, NMISSING 0, 1024^2) -> header celestial
WCS ('A' system) verified against Hipparcos (translation refinement,
>= 50 stars, rms <= 0.6 px) -> 31-px median high-pass -> 2-D
colour-corrected star ZP (quadratic in x, y) -> top-hat
aperture photometry at ICRS-fixed patches with a local-annulus
background.

Flux units: DN/s normalized to the frame's ZP at the patch position,
i.e. F_norm = F_DN * 10^(-0.4 (ZP(x,y) - ZP_REF)) with ZP_REF = 11.10
(the recon frame's B-V = 0.65 zero point), so 1 unit = the DN/s of a
solar-colour V = ZP_REF ... star; magnitudes V_eq = ZP_REF - 2.5 log10 F_norm.
"""

from __future__ import annotations

import hashlib
import json
import time
import urllib.error
import urllib.parse
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
FRAMES = REPO / "runs" / "stereo-hi-crossings" / "frames"
STARCAT = REPO / "runs" / "heliospheric-crossings" / "starcat" / "hip_v95.npz"
TYCHO = REPO / "runs" / "heliospheric-crossings" / "starcat" / "tyc2_v11.npz"

RAL_CGI = "https://www.stereo.rl.ac.uk/cgi-bin/data.py"
RAL_L1 = "lz/L1/a/img/hi_1"
SSC_L2 = "https://stereo-ssc.nascom.nasa.gov/data/ins_data/secchi_hi/L2/a/img/hi_1/"
FRAME_BYTES = 4_219_200

AP_R, ANN_IN, ANN_OUT = 2.5, 5.0, 9.0
HP_SIZE = 31
ZP_REF = 11.10
COLOUR_COEFF = 0.595           # mag/(B-V), recon; re-measured at dev
BV_REF = 0.65
MIN_STARS, MAX_RMS_PX = 50, 0.6

_cat = {k: np.array(v) for k, v in np.load(STARCAT).items()}   # eager: forked workers must not share a lazy npz handle


# ----------------------------------------------------------------- fetch

def _sha(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def fetch_frame(day: str, stem: str, dest_dir: Path = FRAMES,
                prefer: str = "ral", fallback: bool = True) -> tuple[Path | None, str]:
    """Fetch `<stem>_14h1A.fts` for day YYYYMMDD. Returns (path, source).
    Level-1 from the RAL CGI first; the SSC level-2 `24h1A_br01` frame as
    the declared fallback (recorded in the filename)."""
    dest_dir = dest_dir / day
    dest_dir.mkdir(parents=True, exist_ok=True)
    p1 = dest_dir / f"{stem}_14h1A.fts"
    p2 = dest_dir / f"{stem}_24h1A_br01.fts"
    for p, src in ((p1, "ral"), (p2, "ssc")):
        if p.exists() and p.stat().st_size == FRAME_BYTES:
            return p, src + "_cached"
    order = ["ral", "ssc"] if prefer == "ral" else ["ssc", "ral"]
    if not fallback:
        order = order[:1]
    for src in order:
        for attempt in range(3):
            try:
                if src == "ral":
                    data = urllib.parse.urlencode(
                        {"source": "SECCHI", "target": f"{RAL_L1}/{day}/{p1.name}"}).encode()
                    raw = urllib.request.urlopen(urllib.request.Request(RAL_CGI, data=data),
                                                 timeout=180).read()
                    p = p1
                else:
                    raw = urllib.request.urlopen(SSC_L2 + f"{day}/{p2.name}", timeout=180).read()
                    p = p2
                if len(raw) != FRAME_BYTES or not raw.startswith(b"SIMPLE"):
                    raise IOError(f"bad payload {len(raw)} bytes")
                p.write_bytes(raw)
                return p, src
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    break
                time.sleep(3 + 3 * attempt)
            except Exception:
                time.sleep(3 + 3 * attempt)
    return None, "failed"


# ----------------------------------------------------------------- frame

@dataclass
class Frame:
    path: str
    mjd_avg: float
    exptime: float
    n_images: int
    nmissing: float
    level: str                      # "L1" or "L2"
    wcs: WCS
    raw: np.ndarray                 # as delivered (NaN saturation columns)
    img: np.ndarray                 # high-passed
    hpln_c: float = 0.0
    hplt_c: float = 0.0
    translation: np.ndarray = field(default_factory=lambda: np.zeros(2))
    n_match: int = 0
    astrom_rms: float = np.inf
    zp_coef: np.ndarray | None = None   # quadratic in (x, y) scaled to [-1, 1]
    zp_scatter: float = np.inf
    n_calib: int = 0
    calibrators: list = field(default_factory=list)

    def zp_at(self, x, y) -> float:
        if self.zp_coef is None:
            return float("nan")
        u, v = (x - 511.5) / 511.5, (y - 511.5) / 511.5
        return float(np.array([1, u, v, u * u, u * v, v * v]) @ self.zp_coef)


def frame_gate(h: fits.Header) -> str | None:
    if h.get("NAXIS1") != 1024 or h.get("NAXIS2") != 1024:
        return f"shape_{h.get('NAXIS1')}x{h.get('NAXIS2')}"
    if int(h.get("N_IMAGES", 0)) != 30:
        return f"n_images_{h.get('N_IMAGES')}"
    if abs(float(h.get("EXPTIME", 0)) - 1200.0) > 60.0:
        return f"exptime_{h.get('EXPTIME')}"
    if float(h.get("NMISSING", 0) or 0) != 0:
        return f"nmissing_{h.get('NMISSING')}"
    if "CRVAL1A" not in h:
        return "no_celestial_wcs"
    return None


def load_frame(path: Path) -> Frame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        hdu = fits.open(path)[0]
        h = hdu.header
        wcs = WCS(h, key="A")
    raw = np.asarray(hdu.data, np.float64)
    fill = np.nan_to_num(raw, nan=float(np.nanmedian(raw)))
    img = raw - ndimage.median_filter(fill, size=HP_SIZE)
    return Frame(path=str(path), mjd_avg=Time(h["DATE-AVG"], scale="utc").mjd,
                 exptime=float(h["EXPTIME"]), n_images=int(h.get("N_IMAGES", 0)),
                 nmissing=float(h.get("NMISSING", 0) or 0),
                 level="L2" if "24h1A" in Path(path).name else "L1",
                 wcs=wcs, raw=raw, img=img,
                 hpln_c=float(h["CRVAL1"]), hplt_c=float(h["CRVAL2"]))


# ----------------------------------------------------------------- photometry

def aper_flux(img: np.ndarray, x: float, y: float) -> tuple[float, float, float]:
    """(flux, err, valid_fraction): Gaussian-weighted aperture sum with a
    robust local-annulus background; NaNs excluded from the annulus and
    counted against validity in the aperture."""
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
    noise = float(1.4826 * np.median(np.abs(ann - bg))) + 1e-12
    a = rr <= AP_R
    vals = sub[a]
    valid = float(np.mean(np.isfinite(vals)))
    # top-hat aperture (the PSF is undersampled, FWHM 1.7 px: a plain sum
    # is insensitive to sub-pixel position; recon scatter 0.08 mag vs
    # 0.5 mag for a narrow Gaussian weight); NaN pixels are filled with
    # the background and counted against validity
    fin = np.isfinite(vals)
    if fin.sum() == 0:
        return np.nan, np.nan, 0.0
    flux = float(np.sum(vals[fin] - bg))
    err = float(noise * np.sqrt(a.sum()))
    return flux, err, valid


def sky_to_pix(fr: Frame, ra, dec) -> tuple[np.ndarray, np.ndarray]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        px, py = fr.wcs.all_world2pix(np.atleast_1d(ra), np.atleast_1d(dec), 0)
    return px + fr.translation[0], py + fr.translation[1]


def detect(fr: Frame, snr_min: float = 8.0):
    hp = np.where(np.isfinite(fr.img), fr.img, 0.0)
    sig = 1.4826 * np.median(np.abs(hp[np.isfinite(fr.img)])) + 1e-12
    pk = (ndimage.maximum_filter(hp, size=5) == hp) & (hp > snr_min * sig) & np.isfinite(fr.img)
    pk[:6, :] = pk[-6:, :] = pk[:, :6] = pk[:, -6:] = False
    ys, xs = np.nonzero(pk)
    return np.c_[xs, ys].astype(float)


def fit_frame(fr: Frame) -> bool:
    """Verify the header WCS against Hipparcos, refine the translation,
    fit the 2-D colour-corrected ZP. Returns astrometric validity."""
    dyr = (fr.mjd_avg - 48348.5625) / 365.25
    ra = _cat["ra"] + _cat["pmra_masyr"] * dyr / 3.6e6 / np.cos(np.radians(_cat["dec"]))
    dec = _cat["dec"] + _cat["pmde_masyr"] * dyr / 3.6e6
    fr.translation = np.zeros(2)
    px, py = sky_to_pix(fr, ra, dec)
    inside = np.isfinite(px) & (px > 10) & (px < 1013) & (py > 10) & (py < 1013)
    if inside.sum() < MIN_STARS:
        return False
    px, py, vmag, bv = px[inside], py[inside], _cat["vmag"][inside], _cat["bv"][inside]
    det = detect(fr)
    if len(det) < MIN_STARS:
        return False
    tree = cKDTree(det)
    d, i = tree.query(np.c_[px, py], distance_upper_bound=3.0)
    ok = np.isfinite(d)
    if ok.sum() < MIN_STARS:
        return False
    off = np.median(det[i[ok]] - np.c_[px[ok], py[ok]], axis=0)
    fr.translation = off
    d, i = tree.query(np.c_[px + off[0], py + off[1]], distance_upper_bound=2.0)
    ok = np.isfinite(d)
    fr.n_match = int(ok.sum())
    fr.astrom_rms = float(np.sqrt(np.mean(d[ok] ** 2))) if ok.any() else np.inf
    if fr.n_match < MIN_STARS or fr.astrom_rms > MAX_RMS_PX:
        return False
    # 2-D ZP from colour-cut, unsaturated, detected calibrators
    rows = []
    for k in np.nonzero(ok)[0]:
        if not (-0.3 <= bv[k] <= 1.6) or vmag[k] < 4.5:
            continue
        x, y = px[k] + off[0], py[k] + off[1]      # predicted, not detected
        f, e, val = aper_flux(fr.img, x, y)
        if np.isfinite(f) and f > 0 and e > 0 and f / e >= 10 and val >= 0.99:
            # instrumental magnitude of a star of colour B-V is
            # V - c (B-V - 0.65): red stars are brighter in the 630-730 nm band
            zp = vmag[k] - COLOUR_COEFF * (bv[k] - BV_REF) + 2.5 * np.log10(f)
            rows.append((x, y, zp, vmag[k], bv[k], f))
    if len(rows) < 100:
        fr.zp_coef = None
        return True
    R = np.array(rows)
    u, v = (R[:, 0] - 511.5) / 511.5, (R[:, 1] - 511.5) / 511.5
    A = np.c_[np.ones(len(R)), u, v, u * u, u * v, v * v]
    coef, *_ = np.linalg.lstsq(A, R[:, 2], rcond=None)
    resid = R[:, 2] - A @ coef
    mad = 1.4826 * np.median(np.abs(resid))
    keep = np.abs(resid) < 3 * max(mad, 0.05)
    if keep.sum() >= 100:
        coef, *_ = np.linalg.lstsq(A[keep], R[keep, 2], rcond=None)
        resid = R[:, 2] - A @ coef
        mad = 1.4826 * np.median(np.abs(resid[keep]))
    fr.zp_coef, fr.zp_scatter, fr.n_calib = coef, float(mad), int(keep.sum())
    fr.calibrators = [(float(a), float(b), float(c), float(dd), float(ee), float(ff))
                      for (a, b, c, dd, ee, ff), kp in zip(R, keep) if kp]
    return True


def measure_patch(fr: Frame, ra: float, dec: float) -> dict:
    """Photometry at an ICRS position; flux normalized to ZP_REF."""
    x, y = sky_to_pix(fr, ra, dec)
    x, y = float(x[0]), float(y[0])
    f, e, val = aper_flux(fr.img, x, y)
    zp = fr.zp_at(x, y)
    scale = 10 ** (-0.4 * (zp - ZP_REF)) if np.isfinite(zp) else np.nan
    return {"mjd": fr.mjd_avg, "x": round(x, 2), "y": round(y, 2),
            "flux": f * scale if np.isfinite(f) else np.nan,
            "err": e * scale if np.isfinite(e) else np.nan,
            "flux_dn": f, "valid": round(val, 3), "zp": round(zp, 4) if np.isfinite(zp) else None}


def in_footprint(x: float, y: float, margin_px: float = 8.0) -> bool:
    return margin_px <= x <= 1023 - margin_px and margin_px <= y <= 1023 - margin_px
