"""LASCO level-0.5 measurement chain (threshold freeze v1.0).

Per-frame: load -> azimuthal-median radial-profile detrend (0.1 Rsun
bins) -> synthetic celestial WCS (Sun-from-SOHO + Meeus P-angle +
header CROTA, both roll hypotheses) -> Hipparcos star fit (translation
+ roll resolution + star ZP) -> Gaussian-weighted forced photometry at
ICRS positions with local-annulus background.

Conventions validated in the recon star check: with rotation
rot = P + CROTA (+180 for the flipped roll state), image
x = CRPIX1-1 - E'/scale, y = CRPIX2-1 + N'/scale, where (E', N') are
the celestial East/North offsets from Sun center rotated by rot.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.coordinates import get_body_barycentric
from astropy.io import fits
from astropy.time import Time
from scipy.ndimage import maximum_filter, median_filter

REPO = Path(__file__).resolve().parents[3]
SOHO_TAB = REPO / "crossings" / "observers" / "soho_sc_ephemeris.npz"
STARCAT = REPO / "runs" / "heliospheric-crossings" / "starcat" / "hip_v95.npz"

RSUN_KM = 695_700.0
AU_KM = 1.495978707e8
APERTURE_PX = 2.5
ANNULUS_PX = (5.0, 9.0)
GAUSS_SIGMA_PX = 1.2

_soho = np.load(SOHO_TAB)
_cat = np.load(STARCAT)


# ----------------------------------------------------------------- geometry

def soho_xyz_au(mjd: float) -> np.ndarray:
    return np.array([np.interp(mjd, _soho["mjd_utc"], _soho["xyz_au"][:, i])
                     for i in range(3)])


def sun_from_soho(mjd: float) -> tuple[float, float, float, float]:
    """(ra, dec, dist_au, rsun_arcsec) of the Sun as seen from SOHO."""
    t = Time(mjd, format="mjd", scale="utc")
    sun = get_body_barycentric("sun", t).xyz.to_value("AU")
    d = sun - soho_xyz_au(mjd)
    r = float(np.linalg.norm(d))
    ra = float(np.degrees(np.arctan2(d[1], d[0])) % 360.0)
    dec = float(np.degrees(np.arcsin(d[2] / r)))
    rsun_arcsec = float(np.degrees(np.arcsin(RSUN_KM / (r * AU_KM))) * 3600.0)
    return ra, dec, r, rsun_arcsec


def solar_p_angle_deg(mjd: float) -> float:
    jd = Time(mjd, format="mjd", scale="utc").tt.jd
    d = jd - 2451545.0
    eps = np.radians(23.4393 - 3.563e-7 * d)
    incl = np.radians(7.25)
    omega = np.radians(73.6667 + 1.3958333 * (jd - 2396758.0) / 36525.0)
    g = np.radians((357.529 + 0.98560028 * d) % 360.0)
    q = 280.459 + 0.98564736 * d
    lam = np.radians((q + 1.915 * np.sin(g) + 0.020 * np.sin(2 * g)) % 360.0)
    return float(np.degrees(np.arctan(-np.cos(lam) * np.tan(eps))
                            + np.arctan(-np.cos(lam - omega) * np.tan(incl))))


# ----------------------------------------------------------------- frame

@dataclass
class Frame:
    path: str
    camera: str
    mjd_mid: float
    exptime: float
    crpix: np.ndarray          # 0-based Sun-center pixel (header)
    scale_arcsec: float
    crota_deg: float
    img: np.ndarray            # bias-subtracted, radial-profile-detrended
    raw: np.ndarray            # bias-subtracted only
    sun_ra: float = 0.0
    sun_dec: float = 0.0
    rsun_arcsec: float = 954.0
    p_angle: float = 0.0
    # star-fit products
    rot_deg: float | None = None       # resolved rotation (P + CROTA [+180])
    translation: np.ndarray | None = None
    # radius-resolved star ZP (level 0.5 is not flat-fielded: vignetting
    # makes the ZP radial; calibrators are colour-cut 0.2<=B-V<=1.2)
    zp_r: np.ndarray | None = None     # bin centres, Rsun
    zp_v: np.ndarray | None = None     # ZP per bin
    zp_scatter: float = np.inf         # MAD about the radial fit
    n_match: int = 0

    def zp_at(self, r_rsun: float) -> float:
        if self.zp_r is None:
            return float("nan")   # astrometry without ZP (C2 norm, L1)
        return float(np.interp(r_rsun, self.zp_r, self.zp_v))


SYNOPTIC_FILTER = {"c2": "Orange", "c3": "Clear"}


def frame_usable(path: Path) -> str | None:
    """Header-only gate: returns a rejection reason or None (usable).
    Synoptic full-frame unpolarized frames only — polarizer triplets
    attenuate flux and 512x512 binned subframes break the pixel chain."""
    h = fits.getheader(path)
    cam = str(h.get("DETECTOR", "")).strip().lower()
    if h.get("NAXIS1") != 1024 or h.get("NAXIS2") != 1024:
        return f"shape_{h.get('NAXIS1')}x{h.get('NAXIS2')}"
    if str(h.get("FILTER", "")).strip() != SYNOPTIC_FILTER.get(cam):
        return f"filter_{str(h.get('FILTER','')).strip()}"
    if str(h.get("POLAR", "")).strip() != "Clear":
        return f"polar_{str(h.get('POLAR','')).strip()}"
    return None


def load_frame(path: Path) -> Frame:
    hdu = fits.open(path)[0]
    h = hdu.header
    img = np.asarray(hdu.data, float)
    # midpoint time: MID_DATE (MJD) + MID_TIME (s of day); fall back to
    # DATE-OBS/TIME-OBS + EXPTIME/2
    if "MID_DATE" in h and "MID_TIME" in h and h["MID_DATE"] > 0:
        mjd = float(h["MID_DATE"]) + float(h["MID_TIME"]) / 86400.0
    else:
        date = str(h["DATE-OBS"]).replace("/", "-")
        mjd = Time(f"{date}T{h['TIME-OBS']}", scale="utc").mjd \
            + float(h.get("EXPTIME", 0)) / 2 / 86400.0
    img = img - float(h.get("OFFSET", 0.0))
    fr = Frame(
        path=str(path), camera=str(h["DETECTOR"]).strip().lower(),
        mjd_mid=mjd, exptime=float(h.get("EXPTIME", 0)),
        crpix=np.array([float(h["CRPIX1"]) - 1, float(h["CRPIX2"]) - 1]),
        scale_arcsec=float(h.get("PLATESCL") or h["CDELT1"]),
        crota_deg=float(h.get("CROTA1", h.get("CROTA", 0.0))),
        img=img, raw=img.copy(),
    )
    fr.sun_ra, fr.sun_dec, _, fr.rsun_arcsec = sun_from_soho(mjd)
    fr.p_angle = solar_p_angle_deg(mjd)
    fr.img = radial_profile_subtract(fr)
    return fr


def radial_profile_subtract(fr: Frame, bin_rsun: float = 0.1) -> np.ndarray:
    """Azimuthal-median radial profile subtraction (sort-grouped medians
    — a per-bin boolean-mask loop costs ~3 s/frame at 0.1 Rsun bins)."""
    ny, nx = fr.raw.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    rsun_px = fr.rsun_arcsec / fr.scale_arcsec
    r = np.hypot(xx - fr.crpix[0], yy - fr.crpix[1]) / rsun_px
    nb = int(np.ceil(r.max() / bin_rsun)) + 1
    idx = np.minimum((r / bin_rsun).astype(np.int32), nb - 1).ravel()
    flat = fr.raw.ravel()
    order = np.argsort(idx, kind="stable")
    sorted_idx, sorted_val = idx[order], flat[order]
    bounds = np.searchsorted(sorted_idx, np.arange(nb + 1))
    prof = np.zeros(nb)
    for b in range(nb):
        lo, hi = bounds[b], bounds[b + 1]
        if hi > lo:
            prof[b] = np.median(sorted_val[lo:hi])
    return fr.raw - prof[idx.reshape(ny, nx)]


# ----------------------------------------------------------------- star fit

def _sky_to_pix(fr: Frame, ra, dec, rot_deg: float) -> tuple[np.ndarray, np.ndarray]:
    d_e = (np.asarray(ra) - fr.sun_ra) * np.cos(np.radians(np.asarray(dec))) * 3600.0
    d_n = (np.asarray(dec) - fr.sun_dec) * 3600.0
    rot = np.radians(rot_deg)
    e = d_e * np.cos(rot) - d_n * np.sin(rot)
    n = d_e * np.sin(rot) + d_n * np.cos(rot)
    return fr.crpix[0] - e / fr.scale_arcsec, fr.crpix[1] + n / fr.scale_arcsec


def detect_sources(fr: Frame, snr_min: float = 8.0, cap: int = 500):
    hp = fr.img - median_filter(fr.img, size=11)
    sigma = np.median(np.abs(hp)) * 1.4826 + 1e-30
    snr = hp / sigma
    pk = (maximum_filter(snr, size=7) == snr) & (snr > snr_min)
    pk[:16, :] = pk[-16:, :] = pk[:, :16] = pk[:, -16:] = False
    ys, xs = np.nonzero(pk)
    order = np.argsort(snr[ys, xs])[::-1][:cap]
    return np.stack([xs[order], ys[order]], 1).astype(float), snr[ys[order], xs[order]]


def _aper_flux(img: np.ndarray, x: float, y: float) -> tuple[float, float]:
    """Gaussian-weighted aperture flux + error proxy, local-annulus bg."""
    ny, nx = img.shape
    ix, iy = int(round(x)), int(round(y))
    r_out = int(np.ceil(ANNULUS_PX[1]))
    if not (r_out < ix < nx - r_out and r_out < iy < ny - r_out):
        return np.nan, np.nan
    sub = img[iy - r_out: iy + r_out + 1, ix - r_out: ix + r_out + 1]
    yy, xx = np.mgrid[-r_out:r_out + 1, -r_out:r_out + 1]
    rr = np.hypot(xx - (x - ix), yy - (y - iy))
    ann = sub[(rr >= ANNULUS_PX[0]) & (rr <= ANNULUS_PX[1])]
    bg = np.median(ann)
    noise = np.median(np.abs(ann - bg)) * 1.4826 + 1e-30
    w = np.exp(-0.5 * (rr / GAUSS_SIGMA_PX) ** 2) * (rr <= APERTURE_PX)
    flux = float(np.sum(w * (sub - bg)) / np.sum(w ** 2) * np.sum(w))
    err = float(noise * np.sqrt(np.sum(w ** 2)) / np.sum(w ** 2) * np.sum(w))
    return flux, err


def aperture_valid_fraction(fr: Frame, x: float, y: float) -> float:
    ny, nx = fr.raw.shape
    ix, iy = int(round(x)), int(round(y))
    r = int(np.ceil(APERTURE_PX))
    if not (r < ix < nx - r and r < iy < ny - r):
        return 0.0
    sub = fr.raw[iy - r: iy + r + 1, ix - r: ix + r + 1]
    return float(np.mean(np.isfinite(sub) & (sub != 0)))


def fit_frame(fr: Frame, min_match: int | None = None, tol_px: float = 4.0,
              search_px: float = 40.0) -> bool:
    """Resolve roll, translation, and star ZP against Hipparcos.

    Astrometric minimum: 2 stars on C2 (finding L1 — the inner corona
    leaves only the 2-3 brightest catalog stars measurable; a 2-star
    match against a predicted pattern is unambiguous), 5 on C3.
    The ZP is optional (statistics are unit-free z-scores; ZP feeds
    depths later): fr.zp_r stays None when too few calibrators.
    Returns True when the frame passes the frozen photometric gate
    (scatter <= 0.2 with >= 5 calibrators) — astrometric validity is
    fr.rot_deg is not None, recorded separately."""
    if min_match is None:
        min_match = 2 if fr.camera == "c2" else 5
    fov = fr.raw.shape[0] / 2 * fr.scale_arcsec / 3600.0 * 1.5  # deg
    dyr = (fr.mjd_mid - 48348.5625) / 365.25          # from 1991.25
    ra = _cat["ra"] + _cat["pmra_masyr"] * dyr / 3.6e6 / np.cos(np.radians(_cat["dec"]))
    dec = _cat["dec"] + _cat["pmde_masyr"] * dyr / 3.6e6
    sep = np.hypot((ra - fr.sun_ra) * np.cos(np.radians(dec)), dec - fr.sun_dec)
    m = sep < fov
    ra, dec, vmag, bv = ra[m], dec[m], _cat["vmag"][m], _cat["bv"][m]
    if len(ra) < min_match:
        return False
    det, _ = detect_sources(fr)
    if len(det) < min_match:
        return False

    best = None
    for extra in (0.0, 180.0):
        rot = fr.p_angle + fr.crota_deg + extra
        px, py = _sky_to_pix(fr, ra, dec, rot)
        inside = (px > 8) & (px < fr.raw.shape[1] - 8) & (py > 8) & (py < fr.raw.shape[0] - 8)
        if inside.sum() < min_match:
            continue
        pxi, pyi, vi = px[inside], py[inside], vmag[inside]
        order = np.argsort(vi)[:12]                    # brightest anchors
        for a in order[:6]:
            d = det - np.array([pxi[a], pyi[a]])
            for dx, dy in d[(np.abs(d[:, 0]) < search_px) & (np.abs(d[:, 1]) < search_px)]:
                dd = np.hypot(det[:, 0][None, :] - (pxi[:, None] + dx),
                              det[:, 1][None, :] - (pyi[:, None] + dy))
                nearest = dd.min(1)
                hit = nearest < tol_px
                score = (int(hit.sum()), -float(nearest[hit].sum() if hit.any() else 0))
                if best is None or score > best[0]:
                    best = (score, rot, np.array([dx, dy]), inside, hit)
    if best is None or best[0][0] < min_match:
        return False
    _, rot, tr, inside, hit = best
    # refine translation on the matched set (det reused from above)
    px, py = _sky_to_pix(fr, ra[inside], dec[inside], rot)
    dd = np.hypot(det[:, 0][None, :] - (px[:, None] + tr[0]),
                  det[:, 1][None, :] - (py[:, None] + tr[1]))
    j = dd.argmin(1)
    ok = dd[np.arange(len(px)), j] < tol_px
    tr = tr + np.median(det[j[ok]] - np.stack([px[ok] + tr[0], py[ok] + tr[1]], 1), axis=0)
    fr.rot_deg, fr.translation = float(rot), tr
    # radius-resolved star ZP from colour-cut calibrators
    rsun_px = fr.rsun_arcsec / fr.scale_arcsec
    r_lo, r_hi = (2.2, 6.0) if fr.camera == "c2" else (4.4, 29.0)  # adopted annuli
    zps, rads = [], []
    bvi, vmi = bv[inside], vmag[inside]
    for k in np.nonzero(ok)[0]:
        # calibrator gates: colour-cut (reds sit ~+0.6 in the broad
        # Clear band), unsaturated (V < 4.5 bleeds low), in the usable
        # annulus, and DETECTED (S/N >= 10) — sub-detection stars give
        # pure-noise ZPs and edge vignetting gives absurd ones
        if not (-0.3 <= bvi[k] <= 1.2) or vmi[k] < 4.5:
            continue
        x, y = px[k] + tr[0], py[k] + tr[1]
        r_star = np.hypot(x - fr.crpix[0] - tr[0], y - fr.crpix[1] - tr[1]) / rsun_px
        if not (r_lo <= r_star <= r_hi):
            continue
        fx, er = _aper_flux(fr.img, x, y)
        if np.isfinite(fx) and er > 0 and fx / er >= 10.0:
            zps.append(vmi[k] + 2.5 * np.log10(fx))
            rads.append(r_star)
    fr.n_match = int(len(zps))
    if len(zps) < 3:
        return False  # astrometry stands (rot/translation set); no ZP
    zps, rads = np.array(zps), np.array(rads)
    # one MAD-clip pass against the global median (saturation stragglers)
    med = np.median(zps)
    mad = np.median(np.abs(zps - med)) * 1.4826 + 1e-6
    keep = np.abs(zps - med) < 3 * max(mad, 0.15)
    if keep.sum() >= 3:
        zps, rads = zps[keep], rads[keep]
    fr.zp_r, fr.zp_v = _radial_fit(zps, rads, binw=1.5 if fr.camera == "c2" else 6.0)
    resid = zps - np.interp(rads, fr.zp_r, fr.zp_v)
    fr.zp_scatter = float(np.median(np.abs(resid)) * 1.4826)
    fr.n_match = int(len(zps))
    return fr.zp_scatter <= 0.2 and fr.n_match >= 5


def _radial_fit(zps: np.ndarray, rads: np.ndarray, binw: float):
    centres, vals = [], []
    for lo in np.arange(0, rads.max() + binw, binw):
        mm = (rads >= lo) & (rads < lo + binw)
        if mm.sum() >= 3:
            centres.append(lo + binw / 2)
            vals.append(float(np.median(zps[mm])))
    if len(centres) < 2:
        centres = [float(rads.min()), float(rads.max()) + 1e-3]
        vals = [float(np.median(zps))] * 2
    return np.array(centres), np.array(vals)


def forced_photometry(fr: Frame, ra: float, dec: float) -> dict:
    """Photometry at an ICRS position through the fitted transform."""
    if fr.rot_deg is None:
        raise RuntimeError("frame not fitted")
    x, y = _sky_to_pix(fr, ra, dec, fr.rot_deg)
    x, y = float(x + fr.translation[0]), float(y + fr.translation[1])
    valid = aperture_valid_fraction(fr, x, y)
    flux, err = _aper_flux(fr.img, x, y)
    rsun_px = fr.rsun_arcsec / fr.scale_arcsec
    r_rsun = float(np.hypot(x - fr.crpix[0] - fr.translation[0],
                            y - fr.crpix[1] - fr.translation[1]) / rsun_px)
    return {"mjd": fr.mjd_mid, "x": x, "y": y, "flux": flux, "err": err,
            "valid_fraction": valid, "r_rsun": r_rsun,
            "zp": fr.zp_at(r_rsun), "zp_scatter": fr.zp_scatter,
            "n_match": fr.n_match}
