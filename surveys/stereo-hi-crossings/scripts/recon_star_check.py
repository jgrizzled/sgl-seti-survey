"""HI-1A recon star check (the go/no-go astrometry + photometry test).

On the recon level-1 frame (2010-06-15 00:09 UT, `14h1A` DN/s product)
and its level-2 `24h1A_br01` counterpart:

* header celestial WCS ('A' system, RA/DEC-AZP) vs Hipparcos: match
  fraction, residual rms, mean offset;
* per-star zero points V = ZP - 2.5 log10(DN/s) with a B-V colour
  term (HI-1 is a 630-730 nm band);
* PSF width from bright unsaturated stars;
* single-frame noise and 5-sigma point-source depth as a function of
  helioprojective longitude (the in-beam arc is HPLN -6..-4 deg, the
  brightest F-corona in the field);
* L1 vs L2 aperture-flux parity on matched stars.

Output: results/recon_star_check_v1.json
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from scipy import ndimage
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))

REPO = Path(__file__).resolve().parents[3]
RECON = REPO / "runs" / "stereo-hi-crossings" / "recon"
L1 = RECON / "20100615_000901_14h1A.fts"
L2 = RECON / "20100615_000901_24h1A_br01.fts"
CAT = REPO / "runs" / "heliospheric-crossings" / "starcat" / "hip_v95.npz"
OUT = REPO / "surveys" / "stereo-hi-crossings" / "results" / "recon_star_check_v1.json"

AP_R, ANN = 2.5, (5.0, 9.0)


def load(path):
    warnings.simplefilter("ignore")
    h = fits.open(path)[0]
    return h.header, np.asarray(h.data, float)


def highpass(img, size=31):
    fill = np.nan_to_num(img, nan=np.nanmedian(img))
    return img - ndimage.median_filter(fill, size=size)


def aper(img, x, y):
    ix, iy = int(round(x)), int(round(y))
    r = int(np.ceil(ANN[1]))
    if not (r < ix < img.shape[1] - r and r < iy < img.shape[0] - r):
        return np.nan, np.nan, 0.0
    sub = img[iy - r: iy + r + 1, ix - r: ix + r + 1]
    yy, xx = np.mgrid[-r:r + 1, -r:r + 1]
    rr = np.hypot(xx - (x - ix), yy - (y - iy))
    ann = sub[(rr >= ANN[0]) & (rr <= ANN[1])]
    ann = ann[np.isfinite(ann)]
    if len(ann) < 20:
        return np.nan, np.nan, 0.0
    bg = np.median(ann)
    noise = 1.4826 * np.median(np.abs(ann - bg))
    a = rr <= AP_R
    vals = sub[a]
    valid = float(np.mean(np.isfinite(vals)))
    flux = float(np.nansum(vals - bg))
    return flux, float(noise * np.sqrt(a.sum())), valid


def main():
    h1, im1 = load(L1)
    h2, im2 = load(L2)
    wa, wp = WCS(h1, key="A"), WCS(h1)
    cat = np.load(CAT)
    ep = Time(h1["DATE-AVG"]).jyear
    ra = cat["ra"] + cat["pmra_masyr"] / 3.6e6 * (ep - 1991.25) / np.cos(np.radians(cat["dec"]))
    dec = cat["dec"] + cat["pmde_masyr"] / 3.6e6 * (ep - 1991.25)
    px, py = wa.all_world2pix(ra, dec, 0)
    inf = np.isfinite(px) & (px > 10) & (px < 1013) & (py > 10) & (py < 1013)
    vmag, bv = cat["vmag"][inf], cat["bv"][inf]
    px, py = px[inf], py[inf]

    # detections on the high-passed L1 frame
    hp = highpass(im1)
    bad = ~np.isfinite(im1)
    hp[bad] = 0
    sig = 1.4826 * np.nanmedian(np.abs(hp))
    pk = (ndimage.maximum_filter(hp, size=5) == hp) & (hp > 8 * sig) & ~bad
    ys, xs = np.nonzero(pk)
    d, i = cKDTree(np.c_[xs, ys]).query(np.c_[px, py], distance_upper_bound=3.0)
    ok = np.isfinite(d)
    dx, dy = xs[i[ok]] - px[ok], ys[i[ok]] - py[ok]
    astrom = {"hip_in_fov": int(inf.sum()), "matched_3px": int(ok.sum()),
              "median_resid_px": float(np.median(d[ok])),
              "mean_offset_px": [float(dx.mean()), float(dy.mean())],
              "rms_px": [float(dx.std()), float(dy.std())]}

    # photometry on matched stars (centroid = detected peak, refined)
    rows = []
    for k in np.nonzero(ok)[0]:
        x, y = float(xs[i[k]]), float(ys[i[k]])
        f1, e1, v1 = aper(hp, x, y)
        f2, e2, v2 = aper(highpass(im2) if k == 0 else HP2, x, y) if False else (np.nan, np.nan, 0)
        rows.append((vmag[k], bv[k], f1, e1, v1, x, y))
    HP2 = highpass(im2)
    rows2 = [aper(HP2, r[5], r[6]) for r in rows]
    R = np.array([(r[0], r[1], r[2], r[3], r[4], r[5], r[6], q[0], q[1]) for r, q in zip(rows, rows2)])
    good = (R[:, 2] > 0) & (R[:, 4] >= 0.99) & np.isfinite(R[:, 1]) & (R[:, 1] < 90) & (R[:, 0] > 4.5)
    zp = R[good, 0] + 2.5 * np.log10(R[good, 2])
    # colour term: zp = zp0 + c (B-V - 0.65)
    A = np.c_[np.ones(good.sum()), R[good, 1] - 0.65]
    hi = R[good, 2] / R[good, 3] > 10
    coef, *_ = np.linalg.lstsq(A[hi], zp[hi], rcond=None)
    resid = zp[hi] - A[hi] @ coef
    photom = {"n_calibrators_snr10": int(hi.sum()),
              "zp0_at_bv0.65": float(coef[0]), "colour_coeff_mag_per_bv": float(coef[1]),
              "scatter_mad_no_colour": float(1.4826 * np.median(np.abs(zp[hi] - np.median(zp[hi])))),
              "scatter_mad_with_colour": float(1.4826 * np.median(np.abs(resid))),
              "zp_by_vbin": {}}
    for lo in range(5, 10):
        m = hi & (R[good, 0] >= lo) & (R[good, 0] < lo + 1)
        if m.sum() >= 5:
            photom["zp_by_vbin"][f"{lo}-{lo+1}"] = {"n": int(m.sum()), "zp_median": float(np.median(zp[m])),
                                                    "resid_median": float(np.median((zp - A @ coef)[m]))}
    # L1 vs L2 parity
    f1, f2 = R[good, 2], R[good, 7]
    pm = hi & (f2 > 0)
    ratio = f2[pm] / f1[pm]
    parity = {"n": int(pm.sum()), "median_L2_over_L1": float(np.median(ratio)),
              "mad": float(1.4826 * np.median(np.abs(ratio - np.median(ratio))))}

    # PSF width: second moments of bright unsaturated stars (V 5-7)
    fw = []
    for r in R[good & (R[:, 0] < 7.0) & (R[:, 0] > 5.0)]:
        x, y = int(round(r[5])), int(round(r[6]))
        sub = hp[y - 4:y + 5, x - 4:x + 5]
        if sub.shape != (9, 9) or not np.all(np.isfinite(sub)):
            continue
        s = np.clip(sub, 0, None)
        yy, xx = np.mgrid[-4:5, -4:5]
        w = s / s.sum()
        cx, cy = (w * xx).sum(), (w * yy).sum()
        sx = np.sqrt((w * (xx - cx) ** 2).sum())
        sy = np.sqrt((w * (yy - cy) ** 2).sum())
        fw.append(2.355 * 0.5 * (sx + sy))
    psf = {"n_stars": len(fw), "fwhm_px_median": float(np.median(fw)) if fw else None,
           "fwhm_arcsec_median": float(np.median(fw) * 72.0) if fw else None}

    # noise vs helioprojective longitude (columns -> HPLN via primary WCS)
    strips = {}
    for x0 in range(0, 1024, 64):
        hpln = float(wp.all_pix2world([[x0 + 32, 511.5]], 0)[0][0])
        s = hp[:, x0:x0 + 64]
        s = s[np.isfinite(s)]
        sg = 1.4826 * np.median(np.abs(s - np.median(s)))
        bg = float(np.nanmedian(im1[:, x0:x0 + 64]))
        depth = coef[0] - 2.5 * np.log10(5 * sg * np.sqrt(np.pi * AP_R ** 2))
        strips[f"{hpln:.1f}"] = {"bg_dn_s": bg, "hp_sigma_px": float(sg),
                                 "depth_5sig_V_solar_colour": float(depth)}

    out = {"frames": {"l1": L1.name, "l2": L2.name, "date_avg": h1["DATE-AVG"],
                      "exptime_s": float(h1["EXPTIME"]), "n_images": int(h1["N_IMAGES"]),
                      "bunit": h1["BUNIT"], "cdelt_deg": float(h1["CDELT1"]),
                      "nan_pixels_l1": int(bad.sum())},
           "astrometry_header_wcs_vs_hipparcos": astrom,
           "photometry": photom, "l1_l2_parity": parity, "psf": psf,
           "noise_vs_hpln": strips,
           "notes": ["aperture r=2.5 px, annulus 5-9 px on a 31-px median high-pass",
                     "depth = 5 sigma on the aperture noise for a B-V=0.65 source, single 40-min frame",
                     "colour coefficient sign: positive = red stars brighter in DN/s than V predicts"]}
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: out[k] for k in ("astrometry_header_wcs_vs_hipparcos", "photometry",
                                          "l1_l2_parity", "psf")}, indent=1))
    for k, v in strips.items():
        print(k, v)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
