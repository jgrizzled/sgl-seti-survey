"""WISPR-I recon depth vs heliocentric distance (E22 frames), with the
PSF profile, a matched aperture, bright-calibrator zero points, and the
L3 -> L2-unit rescaling ((r/0.2 AU)^2.3, HISTORY of the L3 product).

Per frame (L2 + L3 of the same stamp):
* PSF radial profile from V < 5 stars (peak-normalised, 0.5-px bins);
* ZP (V = ZP - 2.5 log10 F) with a B-V colour term on V 3-7.5 stars,
  S/N > 10, 3-sigma clipped, aperture r 3 px / annulus 5-9 px on the
  rescaled L3 after a 15-px median high-pass; per-V-bin residuals;
* L2 vs rescaled-L3 flux parity;
* single-frame noise vs elongation on the star-masked high-passed
  rescaled L3 -> 5-sigma depth for a solar-colour source.

Output: results/recon_depth_v1.json
"""
from __future__ import annotations
import glob, json, sys, warnings
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import psp_geometry as G

REPO = Path(__file__).resolve().parents[3]
RECON = REPO / "runs" / "wispr-crossings" / "recon"
CAT = REPO / "runs" / "heliospheric-crossings" / "starcat" / "hip_v95.npz"
OUT = REPO / "surveys" / "wispr-crossings" / "results" / "recon_depth_v1.json"
AU_M = G.AU_KM * 1e3
AP_R, ANN = 3.0, (5.0, 9.0)
warnings.simplefilter("ignore")


def aper(img, x, y):
    ix, iy = int(round(x)), int(round(y))
    r = int(np.ceil(ANN[1]))
    if not (r < ix < img.shape[1] - r and r < iy < img.shape[0] - r):
        return np.nan, np.nan
    sub = img[iy - r:iy + r + 1, ix - r:ix + r + 1]
    yy, xx = np.mgrid[-r:r + 1, -r:r + 1]
    rr = np.hypot(xx - (x - ix), yy - (y - iy))
    ann = sub[(rr >= ANN[0]) & (rr <= ANN[1])]
    ann = ann[np.isfinite(ann)]
    if len(ann) < 20:
        return np.nan, np.nan
    bg = np.median(ann)
    a = rr <= AP_R
    return float(np.nansum(sub[a] - bg)), float(1.4826 * np.median(np.abs(ann - bg)) * np.sqrt(a.sum()))


def highpass(img, size=15):
    return img - ndimage.median_filter(np.nan_to_num(img, nan=np.nanmedian(img)), size=size)


def analyse(f2: Path):
    f3 = Path(str(f2).replace("L2", "L3"))
    h = fits.getheader(f2)
    im2, im3 = fits.getdata(f2).astype(float), fits.getdata(f3).astype(float)
    r_au = h["DSUN_OBS"] / AU_M
    scale = (r_au / 0.2) ** 2.3
    im3u = im3 / scale
    wa, wp = WCS(h, key="A"), WCS(h)
    nx, ny = h["NAXIS1"], h["NAXIS2"]
    cat = np.load(CAT)
    ep = Time(h["DATE-AVG"]).jyear
    ra = cat["ra"] + cat["pmra_masyr"] / 3.6e6 * (ep - 1991.25) / np.cos(np.radians(cat["dec"]))
    dec = cat["dec"] + cat["pmde_masyr"] / 3.6e6 * (ep - 1991.25)
    c_ra, c_dec = wa.all_pix2world([[nx / 2, ny / 2]], 0)[0]
    sep = np.degrees(np.arccos(np.clip(G.radec_to_vec(ra, dec) @ G.radec_to_vec(c_ra, c_dec), -1, 1)))
    near = sep < 35
    px, py = wa.all_world2pix(ra[near], dec[near], 0)
    inf = np.isfinite(px) & (px > 12) & (px < nx - 12) & (py > 12) & (py < ny - 12)
    V, BV = cat["vmag"][near][inf], cat["bv"][near][inf]
    px, py = px[inf], py[inf]
    hp3 = highpass(im3u)
    hp2 = highpass(im2)
    # PSF
    prof = []
    for x, y, v in zip(px, py, V):
        if v > 5.0:
            continue
        ix, iy = int(round(x)), int(round(y))
        sub = hp3[iy - 6:iy + 7, ix - 6:ix + 7]
        if sub.shape != (13, 13) or not np.all(np.isfinite(sub)):
            continue
        j = np.unravel_index(np.argmax(sub), sub.shape)
        if abs(j[0] - 6) > 3 or abs(j[1] - 6) > 3:
            continue
        cy, cx = iy + j[0] - 6, ix + j[1] - 6
        sub = hp3[cy - 6:cy + 7, cx - 6:cx + 7]
        if sub.shape != (13, 13):
            continue
        s = np.clip(sub, 0, None)
        yy, xx = np.mgrid[-6:7, -6:7]
        w = s / s.sum()
        ccx, ccy = (w * xx).sum(), (w * yy).sum()
        rr = np.hypot(xx - ccx, yy - ccy)
        prof.append((rr.ravel(), (sub / sub.max()).ravel(), sub.ravel() / sub.sum()))
    rr = np.concatenate([p[0] for p in prof])
    pp = np.concatenate([p[1] for p in prof])
    ff = np.concatenate([p[2] for p in prof])
    radial = [float(np.median(pp[(rr >= a) & (rr < a + 0.5)])) for a in np.arange(0, 6, 0.5)]
    # encircled flux fraction (ensemble, relative to the 13x13 box total)
    enc = {str(R): float(np.sum(ff[rr <= R]) / len(prof)) for R in (1, 2, 3, 4, 5)}
    # photometry
    F = []
    for x, y, v, bv in zip(px, py, V, BV):
        if v > 7.5 or not np.isfinite(bv) or bv > 90:
            continue
        ix, iy = int(round(x)), int(round(y))
        sub = hp3[iy - 3:iy + 4, ix - 3:ix + 4]
        if sub.shape != (7, 7) or not np.all(np.isfinite(sub)):
            continue
        j = np.unravel_index(np.argmax(sub), sub.shape)
        cx, cy = ix + j[1] - 3, iy + j[0] - 3
        fl, er = aper(hp3, cx, cy)
        fl2, _ = aper(hp2, cx, cy)
        hpln, hplt = wp.all_pix2world([[cx, cy]], 0)[0]
        eps = np.degrees(np.arccos(np.cos(np.radians(hpln)) * np.cos(np.radians(hplt))))
        if fl > 0 and fl / er > 10:
            F.append((v, bv, fl, er, fl2, eps))
    F = np.array(F)
    zp = F[:, 0] + 2.5 * np.log10(F[:, 2])
    A = np.c_[np.ones(len(F)), F[:, 1] - 0.65]
    coef, *_ = np.linalg.lstsq(A, zp, rcond=None)
    res = zp - A @ coef
    keep = np.abs(res - np.median(res)) < 3 * 1.4826 * np.median(np.abs(res - np.median(res)))
    coef, *_ = np.linalg.lstsq(A[keep], zp[keep], rcond=None)
    res = zp - A @ coef
    phot = {"n_calibrators": int(keep.sum()), "zp0_at_bv0.65": float(coef[0]),
            "colour_coeff": float(coef[1]),
            "mad": float(1.4826 * np.median(np.abs(res[keep] - np.median(res[keep])))),
            "l2_over_l3u_flux_median": float(np.median(F[keep, 4] / F[keep, 2])),
            "resid_by_vbin": {f"{lo}-{lo+1}": {"n": int((keep & (F[:, 0] >= lo) & (F[:, 0] < lo + 1)).sum()),
                                             "median": float(np.median(res[keep & (F[:, 0] >= lo) & (F[:, 0] < lo + 1)]))}
                              for lo in (3, 4, 5, 6, 7) if (keep & (F[:, 0] >= lo) & (F[:, 0] < lo + 1)).sum() >= 3},
            "resid_by_eps": {f"{int(lo)}-{int(lo)+10}": {"n": int((keep & (F[:, 5] >= lo) & (F[:, 5] < lo + 10)).sum()),
                                                        "median": float(np.median(res[keep & (F[:, 5] >= lo) & (F[:, 5] < lo + 10)]))}
                             for lo in (10, 20, 30, 40, 50) if (keep & (F[:, 5] >= lo) & (F[:, 5] < lo + 10)).sum() >= 3}}
    # noise vs elongation
    mask = np.zeros(hp3.shape, bool)
    for x, y in zip(px, py):
        ix, iy = int(round(x)), int(round(y))
        mask[max(iy - 4, 0):iy + 5, max(ix - 4, 0):ix + 5] = True
    strips = {}
    for x0 in range(0, nx - 63, 96):
        hpln, hplt = wp.all_pix2world([[x0 + 48, ny / 2]], 0)[0]
        eps = float(np.degrees(np.arccos(np.cos(np.radians(hpln)) * np.cos(np.radians(hplt)))))
        ent = {"hpln": round(float(hpln), 2), "l2_bg_msb": float(np.nanmedian(im2[:, x0:x0 + 96]))}
        for name, img in (("l3u", hp3), ("l2", hp2)):
            s = img[:, x0:x0 + 96][~mask[:, x0:x0 + 96]]
            s = s[np.isfinite(s)]
            sg = 1.4826 * np.median(np.abs(s - np.median(s)))
            ent[name + "_px_sigma"] = float(sg)
            ent[name + "_depth_5sig_V"] = float(coef[0] - 2.5 * np.log10(5 * sg * np.sqrt(np.pi * AP_R ** 2)))
        strips[f"{eps:.1f}"] = ent
    return {"frame": f2.name, "date_avg": h["DATE-AVG"], "r_au": round(r_au, 5),
            "xposure_s": float(h["XPOSURE"]), "nsumexp": int(h["NSUMEXP"]), "l3_scale": round(scale, 5),
            "psf_radial_profile_0.5px_bins": [round(x, 3) for x in radial], "psf_n_stars": len(prof),
            "encircled_fraction_by_radius_px": enc, "photometry": phot, "noise_vs_elongation": strips}


def main():
    out = {}
    for f2 in sorted(RECON.glob("psp_L2_wispr_*_1211.fits")):
        r = analyse(f2)
        out[r["frame"]] = r
        print(json.dumps({k: v for k, v in r.items() if k != "noise_vs_elongation"}, indent=None), flush=True)
        for k, v in r["noise_vs_elongation"].items():
            print("  eps", k, {a: (f"{b:.2e}" if isinstance(b, float) and abs(b) < 1e-3 else round(b, 2)) for a, b in v.items()}, flush=True)
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
