"""WISPR recon star check (the go/no-go astrometry + photometry test),
both telescopes, L2 (MSB, calibrated) and L3 (background-subtracted).

Per frame:
* header celestial WCS ('A' system, RA/DEC-ZPN) vs Hipparcos: match
  fraction, residual rms, mean offset;
* per-star zero points V = ZP - 2.5 log10(MSB flux) with a B-V colour
  term (WISPR-I 490-740 nm, WISPR-O 475-725 nm);
* PSF width from bright stars;
* single-frame noise and 5-sigma point-source depth vs elongation
  (column strips, elongation from the primary HPLN/HPLT WCS);
* L2 vs L3 aperture-flux parity;
* the field-model check: the frame's own pointing (CRVAL of the
  helioprojective WCS, the four corners) against the ram frame built
  from the header HCI velocity — is the field on the ram side, what
  are the elongation and orbit-latitude extents.

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
import psp_geometry as G

REPO = Path(__file__).resolve().parents[3]
RECON = REPO / "runs" / "wispr-crossings" / "recon"
CAT = REPO / "runs" / "heliospheric-crossings" / "starcat" / "hip_v95.npz"
OUT = REPO / "surveys" / "wispr-crossings" / "results" / "recon_star_check_v1.json"
FRAMES = {"I": ("psp_L2_wispr_20241224T000015_V1_1211.fits", "psp_L3_wispr_20241224T000015_V1_1211.fits"),
          "O": ("psp_L2_wispr_20241224T000205_V1_2222.fits", "psp_L3_wispr_20241224T000205_V1_2222.fits")}
AP_R, ANN = 2.0, (4.0, 8.0)


def load(path):
    warnings.simplefilter("ignore")
    h = fits.open(path)[0]
    return h.header, np.asarray(h.data, float)


def highpass(img, size=25):
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


def field_model_check(hd):
    """Header pointing vs the psp_geometry ram frame."""
    mjd = Time(hd["DATE-AVG"]).mjd
    wa, wp = WCS(hd, key="A"), WCS(hd)
    nx, ny = hd["NAXIS1"], hd["NAXIS2"]
    pts = [(nx / 2, ny / 2), (0, 0), (nx - 1, 0), (0, ny - 1), (nx - 1, ny - 1),
           (0, ny / 2), (nx - 1, ny / 2), (nx / 2, 0), (nx / 2, ny - 1)]
    names = ["centre", "c00", "cx0", "c0y", "cxy", "left", "right", "bottom", "top"]
    out = {}
    for n, (x, y) in zip(names, pts):
        ra, dec = wa.all_pix2world([[x, y]], 0)[0]
        hpln, hplt = wp.all_pix2world([[x, y]], 0)[0]
        s = G.radec_to_vec(ra, dec)
        eps, lam, bet = G.wispr_coords(s, mjd)
        out[n] = {"hpln": round(float(hpln), 2), "hplt": round(float(hplt), 2),
                  "eps": round(float(eps[0]), 2), "ram_lon": round(float(lam[0]), 2),
                  "orbit_lat": round(float(bet[0]), 2)}
    # header velocity vs our table velocity (HCI frame differs from ICRS,
    # so compare speeds and the Sun-relative geometry only)
    v_hdr = np.hypot(np.hypot(hd["HCIX_VOB"], hd["HCIY_VOB"]), hd["HCIZ_VOB"]) / 1e3
    v_tab = np.linalg.norm(G.psp_vel_au_d(mjd)[0]) * G.AU_KM / 86400
    r_hdr = hd["DSUN_OBS"] / (G.AU_KM * 1e3)
    r_tab = float(G.heliocentric_r_au(mjd)[0])
    out["speed_km_s_header_vs_table"] = [round(float(v_hdr), 2), round(float(v_tab), 2)]
    out["r_au_header_vs_table"] = [round(float(r_hdr), 5), round(r_tab, 5)]
    return out


def main():
    cat = np.load(CAT)
    result = {}
    for tel, (f2, f3) in FRAMES.items():
        h2, im2 = load(RECON / f2)
        h3, im3 = load(RECON / f3)
        wa, wp = WCS(h2, key="A"), WCS(h2)
        ep = Time(h2["DATE-AVG"]).jyear
        ra = cat["ra"] + cat["pmra_masyr"] / 3.6e6 * (ep - 1991.25) / np.cos(np.radians(cat["dec"]))
        dec = cat["dec"] + cat["pmde_masyr"] / 3.6e6 * (ep - 1991.25)
        px, py = wa.all_world2pix(ra, dec, 0)
        nx, ny = h2["NAXIS1"], h2["NAXIS2"]
        inf = np.isfinite(px) & (px > 10) & (px < nx - 10) & (py > 10) & (py < ny - 10)
        # ZPN projections can wrap far-off-axis stars back into the frame: keep only
        # stars whose true angular distance from the frame centre is < 40 deg
        c_ra, c_dec = wa.all_pix2world([[nx / 2, ny / 2]], 0)[0]
        cv = G.radec_to_vec(c_ra, c_dec)
        sep = np.degrees(np.arccos(np.clip(G.radec_to_vec(ra, dec) @ cv, -1, 1)))
        inf &= sep < 45
        vmag, bv = cat["vmag"][inf], cat["bv"][inf]
        px, py = px[inf], py[inf]

        hp = highpass(im2)
        bad = ~np.isfinite(im2)
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
                  "rms_px": [float(dx.std()), float(dy.std())],
                  "pixel_arcsec": float(h2["CDELT1"]) * 3600}

        HP3 = highpass(im3)
        HP3[~np.isfinite(im3)] = 0
        rows = []
        for k in np.nonzero(ok)[0]:
            x, y = float(xs[i[k]]), float(ys[i[k]])
            f2_, e2_, v2_ = aper(hp, x, y)
            f3_, e3_, v3_ = aper(HP3, x, y)
            hpln, hplt = wp.all_pix2world([[x, y]], 0)[0]
            rows.append((vmag[k], bv[k], f2_, e2_, v2_, x, y, f3_, e3_, hpln, hplt))
        R = np.array(rows)
        good = (R[:, 2] > 0) & (R[:, 4] >= 0.99) & np.isfinite(R[:, 1]) & (R[:, 1] < 90) & (R[:, 0] > 3.0)
        zp = R[good, 0] + 2.5 * np.log10(R[good, 2])
        A = np.c_[np.ones(good.sum()), R[good, 1] - 0.65]
        hi = R[good, 2] / R[good, 3] > 10
        coef, *_ = np.linalg.lstsq(A[hi], zp[hi], rcond=None)
        resid = zp[hi] - A[hi] @ coef
        photom = {"n_calibrators_snr10": int(hi.sum()),
                  "zp0_at_bv0.65": float(coef[0]), "colour_coeff_mag_per_bv": float(coef[1]),
                  "scatter_mad_no_colour": float(1.4826 * np.median(np.abs(zp[hi] - np.median(zp[hi])))),
                  "scatter_mad_with_colour": float(1.4826 * np.median(np.abs(resid))),
                  "zp_by_vbin": {}}
        for lo in range(3, 10):
            m = hi & (R[good, 0] >= lo) & (R[good, 0] < lo + 1)
            if m.sum() >= 5:
                photom["zp_by_vbin"][f"{lo}-{lo+1}"] = {
                    "n": int(m.sum()), "zp_median": float(np.median(zp[m])),
                    "resid_median": float(np.median((zp - A @ coef)[m]))}
        # residual vs elongation (flat-field / vignetting check)
        eps_star = np.degrees(np.arccos(np.cos(np.radians(R[good, 9])) * np.cos(np.radians(R[good, 10]))))
        photom["resid_by_elongation"] = {}
        for lo in np.arange(np.floor(eps_star.min() / 10) * 10, eps_star.max(), 10):
            m = hi & (eps_star >= lo) & (eps_star < lo + 10)
            if m.sum() >= 5:
                photom["resid_by_elongation"][f"{int(lo)}-{int(lo)+10}"] = {
                    "n": int(m.sum()), "resid_median": float(np.median((zp - A @ coef)[m])),
                    "mad": float(1.4826 * np.median(np.abs((zp - A @ coef)[m] - np.median((zp - A @ coef)[m]))))}

        f2v, f3v = R[good, 2], R[good, 7]
        pm = hi & (f3v > 0)
        ratio = f3v[pm] / f2v[pm]
        parity = {"n": int(pm.sum()), "median_L3_over_L2": float(np.median(ratio)),
                  "mad": float(1.4826 * np.median(np.abs(ratio - np.median(ratio))))}

        fw = []
        for r in R[good & (R[:, 0] < 6.0) & (R[:, 0] > 3.0)]:
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
               "fwhm_arcsec_median": float(np.median(fw) * h2["CDELT1"] * 3600) if fw else None}

        strips = {}
        for x0 in range(0, nx - 63, 64):
            hpln, hplt = wp.all_pix2world([[x0 + 32, ny / 2]], 0)[0]
            eps = float(np.degrees(np.arccos(np.cos(np.radians(hpln)) * np.cos(np.radians(hplt)))))
            for lvl, img, hpi in (("L2", im2, hp), ("L3", im3, HP3)):
                s = hpi[:, x0:x0 + 64]
                s = s[np.isfinite(s) & (s != 0)]
                sg = 1.4826 * np.median(np.abs(s - np.median(s)))
                bg = float(np.nanmedian(img[:, x0:x0 + 64]))
                depth = coef[0] - 2.5 * np.log10(5 * sg * np.sqrt(np.pi * AP_R ** 2))
                strips.setdefault(f"{eps:.1f}", {"hpln": round(float(hpln), 2)})[lvl] = {
                    "bg_msb": bg, "hp_sigma_px": float(sg), "depth_5sig_V_solar_colour": float(depth)}

        result[f"WISPR-{tel}"] = {
            "frames": {"l2": f2, "l3": f3, "date_avg": h2["DATE-AVG"], "xposure_s": float(h2["XPOSURE"]),
                       "nsumexp": int(h2["NSUMEXP"]), "bunit": h2["BUNIT"], "cdelt_deg": float(h2["CDELT1"]),
                       "naxis": [int(nx), int(ny)], "nan_pixels": int(bad.sum()),
                       "dsun_au": float(h2["DSUN_OBS"]) / (G.AU_KM * 1e3)},
            "field_model_check": field_model_check(h2),
            "astrometry_header_wcs_vs_hipparcos": astrom,
            "photometry": photom, "l2_l3_parity": parity, "psf": psf,
            "noise_vs_elongation": strips}
    result["notes"] = ["aperture r=2 px, annulus 4-8 px on a 25-px median high-pass",
                       "depth = 5 sigma on the aperture noise for a B-V=0.65 source, one L2 frame",
                       "ZP: V = ZP - 2.5 log10(sum MSB in aperture); colour coefficient positive = red stars brighter",
                       "field_model_check: ram_lon > 0 means the ram side per psp_geometry (t = prograde tangential)"]
    OUT.write_text(json.dumps(result, indent=1) + "\n")
    for tel in ("WISPR-I", "WISPR-O"):
        r = result[tel]
        print("==", tel, r["frames"])
        print(json.dumps(r["field_model_check"], indent=None))
        print(json.dumps(r["astrometry_header_wcs_vs_hipparcos"], indent=None))
        print(json.dumps({k: v for k, v in r["photometry"].items() if k != "zp_by_vbin"}, indent=None))
        print(json.dumps(r["photometry"]["zp_by_vbin"], indent=None))
        print(json.dumps(r["l2_l3_parity"]), json.dumps(r["psf"]))
        for k, v in r["noise_vs_elongation"].items():
            print(" eps", k, v)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
