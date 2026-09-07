"""SoloHI recon star check (go/no-go astrometry + photometry), L2 (MSB,
calibrated; SOAR serves no L3 for SoloHI). Adapted from the WISPR
`recon_star_check.py`.

Per frame (every tile, four epochs at r 0.30 / 0.39 / 0.62 / 1.01 AU):
* header celestial WCS ('A', RA/DEC-ZPN) vs Hipparcos: match fraction,
  residual rms, mean offset;
* per-star zero points V = ZP - 2.5 log10(MSB flux) with a B-V colour
  term; ZP vs V bin and vs elongation (flat-field / vignetting);
* PSF width from bright stars;
* single-frame noise and 5-sigma point-source depth vs elongation
  (column strips, elongation from the primary HPLN/HPLT WCS);
* the field-model check: the frame's pointing mapped into the
  solo_geometry orbit-plane frame.

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
import solo_geometry as G

REPO = Path(__file__).resolve().parents[3]
RECON = REPO / "runs" / "solohi-crossings" / "recon"
CAT = REPO / "runs" / "heliospheric-crossings" / "starcat" / "hip_v95.npz"
OUT = REPO / "surveys" / "solohi-crossings" / "results" / "recon_star_check_v1.json"
FRAMES = [
    "solo_L2_solohi-1ft_20221010T000035_V02.fits", "solo_L2_solohi-2ft_20221010T000236_V02.fits",
    "solo_L2_solohi-3ft_20221010T000524_V02.fits",
    "solo_L2_solohi-1ft_20250320T000259_V01.fits", "solo_L2_solohi-2ft_20250320T000459_V01.fits",
    "solo_L2_solohi-3fg_20250320T000920_V01.fits", "solo_L2_solohi-4fg_20250320T002120_V01.fits",
    "solo_L2_solohi-1ft_20260201T002215_V01.fits", "solo_L2_solohi-2ft_20260201T003015_V01.fits",
    "solo_L2_solohi-3ft_20260201T004005_V01.fits",
    "solo_L2_solohi-1ft_20211227T002044_V02.fits", "solo_L2_solohi-2ft_20211227T031444_V02.fits",
    "solo_L2_solohi-3fg_20211227T033213_V02.fits", "solo_L2_solohi-4fg_20211227T001213_V02.fits",
]
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
        eps, lam, bet = G.solohi_coords(s, mjd)
        out[n] = {"hpln": round(float((hpln + 180) % 360 - 180), 2), "hplt": round(float(hplt), 2),
                  "eps": round(float(eps[0]), 2), "ram_lon": round(float(lam[0]), 2),
                  "orbit_lat": round(float(bet[0]), 2)}
    v_hdr = np.hypot(np.hypot(hd["HCIX_VOB"], hd["HCIY_VOB"]), hd["HCIZ_VOB"]) / 1e3
    v_tab = np.linalg.norm(G.solo_vel_au_d(mjd)[0]) * G.AU_KM / 86400
    r_hdr = hd["DSUN_OBS"] / (G.AU_KM * 1e3)
    r_tab = float(G.heliocentric_r_au(mjd)[0])
    out["speed_km_s_header_vs_table"] = [round(float(v_hdr), 2), round(float(v_tab), 2)]
    out["r_au_header_vs_table"] = [round(float(r_hdr), 5), round(r_tab, 5)]
    out["sc_roll_deg"] = float(hd["SC_ROLL"])
    out["hglt_obs_deg"] = float(hd["HGLT_OBS"])
    return out


def check_frame(cat, f2):
    h2, im2 = load(RECON / f2)
    wa, wp = WCS(h2, key="A"), WCS(h2)
    ep = Time(h2["DATE-AVG"]).jyear
    ra = cat["ra"] + cat["pmra_masyr"] / 3.6e6 * (ep - 1991.25) / np.cos(np.radians(cat["dec"]))
    dec = cat["dec"] + cat["pmde_masyr"] / 3.6e6 * (ep - 1991.25)
    px, py = wa.all_world2pix(ra, dec, 0)
    nx, ny = h2["NAXIS1"], h2["NAXIS2"]
    inf = np.isfinite(px) & (px > 10) & (px < nx - 10) & (py > 10) & (py < ny - 10)
    c_ra, c_dec = wa.all_pix2world([[nx / 2, ny / 2]], 0)[0]
    cv = G.radec_to_vec(c_ra, c_dec)
    sep = np.degrees(np.arccos(np.clip(G.radec_to_vec(ra, dec) @ cv, -1, 1)))
    inf &= sep < 30
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
              "median_resid_px": float(np.median(d[ok])) if ok.any() else None,
              "mean_offset_px": [float(dx.mean()), float(dy.mean())] if ok.any() else None,
              "rms_px": [float(dx.std()), float(dy.std())] if ok.any() else None,
              "pixel_arcsec": float(h2["CDELT1"]) * 3600}

    rows = []
    for k in np.nonzero(ok)[0]:
        x, y = float(xs[i[k]]), float(ys[i[k]])
        f2_, e2_, v2_ = aper(hp, x, y)
        hpln, hplt = wp.all_pix2world([[x, y]], 0)[0]
        rows.append((vmag[k], bv[k], f2_, e2_, v2_, x, y, (hpln + 180) % 360 - 180, hplt))
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
              "zp_by_vbin": {}, "resid_by_elongation": {}}
    for lo in range(3, 10):
        m = hi & (R[good, 0] >= lo) & (R[good, 0] < lo + 1)
        if m.sum() >= 5:
            photom["zp_by_vbin"][f"{lo}-{lo+1}"] = {
                "n": int(m.sum()), "zp_median": float(np.median(zp[m])),
                "resid_median": float(np.median((zp - A @ coef)[m]))}
    eps_star = np.degrees(np.arccos(np.cos(np.radians(R[good, 7])) * np.cos(np.radians(R[good, 8]))))
    for lo in np.arange(np.floor(eps_star.min() / 5) * 5, eps_star.max(), 5):
        m = hi & (eps_star >= lo) & (eps_star < lo + 5)
        if m.sum() >= 5:
            rr = (zp - A @ coef)[m]
            photom["resid_by_elongation"][f"{int(lo)}-{int(lo)+5}"] = {
                "n": int(m.sum()), "resid_median": float(np.median(rr)),
                "mad": float(1.4826 * np.median(np.abs(rr - np.median(rr))))}

    fw = []
    for r in R[good & (R[:, 0] < 6.5) & (R[:, 0] > 3.0)]:
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
        s = hp[:, x0:x0 + 64]
        s = s[np.isfinite(s) & (s != 0)]
        sg = 1.4826 * np.median(np.abs(s - np.median(s)))
        bg = float(np.nanmedian(im2[:, x0:x0 + 64]))
        depth = coef[0] - 2.5 * np.log10(5 * sg * np.sqrt(np.pi * AP_R ** 2))
        strips[f"{eps:.1f}"] = {"hpln": round(float((hpln + 180) % 360 - 180), 2), "bg_msb": bg,
                                "hp_sigma_px": float(sg), "depth_5sig_V_solar_colour": float(depth)}
    return {
        "frame": {"l2": f2, "date_avg": h2["DATE-AVG"], "xposure_s": float(h2["XPOSURE"]),
                  "nsumexp": int(h2["NSUMEXP"]), "bunit": h2["BUNIT"], "cdelt_deg": float(h2["CDELT1"]),
                  "naxis": [int(nx), int(ny)], "nan_pixels": int(bad.sum()),
                  "dsun_au": float(h2["DSUN_OBS"]) / (G.AU_KM * 1e3), "obs_mode": h2.get("OBS_MODE"),
                  "gainmode": h2.get("GAINMODE"), "detector": h2.get("DETECTOR"),
                  "version": h2.get("VERSION"), "vers_cal": h2.get("VERS_CAL"),
                  "datasat": h2.get("DATASAT"), "datamax": float(h2.get("DATAMAX", np.nan))},
        "field_model_check": field_model_check(h2),
        "astrometry_header_wcs_vs_hipparcos": astrom,
        "photometry": photom, "psf": psf, "noise_vs_elongation": strips}


def main():
    cat = np.load(CAT)
    result = {"frames": {}}
    for f2 in FRAMES:
        r = check_frame(cat, f2)
        result["frames"][f2] = r
        fr = r["frame"]
        print(f"== {f2}  r={fr['dsun_au']:.3f} AU  xposure {fr['xposure_s']} s  mode {fr['obs_mode']} gain {fr['gainmode']}")
        fm = r["field_model_check"]
        print("   centre", fm["centre"], " roll", fm["sc_roll_deg"], " speed", fm["speed_km_s_header_vs_table"])
        print("  ", json.dumps(r["astrometry_header_wcs_vs_hipparcos"]))
        print("  ", json.dumps({k: v for k, v in r["photometry"].items() if k not in ("zp_by_vbin", "resid_by_elongation")}))
        print("   zp_by_vbin", json.dumps(r["photometry"]["zp_by_vbin"]))
        print("   resid_by_eps", json.dumps(r["photometry"]["resid_by_elongation"]))
        print("   psf", json.dumps(r["psf"]))
        for k, v in r["noise_vs_elongation"].items():
            print(f"    eps {k:>5s}  bg {v['bg_msb']:.2e}  sig {v['hp_sigma_px']:.2e}  depth V {v['depth_5sig_V_solar_colour']:.2f}")
    result["notes"] = ["aperture r=2 px, annulus 4-8 px on a 25-px median high-pass",
                       "depth = 5 sigma on the aperture noise for a B-V=0.65 source, one L2 frame",
                       "ZP: V = ZP - 2.5 log10(sum MSB in aperture); colour coefficient positive = red stars brighter",
                       "field_model_check: ram_lon < 0 means the anti-ram side per solo_geometry (t = prograde tangential)"]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=1) + "\n")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
