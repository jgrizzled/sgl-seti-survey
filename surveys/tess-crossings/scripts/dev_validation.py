"""Dev-stage machinery validation (threshold freeze v1.1).

No searchable dev unit exists (amendment v1.1: ross-128 is off the
science array), so the frozen machinery is validated end-to-end on
the ZERO-TRIAL channel-A reference rows — teegarden s71 (200 s) and
van-maanen s43 (600 s) — before any confirmatory B unit is touched:

* per-cadence matched-filter forced photometry through the identical
  kernel machinery (flux_map_from_arrays, Gaussian kernel at the
  per-cutout PSF FWHM fitted from the star stack);
* off-window baseline (robust mean, 3x3sigma clip), empirical
  variance rescale k, WEIGHT_CAP chord statistic S_c and per-cadence
  pulse maximum S_p — reference values, threshold-free by freeze;
* flare screening annotation (in-window cadences > 5 sigma);
* the TIC flux-scale gate: field stars measured on a median stack
  through the same kernel; PASS requires scatter <= 0.2 mag
  (SPHEREx-gate rule). Run on both A cubes and the wolf-359 B cube
  (the B fields carry their own calibrators).

Writes results/dev_validation_v1.json (+ .md summary by hand).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import sglsurvey.photometry as phot
from sglsurvey.photometry import _gaussian_kernel, flux_map_from_arrays
from tesscut_lib import PRODUCTS, load_cube, tic_cone

# register the TESS band in the kernel table (the explicit kernel
# argument overrides it, but the lookup runs unconditionally)
phot.PSF_FWHM_ARCSEC.setdefault("T", 1.5 * 21.0)

OUT = REPO / "surveys" / "tess-crossings" / "results"
WEIGHT_CAP = 20.0
ZP_NOMINAL = 20.44
PIX_ARCSEC = 21.0

FREEZE = json.loads((REPO / "surveys" / "tess-crossings" / "configs"
                     / "threshold_freeze_v1.json").read_text())
A_ROWS = FREEZE["a_reference_rows"]["rows"]


def cube_for(prefix):
    for d in sorted(PRODUCTS.glob(f"{prefix}*")):
        p = sorted(d.glob("*.fits"))
        if p:
            return load_cube(p[0], d.name)
    raise FileNotFoundError(prefix)


def fit_fwhm(stack, cx, cy):
    """Crude radial Gaussian FWHM fit around (cx, cy)."""
    ny, nx = stack.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    r2 = (xx - cx) ** 2 + (yy - cy) ** 2
    m = (r2 <= 25) & np.isfinite(stack)
    z = stack[m] - np.nanmedian(stack)
    z = np.clip(z, 1e-3, None)
    w = z / z.sum()
    sig2 = float(np.sum(w * r2[m]) / 2.0)
    return float(np.clip(2.355 * np.sqrt(max(sig2, 0.2)), 1.0, 3.0))


def per_cadence_series(cube, x, y, fwhm_pix):
    """Matched-filter (flux, var, good_frac) at a fixed pixel position
    for every primary cadence, via the identical kernel machinery."""
    kern = _gaussian_kernel(fwhm_pix, max(3, int(np.ceil(2 * fwhm_pix))))
    ok = cube.quality == 0
    f_out, v_out = [], []
    idx = np.where(ok)[0]
    for i in idx:
        img = np.asarray(cube.flux[i], float)
        err = np.asarray(cube.flux_err[i], float)
        good = np.isfinite(img) & np.isfinite(err) & (err > 0)
        var = np.where(good, err ** 2, 1e30)
        fm = flux_map_from_arrays(img, var, good, cube.wcs, "T",
                                  float(cube.mjd_utc[i]), kernel=kern,
                                  pix_arcsec=PIX_ARCSEC)
        xi, yi = int(round(x)), int(round(y))
        f_out.append(float(fm.flux[yi, xi]))
        v_out.append(float(fm.var[yi, xi]))
    return (cube.mjd_utc[idx], np.asarray(f_out), np.asarray(v_out))


def chord_pulse(mjd, f, v, lo, hi):
    inw = (mjd >= lo) & (mjd <= hi)
    off = ~inw
    fo = f[off]
    for _ in range(3):
        med = np.median(fo)
        sd = 1.4826 * np.median(np.abs(fo - med)) or 1.0
        fo = fo[np.abs(fo - med) <= 3 * sd]
    base = float(np.median(fo))
    vo = v[off][np.isin(f[off], fo)]
    k = max(1.0, float(np.median((fo - base) ** 2
                                 / np.maximum(vo[:len(fo)], 1e-30))
                       / 0.4549)) if len(fo) >= 6 else 1.0
    w = 1.0 / (v[inw] * k)
    w = np.minimum(w, WEIGHT_CAP * np.median(w))
    S_c = float(np.sum(w * (f[inw] - base)) / np.sqrt(np.sum(w)))
    snr = (f[inw] - base) / np.sqrt(v[inw] * k)
    S_p = float(np.max(snr))
    n_flare = int(np.sum(snr > 5.0))
    return base, k, S_c, S_p, n_flare, int(inw.sum()), int(off.sum())


def zp_check(cube, exclude_xy=None):
    """TIC stars measured on a median stack through the same kernel."""
    ok = cube.quality == 0
    stack = np.nanmedian(cube.flux[ok][::10], axis=0)
    err = np.nanmedian(cube.flux_err[ok][::10], axis=0)
    good = np.isfinite(stack) & np.isfinite(err) & (err > 0)
    ny, nx = stack.shape
    ctr = cube.wcs.wcs_pix2world([[nx / 2, ny / 2]], 0)[0]
    rows = tic_cone(float(ctr[0]), float(ctr[1]),
                    max(nx, ny) * PIX_ARCSEC / 3600.0 / 1.4)
    kern = _gaussian_kernel(1.5, 5)
    fm = flux_map_from_arrays(stack, np.where(good, err ** 2, 1e30),
                              good, cube.wcs, "T", 0.0, kernel=kern,
                              pix_arcsec=PIX_ARCSEC)
    zps = []
    for r in rows:
        t = r.get("Tmag")
        if t is None or not (9.0 <= float(t) <= 15.0):
            continue
        pix = cube.wcs.wcs_world2pix([[float(r["ra"]),
                                       float(r["dec"])]], 0)[0]
        x, y = pix
        if not (2 <= x <= nx - 3 and 2 <= y <= ny - 3):
            continue
        if exclude_xy is not None and np.hypot(
                x - exclude_xy[0], y - exclude_xy[1]) < 2.5:
            continue
        fl = float(fm.flux[int(round(y)), int(round(x))])
        if fl > 5 * float(np.sqrt(fm.var[int(round(y)),
                                         int(round(x))])):
            zps.append(float(t) + 2.5 * np.log10(fl))
    zps = np.asarray(zps)
    if len(zps) < 3:
        return {"n_stars": int(len(zps)), "status": "sparse"}
    return {"n_stars": int(len(zps)),
            "zp_median": round(float(np.median(zps)), 3),
            "zp_scatter_mad": round(float(
                1.4826 * np.median(np.abs(zps - np.median(zps)))), 3),
            "offset_vs_nominal": round(float(
                np.median(zps) - ZP_NOMINAL), 3),
            "status": "pass" if 1.4826 * np.median(
                np.abs(zps - np.median(zps))) <= 0.2 else "FAIL"}


def main():
    results = {"a_reference": [], "zp_checks": {}}
    for row in A_ROWS:
        tid, sec = row["target_id"], row["sector"]
        cube = cube_for(f"A-{tid}-s{sec:04d}")
        ok = cube.quality == 0
        stack = np.nanmedian(cube.flux[ok][::10], axis=0)
        ny, nx = stack.shape
        cx, cy = nx // 2, ny // 2
        fwhm = fit_fwhm(stack, cx, cy)
        mjd, f, v = per_cadence_series(cube, cx, cy, fwhm)
        hd = row["window_days"] / 2.0
        lo, hi = row["t_ca_mjd"] - hd, row["t_ca_mjd"] + hd
        base, k, S_c, S_p, n_flare, n_in, n_off = chord_pulse(
            mjd, f, v, lo, hi)
        tmag_implied = ZP_NOMINAL - 2.5 * np.log10(max(base, 1e-3))
        rec = {"target_id": tid, "sector": sec,
               "fwhm_pix": round(fwhm, 2),
               "n_in": n_in, "n_off": n_off,
               "star_flux_e_s": round(base, 1),
               "implied_T_nominal_zp": round(float(tmag_implied), 2),
               "k_rescale": round(k, 2),
               "S_c_reference": round(S_c, 2),
               "S_p_reference": round(S_p, 2),
               "n_inwindow_gt5sigma": n_flare}
        results["a_reference"].append(rec)
        print(rec, flush=True)
        results["zp_checks"][f"A-{tid}-s{sec}"] = zp_check(
            cube, exclude_xy=(cx, cy))
        print("  zp:", results["zp_checks"][f"A-{tid}-s{sec}"],
              flush=True)
    bcube = cube_for("B-wolf-359-0714db-s0042")
    results["zp_checks"]["B-wolf-359-s42"] = zp_check(bcube)
    print("  zp B:", results["zp_checks"]["B-wolf-359-s42"], flush=True)
    (OUT / "dev_validation_v1.json").write_text(
        json.dumps(results, indent=1) + "\n")


if __name__ == "__main__":
    main()
