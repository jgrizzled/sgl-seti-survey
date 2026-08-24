"""Confirmatory search (threshold freeze v1.1).

The 6 channel-B units (wolf-359 s42 and teegarden s91, three nested
rungs each), exactly as frozen: per-cadence matched-filter forced
photometry on the usable 31x31 TESScut cubes (Gaussian kernel at the
per-cutout PSF FWHM fitted from the brightest field star), track
positions per cadence for the nested z family (550/1000/2500/5500/
10000 AU — the programme's nested-family convention: unit statistic =
max over z, controls take the same max), per-(z, trajectory)
off-window baseline with 3x3sigma clip and empirical variance rescale
k (floor 1), WEIGHT_CAP 20x; two frozen statistics per unit — the
chord amplitude S_c and the per-cadence pulse maximum S_p — against
the 8 designated ring-control trajectories; exceedance S > max(T, 0).
Primary quality mask QUALITY == 0. Annotations recorded per unit:
straylight fraction of all in-window cadences, argmax cadence and sky
position of S_p (census input), per-z values.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.coordinates import get_body_barycentric
from astropy.table import Table
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import sglsurvey.photometry as phot
from sglsurvey.photometry import _gaussian_kernel, _matched_filter_full
from tesscut_lib import PRODUCTS, load_cube

phot.PSF_FWHM_ARCSEC.setdefault("T", 1.5 * 21.0)

D = REPO / "surveys" / "tess-crossings"
OUT = D / "results"
FREEZE = json.loads((D / "configs"
                     / "threshold_freeze_v1.json").read_text())
Z_GRID = (550.0, 1000.0, 2500.0, 5500.0, 10000.0)
RINGS = [tuple(o) for o in FREEZE["controls"]["offsets_arcsec"]]
WEIGHT_CAP = 20.0
GOOD_MIN = 0.7
PIX_ARCSEC = 21.0

CUBE_PREFIX = {"wolf-359": "B-wolf-359-0714db-s0042",
               "teegarden": "B-teegarden-fa8d3d-s0091"}


def event_row(eid):
    t = Table.read(REPO / "crossings" / "tess_v1" / "events.ecsv")
    m = [str(x) == eid for x in t["event_id"]]
    return t[np.asarray(m)][0]


def relay_tracks(ev, mjd):
    """(nz, n, 2) apparent relay ra/dec for the z grid, vectorized."""
    tt = Time(mjd, format="mjd", scale="utc")
    sun = get_body_barycentric("sun", tt).xyz.to_value("AU").T
    earth = get_body_barycentric("earth", tt).xyz.to_value("AU").T
    ra = np.radians(float(ev["star_icrs_ra_deg"]))
    de = np.radians(float(ev["star_icrs_dec_deg"]))
    u = np.array([np.cos(de) * np.cos(ra), np.cos(de) * np.sin(ra),
                  np.sin(de)])
    out = np.empty((len(Z_GRID), len(mjd), 2))
    for zi, z in enumerate(Z_GRID):
        v = (sun - z * u[None, :]) - earth
        v /= np.linalg.norm(v, axis=1)[:, None]
        out[zi, :, 0] = np.degrees(np.arctan2(v[:, 1], v[:, 0])) % 360
        out[zi, :, 1] = np.degrees(np.arcsin(np.clip(v[:, 2], -1, 1)))
    return out


def fit_fwhm_field(stack):
    iy, ix = np.unravel_index(np.nanargmax(stack), stack.shape)
    ny, nx = stack.shape
    if not (3 <= ix <= nx - 4 and 3 <= iy <= ny - 4):
        return 1.5, "edge_star_fallback"
    yy, xx = np.mgrid[0:ny, 0:nx]
    r2 = (xx - ix) ** 2 + (yy - iy) ** 2
    m = (r2 <= 25) & np.isfinite(stack)
    z = np.clip(stack[m] - np.nanmedian(stack), 1e-3, None)
    w = z / z.sum()
    sig2 = float(np.sum(w * r2[m]) / 2.0)
    return float(np.clip(2.355 * np.sqrt(max(sig2, 0.2)), 1.0, 3.0)), \
        f"field_star_at_{ix}_{iy}"


def build_maps(cube, fwhm_pix):
    """Per-primary-cadence matched-filter stacks."""
    kern = _gaussian_kernel(fwhm_pix,
                            max(3, int(np.ceil(2 * fwhm_pix))))
    idx = np.where((cube.quality == 0)
                   & np.isfinite(cube.mjd_utc))[0]
    n = len(idx)
    ny, nx = cube.flux.shape[1:]
    F = np.empty((n, ny, nx), np.float32)
    V = np.empty((n, ny, nx), np.float32)
    G = np.empty((n, ny, nx), np.float32)
    for j, i in enumerate(idx):
        img = np.asarray(cube.flux[i], float)
        err = np.asarray(cube.flux_err[i], float)
        good = np.isfinite(img) & np.isfinite(err) & (err > 0)
        var = np.where(good, err ** 2, 1e30)
        f, v, g, *_ = _matched_filter_full(img, var, good, fwhm_pix,
                                           True, kern)
        F[j], V[j], G[j] = f, v, g
    return idx, F, V, G


def bilinear(A, x, y):
    """Sample A[i] at (x[i], y[i]) bilinearly; NaN outside."""
    n, ny, nx = A.shape
    out = np.full(n, np.nan)
    ok = (x >= 0) & (x <= nx - 1.001) & (y >= 0) & (y <= ny - 1.001) \
        & np.isfinite(x) & np.isfinite(y)
    xi = np.floor(x[ok]).astype(int)
    yi = np.floor(y[ok]).astype(int)
    fx, fy = x[ok] - xi, y[ok] - yi
    ii = np.where(ok)[0]
    out[ok] = (A[ii, yi, xi] * (1 - fx) * (1 - fy)
               + A[ii, yi, xi + 1] * fx * (1 - fy)
               + A[ii, yi + 1, xi] * (1 - fx) * fy
               + A[ii, yi + 1, xi + 1] * fx * fy)
    return out


def stats_for(mjd, f, v, g, lo, hi):
    """Frozen chord + pulse statistics for one trajectory sample set."""
    ok = np.isfinite(f) & np.isfinite(v) & (v > 0) & np.isfinite(g) \
        & (g >= GOOD_MIN)
    mjd, f, v = mjd[ok], f[ok], v[ok]
    inw = (mjd >= lo) & (mjd <= hi)
    if inw.sum() == 0 or (~inw).sum() < 20:
        return None
    fo, vo = f[~inw], v[~inw]
    keep = np.ones(len(fo), bool)
    for _ in range(3):
        med = np.median(fo[keep])
        sd = 1.4826 * np.median(np.abs(fo[keep] - med)) or 1.0
        keep &= np.abs(fo - med) <= 3 * sd
    base = float(np.median(fo[keep]))
    k = 1.0
    if keep.sum() >= 6:
        k = max(1.0, float(np.median(
            (fo[keep] - base) ** 2 / vo[keep]) / 0.4549))
    w = 1.0 / (v[inw] * k)
    w = np.minimum(w, WEIGHT_CAP * np.median(w))
    S_c = float(np.sum(w * (f[inw] - base)) / np.sqrt(np.sum(w)))
    snr = (f[inw] - base) / np.sqrt(v[inw] * k)
    j = int(np.argmax(snr))
    return {"S_c": S_c, "S_p": float(snr[j]),
            "S_p_mjd": float(mjd[inw][j]), "k": k,
            "n_in": int(inw.sum()), "n_off": int(keep.sum())}


def main():
    results = []
    for tid, prefix in CUBE_PREFIX.items():
        units = [u for u in FREEZE["search_units"]["units"]
                 if u["target_id"] == tid]
        cube_dir = sorted(PRODUCTS.glob(f"{prefix}*"))[0]
        cube = load_cube(sorted(cube_dir.glob("*.fits"))[0], tid)
        idx0 = np.where((cube.quality == 0)
                        & np.isfinite(cube.mjd_utc))[0]
        stack = np.nanmedian(cube.flux[idx0][::10], axis=0)
        fwhm, fwhm_src = fit_fwhm_field(stack)
        print(f"[{tid}] fwhm={fwhm:.2f}px ({fwhm_src}); building "
              f"{len(idx0)} maps", flush=True)
        idx, F, V, G = build_maps(cube, fwhm)
        mjd = cube.mjd_utc[idx]
        ev = event_row(units[0]["event_id"])
        tracks = relay_tracks(ev, mjd)          # (nz, n, 2)
        cosd = np.cos(np.radians(tracks[..., 1]))

        for u in units:
            hd = u["window_days"] / 2.0
            lo, hi = u["t_ca_mjd"] - hd, u["t_ca_mjd"] + hd
            per_traj = []
            for dx, dy in [(0.0, 0.0)] + list(RINGS):
                best = None
                zvals = []
                for zi in range(len(Z_GRID)):
                    ra = tracks[zi, :, 0] + dx / 3600.0 \
                        / np.maximum(cosd[zi], 0.05)
                    de = tracks[zi, :, 1] + dy / 3600.0
                    pix = cube.wcs.wcs_world2pix(
                        np.column_stack([ra, de]), 0)
                    f = bilinear(F, pix[:, 0], pix[:, 1])
                    v = bilinear(V, pix[:, 0], pix[:, 1])
                    g = bilinear(G, pix[:, 0], pix[:, 1])
                    st = stats_for(mjd, f, v, g, lo, hi)
                    zvals.append(st)
                    if st and (best is None
                               or st["S_c"] > best["S_c"]):
                        best = {**st, "z": Z_GRID[zi]}
                sp = max((z["S_p"] for z in zvals if z), default=np.nan)
                sp_meta = max((z for z in zvals if z),
                              key=lambda z: z["S_p"], default=None)
                per_traj.append({"dx": dx, "dy": dy,
                                 "S_c": None if best is None
                                 else best["S_c"],
                                 "S_c_z": None if best is None
                                 else best["z"],
                                 "k": None if best is None
                                 else best["k"],
                                 "S_p": None if not np.isfinite(sp)
                                 else sp,
                                 "S_p_mjd": None if sp_meta is None
                                 else sp_meta["S_p_mjd"]})
            real, ctrls = per_traj[0], per_traj[1:]
            Tc = max((c["S_c"] for c in ctrls
                      if c["S_c"] is not None), default=np.nan)
            Tp = max((c["S_p"] for c in ctrls
                      if c["S_p"] is not None), default=np.nan)
            allw = (cube.mjd_utc >= lo) & (cube.mjd_utc <= hi) \
                & np.isfinite(cube.mjd_utc)
            stray = float(np.mean(cube.quality[allw] != 0)) \
                if allw.any() else np.nan
            exc_c = bool(real["S_c"] is not None
                         and real["S_c"] > max(Tc, 0.0))
            exc_p = bool(real["S_p"] is not None
                         and real["S_p"] > max(Tp, 0.0))
            row = {
                "target_id": tid, "event_id": u["event_id"],
                "radius_au": u["radius_au"], "b_rsun": u["b_rsun"],
                "sector": u["sector"], "window_days":
                u["window_days"], "fwhm_pix": round(fwhm, 2),
                "S_c": round(real["S_c"], 3), "S_c_z": real["S_c_z"],
                "T_c": round(Tc, 3),
                "margin_c": round(real["S_c"] - max(Tc, 0.0), 3),
                "exceedance_c": exc_c,
                "S_p": round(real["S_p"], 3),
                "S_p_mjd": real["S_p_mjd"],
                "T_p": round(Tp, 3),
                "margin_p": round(real["S_p"] - max(Tp, 0.0), 3),
                "exceedance_p": exc_p,
                "k": round(real["k"], 3),
                "inwindow_nonzero_quality_fraction":
                    round(stray, 3),
                "n_controls_c": sum(1 for c in ctrls
                                    if c["S_c"] is not None),
                "controls_S_c": [None if c["S_c"] is None else
                                 round(c["S_c"], 3) for c in ctrls],
                "controls_S_p": [None if c["S_p"] is None else
                                 round(c["S_p"], 3) for c in ctrls],
            }
            results.append(row)
            print(json.dumps({k: row[k] for k in
                              ("target_id", "radius_au", "S_c", "T_c",
                               "exceedance_c", "S_p", "T_p",
                               "exceedance_p", "k")}), flush=True)

    (OUT / "confirmatory_v1.json").write_text(
        json.dumps({"freeze": FREEZE["freeze_content_hash"],
                    "z_grid_au": list(Z_GRID),
                    "note": "nested z family: unit statistic = max "
                            "over z, controls take the same max (the "
                            "programme convention, made explicit "
                            "here)",
                    "units": results}, indent=1) + "\n")
    n_exc = sum(r["exceedance_c"] + r["exceedance_p"]
                for r in results)
    print(f"\n{len(results)} units x 2 statistics; "
          f"{n_exc} exceedances "
          f"(expected {FREEZE['threshold_rule']['expected_control_crossings']})")


if __name__ == "__main__":
    main()
