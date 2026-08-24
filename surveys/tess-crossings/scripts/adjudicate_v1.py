"""Adjudication of the confirmatory chord exceedances (freeze v1.1).

Three chord-statistic exceedances (wolf-359 2.5 Rsun margin +0.46;
teegarden 2.5 Rsun +23.4; teegarden 0.1 AU +27.6) against 1.3
expected; all pulse statistics are null. The frozen ladder for chord
exceedances:

* strict-mask re-run: drop cadences within 0.05 d of any
  nonzero-QUALITY cadence and cadences with |POS_CORR| > 0.5 px;
* chord-shape/timing test: a beam chord is flat-topped across the
  window — split-half amplitudes (first vs second half of the
  in-window span) must agree; a scattered-light ramp is monotonic
  and fails;
* background-correlation annotation: Pearson r between the track
  sample series and the per-cadence field background (cutout median),
  in-window.

Census is not applicable (chord exceedances are not point events; the
pulse statistic, which census guards, is null everywhere).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from confirmatory_search import (CUBE_PREFIX, FREEZE, RINGS, Z_GRID,
                                 bilinear, build_maps, event_row,
                                 fit_fwhm_field, relay_tracks)
from tesscut_lib import PRODUCTS, load_cube

OUT = REPO / "surveys" / "tess-crossings" / "results"
WEIGHT_CAP = 20.0
GOOD_MIN = 0.7
STRICT_MARGIN_D = 0.05
POSCORR_MAX_PX = 0.5

CONF = json.loads((OUT / "confirmatory_v1.json").read_text())
EXC = [u for u in CONF["units"] if u["exceedance_c"]]


def chord_amp(mjd, f, v, g, lo, hi, keep=None):
    ok = np.isfinite(f) & np.isfinite(v) & (v > 0) & np.isfinite(g) \
        & (g >= GOOD_MIN)
    if keep is not None:
        ok &= keep
    mjd, f, v = mjd[ok], f[ok], v[ok]
    inw = (mjd >= lo) & (mjd <= hi)
    if inw.sum() < 5 or (~inw).sum() < 20:
        return None
    fo, vo = f[~inw], v[~inw]
    kmask = np.ones(len(fo), bool)
    for _ in range(3):
        med = np.median(fo[kmask])
        sd = 1.4826 * np.median(np.abs(fo[kmask] - med)) or 1.0
        kmask &= np.abs(fo - med) <= 3 * sd
    base = float(np.median(fo[kmask]))
    k = max(1.0, float(np.median((fo[kmask] - base) ** 2
                                 / vo[kmask]) / 0.4549)) \
        if kmask.sum() >= 6 else 1.0

    def amp(sel):
        if sel.sum() < 3:
            return None
        w = 1.0 / (v[sel] * k)
        w = np.minimum(w, WEIGHT_CAP * np.median(w))
        return float(np.sum(w * (f[sel] - base)) / np.sqrt(np.sum(w)))

    mid = 0.5 * (max(lo, mjd[inw].min()) + min(hi, mjd[inw].max()))
    return {"S_c": amp(inw), "k": k,
            "S_first": amp(inw & (mjd <= mid)),
            "S_second": amp(inw & (mjd > mid)),
            "n_in": int(inw.sum())}


def main():
    results = []
    for tid, prefix in CUBE_PREFIX.items():
        exc_units = [u for u in EXC if u["target_id"] == tid]
        if not exc_units:
            continue
        cube_dir = sorted(PRODUCTS.glob(f"{prefix}*"))[0]
        fpath = sorted(cube_dir.glob("*.fits"))[0]
        cube = load_cube(fpath, tid)
        with fits.open(fpath) as f:
            pc1 = np.asarray(f[1].data["POS_CORR1"], float)
            pc2 = np.asarray(f[1].data["POS_CORR2"], float)
        idx0 = np.where((cube.quality == 0)
                        & np.isfinite(cube.mjd_utc))[0]
        stack = np.nanmedian(cube.flux[idx0][::10], axis=0)
        fwhm, _ = fit_fwhm_field(stack)
        idx, F, V, G = build_maps(cube, fwhm)
        mjd = cube.mjd_utc[idx]
        # strict mask: distance to nearest nonzero-quality cadence
        bad_t = cube.mjd_utc[(cube.quality != 0)
                             & np.isfinite(cube.mjd_utc)]
        if len(bad_t):
            d = np.min(np.abs(mjd[:, None] - bad_t[None, :]), axis=1)
        else:
            d = np.full(len(mjd), np.inf)
        strict = (d >= STRICT_MARGIN_D) \
            & (np.abs(np.nan_to_num(pc1[idx])) <= POSCORR_MAX_PX) \
            & (np.abs(np.nan_to_num(pc2[idx])) <= POSCORR_MAX_PX)
        # per-cadence field background (cutout median)
        bg = np.nanmedian(cube.flux[idx].reshape(len(idx), -1), axis=1)
        ev = event_row(exc_units[0]["event_id"])
        tracks = relay_tracks(ev, mjd)
        cosd = np.cos(np.radians(tracks[..., 1]))

        fz = {(x["target_id"], x["radius_au"]): x
              for x in FREEZE["search_units"]["units"]}
        for u in exc_units:
            t_ca = fz[(u["target_id"], u["radius_au"])]["t_ca_mjd"]
            hd = u["window_days"] / 2.0
            lo, hi = t_ca - hd, t_ca + hd
            zi = list(Z_GRID).index(u["S_c_z"])

            def samples(dx, dy):
                ra = tracks[zi, :, 0] + dx / 3600.0 \
                    / np.maximum(cosd[zi], 0.05)
                de = tracks[zi, :, 1] + dy / 3600.0
                pix = cube.wcs.wcs_world2pix(
                    np.column_stack([ra, de]), 0)
                return (bilinear(F, pix[:, 0], pix[:, 1]),
                        bilinear(V, pix[:, 0], pix[:, 1]),
                        bilinear(G, pix[:, 0], pix[:, 1]))

            f0, v0, g0 = samples(0.0, 0.0)
            primary = chord_amp(mjd, f0, v0, g0, lo, hi)
            strict_res = chord_amp(mjd, f0, v0, g0, lo, hi, keep=strict)
            strict_ctrl = []
            for dx, dy in RINGS:
                fc, vc, gc = samples(dx, dy)
                st = chord_amp(mjd, fc, vc, gc, lo, hi, keep=strict)
                if st and st["S_c"] is not None:
                    strict_ctrl.append(st["S_c"])
            T_strict = max(strict_ctrl) if strict_ctrl else np.nan
            inw = (mjd >= lo) & (mjd <= hi) & np.isfinite(f0)
            r_bg = float(np.corrcoef(f0[inw], bg[inw])[0, 1]) \
                if inw.sum() > 5 else np.nan
            sf, ss = primary["S_first"], primary["S_second"]
            split_ratio = None
            if sf is not None and ss is not None and \
                    max(abs(sf), abs(ss)) > 0:
                split_ratio = round(min(sf, ss) / max(sf, ss), 3)
            row = {
                "target_id": tid, "radius_au": u["radius_au"],
                "S_c": u["S_c"], "T_c": u["T_c"],
                "margin_c": u["margin_c"], "S_c_z": u["S_c_z"],
                "split_half": {"first": None if sf is None
                               else round(sf, 2),
                               "second": None if ss is None
                               else round(ss, 2),
                               "ratio": split_ratio},
                "background_corr_inwindow": round(r_bg, 3),
                "strict": {"S_c": None if strict_res is None
                           else round(strict_res["S_c"], 3),
                           "T_c": round(T_strict, 3),
                           "exceedance": bool(
                               strict_res is not None
                               and strict_res["S_c"] is not None
                               and np.isfinite(T_strict)
                               and strict_res["S_c"]
                               > max(T_strict, 0.0)),
                           "n_in": None if strict_res is None
                           else strict_res["n_in"]},
            }
            results.append(row)
            print(json.dumps(row, indent=1), flush=True)

    (OUT / "adjudication_v1.json").write_text(
        json.dumps(results, indent=1) + "\n")


if __name__ == "__main__":
    main()
