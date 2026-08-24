"""Channel-A dev re-run under freeze amendment v1.1.

Same-rung-only pseudo-window exclusion, parallax-factor systematics
template fitted on off-window epochs, exceedance rule S > max(T, 0).
Reuses runs/ztf-crossings/products/dev cutouts; fetches only epochs the
v1.0 run did not need (new valid offsets + template fit subsample).
Stages: --plan, --fetch, --stats.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import yaml
from astropy.coordinates import get_body_barycentric
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from dev_search import (BANDCODE, CUT_A, D, ERA, EVENTS, FREEZE, COV, RUN,
                        _epoch_dest, _fluxmap, _snr_stack, target_windows)
from sglsurvey.adapters.base import ConeRegion, CutoutSpec
from sglsurvey.adapters.irsa_ztf import ZtfNominalFootprint, ZtfSciAdapter
from sglsurvey.snapshots import SnapshotStore

DEV_A = FREEZE["split"]["A"]["dev"]
OFFSETS = [-97, -71, -47, -23, 23, 47, 71, 97]
REDRAWS = [-127, -113, 113, 127]
FIT_MAX, FIT_MIN = 150, 30
REG = yaml.safe_load(open(REPO / "registries" / "pilot_wise_2026.yaml"))["targets"]
PLAN_PATH = D / "configs" / "dev_plan_v11.json"


def valid_offsets_same_rung(tid, rad):
    rung = [w for w in target_windows(tid) if w[2] == rad]
    chosen = []
    for off in OFFSETS + REDRAWS:
        ok = True
        for lo, hi, _, _ in rung:
            slo, shi = lo + off, hi + off
            if shi < ERA.start_mjd_utc or slo > ERA.stop_mjd_utc:
                continue
            if any(slo < whi and shi > wlo for wlo, whi, _, _ in rung):
                ok = False
                break
        if ok:
            chosen.append(off)
        if len(chosen) == 8:
            break
    return chosen


def parallax_factors(tid, mjd):
    st = REG[tid]["state"]["astrometry"]
    plx = st["parallax_mas"] / 1000.0            # arcsec
    ra, de = np.radians(st["ra_deg"]), np.radians(st["dec_deg"])
    t = Time(np.atleast_1d(mjd), format="mjd", scale="utc")
    xyz = get_body_barycentric("earth", t).xyz.to_value("AU")
    X, Y, Z = xyz[0], xyz[1], xyz[2]
    p_ra = plx * (X * np.sin(ra) - Y * np.cos(ra))
    p_de = plx * (X * np.cos(ra) * np.sin(de) + Y * np.sin(ra) * np.sin(de)
                  - Z * np.cos(de))
    return p_ra, p_de


def unit_epochs(tid, band, ra, dec, omap):
    eps = []
    for k, o in omap.items():
        if o.band != BANDCODE[band]:
            continue
        q = o.quality_flags
        if q.get("bad_quality") or (q.get("seeing") or 99) > 4.0:
            continue
        if not ZtfNominalFootprint(o).contains(ra, dec):
            continue
        eps.append(o)
    return sorted(eps, key=lambda o: o.t_mid_mjd_utc)


def build_plan(discover=True):
    adapter, store = ZtfSciAdapter(), SnapshotStore(RUN)
    v10 = json.loads((D / "configs" / "dev_plan_v1.json").read_text())
    units = [u for u in FREEZE["search_units"]["A"]["units"]
             if u["target_id"] in DEV_A]
    pos = {p["target_id"]: (p["ra"], p["dec"]) for p in v10["A"]}
    plan = {"A": []}
    obs_cache = {}
    for u in units:
        tid, band, rad = u["target_id"], u["band"], u["radius_au"]
        ra, dec = pos[tid]
        if tid not in obs_cache:
            obs_cache[tid] = {json.dumps(o.native_key, sort_keys=True): o
                              for o in adapter.discover(
                                  ConeRegion(ra, dec, 0.05), ERA, store)}
        omap = obs_cache[tid]
        offs = valid_offsets_same_rung(tid, rad)
        all_wins = [(lo, hi) for lo, hi, _, _ in target_windows(tid)]
        rung_wins = [(lo, hi) for lo, hi, r_, _ in target_windows(tid)
                     if r_ == rad]
        eps = unit_epochs(tid, band, ra, dec, omap)
        off_window = [o for o in eps
                      if not any(lo <= o.t_mid_mjd_utc <= hi
                                 for lo, hi in all_wins)]
        step = max(1, len(off_window) // FIT_MAX)
        fit_eps = off_window[::step][:FIT_MAX]
        spans = list(rung_wins)
        spans += [(lo + o, hi + o) for o in offs for lo, hi in rung_wins]
        stat_eps = [o for o in eps
                    if any(lo <= o.t_mid_mjd_utc <= hi for lo, hi in spans)]
        keys = {json.dumps(o.native_key, sort_keys=True)
                for o in fit_eps + stat_eps}
        plan["A"].append({
            "target_id": tid, "band": band, "radius_au": rad,
            "ra": ra, "dec": dec, "valid_offsets": offs,
            "control_defined": len(offs) == 8,
            "rung_windows": [[lo, hi] for lo, hi in rung_wins],
            "all_windows": [[lo, hi] for lo, hi in all_wins],
            "n_fit_epochs": len(fit_eps), "n_stat_epochs": len(stat_eps),
            "epoch_keys": sorted(keys)})
    PLAN_PATH.write_text(json.dumps(plan, indent=1) + "\n")
    need_fetch = set()
    for p in plan["A"]:
        omap = obs_cache[p["target_id"]]
        for k in p["epoch_keys"]:
            o = omap[k]
            if not _epoch_dest(p["target_id"], "A", o).exists():
                need_fetch.add((p["target_id"], k))
    print(json.dumps({
        "units": [{k: p[k] for k in ("target_id", "band", "radius_au",
                                     "valid_offsets", "control_defined",
                                     "n_fit_epochs", "n_stat_epochs")}
                  for p in plan["A"]],
        "new_epochs_to_fetch": len(need_fetch)}, indent=1))
    return plan, obs_cache, adapter, sorted(need_fetch)


def fetch():
    plan, obs_cache, adapter, need = build_plan()
    pos = {p["target_id"]: (p["ra"], p["dec"]) for p in plan["A"]}
    print(f"{len(need)} new epochs", flush=True)
    ok = fail = 0
    for i, (tid, k) in enumerate(need):
        o = obs_cache[tid][k]
        ra, dec = pos[tid]
        dest = _epoch_dest(tid, "A", o)
        cut = CutoutSpec(ra_deg=ra, dec_deg=dec, size_pix=CUT_A)
        good = True
        for kind in ("sci", "msk", "diff"):
            try:
                adapter.fetch(o, (kind,), dest, cutout=cut)
            except FileNotFoundError:
                if kind != "diff":
                    good = False
                    break
            except Exception as e:
                print(f"  ! {dest.name}/{kind}: {e}", flush=True)
                good = False
                fail += 1
                break
        ok += bool(good)
        if (i + 1) % 200 == 0:
            print(f"  {i+1}/{len(need)} ok={ok} fail={fail}", flush=True)
    print(f"done {ok}/{len(need)} (fail={fail})", flush=True)


def stats():
    plan, obs_cache, adapter, _ = build_plan()
    results = []
    for p in plan["A"]:
        k_scale = 1.0
        tid, band = p["target_id"], p["band"]
        omap = obs_cache[tid]
        ra, dec = p["ra"], p["dec"]
        # per-epoch samples for every planned epoch
        samples = {}
        for k in p["epoch_keys"]:
            o = omap[k]
            fm = _fluxmap(tid, "A", o)
            if fm is None:
                continue
            f, v, g = fm.sample(ra, dec)
            if np.isfinite(g[0]) and g[0] >= 0.7 and np.isfinite(f[0]):
                samples[k] = (o.t_mid_mjd_utc, float(f[0]), float(v[0]),
                              float(o.quality_flags.get("seeing") or 2.0))
        all_wins = p["all_windows"]
        off_items = [(k, s) for k, s in samples.items()
                     if not any(lo <= s[0] <= hi for lo, hi in all_wins)]
        # template fit (3 x 3-sigma clip) on off-window epochs
        t_arr = np.array([s[0] for _, s in off_items])
        f_arr = np.array([s[1] for _, s in off_items])
        see = np.array([s[3] for _, s in off_items])
        gate = {"n_fit": len(off_items)}
        coef, med_see = None, None
        if len(off_items) >= FIT_MIN:
            pra, pde = parallax_factors(tid, t_arr)
            med_see = float(np.median(see))
            B = np.column_stack([np.ones_like(t_arr),
                                 (t_arr - 59800.0) / 365.25, pra, pde,
                                 see - med_see])
            keep = np.ones(len(t_arr), bool)
            for _ in range(3):
                coef, *_ = np.linalg.lstsq(B[keep], f_arr[keep], rcond=None)
                res = f_arr - B @ coef
                sd = 1.4826 * np.median(np.abs(res[keep]
                                               - np.median(res[keep])))
                keep = np.abs(res - np.median(res[keep])) <= 3 * sd
            gate["n_fit_kept"] = int(keep.sum())
            gate["fit_rms_before"] = float(np.std(f_arr))
            gate["fit_rms_after"] = float(np.std((f_arr - B @ coef)[keep]))
            if V12:
                v_arr = np.array([s_[2] for _, s_ in off_items])
                r2v = ((f_arr - B @ coef) ** 2 / v_arr)[keep]
                k_scale = max(1.0, float(np.median(r2v) / 0.4549))
                gate["k_scale"] = round(k_scale, 1)
        else:
            gate["n_fit_kept"] = len(off_items)

        def resid(k):
            t, f, v, s = samples[k]
            if coef is None:
                return f, v
            pra, pde = parallax_factors(tid, t)
            b = np.array([1.0, (t - 59800.0) / 365.25,
                          float(pra[0]), float(pde[0]), s - med_see])
            return f - float(b @ coef), v * k_scale

        def span_S(spans):
            per_win = []
            for lo, hi in spans:
                sm = [resid(k) for k, s in samples.items()
                      if lo <= s[0] <= hi]
                if sm:
                    per_win.append(_snr_stack(sm))
            if not per_win:
                return np.nan, 0
            Sws = np.array([s for s, _ in per_win])
            return float(np.sum(Sws) / np.sqrt(len(Sws))), len(per_win)

        rung = [tuple(w) for w in p["rung_windows"]]
        S, nw = span_S(rung)
        controls = []
        for off in p["valid_offsets"]:
            Sc, _ = span_S([(lo + off, hi + off) for lo, hi in rung])
            if np.isfinite(Sc):
                controls.append(round(Sc, 3))
        T = max(controls) if controls else np.nan
        med_c = float(np.median(controls)) if controls else np.nan
        searchable = (p["control_defined"]
                      and gate["n_fit_kept"] >= FIT_MIN
                      and np.isfinite(med_c) and abs(med_c) <= 1.0)
        exceed = bool(searchable and np.isfinite(S)
                      and S > max(T if np.isfinite(T) else 0.0, 0.0))
        results.append({
            "target_id": tid, "band": band, "radius_au": p["radius_au"],
            "n_valid_offsets": len(p["valid_offsets"]),
            "control_defined": p["control_defined"],
            "gate": gate, "median_control": None if np.isnan(med_c)
            else round(med_c, 3),
            "S": None if np.isnan(S) else round(S, 3),
            "n_windows_with_data": nw, "controls": controls,
            "T": None if np.isnan(T) else round(T, 3),
            "margin": (None if np.isnan(S) or np.isnan(T)
                       else round(S - max(T, 0.0), 3)),
            "status": ("searchable" if searchable else "constraint_only"),
            "exceedance": exceed})
    out = D / "results" / ("dev_search_v12.json" if V12
                           else "dev_search_v11.json")
    out.write_text(json.dumps({"A": results}, indent=1) + "\n")
    print(json.dumps({"A": results}, indent=1))


V12 = "--v12" in sys.argv

if __name__ == "__main__":
    if "--plan" in sys.argv:
        build_plan()
    elif "--fetch" in sys.argv:
        fetch()
    elif "--stats" in sys.argv:
        stats()
