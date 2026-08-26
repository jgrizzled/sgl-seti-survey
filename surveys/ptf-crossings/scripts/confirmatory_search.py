"""Confirmatory search (threshold freeze v1.0, confirmatory split).

Runs the six confirmatory channel-B units under the frozen
constructions exactly as validated by the dev run
(results/dev_search_v1.md): scie-direct track statistics with ring
controls, dwarf-locus calibrator selection, 256-px cutouts, k rescale
per (unit, band) from off-window antipode samples. Blind discipline:
fetch, measure, reduce once; any exceedance is adjudicated by the
frozen veto ladder, never silently dropped.

Units (freeze v1.0 `search_units.B`, split=confirmatory):
  van-maanen B 2.5Rsun g / B 0.1AU g / B 0.1AU R;
  ross-128  B 2.5Rsun R / B 0.1AU R / B 0.1AU g — 10 trials
  (S_event everywhere except the two single-multi-event units;
  S_stack where >= 2 covered events).

Statistics per unit (freeze D4):
  S_event = max over the unit's multi-epoch (>= 2 usable in-window
  epochs at measure time) events of max_z of the weighted in-window
  track stack; controls take the identical double max on the 8 ring
  trajectories. S_stack = max_z of the weighted stack pooling every
  covered event's in-window samples at that z (one shared relay
  distance across recurrences); ring controls identical. Events that
  fall below 2 usable epochs at measure time are reported; if no
  multi-epoch event survives, S_event is constraint-only (frozen
  rule).

Off-window sample pool per (target, band): every usable era epoch at
the per-era antipode point outside every 0.1 AU window of that
target (the widest rung, so the pool is signal-free under every
rung hypothesis).

Stages: --plan, --fetch, --stats. Products under
runs/ptf-crossings/products/conf/.

Amendment v1.1 (2026-08-26, frozen before re-reduction): confirmatory
cutouts 256 -> 384 px. The 256-px realization was inoperable at the
ross-128 antipode — 142/144 frames held < 5 in-frame calibrators
(surface-density, not selection: at 384 px the UNCHANGED dev-fixed
calibrator rule yields 9-13 calibrators at 0.04-0.11 mag scatter,
tested on in- and off-window frames). No statistic, threshold, gate
or selection element changed; all units re-fetched and re-measured
uniformly under the amended support. The superseded 256-px pass
(results/confirmatory_v1_cut256_superseded.json) had measured the
van-maanen units only: S 0.8-2.0, every margin negative — nothing
about the amendment was conditioned on those values.
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import requests

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from sglsurvey.adapters.base import CutoutSpec
from sglsurvey.adapters.irsa_ptf import PtfLevel1Adapter
from sglsurvey.adapters.mast_ps1 import MEAN_COLS, catalog_cone
from sglsurvey.snapshots import SnapshotStore

import dev_search as ds
from coverage_intersect import load_events, relay_apparent, window  # noqa: E402

D = ds.D
RUN = ds.RUN
CONF = RUN / "products" / "conf"
ds.PRODUCTS = CONF                # dev helpers now point at conf tree
Z_GRID = ds.Z_GRID
RING = ds.RING
CUT = 384          # amendment v1.1 (was 256)
GOOD_FRAC_MIN = ds.GOOD_FRAC_MIN

UNITS = [u for u in ds.FREEZE["search_units"]["B"]
         if u["split"] == "confirmatory" and u["class"] == "searched"]


def target_events(tid):
    tab = load_events()["B"]
    return tab[np.asarray([str(x) == tid for x in tab["target_id"]])]


def anti_of(ev):
    return (float(ev["relay_icrs_ra_deg"]), float(ev["relay_icrs_dec_deg"]))


def cal_name(tid):
    return f"conf-{tid}-anti"


def plan():
    adapter = PtfLevel1Adapter()
    store = SnapshotStore(RUN)
    out = {"units": []}
    for u in UNITS:
        tid, band, rad = u["target_id"], u["band"], float(u["radius_au"])
        tab = target_events(tid)
        # per-era representative antipode = median over unit events
        eids = sorted(u["epochs_per_event"])
        evs = {eid: tab[np.asarray([str(x) == eid
                                    for x in tab["event_id"]])][0]
               for eid in eids}
        anti0 = (float(np.median([anti_of(e)[0] for e in evs.values()])),
                 float(np.median([anti_of(e)[1] for e in evs.values()])))
        obs = ds.discover_at(adapter, store, f"conf-{tid}", *anti0)
        # all 0.1 AU windows of the target (off-window exclusion)
        all01 = [window(e, 0.1) for e in tab]
        events = []
        for eid, ev in evs.items():
            lo, hi = window(ev, rad)
            ins = []
            for o in obs:
                if o.band != band:
                    continue
                t = o.t_mid_mjd_utc
                if lo <= t <= hi:
                    pos = {z: relay_apparent(ev, z, t) for z in Z_GRID}
                    fp = adapter.nominal_footprint(o)
                    if any(fp.contains(r_, d_) for r_, d_ in pos.values()):
                        ins.append({"key": ds.nkey(o), "t": t,
                                    "pos": {str(z): list(p)
                                            for z, p in pos.items()}})
            events.append({"event_id": eid,
                           "t_ca_mjd": float(ev["t_ca_mjd"]),
                           "b_min_rsun":
                               float(ev["b_min_au"]) / 0.00465047,
                           "window": [lo, hi], "n_in": len(ins),
                           "in_epochs": ins})
        offs = []
        for o in obs:
            if o.band != band:
                continue
            t = o.t_mid_mjd_utc
            if any(lo <= t <= hi for lo, hi in all01):
                continue
            if adapter.nominal_footprint(o).contains(*anti0):
                offs.append({"key": ds.nkey(o), "t": t})
        out["units"].append({
            "target_id": tid, "band": band, "radius_au": rad,
            "statistics": u["statistics"], "anti": list(anti0),
            "events": events, "n_off": len(offs), "off_epochs": offs})
        print(f"{tid} {rad} {band}: "
              f"{[e['n_in'] for e in events]} in-window, "
              f"{len(offs)} off-window", flush=True)
    (D / "configs" / "conf_plan_v1.json").write_text(
        json.dumps(out, indent=1) + "\n")


def _plan():
    return json.loads((D / "configs" / "conf_plan_v1.json").read_text())


def fetch():
    adapter = PtfLevel1Adapter()
    store = SnapshotStore(RUN)
    p = _plan()
    jobs = {}
    for u in p["units"]:
        tid = u["target_id"]
        omap = {ds.nkey(o): o for o in ds.discover_at(
            adapter, store, f"conf-{tid}", *u["anti"])}
        grp = f"conf-{tid}"
        for e in u["events"]:
            for ep in e["in_epochs"]:
                ra, dec = ep["pos"]["2500.0"]
                jobs[(grp, ep["key"])] = (omap[ep["key"]], ra, dec)
        for ep in u["off_epochs"]:
            jobs.setdefault((grp, ep["key"]),
                            (omap[ep["key"]], u["anti"][0], u["anti"][1]))
    print(f"{len(jobs)} epoch-cutouts to fetch", flush=True)
    ok = fail = 0
    for i, ((grp, key), (o, ra, dec)) in enumerate(sorted(jobs.items())):
        dest = ds._dest(grp, key)
        try:
            adapter.fetch(o, ("sci", "msk"), dest,
                          cutout=CutoutSpec(ra, dec, CUT))
            ok += 1
        except FileNotFoundError:
            print(f"  404 {dest.name}", flush=True)
            fail += 1
        except Exception as exc:
            print(f"  ! {dest.name}: {exc}", flush=True)
            fail += 1
        if (i + 1) % 50 == 0:
            print(f"  {i+1}/{len(jobs)} (ok={ok} fail={fail})", flush=True)
    print(f"cutouts done: {ok}/{len(jobs)} ok, {fail} failed", flush=True)
    for tid, anti in {(u["target_id"], tuple(u["anti"]))
                      for u in p["units"]}:
        dest = RUN / "calibrators" / f"{cal_name(tid)}.json"
        if dest.exists():
            continue
        rows, snaps = catalog_cone("mean", anti[0], anti[1],
                                   ds.CAL_CONE_DEG, MEAN_COLS, store)
        dest.write_text(json.dumps(
            {"center": list(anti), "radius_deg": ds.CAL_CONE_DEG,
             "snapshots": snaps, "rows": rows}) + "\n")
        print(f"calibrators {tid}: {len(rows)}", flush=True)


def unit_stats(adapter, store, u):
    tid, band, rad = u["target_id"], u["band"], u["radius_au"]
    grp = f"conf-{tid}"
    omap = {ds.nkey(o): o for o in ds.discover_at(
        adapter, store, f"conf-{tid}", *u["anti"])}
    gate = defaultdict(int)

    # off-window pool -> k
    off = []
    for ep in u["off_epochs"]:
        fm = ds.fluxmap(grp, ep["key"], omap[ep["key"]], band,
                        cal_name(tid), CUT)
        if not hasattr(fm, "sample"):
            gate[fm[0] if fm else "missing"] += 1
            continue
        f, v, g = fm.sample(u["anti"][0], u["anti"][1])
        if np.isfinite(g[0]) and g[0] >= GOOD_FRAC_MIN:
            off.append((f[0], v[0]))
    k = 1.0
    k_status = f"undefined_n{len(off)}"
    if len(off) >= 6:
        fo = np.array([f for f, _ in off])
        vo = np.array([v for _, v in off])
        for _ in range(3):
            med = np.median(fo)
            sd = 1.4826 * np.median(np.abs(fo - med)) or 1.0
            m = np.abs(fo - med) <= 3 * sd
            fo, vo = fo[m], vo[m]
        if len(fo) >= 6:
            med = np.median(fo)
            k = max(1.0, float(np.median((fo - med) ** 2 / vo) / 0.4549))
            k_status = f"ok_n{len(fo)}"

    # per-event usable in-window flux maps
    ev_fms = {}
    tab = target_events(tid)
    for e in u["events"]:
        fms = []
        for ep in e["in_epochs"]:
            fm = ds.fluxmap(grp, ep["key"], omap[ep["key"]], band,
                            cal_name(tid), CUT)
            if hasattr(fm, "sample"):
                fms.append((omap[ep["key"]], fm))
            else:
                gate[fm[0] if fm else "missing"] += 1
        row = tab[np.asarray([str(x) == e["event_id"]
                              for x in tab["event_id"]])][0]
        ev_fms[e["event_id"]] = (row, fms)

    def ev_samples(eid, z, dx, dy):
        row, fms = ev_fms[eid]
        sm = []
        for o, fm in fms:
            ra, dec = relay_apparent(row, z, o.t_mid_mjd_utc)
            ra += dx / 3600.0 / max(np.cos(np.radians(dec)), .05)
            dec += dy / 3600.0
            f, v, g = fm.sample(ra, dec)
            if np.isfinite(g[0]) and g[0] >= GOOD_FRAC_MIN:
                sm.append((f[0], v[0] * k))
        return sm

    multi = [e["event_id"] for e in u["events"]
             if len(ev_fms[e["event_id"]][1]) >= 2]
    covered = [e["event_id"] for e in u["events"]
               if len(ev_fms[e["event_id"]][1]) >= 1]

    def S_event(dx, dy):
        best = -np.inf
        for eid in multi:
            for z in Z_GRID:
                sm = ev_samples(eid, z, dx, dy)
                if len(sm) >= 2:
                    best = max(best, ds._wstack(sm))
        return best if np.isfinite(best) else np.nan

    def S_stack(dx, dy):
        best = -np.inf
        for z in Z_GRID:
            sm = []
            for eid in covered:
                sm += ev_samples(eid, z, dx, dy)
            if sm:
                best = max(best, ds._wstack(sm))
        return best if np.isfinite(best) else np.nan

    res = {"target_id": tid, "band": band, "radius_au": rad,
           "n_off_usable": len(off), "k_rescale": round(k, 3),
           "k_status": k_status,
           "calibration_gate_failures": dict(gate),
           "events_usable_epochs": {e["event_id"]:
                                    len(ev_fms[e["event_id"]][1])
                                    for e in u["events"]},
           "multi_epoch_events": multi, "statistics": {}}
    for name, fn, applicable in (
            ("S_event", S_event, len(multi) >= 1),
            ("S_stack", S_stack, len(covered) >= 2)):
        if name not in u["statistics"]:
            continue
        if not applicable:
            res["statistics"][name] = {
                "status": "constraint_only_gate_failed_at_measure"}
            continue
        S = fn(0.0, 0.0)
        ctrl = [c for c in (fn(dx, dy) for dx, dy in RING)
                if np.isfinite(c)]
        T = max(ctrl) if ctrl else np.nan
        res["statistics"][name] = {
            "status": ("ok" if len(ctrl) == 8
                       else f"controls_partial_{len(ctrl)}of8"),
            "S": None if np.isnan(S) else round(S, 3),
            "n_controls": len(ctrl),
            "controls": [round(c, 3) for c in ctrl],
            "T": None if not ctrl else round(T, 3),
            "margin": (None if not ctrl or np.isnan(S)
                       else round(S - max(T, 0.0), 3)),
            "exceedance": bool(np.isfinite(S) and ctrl
                               and S > max(T, 0.0)),
        }

    # static annotation
    try:
        from datetime import datetime, timezone
        url = "https://irsa.ipac.caltech.edu/TAP/sync"
        q = ("SELECT ra,dec,ngoodobs FROM ptf_objects WHERE "
             f"CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',"
             f"{u['anti'][0]:.5f},{u['anti'][1]:.5f},0.03))=1")
        resp = requests.get(url, params={"QUERY": q, "FORMAT": "csv"},
                            timeout=300)
        store.store(service_url=url, query=q,
                    request_utc=datetime.now(timezone.utc).isoformat(),
                    response_bytes=resp.content, row_count=None,
                    http_status=resp.status_code)
        import csv as _csv
        import io as _io
        rows = [r for r in _csv.DictReader(_io.StringIO(resp.text))
                if float(r["ngoodobs"] or 0) >= 3]
        dmin = np.inf
        for e in u["events"]:
            for ep in e["in_epochs"]:
                for z in Z_GRID:
                    ra, dec = ep["pos"][str(z)]
                    for r in rows:
                        dmin = min(dmin, np.hypot(
                            (float(r["ra"]) - ra)
                            * np.cos(np.radians(dec)),
                            float(r["dec"]) - dec) * 3600.0)
        res["static_min_arcsec"] = (None if not np.isfinite(dmin)
                                    else round(float(dmin), 2))
    except Exception as exc:
        res["static_min_arcsec"] = f"query_failed: {exc}"
    return res


def stats():
    adapter = PtfLevel1Adapter()
    store = SnapshotStore(RUN)
    p = _plan()
    results = []
    for u in p["units"]:
        print(f"measuring {u['target_id']} {u['radius_au']} {u['band']}",
              flush=True)
        results.append(unit_stats(adapter, store, u))
    n_trials = n_exc = 0
    for r in results:
        for s in r["statistics"].values():
            if "S" in s:
                n_trials += 1
                n_exc += bool(s["exceedance"])
    out = {"n_trials_run": n_trials, "n_exceedances": n_exc,
           "expected_control_crossings": round(n_trials / 9.0, 2),
           "units": results}
    (D / "results" / "confirmatory_v1.json").write_text(
        json.dumps(out, indent=1) + "\n")
    print(json.dumps({"n_trials_run": n_trials, "n_exceedances": n_exc},
                     indent=1))
    for r in results:
        for name, s in r["statistics"].items():
            print(f"  {r['target_id']} {r['radius_au']} {r['band']} "
                  f"{name}: S={s.get('S')} T={s.get('T')} "
                  f"margin={s.get('margin')} "
                  f"{'EXCEEDANCE' if s.get('exceedance') else ''}")


if __name__ == "__main__":
    if "--plan" in sys.argv:
        plan()
    elif "--fetch" in sys.argv:
        fetch()
    elif "--stats" in sys.argv:
        stats()
    else:
        print("usage: confirmatory_search.py --plan|--fetch|--stats")
