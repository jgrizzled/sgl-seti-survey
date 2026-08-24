"""Confirmatory search (threshold freeze v1.0, confirmatory split).

Runs every non-dev search unit under the frozen constructions exactly
as validated by the dev run (results/dev_search_v1.md): channel B
warp-direct track statistic with ring controls, empirical variance
rescale from off-window antipode samples (dev finding P3), skycell
dedup, `track_masked` classification (P1); channel A window-locked
excess vs the same star's off-window baseline with temporal
pseudo-window controls under the frozen gates (>= 2 covered windows,
>= 8 off-window primary epochs, >= 8 valid same-rung offsets), on
per-epoch propagated star positions (linear PM fit through the era
events' tabulated star positions; teegarden moves 29" over the era).
A single-window/single-epoch units are computed and reported
constraint-only, never candidates.

Stages: --plan, --fetch, --stats. Products under
runs/ps1-crossings/products/conf/.
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from sglsurvey.adapters.base import ConeRegion, CutoutSpec
from sglsurvey.adapters.mast_ps1 import Ps1WarpAdapter
from sglsurvey.snapshots import SnapshotStore
from sglsurvey.vetting import load_ps1_mean

import dev_search as ds
from coverage_intersect import load_events, relay_apparent, window  # noqa: E402

D = ds.D
RUN = ds.RUN
CONF_PRODUCTS = RUN / "products" / "conf"
ds.PRODUCTS = CONF_PRODUCTS          # dev_search helpers now point here
FREEZE = ds.FREEZE
COV = ds.COV
Z_GRID = ds.Z_GRID
RING = ds.RING
CUT_IN, CUT_OFF, CUT_A = ds.CUT_IN, ds.CUT_OFF, 128
OFFSETS = [-97, -71, -47, -23, 23, 47, 71, 97]
REDRAWS = [-127, -113, 113, 127]
DEV_B = set(FREEZE["split"]["B_wide"]["dev"])
PLAN_PATH = D / "configs" / "conf_plan_v1.json"

_eventsA = None


def events_a():
    global _eventsA
    if _eventsA is None:
        _eventsA = load_events()["A"]
    return _eventsA


def conf_units_b():
    return [u for u in FREEZE["search_units"]["B"]["units"]
            if u["target_id"] not in DEV_B or u["radius_au"] < 0.1]


def singles_b():
    return [u for u in FREEZE["search_units"]["B"]["single_epoch"]
            if u["target_id"] not in DEV_B]


def units_a():
    return FREEZE["search_units"]["A"]["units"]


def singles_a():
    return FREEZE["search_units"]["A"]["single_window"]


def star_track(tid):
    """Linear (PM) fit through the era events' tabulated star positions:
    returns callable t_mjd -> (ra, dec)."""
    sub = events_a()[np.asarray([str(x) == tid
                                 for x in events_a()["target_id"]])]
    t = np.asarray(sub["t_ca_mjd"], float)
    ra = np.asarray(sub["star_icrs_ra_deg"], float)
    de = np.asarray(sub["star_icrs_dec_deg"], float)
    cosd = np.cos(np.radians(de.mean()))
    px = np.polyfit(t, ra * cosd, 1)
    py = np.polyfit(t, de, 1)

    def at(tm):
        return (float(np.polyval(px, tm)) / cosd, float(np.polyval(py, tm)))
    return at


def a_windows_band(tid, band, rad=0.1):
    sub = COV[(COV["channel"] == "A")
              & np.asarray([str(x) == tid for x in COV["target_id"]])]
    sub = sub[np.asarray(sub["radius_au"]) == rad]
    return [{"event_id": str(r["event_id"]),
             "lo": float(r["t_ca_mjd"] - r["window_days"] / 2),
             "hi": float(r["t_ca_mjd"] + r["window_days"] / 2),
             "n_pri": int(r[f"n_{band}_pri"])} for r in sub]


def valid_offsets(wins):
    """Frozen same-rung rule with epoch support: an offset is valid iff
    no shifted window overlaps any real same-rung window AND the
    pseudo-window set contains >= 2 windows with >= 1 primary epoch
    (epoch support resolved at stats time; here geometric part only)."""
    chosen = []
    for off in OFFSETS + REDRAWS:
        ok = True
        for w in wins:
            slo, shi = w["lo"] + off, w["hi"] + off
            for w2 in wins:
                if slo < w2["hi"] and shi > w2["lo"]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            chosen.append(off)
        if len(chosen) == 8:
            break
    return chosen


_a_obs_cache = {}


def discover_a(adapter, store, tid):
    if tid not in _a_obs_cache:
        at = star_track(tid)
        ra0, de0 = at(56000.0)
        _a_obs_cache[tid] = list(adapter.discover(
            ConeRegion(ra0, de0, 0.13), ds.ERA, store))
    return _a_obs_cache[tid]


def plan():
    adapter = Ps1WarpAdapter(skycell_cache=ds.SKYCELL_CACHE)
    store = SnapshotStore(RUN)
    out = {"B": [], "A": []}

    # ---- channel B (identical construction to dev)
    for u in conf_units_b() + singles_b():
        tid, band, eid = u["target_id"], u["band"], u["event_id"]
        _, sub = ds.target_cone(tid)
        ev = ds.event_row(sub, eid)
        rad = u["radius_au"]
        lo, hi = window(ev, rad)
        anti = (float(ev["relay_icrs_ra_deg"]),
                float(ev["relay_icrs_dec_deg"]))
        obs = ds.discover(adapter, store, tid)
        ins, offs = [], []
        for o in obs:
            if o.band != band or not ds.okq(o):
                continue
            if o.extra["skycell"] not in adapter._skycells:
                continue
            fp = adapter.nominal_footprint(o)
            t = o.t_mid_mjd_utc
            if lo <= t <= hi:
                pos = {z: relay_apparent(ev, z, t) for z in Z_GRID}
                if any(fp.contains(r_, d_) for r_, d_ in pos.values()):
                    ins.append({"key": ds.nkey(o), "t": t,
                                "pos": {str(z): list(p)
                                        for z, p in pos.items()}})
            elif fp.contains(*anti):
                offs.append({"key": ds.nkey(o), "t": t})
        out["B"].append({
            "target_id": tid, "band": band, "event_id": eid,
            "radius_au": rad, "t_ca_mjd": float(ev["t_ca_mjd"]),
            "window": [lo, hi], "anti": list(anti),
            "n_in": len(ins), "n_off": len(offs),
            "in_epochs": ins, "off_epochs": offs,
            "class": ("search_unit" if u.get("n_epochs_max_z", 1) >= 2
                      else "single_epoch")})

    # ---- channel A
    for u in units_a() + singles_a():
        tid, band = u["target_id"], u["band"]
        at = star_track(tid)
        wins = a_windows_band(tid, band)
        offs_geo = valid_offsets(wins)
        obs = discover_a(adapter, store, tid)
        eps = []
        for o in obs:
            if o.band != band or not ds.okq(o):
                continue
            if o.extra["skycell"] not in adapter._skycells:
                continue
            ra, dec = at(o.t_mid_mjd_utc)
            if not adapter.nominal_footprint(o).contains(ra, dec):
                continue
            eps.append({"key": ds.nkey(o), "t": o.t_mid_mjd_utc,
                        "pos": [ra, dec]})
        n_in = sum(1 for e in eps
                   if any(w["lo"] <= e["t"] <= w["hi"] for w in wins))
        out["A"].append({
            "target_id": tid, "band": band, "radius_au": 0.1,
            "windows": wins, "valid_offsets_geometric": offs_geo,
            "n_epochs": len(eps), "n_in": n_in,
            "n_off": len(eps) - n_in, "epochs": eps,
            "class": ("search_unit" if u in units_a()
                      else "single_window")})

    PLAN_PATH.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(
        {"B": [{k: p[k] for k in ("target_id", "band", "event_id",
                                  "radius_au", "n_in", "n_off")}
               for p in out["B"]],
         "A": [{k: p[k] for k in ("target_id", "band", "n_epochs",
                                  "n_in", "n_off")}
               | {"n_offsets_geo": len(p["valid_offsets_geometric"])}
               for p in out["A"]]}, indent=2))


def fetch():
    adapter = Ps1WarpAdapter(skycell_cache=ds.SKYCELL_CACHE)
    store = SnapshotStore(RUN)
    plan_doc = json.loads(PLAN_PATH.read_text())
    jobs = {}
    for p in plan_doc["B"]:
        tid = p["target_id"]
        omap = {ds.nkey(o): o for o in ds.discover(adapter, store, tid)}
        for e in p["in_epochs"]:
            jobs.setdefault((tid, e["key"], CUT_IN),
                            (omap[e["key"]], *e["pos"]["2500.0"]))
        for e in p["off_epochs"]:
            jobs.setdefault((tid, e["key"], CUT_OFF),
                            (omap[e["key"]], *p["anti"]))
    for p in plan_doc["A"]:
        tid = p["target_id"]
        omap = {ds.nkey(o): o
                for o in discover_a(adapter, store, tid)}
        for e in p["epochs"]:
            jobs.setdefault((tid, e["key"], CUT_A),
                            (omap[e["key"]], *e["pos"]))
    print(f"{len(jobs)} epoch-cutouts to fetch", flush=True)
    ok = fail = 0
    for i, ((tid, key, size), (o, ra, dec)) in enumerate(sorted(jobs.items())):
        dest = ds._dest(tid, key)
        cut = CutoutSpec(ra_deg=ra, dec_deg=dec, size_pix=size)
        try:
            for kind in ("img", "wt", "msk"):
                adapter.fetch(o, (kind,), dest, cutout=cut)
            ok += 1
        except FileNotFoundError:
            print(f"  404 {dest.name}", flush=True)
            fail += 1
        except Exception as exc:
            print(f"  ! {dest.name}: {exc}", flush=True)
            fail += 1
        if (i + 1) % 50 == 0:
            print(f"  {i+1}/{len(jobs)} (ok={ok} fail={fail})", flush=True)
    print(f"done: {ok}/{len(jobs)} ok, {fail} failed", flush=True)


def robust_baseline(f):
    f = np.asarray(f, float)
    f = f[np.isfinite(f)]
    for _ in range(3):
        med = np.median(f)
        sd = 1.4826 * np.median(np.abs(f - med)) or 1.0
        f = f[np.abs(f - med) <= 3 * sd]
    return float(np.median(f)), f


def stats():
    adapter = Ps1WarpAdapter(skycell_cache=ds.SKYCELL_CACHE)
    store = SnapshotStore(RUN)
    plan_doc = json.loads(PLAN_PATH.read_text())
    results = {"B": [], "A": []}
    cat_cache = {}

    # ---------------- channel B (dev-validated construction)
    for p in plan_doc["B"]:
        tid, band, eid = p["target_id"], p["band"], p["event_id"]
        _, sub = ds.target_cone(tid)
        ev = ds.event_row(sub, eid)
        omap = {ds.nkey(o): o for o in ds.discover(adapter, store, tid)}
        off = []
        zp_src = defaultdict(int)
        for e in p["off_epochs"]:
            o = omap.get(e["key"])
            if o is None:
                continue
            fm = ds._fluxmap(tid, o, CUT_OFF)
            if fm is None:
                continue
            zp_src[fm.zp_source] += 1
            f, v, g = fm.sample(p["anti"][0], p["anti"][1])
            if np.isfinite(g[0]) and g[0] >= ds.GOOD_FRAC_MIN:
                off.append((round(o.t_mid_mjd_utc / ds.DUP_BIN_DAYS),
                            g[0], f[0], v[0]))
        offd = ds._dedup_best(off)
        k, k_status = 1.0, f"undefined_n{len(offd)}"
        if len(offd) >= 6:
            fo = np.array([f for f, _ in offd])
            vo = np.array([v for _, v in offd])
            keep = np.isfinite(fo) & np.isfinite(vo) & (vo > 0)
            fo, vo = fo[keep], vo[keep]
            for _ in range(3):
                med = np.median(fo)
                sd = 1.4826 * np.median(np.abs(fo - med)) or 1.0
                m = np.abs(fo - med) <= 3 * sd
                fo, vo = fo[m], vo[m]
            if len(fo) >= 6:
                med = np.median(fo)
                k = max(1.0, float(np.median((fo - med) ** 2 / vo)
                                   / 0.4549))
                k_status = f"ok_n{len(fo)}"
        infms = []
        for e in p["in_epochs"]:
            o = omap.get(e["key"])
            if o is None:
                continue
            fm = ds._fluxmap(tid, o, CUT_IN)
            if fm is not None:
                zp_src[fm.zp_source] += 1
                infms.append((o, fm))

        def track_S(dx, dy):
            best, best_z = -np.inf, None
            for z in Z_GRID:
                sm = []
                for o, fm in infms:
                    ra, dec = relay_apparent(ev, z, o.t_mid_mjd_utc)
                    ra += dx / 3600.0 / max(np.cos(np.radians(dec)), .05)
                    dec += dy / 3600.0
                    f, v, g = fm.sample(ra, dec)
                    if np.isfinite(g[0]) and g[0] >= ds.GOOD_FRAC_MIN:
                        sm.append((round(o.t_mid_mjd_utc / ds.DUP_BIN_DAYS),
                                   g[0], f[0], v[0] * k))
                smd = ds._dedup_best(sm)
                if smd:
                    S, _ = ds._snr_stack(smd)
                    if S > best:
                        best, best_z = S, z
            return (best if np.isfinite(best) else np.nan), best_z

        S, S_z = track_S(0.0, 0.0)
        controls = [c for c, _ in (track_S(dx, dy) for dx, dy in RING)
                    if np.isfinite(c)]
        T = max(controls) if controls else np.nan
        corr = ds.corridor_name(tid)
        if corr not in cat_cache:
            try:
                cat_cache[corr] = load_ps1_mean(ds.SCREEN, corr, band="r")
            except Exception:
                cat_cache[corr] = None
        cat = cat_cache[corr]
        static_min = None
        if cat is not None:
            dmin = np.inf
            for e in p["in_epochs"]:
                for z in Z_GRID:
                    ra, dec = e["pos"][str(z)]
                    d = np.hypot((cat.ra - ra) * np.cos(np.radians(dec)),
                                 cat.dec - dec) * 3600.0
                    m = cat.ndet >= 3
                    if m.any():
                        dmin = min(dmin, float(np.min(d[m])))
            static_min = None if not np.isfinite(dmin) else round(dmin, 2)
        n_dedup = len(ds._dedup_best(
            [(round(o.t_mid_mjd_utc / ds.DUP_BIN_DAYS), 1.0, 0.0, 1.0)
             for o, _ in infms]))
        status = ("track_masked" if np.isnan(S)
                  else "ok" if len(controls) == 8
                  else f"controls_partial_{len(controls)}of8")
        exceed = bool(np.isfinite(S) and S > max(
            T if np.isfinite(T) else -np.inf, 0.0))
        results["B"].append({
            "target_id": tid, "band": band, "event_id": eid,
            "radius_au": p["radius_au"], "t_ca_mjd": p["t_ca_mjd"],
            "status": status,
            "n_in_epochs_raw": len(p["in_epochs"]),
            "n_in_exposures_dedup": n_dedup,
            "n_off_samples": len(offd),
            "k_rescale": round(k, 3), "k_status": k_status,
            "zp_sources": dict(zp_src),
            "S": None if np.isnan(S) else round(S, 3),
            "S_best_z": S_z,
            "n_controls": len(controls),
            "controls": [round(c, 3) for c in controls],
            "T": None if not controls else round(T, 3),
            "margin": (None if not controls or np.isnan(S)
                       else round(S - max(T, 0.0), 3)),
            "static_min_arcsec": static_min,
            "static_annotated": (static_min is not None
                                 and static_min < ds.STATIC_ARCSEC),
            "class": ("search_unit" if n_dedup >= 2 else "single_epoch"),
            "exceedance": exceed})

    # ---------------- channel A
    for p in plan_doc["A"]:
        tid, band = p["target_id"], p["band"]
        omap = {ds.nkey(o): o
                for o in discover_a(adapter, store, tid)}
        wins = p["windows"]
        samples = []       # (t, f, v) scaled
        zp_src = defaultdict(int)
        for e in p["epochs"]:
            o = omap.get(e["key"])
            if o is None:
                continue
            fm = ds._fluxmap(tid, o, CUT_A)
            if fm is None:
                continue
            zp_src[fm.zp_source] += 1
            f, v, g = fm.sample(e["pos"][0], e["pos"][1])
            if np.isfinite(g[0]) and g[0] >= ds.GOOD_FRAC_MIN:
                samples.append((o.t_mid_mjd_utc, f[0], v[0]))
        in_any = [s for s in samples
                  if any(w["lo"] <= s[0] <= w["hi"] for w in wins)]
        off_s = [s for s in samples
                 if not any(w["lo"] <= s[0] <= w["hi"] for w in wins)]
        base, kept = robust_baseline([f for _, f, _ in off_s])
        koff = [(f, v) for t, f, v in off_s
                if np.isfinite(f) and f in kept] or \
               [(f, v) for t, f, v in off_s]
        k = 1.0
        if len(koff) >= 6:
            fo = np.array([f for f, _ in koff])
            vo = np.array([v for _, v in koff])
            k = max(1.0, float(np.median((fo - base) ** 2 / vo) / 0.4549))

        def span_S(offset):
            per_win = []
            for w in wins:
                sm = [(f - base, v * k) for t, f, v in samples
                      if w["lo"] + offset <= t <= w["hi"] + offset]
                if sm:
                    per_win.append(ds._snr_stack(sm))
            if not per_win:
                return np.nan, 0
            Sws = np.array([s for s, _ in per_win])
            return float(np.sum(Sws) / np.sqrt(len(Sws))), len(per_win)

        S, nw = span_S(0.0)
        valid, controls = [], []
        for off in p["valid_offsets_geometric"]:
            Sc, nc = span_S(off)
            if np.isfinite(Sc) and nc >= 2:
                valid.append(off)
                controls.append(Sc)
        T = max(controls) if controls else np.nan
        gates = {"n_windows_with_data": nw,
                 "n_off_epochs": len(off_s),
                 "n_valid_offsets": len(valid)}
        searchable = (p["class"] == "search_unit" and nw >= 2
                      and len(off_s) >= 8 and len(valid) >= 8)
        exceed = bool(searchable and np.isfinite(S)
                      and S > max(T if np.isfinite(T) else -np.inf, 0.0))
        results["A"].append({
            "target_id": tid, "band": band, "class": p["class"],
            "gates": gates,
            "status": ("searchable" if searchable else "constraint_only"),
            "baseline_zp25": round(base, 2) if np.isfinite(base) else None,
            "k_rescale": round(k, 3),
            "zp_sources": dict(zp_src),
            "S": None if np.isnan(S) else round(S, 3),
            "controls": [round(c, 3) for c in controls],
            "T": None if not controls else round(T, 3),
            "margin": (None if not controls or np.isnan(S)
                       else round(S - max(T, 0.0), 3)),
            "exceedance": exceed})

    out = D / "results" / "confirmatory_v1.json"
    out.write_text(json.dumps(results, indent=1) + "\n")
    nB = sum(1 for r in results["B"] if r["class"] == "search_unit"
             and r["status"] != "track_masked")
    nA = sum(1 for r in results["A"] if r["status"] == "searchable")
    print(json.dumps({
        "B_rows": len(results["B"]), "B_searchable": nB,
        "A_rows": len(results["A"]), "A_searchable": nA,
        "exceedances": [
            f"{r['target_id']}/{r['band']}/{r.get('event_id', 'A')}"
            for ch in ("B", "A") for r in results[ch]
            if r["exceedance"]],
        "track_masked": [
            f"{r['target_id']}/{r['band']}/{r['event_id']}"
            for r in results["B"] if r["status"] == "track_masked"],
    }, indent=2))


if __name__ == "__main__":
    if "--plan" in sys.argv:
        plan()
    elif "--fetch" in sys.argv:
        fetch()
    elif "--stats" in sys.argv:
        stats()
    else:
        print("usage: confirmatory_search.py --plan | --fetch | --stats")
