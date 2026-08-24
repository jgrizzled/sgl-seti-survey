"""Dev-target search (threshold freeze v1.0, dev split only).

Stages: --plan (enumerate windows, pseudo-window validity, needed
epochs; write dev_plan_v1.json), --fetch (resumable sci/diff/msk cutout
downloads), --stats (matched-filter photometry, frozen statistics,
dev_search_v1.json/md). The frozen rules are applied literally; where a
construction fails (e.g. no valid pseudo-window offsets) the unit is
reported control_undefined, never silently patched.
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from sglsurvey.adapters.base import ConeRegion, CutoutSpec, MjdRange
from sglsurvey.adapters.irsa_ztf import (MASK_FATAL_TEMPLATE,
                                         ZtfNominalFootprint, ZtfSciAdapter)
from sglsurvey.photometry import build_flux_map_ztf
from sglsurvey.snapshots import SnapshotStore

D = REPO / "surveys" / "ztf-crossings"
FREEZE = json.loads((D / "configs" / "threshold_freeze_v1.json").read_text())
COV = Table.read(D / "results" / "coverage_v1_events.ecsv")
EVENTS = Table.read(REPO / "crossings" / "universal_v1" / "events.ecsv")
RUN = REPO / "runs" / "ztf-crossings"
ERA = MjdRange(58178.0, 61275.0)
KM_PER_AU = 1.495978707e8
OFFSETS = [-97, -71, -47, -23, 23, 47, 71, 97]
REDRAWS = [-127, -113, 113, 127]
BANDCODE = {"g": "zg", "r": "zr", "i": "zi"}
Z_GRID = (550.0, 1000.0, 2500.0, 5500.0, 10000.0)
CUT_A, CUT_B = 96, 160          # cutout size (pix)
WEIGHT_CAP = 20.0

sys.path.insert(0, str(D / "scripts"))
from coverage_intersect import relay_apparent  # noqa: E402  (same repo)

DEV_A = FREEZE["split"]["A"]["dev"]
DEV_B = FREEZE["split"]["B_wide"]["dev"]


def half_days(b, r, vperp):
    return float(np.sqrt(r * r - b * b) * KM_PER_AU / vperp / 86400.0)


def target_windows(tid):
    """All real windows (both rungs) for one channel-A target: list of
    (lo, hi, radius, event row)."""
    sub = COV[(COV["channel"] == "A")
              & np.asarray([str(x) == tid for x in COV["target_id"]])]
    wins = []
    for row in sub:
        wins.append((row["t_ca_mjd"] - row["window_days"] / 2,
                     row["t_ca_mjd"] + row["window_days"] / 2,
                     float(row["radius_au"]), str(row["event_id"])))
    return wins


def valid_offsets(tid, rad):
    """Frozen pseudo-window rule: an offset is valid if NO shifted real
    window of this rung overlaps ANY real window of either rung."""
    all_wins = target_windows(tid)
    rung = [w for w in all_wins if w[2] == rad]
    chosen, pool = [], OFFSETS + REDRAWS
    for off in pool:
        ok = True
        for lo, hi, _, _ in rung:
            slo, shi = lo + off, hi + off
            if shi < ERA.start_mjd_utc or slo > ERA.stop_mjd_utc:
                continue
            for wlo, whi, _, _ in all_wins:
                if slo < whi and shi > wlo:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            chosen.append(off)
        if len(chosen) == 8:
            break
    return chosen


def plan():
    adapter, store = ZtfSciAdapter(), SnapshotStore(RUN)
    units_A = [u for u in FREEZE["search_units"]["A"]["units"]
               if u["target_id"] in DEV_A]
    units_B = [u for u in FREEZE["search_units"]["B"]["units"]
               if u["target_id"] in DEV_B and u["radius_au"] == 0.1]
    ev_by_id = {str(e["event_id"]): e for e in EVENTS}

    plan = {"A": [], "B": [], "epochs": {}}
    need = defaultdict(set)     # (tid, channel) -> set of native pid keys
    obs_cache = {}

    def cone_obs(tid, ch, ra, dec, radius):
        key = (tid, ch)
        if key not in obs_cache:
            obs_cache[key] = list(adapter.discover(
                ConeRegion(ra, dec, radius), ERA, store))
        return obs_cache[key]

    for u in units_A:
        tid, band, rad = u["target_id"], u["band"], u["radius_au"]
        sub = COV[(COV["channel"] == "A")
                  & np.asarray([str(x) == tid for x in COV["target_id"]])]
        sub = sub[np.asarray(sub["radius_au"]) == rad]
        ev0 = ev_by_id[str(sub[0]["event_id"])]
        ra, dec = float(ev0["star_icrs_ra_deg"]), float(ev0["star_icrs_dec_deg"])
        offs = valid_offsets(tid, rad)
        wins = [(row["t_ca_mjd"] - row["window_days"] / 2,
                 row["t_ca_mjd"] + row["window_days"] / 2,
                 str(row["event_id"])) for row in sub]
        spans = [(lo, hi) for lo, hi, _ in wins]
        spans += [(lo + o, hi + o) for o in offs for lo, hi, _ in wins]
        obs = cone_obs(tid, "A", ra, dec, 0.05)
        eps = []
        for o in obs:
            if o.band != BANDCODE[band]:
                continue
            q = o.quality_flags
            if q.get("bad_quality") or (q.get("seeing") or 99) > 4.0:
                continue
            t = o.t_mid_mjd_utc
            if not any(lo <= t <= hi for lo, hi in spans):
                continue
            if not ZtfNominalFootprint(o).contains(ra, dec):
                continue
            eps.append(o)
            need[(tid, "A")].add(json.dumps(o.native_key, sort_keys=True))
        plan["A"].append({
            "target_id": tid, "band": band, "radius_au": rad,
            "ra": ra, "dec": dec, "n_windows": len(wins),
            "valid_offsets": offs, "control_defined": len(offs) == 8,
            "windows": [[lo, hi, eid] for lo, hi, eid in wins],
            "n_epochs": len(eps)})

    for u in units_B + [u for u in FREEZE["search_units"]["B"]["single_epoch"]
                        if u["target_id"] in DEV_B and u["radius_au"] == 0.1]:
        tid, band = u["target_id"], u["band"]
        eid = u["event_id"]
        row = COV[(COV["channel"] == "B")
                  & np.asarray([str(x) == eid for x in COV["event_id"]])][0]
        ev = ev_by_id[eid]
        lo = row["t_ca_mjd"] - row["window_days"] / 2
        hi = row["t_ca_mjd"] + row["window_days"] / 2
        ra0, de0 = float(ev["relay_icrs_ra_deg"]), float(ev["relay_icrs_dec_deg"])
        obs = cone_obs(tid, "B", ra0, de0, 0.15)
        eps = []
        for o in obs:
            if o.band != BANDCODE[band]:
                continue
            q = o.quality_flags
            if q.get("bad_quality") or (q.get("seeing") or 99) > 4.0:
                continue
            if not (lo <= o.t_mid_mjd_utc <= hi):
                continue
            pos = [relay_apparent(ev, z, o.t_mid_mjd_utc) for z in Z_GRID]
            if not any(ZtfNominalFootprint(o).contains(r_, d_)
                       for r_, d_ in pos):
                continue
            eps.append((o, pos))
            need[(tid, "B")].add(json.dumps(o.native_key, sort_keys=True))
        plan["B"].append({
            "target_id": tid, "band": band, "event_id": eid,
            "t_ca_mjd": float(row["t_ca_mjd"]), "window": [float(lo), float(hi)],
            "n_epochs": len(eps),
            "positions": {json.dumps(o.native_key, sort_keys=True):
                          pos for o, pos in eps}})

    plan["fetch_totals"] = {f"{k[0]}/{k[1]}": len(v) for k, v in need.items()}
    plan["n_products"] = 3 * sum(len(v) for v in need.values())
    (D / "configs" / "dev_plan_v1.json").write_text(
        json.dumps(plan, indent=1) + "\n")
    print(json.dumps({"A_units": len(plan["A"]),
                      "B_rows": len(plan["B"]),
                      "control_defined_A": sum(p["control_defined"]
                                               for p in plan["A"]),
                      "valid_offset_counts": {p["target_id"] + "/" + p["band"] +
                                              f"/{p['radius_au']}":
                                              len(p["valid_offsets"])
                                              for p in plan["A"]},
                      "fetch_totals": plan["fetch_totals"],
                      "n_products": plan["n_products"]}, indent=2))


if __name__ == "__main__":
    if "--plan" in sys.argv:
        plan()


# ---------------------------------------------------------------- fetch/stats

PRODUCTS_DIR = RUN / "products" / "dev"


def _discover_all():
    adapter, store = ZtfSciAdapter(), SnapshotStore(RUN)
    plan_doc = json.loads((D / "configs" / "dev_plan_v1.json").read_text())
    ev_by_id = {str(e["event_id"]): e for e in EVENTS}
    cones = {}
    for p in plan_doc["A"]:
        cones[(p["target_id"], "A")] = (p["ra"], p["dec"], 0.05)
    for p in plan_doc["B"]:
        ev = ev_by_id[p["event_id"]]
        cones.setdefault((p["target_id"], "B"),
                         (float(ev["relay_icrs_ra_deg"]),
                          float(ev["relay_icrs_dec_deg"]), 0.15))
    obs_by_key = {}
    for (tid, ch), (ra, dec, rad) in cones.items():
        obs = list(adapter.discover(ConeRegion(ra, dec, rad), ERA, store))
        obs_by_key[(tid, ch)] = {json.dumps(o.native_key, sort_keys=True): o
                                 for o in obs}
    return adapter, plan_doc, obs_by_key


def _epoch_dest(tid, ch, o):
    return PRODUCTS_DIR / f"{ch}-{tid}" / str(o.native_key["pid"])


def _needed(plan_doc, obs_by_key):
    out = []
    for p in plan_doc["A"]:
        tid = p["target_id"]
        omap = obs_by_key[(tid, "A")]
        spans = [(lo, hi) for lo, hi, _ in p["windows"]]
        spans += [(lo + o, hi + o) for o in p["valid_offsets"]
                  for lo, hi, _ in p["windows"]]
        for k, o in omap.items():
            if o.band != BANDCODE[p["band"]]:
                continue
            q = o.quality_flags
            if q.get("bad_quality") or (q.get("seeing") or 99) > 4.0:
                continue
            if not any(lo <= o.t_mid_mjd_utc <= hi for lo, hi in spans):
                continue
            if not ZtfNominalFootprint(o).contains(p["ra"], p["dec"]):
                continue
            out.append((tid, "A", o, p["ra"], p["dec"], CUT_A))
    for p in plan_doc["B"]:
        tid = p["target_id"]
        omap = obs_by_key[(tid, "B")]
        for k, pos in p["positions"].items():
            o = omap.get(k)
            if o is None or o.band != BANDCODE[p["band"]]:
                continue
            ra, dec = pos[2]        # z=2500 center
            out.append((tid, "B", o, ra, dec, CUT_B))
    # dedup per (dest, size): same epoch may serve several units
    seen, ded = set(), []
    for tid, ch, o, ra, dec, size in out:
        key = (str(_epoch_dest(tid, ch, o)), size)
        if key not in seen:
            seen.add(key)
            ded.append((tid, ch, o, ra, dec, size))
    return ded


def fetch():
    adapter, plan_doc, obs_by_key = _discover_all()
    need = _needed(plan_doc, obs_by_key)
    print(f"{len(need)} epochs to fetch", flush=True)
    ok = missing_diff = failed = 0
    for i, (tid, ch, o, ra, dec, size) in enumerate(need):
        dest = _epoch_dest(tid, ch, o)
        cut = CutoutSpec(ra_deg=ra, dec_deg=dec, size_pix=size)
        got = {}
        for kind in ("sci", "msk", "diff"):
            try:
                ps = adapter.fetch(o, (kind,), dest, cutout=cut)
                got[kind] = ps.products[0].path
            except FileNotFoundError:
                if kind == "diff":
                    missing_diff += 1
                else:
                    got = None
                    break
            except Exception as e:      # transient; retry next run
                print(f"  ! {dest.name}/{kind}: {e}", flush=True)
                got = None
                failed += 1
                break
        if got:
            ok += 1
        if (i + 1) % 200 == 0:
            print(f"  {i+1}/{len(need)} (ok={ok} nodiff={missing_diff} "
                  f"fail={failed})", flush=True)
    print(f"done: {ok}/{len(need)} ok, {missing_diff} without diff, "
          f"{failed} failed", flush=True)


def _fluxmap(tid, ch, o):
    dest = _epoch_dest(tid, ch, o)
    if not dest.exists():
        return None
    paths = {p.name.split("-", 1)[0].replace(f"cut{CUT_A}", "").replace(
        f"cut{CUT_B}", ""): p for p in dest.iterdir()}
    sci = next((p for p in dest.iterdir() if p.name.endswith("sciimg.fits")), None)
    msk = next((p for p in dest.iterdir() if p.name.endswith("mskimg.fits")), None)
    diff = next((p for p in dest.iterdir()
                 if p.name.endswith("scimrefdiffimg.fits.fz")), None)
    if sci is None or msk is None:
        return None
    try:
        return build_flux_map_ztf(sci, diff, msk, MASK_FATAL_TEMPLATE,
                                  o.band, o.t_mid_mjd_utc)
    except Exception:
        return None


def _snr_stack(samples):
    """samples: list of (flux, var). Returns (S, n_eff) with WEIGHT_CAP."""
    f = np.array([s[0] for s in samples], float)
    v = np.array([s[1] for s in samples], float)
    okm = np.isfinite(f) & np.isfinite(v) & (v > 0)
    if not okm.any():
        return np.nan, 0.0
    f, v = f[okm], v[okm]
    w = 1.0 / v
    w = np.minimum(w, WEIGHT_CAP * np.median(w))
    S = float(np.sum(w * f) / np.sqrt(np.sum(w)))
    n_eff = float(np.sum(w) ** 2 / np.sum(w ** 2))
    return S, n_eff


def stats():
    adapter, plan_doc, obs_by_key = _discover_all()
    ev_by_id = {str(e["event_id"]): e for e in EVENTS}
    results = {"A": [], "B": []}

    for p in plan_doc["A"]:
        tid, band = p["target_id"], p["band"]
        omap = obs_by_key[(tid, "A")]
        eps = []
        for k, o in omap.items():
            if o.band != BANDCODE[band]:
                continue
            q = o.quality_flags
            if q.get("bad_quality") or (q.get("seeing") or 99) > 4.0:
                continue
            if not ZtfNominalFootprint(o).contains(p["ra"], p["dec"]):
                continue
            eps.append(o)

        def span_S(spans):
            per_win = []
            for lo, hi in spans:
                sm = []
                for o in eps:
                    if lo <= o.t_mid_mjd_utc <= hi:
                        fm = _fluxmap(tid, "A", o)
                        if fm is None:
                            continue
                        f, v, g = fm.sample(p["ra"], p["dec"])
                        if np.isfinite(g[0]) and g[0] >= 0.7:
                            sm.append((f[0], v[0]))
                if sm:
                    per_win.append(_snr_stack(sm))
            if not per_win:
                return np.nan, 0
            Sws = np.array([s for s, _ in per_win])
            return float(np.sum(Sws) / np.sqrt(len(Sws))), len(per_win)

        real = [(lo, hi) for lo, hi, _ in p["windows"]]
        S, nw = span_S(real)
        controls = []
        for off in p["valid_offsets"]:
            Sc, nc = span_S([(lo + off, hi + off) for lo, hi in real])
            if np.isfinite(Sc):
                controls.append(Sc)
        T = max(controls) if controls else np.nan
        results["A"].append({
            "target_id": tid, "band": band, "radius_au": p["radius_au"],
            "control_defined": p["control_defined"],
            "n_valid_offsets": len(p["valid_offsets"]),
            "S": None if np.isnan(S) else round(S, 3),
            "n_windows_with_data": nw,
            "controls_partial": [round(c, 3) for c in controls],
            "T_partial": None if not controls else round(T, 3),
            "R_partial": (None if not controls or np.isnan(S) or T <= 0
                          else round(S / T, 3)),
            "status": "control_undefined (diagnostic only)"})

    RING = [(20.0, 0.0), (-20.0, 0.0), (30.0, 0.0), (-30.0, 0.0),
            (40.0, 0.0), (-40.0, 0.0), (0.0, 30.0), (0.0, -30.0)]
    for p in plan_doc["B"]:
        tid, band, eid = p["target_id"], p["band"], p["event_id"]
        ev = ev_by_id[eid]
        omap = obs_by_key[(tid, "B")]
        fms = []
        for k in p["positions"]:
            o = omap.get(k)
            if o is None or o.band != BANDCODE[band]:
                continue
            fm = _fluxmap(tid, "B", o)
            if fm is not None:
                fms.append((o, fm))

        def track_S(dx_arcsec, dy_arcsec):
            best = -np.inf
            for z in Z_GRID:
                sm = []
                for o, fm in fms:
                    ra, dec = relay_apparent(ev, z, o.t_mid_mjd_utc)
                    ra += dx_arcsec / 3600.0 / max(np.cos(np.radians(dec)), .05)
                    dec += dy_arcsec / 3600.0
                    f, v, g = fm.sample(ra, dec)
                    if np.isfinite(g[0]) and g[0] >= 0.7:
                        sm.append((f[0], v[0]))
                if sm:
                    S, _ = _snr_stack(sm)
                    best = max(best, S)
            return best if np.isfinite(best) else np.nan

        S = track_S(0.0, 0.0)
        controls = [track_S(dx, dy) for dx, dy in RING]
        controls = [c for c in controls if np.isfinite(c)]
        T = max(controls) if controls else np.nan
        n_ep = len(fms)
        results["B"].append({
            "target_id": tid, "band": band, "event_id": eid,
            "t_ca_mjd": p["t_ca_mjd"], "n_epochs": n_ep,
            "S": None if np.isnan(S) else round(S, 3),
            "n_controls": len(controls),
            "T": None if np.isnan(T) else round(T, 3),
            "R": (None if np.isnan(S) or np.isnan(T) or T <= 0
                  else round(S / T, 3)),
            "class": ("search_unit" if n_ep >= 2 else "single_epoch"),
            "exceedance": bool(np.isfinite(S) and np.isfinite(T)
                               and T > 0 and S / T > 1)})

    out = D / "results" / "dev_search_v1.json"
    out.write_text(json.dumps(results, indent=1) + "\n")
    print(json.dumps(results, indent=1))


if __name__ == "__main__":
    if "--fetch" in sys.argv:
        fetch()
    elif "--stats" in sys.argv:
        stats()
