"""PS1 crossings coverage intersection v1 (hypotheses.md freeze v1.0).

Port of surveys/ztf-crossings/scripts/coverage_intersect.py onto the
MAST PS1 DR2 warp adapter. Intersects the universal Earth-center
crossing list with actual warp coverage: one hex-sampled ps1filenames
discovery cone per (channel, target) over the full era, snapshotted
under runs/ps1-crossings; then a local intersection: a warp covers an
event at ladder radius r if its t_mid lies in t_ca +/- sqrt(r^2-b^2)/
v_perp and its skycell footprint (WCS from one cached mask per
skycell; cache seeded from the Pipeline-A corridor cache) contains the
predicted source position (channel A: the star; channel B: the
apparent relay position for each z on the 5-point grid at that epoch).
Coverage counting only - no pixels are touched and no signal statistic
is formed. Skycell footprints are nominal: Pipeline-A experience is
~75% of nominal hits usable under the exact warp mask, and that
attrition is applied at the search stage, not here.

Extra columns vs the ZTF version: per-band counts for all of grizy,
and per-band TTI-pair counts (in-window epochs with a same-band
neighbor within 0.05 d) - the raw material of the pair-motion veto.
"""

from __future__ import annotations

import json
import shutil
import sys
import time as _time
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.coordinates import get_body_barycentric
from astropy.table import Table
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from sglsurvey.adapters.base import ConeRegion, MjdRange
from sglsurvey.adapters.mast_ps1 import Ps1WarpAdapter
from sglsurvey.snapshots import SnapshotStore

ERA = MjdRange(54900.0, 57300.0)          # DR2 warps actually 54985..57067
KM_PER_AU = 1.495978707e8
Z_GRID_AU = (550.0, 1000.0, 2500.0, 5500.0, 10000.0)
LADDER = {"A": (0.1, 1.0), "B": (0.0056, 0.0116, 0.1)}
CONE_MARGIN_DEG = {"A": 0.13, "B": 0.15}   # A: >= one hex ring so
# overlapping neighbor skycells are listed; B: in-window z-track stays
# within ~40 arcsec of the tabulated relay direction, margin dominates
BANDS = ("g", "r", "i", "z", "y")
PAIR_DT_DAYS = 0.05
OUT = REPO / "surveys" / "ps1-crossings" / "results"
RUN = REPO / "runs" / "ps1-crossings"
MSK_DIR = RUN / "products" / "msk"
SKYCELL_CACHE = RUN / "skycell_wcs.json"
SEED_CACHE = REPO / "runs" / "panstarrs" / "skycell_wcs.json"


def load_events():
    t = Table.read(REPO / "crossings" / "universal_v1" / "events.ecsv")
    mjd = Time(list(t["t_ca_utc"]), format="isot", scale="utc").mjd
    t["t_ca_mjd"] = mjd
    era = (mjd >= ERA.start_mjd_utc) & (mjd <= ERA.stop_mjd_utc)
    side = np.asarray(t["axis_distance_au"])
    tA = t[era & (t["link_direction"] == "inbound") & (side > 0)]
    tB = t[era & (t["link_direction"] == "outbound") & (side < 0)]
    tA = tA[np.asarray(tA["b_min_au"]) < max(LADDER["A"])]
    tB = tB[np.asarray(tB["b_min_au"]) < max(LADDER["B"])]
    return {"A": tA, "B": tB}


def source_radec(ch, ev):
    if ch == "A":
        return float(ev["star_icrs_ra_deg"]), float(ev["star_icrs_dec_deg"])
    return float(ev["relay_icrs_ra_deg"]), float(ev["relay_icrs_dec_deg"])


def relay_apparent(ev, z_au, t_mjd):
    """Apparent ICRS ra/dec of a relay at z_au on the anti-star axis."""
    t = Time(t_mjd, format="mjd", scale="utc")
    sun = get_body_barycentric("sun", t).xyz.to_value("AU")
    earth = get_body_barycentric("earth", t).xyz.to_value("AU")
    ra = np.radians(float(ev["star_icrs_ra_deg"]))
    de = np.radians(float(ev["star_icrs_dec_deg"]))
    u_star = np.array([np.cos(de) * np.cos(ra), np.cos(de) * np.sin(ra),
                       np.sin(de)])
    v = (sun - z_au * u_star) - earth
    v /= np.linalg.norm(v)
    return (float(np.degrees(np.arctan2(v[1], v[0])) % 360.0),
            float(np.degrees(np.arcsin(np.clip(v[2], -1, 1)))))


def window(ev, r):
    b = float(ev["b_min_au"])
    half_d = (np.sqrt(r * r - b * b) * KM_PER_AU
              / float(ev["v_perp_km_s"]) / 86400.0)
    return float(ev["t_ca_mjd"]) - half_d, float(ev["t_ca_mjd"]) + half_d


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if SEED_CACHE.exists() and not SKYCELL_CACHE.exists():
        SKYCELL_CACHE.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(SEED_CACHE, SKYCELL_CACHE)
        print(f"seeded skycell WCS cache from {SEED_CACHE} "
              f"({len(json.loads(SKYCELL_CACHE.read_text()))} skycells)")
    events = load_events()
    adapter = Ps1WarpAdapter(skycell_cache=SKYCELL_CACHE)
    store = SnapshotStore(RUN)

    # one discovery cone per (channel, target)
    obs_by_key = {}
    for ch, tab in events.items():
        for tid in np.unique(tab["target_id"]):
            sub = tab[tab["target_id"] == tid]
            ras = np.array([source_radec(ch, e)[0] for e in sub])
            des = np.array([source_radec(ch, e)[1] for e in sub])
            ra0, de0 = float(np.median(ras)), float(np.median(des))
            cosd = max(np.cos(np.radians(de0)), 0.05)
            spread = max(np.ptp(ras) * cosd, np.ptp(des)) / 2.0
            radius = spread + CONE_MARGIN_DEG[ch]
            region = ConeRegion(ra0, de0, radius)
            t0 = _time.monotonic()
            obs = list(adapter.discover(region, ERA, store))
            print(f"[{ch}] {tid}: cone ({ra0:.4f},{de0:+.4f}) "
                  f"r={radius:.3f} deg ({len(sub)} events) -> "
                  f"{len(obs)} warps ({_time.monotonic() - t0:.0f}s)",
                  flush=True)
            obs_by_key[(ch, str(tid))] = obs

    # skycell WCS for every warp inside any window of its target
    for (ch, tid), obs in obs_by_key.items():
        tab = events[ch]
        sub = tab[np.asarray([str(x) == tid for x in tab["target_id"]])]
        need = {}
        for ev in sub:
            for r in LADDER[ch]:
                if float(ev["b_min_au"]) >= r:
                    continue
                lo, hi = window(ev, r)
                for o in obs:
                    if lo <= o.t_mid_mjd_utc <= hi:
                        need.setdefault(o.extra["skycell"], []).append(o)
        for sc, cands in sorted(need.items()):
            if sc in adapter._skycells:
                continue
            for o in cands[:10]:   # ~1% of masks 404
                try:
                    adapter.skycell_wcs(o, MSK_DIR)
                    print(f"      wcs {sc} <- mask", flush=True)
                    break
                except FileNotFoundError:
                    continue
            else:
                raise RuntimeError(f"no served mask for skycell {sc} "
                                   f"({ch} {tid})")

    rows = []
    for ch, tab in events.items():
        for ev in tab:
            tid = str(ev["target_id"])
            obs = obs_by_key[(ch, tid)]
            b = float(ev["b_min_au"])
            for r in LADDER[ch]:
                if b >= r:
                    continue
                lo, hi = window(ev, r)
                in_t = [o for o in obs if lo <= o.t_mid_mjd_utc <= hi]
                z_list = Z_GRID_AU if ch == "B" else (None,)
                for z in z_list:
                    n_tot = defaultdict(int)
                    n_pri = defaultdict(int)
                    mjd_by_band = defaultdict(list)
                    mjds = []
                    for o in in_t:
                        if o.extra["skycell"] not in adapter._skycells:
                            continue
                        if ch == "A":
                            ra, de = source_radec(ch, ev)
                        else:
                            ra, de = relay_apparent(ev, z, o.t_mid_mjd_utc)
                        if not adapter.nominal_footprint(o).contains(ra, de):
                            continue
                        n_tot[o.band] += 1
                        if o.quality_flags.get("badflag", 0) == 0:
                            n_pri[o.band] += 1
                        mjd_by_band[o.band].append(o.t_mid_mjd_utc)
                        mjds.append(o.t_mid_mjd_utc)
                    pairs = {}
                    for band in BANDS:
                        ts = np.sort(mjd_by_band[band])
                        pairs[band] = int(np.sum(
                            np.diff(ts) <= PAIR_DT_DAYS)) if len(ts) > 1 else 0
                    row = {
                        "channel": ch, "event_id": str(ev["event_id"]),
                        "target_id": tid, "t_ca_mjd": float(ev["t_ca_mjd"]),
                        "b_min_au": b,
                        "v_perp_km_s": float(ev["v_perp_km_s"]),
                        "radius_au": r, "z_au": z if z else np.nan,
                        "window_days": hi - lo,
                    }
                    for band in BANDS:
                        row[f"n_{band}"] = n_tot[band]
                        row[f"n_{band}_pri"] = n_pri[band]
                        row[f"pair_{band}"] = pairs[band]
                    row.update({
                        "n_total": sum(n_tot.values()),
                        "n_primary": sum(n_pri.values()),
                        "n_pairs": sum(pairs.values()),
                        "mjd_first": min(mjds) if mjds else np.nan,
                        "mjd_last": max(mjds) if mjds else np.nan,
                    })
                    rows.append(row)

    out = Table(rows=rows)
    out.write(OUT / "coverage_v1_events.ecsv", format="ascii.ecsv",
              overwrite=True)

    summary = {}
    for ch in events:
        s = out[out["channel"] == ch]
        summary[ch] = {}
        for r in LADDER[ch]:
            sr = s[np.asarray(s["radius_au"]) == r]
            if ch == "B":   # best z per event
                per_event = {}
                for row in sr:
                    k = row["event_id"]
                    per_event[k] = max(per_event.get(k, 0),
                                       int(row["n_primary"]))
                covered = sum(1 for v in per_event.values() if v > 0)
                n_events = len(per_event)
            else:
                covered = int((np.asarray(sr["n_primary"]) > 0).sum())
                n_events = len(sr)
            summary[ch][str(r)] = {
                "n_events": n_events, "n_covered_primary": covered,
                "median_primary_epochs":
                    float(np.median(np.asarray(sr["n_primary"])))
                    if len(sr) else 0.0,
            }
    (OUT / "coverage_v1_summary.json").write_text(
        json.dumps({"era_mjd": [ERA.start_mjd_utc, ERA.stop_mjd_utc],
                    "z_grid_au": list(Z_GRID_AU), "ladder": LADDER,
                    "pair_dt_days": PAIR_DT_DAYS,
                    "summary": summary}, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
