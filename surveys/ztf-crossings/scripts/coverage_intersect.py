"""ZTF crossings coverage intersection v1 (hypotheses.md freeze v1.0).

Intersects the universal Earth-center crossing list with actual ZTF
public single-epoch coverage. One IBE discovery cone per (channel,
target) over the full era, snapshotted under runs/ztf-crossings; then a
local intersection: an epoch covers an event at ladder radius r if its
t_mid lies in t_ca +/- sqrt(r^2-b^2)/v_perp and its nominal footprint
contains the predicted source position (channel A: the star; channel B:
the apparent relay position for each z on a 5-point grid, computed from
the astropy ephemeris at that epoch). Coverage counting only - no
pixels are touched and no signal statistic is formed.
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.coordinates import get_body_barycentric
from astropy.table import Table
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from sglsurvey.adapters.base import ConeRegion, MjdRange
from sglsurvey.adapters.irsa_ztf import ZtfNominalFootprint, ZtfSciAdapter
from sglsurvey.snapshots import SnapshotStore

ERA = MjdRange(58178.0, 61275.0)          # 2018-03-01 .. 2026-08-23
KM_PER_AU = 1.495978707e8
Z_GRID_AU = (550.0, 1000.0, 2500.0, 5500.0, 10000.0)
LADDER = {"A": (0.1, 1.0), "B": (0.0056, 0.0116, 0.1)}
CONE_MARGIN_DEG = {"A": 0.05, "B": 0.15}
OUT = REPO / "surveys" / "ztf-crossings" / "results"
RUN = REPO / "runs" / "ztf-crossings"


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


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    events = load_events()
    adapter = ZtfSciAdapter()
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
            print(f"[{ch}] {tid}: cone ({ra0:.4f},{de0:.4f}) r={radius:.3f} deg "
                  f"({len(sub)} events)", flush=True)
            obs = list(adapter.discover(region, ERA, store))
            print(f"      {len(obs)} quadrant-exposures", flush=True)
            obs_by_key[(ch, str(tid))] = obs

    rows = []
    for ch, tab in events.items():
        for ev in tab:
            tid = str(ev["target_id"])
            obs = obs_by_key[(ch, tid)]
            b = float(ev["b_min_au"])
            vperp = float(ev["v_perp_km_s"])
            for r in LADDER[ch]:
                if b >= r:
                    continue
                half_d = np.sqrt(r * r - b * b) * KM_PER_AU / vperp / 86400.0
                lo, hi = ev["t_ca_mjd"] - half_d, ev["t_ca_mjd"] + half_d
                in_t = [o for o in obs
                        if lo <= o.t_mid_mjd_utc <= hi
                        and o.t_mid_mjd_utc >= ERA.start_mjd_utc]
                z_list = Z_GRID_AU if ch == "B" else (None,)
                for z in z_list:
                    n_tot = defaultdict(int)
                    n_pri = defaultdict(int)
                    mjds = []
                    for o in in_t:
                        if ch == "A":
                            ra, de = source_radec(ch, ev)
                        else:
                            ra, de = relay_apparent(ev, z, o.t_mid_mjd_utc)
                        if not ZtfNominalFootprint(o).contains(ra, de):
                            continue
                        q = o.quality_flags
                        n_tot[o.band] += 1
                        if (not q.get("bad_quality")
                                and (q.get("seeing") or 99) <= 4.0):
                            n_pri[o.band] += 1
                        mjds.append(o.t_mid_mjd_utc)
                    rows.append({
                        "channel": ch, "event_id": str(ev["event_id"]),
                        "target_id": tid, "t_ca_mjd": float(ev["t_ca_mjd"]),
                        "b_min_au": b, "v_perp_km_s": vperp,
                        "radius_au": r, "z_au": z if z else np.nan,
                        "window_days": 2 * half_d,
                        "n_g": n_tot["zg"], "n_r": n_tot["zr"],
                        "n_i": n_tot["zi"],
                        "n_g_pri": n_pri["zg"], "n_r_pri": n_pri["zr"],
                        "n_i_pri": n_pri["zi"],
                        "n_total": sum(n_tot.values()),
                        "n_primary": sum(n_pri.values()),
                        "mjd_first": min(mjds) if mjds else np.nan,
                        "mjd_last": max(mjds) if mjds else np.nan,
                    })

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
                    float(np.median(np.asarray(sr["n_primary"]))) if len(sr) else 0.0,
            }
    (OUT / "coverage_v1_summary.json").write_text(
        json.dumps({"era_mjd": [ERA.start_mjd_utc, ERA.stop_mjd_utc],
                    "z_grid_au": list(Z_GRID_AU), "ladder": LADDER,
                    "summary": summary}, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
