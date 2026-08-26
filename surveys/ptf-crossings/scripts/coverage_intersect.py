"""PTF crossings coverage intersection v1 (hypotheses.md freeze v1.0).

Port of surveys/ps1-crossings/scripts/coverage_intersect.py onto the
IRSA PTF level-1 adapter. Intersects the universal Earth-center
crossing list with actual PTF epochal coverage: one IBE discovery box
per (channel, target) over the full era, snapshotted under
runs/ptf-crossings; then a local intersection: an exposure covers an
event at ladder radius r if its t_mid lies in t_ca +/- sqrt(r^2-b^2)/
v_perp and its nominal CCD footprint (linear WCS from metadata; PV
distortion absorbed by the pad) contains the predicted source position
(channel A: the per-event star position; channel B: the apparent relay
position for each z on the 5-point grid at that epoch). Coverage
counting only - no pixels are touched and no signal statistic is
formed. Per freeze D5 there is no listing-level quality gate: every
imgtype='object', fid<=2 row counts; dmask and calibration-gate
attrition happens at the search stage, after coverage freezes.

Columns: per-band (g, R) epoch counts, per-band photcalflag=1 counts
(informational - the calibration chain is field-star primary), and
per-band same-night pair counts (epochs with a same-band neighbor
within 0.05 d) - the raw material of the same-night repeat veto.
"""

from __future__ import annotations

import json
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
from sglsurvey.adapters.irsa_ptf import PtfLevel1Adapter
from sglsurvey.snapshots import SnapshotStore

ERA = MjdRange(54891.0, 57051.0)          # measured ptf_procimg span
KM_PER_AU = 1.495978707e8
Z_GRID_AU = (550.0, 1000.0, 2500.0, 5500.0, 10000.0)
LADDER = {"A": (0.1, 1.0), "B": (0.0056, 0.0116, 0.1)}
#: box half-size margin beyond the per-event position spread: covers
#: the in-window B z-track drift (<= ~40 arcsec at z=550) plus
#: TAN-vs-PV slop; the CCD (0.57 x 1.15 deg) dwarfs it.
CONE_MARGIN_DEG = {"A": 0.05, "B": 0.05}
#: nominal-footprint pad: 0 (coverage counts nominal hits; the exact
#: dmask WCS decides usability at the search stage)
PAD_ARCSEC = 0.0
BANDS = ("g", "R")
PAIR_DT_DAYS = 0.05
OUT = REPO / "surveys" / "ptf-crossings" / "results"
RUN = REPO / "runs" / "ptf-crossings"


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
    events = load_events()
    adapter = PtfLevel1Adapter()
    store = SnapshotStore(RUN)

    # one discovery box per (channel, target)
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
            keys = [(o.t_start_mjd_utc, o.native_key["ccdid"],
                     o.native_key["ptffield"]) for o in obs]
            assert len(keys) == len(set(keys)), \
                f"duplicate exposures in {ch} {tid}"
            print(f"[{ch}] {tid}: box ({ra0:.4f},{de0:+.4f}) "
                  f"r={radius:.3f} deg ({len(sub)} events) -> "
                  f"{len(obs)} exposures ({_time.monotonic() - t0:.0f}s)",
                  flush=True)
            obs_by_key[(ch, str(tid))] = obs

    fp_cache = {}

    def footprint(o):
        if o.observation_id not in fp_cache:
            fp_cache[o.observation_id] = adapter.nominal_footprint(
                o, PAD_ARCSEC)
        return fp_cache[o.observation_id]

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
                    n_cal = defaultdict(int)
                    mjd_by_band = defaultdict(list)
                    mjds = []
                    for o in in_t:
                        if ch == "A":
                            ra, de = source_radec(ch, ev)
                        else:
                            ra, de = relay_apparent(ev, z, o.t_mid_mjd_utc)
                        if not footprint(o).contains(ra, de):
                            continue
                        n_tot[o.band] += 1
                        if o.quality_flags.get("photcalflag", 0) == 1:
                            n_cal[o.band] += 1
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
                        row[f"n_{band}_cal"] = n_cal[band]
                        row[f"pair_{band}"] = pairs[band]
                    row.update({
                        "n_total": sum(n_tot.values()),
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
                                       int(row["n_total"]))
                covered = sum(1 for v in per_event.values() if v > 0)
                n_events = len(per_event)
            else:
                covered = int((np.asarray(sr["n_total"]) > 0).sum())
                n_events = len(sr)
            summary[ch][str(r)] = {
                "n_events": n_events, "n_covered": covered,
                "median_epochs":
                    float(np.median(np.asarray(sr["n_total"])))
                    if len(sr) else 0.0,
            }
    (OUT / "coverage_v1_summary.json").write_text(
        json.dumps({"era_mjd": [ERA.start_mjd_utc, ERA.stop_mjd_utc],
                    "z_grid_au": list(Z_GRID_AU), "ladder": LADDER,
                    "pair_dt_days": PAIR_DT_DAYS,
                    "listing_gate": "none (freeze D5): imgtype='object' "
                                    "and fid<=2 only",
                    "summary": summary}, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
