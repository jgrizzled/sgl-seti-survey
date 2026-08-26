"""GALEX crossings coverage intersection v1 (hypotheses.md freeze v1.0).

Intersects the universal Earth-center crossing list with actual GALEX
aspect coverage from the gPhoton database: one era-wide
fGetNearbyAspectEq discovery query per (channel, target) at the median
event position (snapshotted under runs/galex-crossings), then a local
intersection with the flat-chord windows; for every candidate
in-window visit the full aspect rows (boresight ra0/dec0, band, flag)
are pulled and the freeze gates applied exactly: per-second angular
distance from the boresight to the per-event predicted position
(channel A: the event star position; channel B: the apparent relay
position for each z on the 5-point grid at the visit epoch) <= 33
arcmin, aspect flag % 2 == 0 (amendment v1.1: the gPhoton
PhotonPipe convention — only bit 0 marks a failed aspect solution),
per-band live seconds by band substring ('FUV/NUV' = both on).
Coverage counting only — no photon is touched
and no signal statistic is formed. Grazing rungs (B 1.2/2.5 Rsun) are
computed for the coverage-without-statistic ledger; A 1.0 AU is out of
scope (freeze section 2).
"""

from __future__ import annotations

import json
import sys
import time as _time
from pathlib import Path

import numpy as np
from astropy.coordinates import get_body_barycentric
from astropy.table import Table
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from sglsurvey.adapters.mast_gphoton import (
    GphotonClient, galex_ms_to_mjd, group_visits, mjd_to_galex_ms)
from sglsurvey.snapshots import SnapshotStore

#: frozen era = measured aspect-table span (freeze D2)
ERA_UTC = ("2003-06-07T05:02:17.995", "2013-05-01T18:02:22.995")
KM_PER_AU = 1.495978707e8
RSUN_KM = 695_700.0
Z_GRID_AU = (550.0, 1000.0, 2500.0, 5500.0, 10000.0)
LADDER = {"A": (0.1,), "B": (0.0056, 0.0116, 0.1)}
#: discovery radius margin (arcmin) beyond the per-target event-position
#: spread: covers the in-window B z-track offset (<= ~40 arcsec) and
#: records rim-adjacent visits well beyond the 33-arcmin usable gate.
DISCOVERY_MARGIN_ARCMIN = 40.0
USABLE_ARCMIN = 33.0          # freeze section 7
BANDS = ("NUV", "FUV")
OUT = REPO / "surveys" / "galex-crossings" / "results"
RUN = REPO / "runs" / "galex-crossings"


def load_events():
    t = Table.read(REPO / "crossings" / "universal_v1" / "events.ecsv")
    mjd = Time(list(t["t_ca_utc"]), format="isot", scale="utc").mjd
    t["t_ca_mjd"] = mjd
    era = Time(list(ERA_UTC), format="isot", scale="utc").mjd
    inera = (mjd >= era[0]) & (mjd <= era[1]) & (t["validity"] == "valid")
    side = np.asarray(t["axis_distance_au"])
    tA = t[inera & (t["link_direction"] == "inbound") & (side > 0)]
    tB = t[inera & (t["link_direction"] == "outbound") & (side < 0)]
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


def window_mjd(ev, r_au):
    b = float(ev["b_min_au"])
    half_d = (np.sqrt(r_au * r_au - b * b) * KM_PER_AU
              / float(ev["v_perp_km_s"]) / 86400.0)
    return float(ev["t_ca_mjd"]) - half_d, float(ev["t_ca_mjd"]) + half_d


def ang_dist_arcmin(ra1, de1, ra2, de2):
    r1, d1, r2, d2 = map(np.radians, (ra1, de1, ra2, de2))
    c = (np.sin(d1) * np.sin(d2)
         + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2))
    return float(np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))) * 60.0


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    events = load_events()
    client = GphotonClient()
    store = SnapshotStore(RUN)
    era_ms = [mjd_to_galex_ms(Time(u, format="isot", scale="utc").mjd)
              for u in ERA_UTC]

    for ch, tab in events.items():
        print(f"[{ch}] {len(tab)} in-era events, "
              f"{len(set(map(str, tab['target_id'])))} targets", flush=True)

    # one era-wide discovery query per (channel, target)
    disco = {}
    for ch, tab in events.items():
        for tid in sorted(set(map(str, tab["target_id"]))):
            sub = tab[tab["target_id"] == tid]
            ras = np.array([source_radec(ch, e)[0] for e in sub])
            des = np.array([source_radec(ch, e)[1] for e in sub])
            ra0, de0 = float(np.median(ras)), float(np.median(des))
            cosd = max(np.cos(np.radians(de0)), 0.05)
            spread_am = max(np.ptp(ras) * cosd, np.ptp(des)) * 60.0 / 2.0
            radius = spread_am + DISCOVERY_MARGIN_ARCMIN
            t0 = _time.monotonic()
            rows = client.aspect_near(ra0, de0, radius,
                                      era_ms[0], era_ms[1], store)
            times = sorted(set(t for t, _, _ in rows))
            print(f"[{ch}] {tid}: ({ra0:.4f},{de0:+.4f}) r={radius:.1f}' "
                  f"({len(sub)} ev) -> {len(times)} aspect s, "
                  f"{len(group_visits(times))} visits "
                  f"({_time.monotonic() - t0:.0f}s)", flush=True)
            disco[(ch, tid)] = times

    range_cache = {}

    def aspect_exact(t0_ms, t1_ms):
        key = (t0_ms, t1_ms)
        if key not in range_cache:
            range_cache[key] = client.aspect_range(t0_ms, t1_ms + 1000,
                                                   store)
        return range_cache[key]

    rows_out = []
    for ch, tab in events.items():
        for ev in tab:
            tid = str(ev["target_id"])
            times = disco[(ch, tid)]
            b = float(ev["b_min_au"])
            for r_au in LADDER[ch]:
                if b >= r_au:
                    continue
                lo, hi = window_mjd(ev, r_au)
                w0, w1 = mjd_to_galex_ms(lo), mjd_to_galex_ms(hi)
                cand = [t for t in times if w0 <= t < w1]
                visits = group_visits(cand)
                for z in (Z_GRID_AU if ch == "B" else (None,)):
                    live = {bd: 0 for bd in BANDS}
                    n_flagged = 0
                    n_rim = 0
                    dmin = np.inf
                    n_visits_used = 0
                    for v0, v1 in visits:
                        if ch == "A":
                            ra, de = source_radec(ch, ev)
                        else:
                            ra, de = relay_apparent(
                                ev, z, galex_ms_to_mjd((v0 + v1) / 2.0))
                        seen = set()
                        used = False
                        for t_ms, ra0, de0, band, flag in \
                                aspect_exact(v0, v1):
                            if not (w0 <= t_ms < w1) or t_ms in seen:
                                continue
                            seen.add(t_ms)
                            d = ang_dist_arcmin(ra0, de0, ra, de)
                            dmin = min(dmin, d)
                            if d > USABLE_ARCMIN:
                                n_rim += 1
                                continue
                            if flag % 2 != 0:
                                n_flagged += 1
                                continue
                            used = True
                            for bd in BANDS:
                                if bd in band:
                                    live[bd] += 1
                        n_visits_used += int(used)
                    rows_out.append({
                        "channel": ch, "event_id": str(ev["event_id"]),
                        "target_id": tid,
                        "t_ca_mjd": float(ev["t_ca_mjd"]),
                        "t_ca_utc": str(ev["t_ca_utc"]),
                        "b_min_au": b,
                        "b_min_rsun": b * KM_PER_AU / RSUN_KM,
                        "v_perp_km_s": float(ev["v_perp_km_s"]),
                        "radius_au": r_au,
                        "z_au": z if z else np.nan,
                        "window_days": hi - lo,
                        "n_cand_visits": len(visits),
                        "n_visits_used": n_visits_used,
                        "nuv_live_s": live["NUV"],
                        "fuv_live_s": live["FUV"],
                        "n_rim_s": n_rim,
                        "n_flagged_s": n_flagged,
                        "dmin_arcmin": (float(dmin) if np.isfinite(dmin)
                                        else np.nan),
                        "covered": int(live["NUV"] + live["FUV"] > 0),
                    })

    out = Table(rows=rows_out)
    out.write(OUT / "coverage_v1_events.ecsv", format="ascii.ecsv",
              overwrite=True)

    summary = {}
    for ch in events:
        s = out[out["channel"] == ch]
        summary[ch] = {}
        for r_au in LADDER[ch]:
            sr = s[np.asarray(s["radius_au"]) == r_au]
            per_event = {}
            for row in sr:
                k = str(row["event_id"])
                per_event[k] = max(per_event.get(k, 0), int(row["covered"]))
            covered_rows = sr[np.asarray(sr["covered"]) == 1]
            units = sorted(set(
                (str(r["target_id"]), str(r["t_ca_utc"])[:10])
                for r in covered_rows))
            summary[ch][str(r_au)] = {
                "n_events": len(per_event),
                "n_covered": sum(per_event.values()),
                "covered_units": [" ".join(u) for u in units],
            }
    payload = {
        "era_utc": list(ERA_UTC),
        "z_grid_au": list(Z_GRID_AU),
        "ladder": {k: list(v) for k, v in LADDER.items()},
        "usable_arcmin": USABLE_ARCMIN,
        "discovery_margin_arcmin": DISCOVERY_MARGIN_ARCMIN,
        "gates": "boresight <= 33 arcmin, aspect flag % 2 == 0 "
                 "(amendment v1.1), per-band live seconds by band "
                 "substring; freeze v1.0 sections 2/7 + v1.1",
        "summary": summary,
    }
    (OUT / "coverage_v1_summary.json").write_text(
        json.dumps(payload, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
