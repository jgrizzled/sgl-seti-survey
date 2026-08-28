"""Rubin DP2 crossings coverage intersection v1 (hypotheses.md v1.0).

Intersects the universal Earth-center crossing list with DP2 visit
coverage, applying the freeze gates exactly: era re-measured from
dp2.Visit (D2), validity == 'valid', side-of-axis cuts (A: inbound,
axis > 0; B: outbound, axis < 0), flat-chord rung windows, and the
section-7 detector-footprint gate — every predicted position (channel
A: star; channel B: apparent relay per z on the frozen 5-point grid at
the visit epoch) is tested point-in-quadrilateral against the
dp2.VisitDetector corner coordinates with non-null magLim. Coverage
counting only — no DiaSource row is touched and no signal statistic is
formed. Grazing rungs (B 1.2/2.5 Rsun) are computed for the
coverage-without-statistic ledger; A 1.0 AU is geometry-only counts
(freeze section 2). All TAP queries snapshotted under
runs/rubin-crossings/coverage_v1/.
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

from sglsurvey.adapters.rsp_tap import RspTapClient
from sglsurvey.snapshots import SnapshotStore

KM_PER_AU = 1.495978707e8
RSUN_KM = 695_700.0
Z_GRID_AU = (550.0, 1000.0, 2500.0, 5500.0, 10000.0)
LADDER = {"A": (0.1,), "B": (0.0056, 0.0116, 0.1)}
DISCOVERY_DEG = 2.0           # boresight discovery cone; footprint decides
OUT = REPO / "surveys" / "rubin-crossings" / "results"
RUN = REPO / "runs" / "rubin-crossings" / "coverage_v1"


def relay_apparent(ev, z_au, t_mjd):
    """Apparent ICRS ra/dec of a relay at z_au on the anti-star axis
    (the GALEX coverage construction, freeze section 5)."""
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
    t_ca = float(ev["t_ca_tdb_jd"]) - 2400000.5
    return t_ca - half_d, t_ca + half_d


def sep_deg(ra1, de1, ra2, de2):
    r1, d1, r2, d2 = map(np.radians, (ra1, de1, ra2, de2))
    c = (np.sin(d1) * np.sin(d2)
         + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2))
    return float(np.degrees(np.arccos(np.clip(c, -1.0, 1.0))))


def point_in_quad(ra, dec, corners):
    """Winding test in a local tangent frame centred on the point.

    corners: [(ra, dec) x 4] in detector-corner order (llc, lrc, urc,
    ulc). Detector scale ~0.25 deg: planar approximation with
    cos(dec) scaling is exact to << 1 arcsec here.
    """
    cosd = np.cos(np.radians(dec))
    dra = [((c[0] - ra + 180.0) % 360.0 - 180.0) * cosd for c in corners]
    dde = [c[1] - dec for c in corners]
    inside = False
    n = len(corners)
    for i in range(n):
        x1, y1 = dra[i], dde[i]
        x2, y2 = dra[(i + 1) % n], dde[(i + 1) % n]
        if (y1 > 0) != (y2 > 0):
            x_int = x1 + (0 - y1) * (x2 - x1) / (y2 - y1)
            if x_int > 0:
                inside = not inside
    return inside


def load_events(era):
    t = Table.read(REPO / "crossings" / "universal_v1" / "events.ecsv")
    t["t_ca_mjd"] = np.array(t["t_ca_tdb_jd"]) - 2400000.5
    valid = np.array(t["validity"]) == "valid"
    side = np.asarray(t["axis_distance_au"])
    tA = t[valid & (t["link_direction"] == "inbound") & (side > 0)]
    tB = t[valid & (t["link_direction"] == "outbound") & (side < 0)]
    out = {}
    for ch, tab, wide in (("A", tA, False), ("B", tB, False),
                          ("A_1AU", tA, True)):
        lim = 1.0 if wide else max(LADDER[ch[0]])
        keep = []
        for ev in tab:
            if float(ev["b_min_au"]) >= lim:
                continue
            lo, hi = window_mjd(ev, lim)
            if hi >= era[0] and lo <= era[1]:
                keep.append(ev)
        out[ch] = keep
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    client = RspTapClient()
    store = SnapshotStore(RUN)

    # D2: re-measure the era
    cols, rows = client.query(
        "SELECT COUNT(*) AS n, MIN(expMidptMJD) AS lo, "
        "MAX(expMidptMJD) AS hi FROM dp2.Visit", store)
    n_visits, era = int(rows[0][0]), (float(rows[0][1]), float(rows[0][2]))
    print(f"era re-measured: {n_visits} visits, MJD {era[0]:.3f}-{era[1]:.3f}")

    events = load_events(era)
    for ch in ("A", "B"):
        print(f"[{ch}] {len(events[ch])} in-era events, "
              f"{len(set(str(e['target_id']) for e in events[ch]))} targets")

    rows_out = []
    vd_cache = {}

    def visit_detectors(visit_id):
        if visit_id not in vd_cache:
            cols, rows = client.query(
                "SELECT detector, llcra, llcdec, lrcra, lrcdec, urcra, "
                "urcdec, ulcra, ulcdec, magLim, seeing FROM "
                f"dp2.VisitDetector WHERE visitId = {visit_id}", store)
            vd_cache[visit_id] = [
                {"detector": int(r[0]),
                 "corners": [(float(r[1]), float(r[2])),
                             (float(r[3]), float(r[4])),
                             (float(r[5]), float(r[6])),
                             (float(r[7]), float(r[8]))],
                 "magLim": (float(r[9]) if r[9] not in ("", None) else None),
                 "seeing": (float(r[10]) if r[10] not in ("", None)
                            else None)}
                for r in rows
                # rows with masked corner coordinates carry no usable
                # footprint (failed per-detector processing) — excluded
                # by the section-7 gate by construction
                if all(x not in ("", None) for x in r[1:9])]
        return vd_cache[visit_id]

    for ch in ("A", "B"):
        for ev in events[ch]:
            tid = str(ev["target_id"])
            b = float(ev["b_min_au"])
            t_ca = float(ev["t_ca_mjd"])
            if ch == "A":
                pos_ra = float(ev["star_icrs_ra_deg"])
                pos_de = float(ev["star_icrs_dec_deg"])
            else:
                pos_ra = float(ev["relay_icrs_ra_deg"])
                pos_de = float(ev["relay_icrs_dec_deg"])
            for r_au in LADDER[ch]:
                if b >= r_au:
                    continue
                lo, hi = window_mjd(ev, r_au)
                cols, vrows = client.query(
                    "SELECT visit, ra, dec, band, expMidptMJD, expTime "
                    "FROM dp2.Visit WHERE expMidptMJD BETWEEN "
                    f"{lo:.6f} AND {hi:.6f} AND CONTAINS(POINT('ICRS', "
                    f"ra, dec), CIRCLE('ICRS', {pos_ra:.6f}, {pos_de:.6f}, "
                    f"{DISCOVERY_DEG})) = 1", store)
                for vr in vrows:
                    visit_id = int(vr[0])
                    v_mjd = float(vr[4])
                    dt = v_mjd - t_ca
                    b_at = np.hypot(b, float(ev["v_perp_km_s"]) * abs(dt)
                                    * 86400.0 / KM_PER_AU)
                    dets = visit_detectors(visit_id)
                    for z in (Z_GRID_AU if ch == "B" else (None,)):
                        if ch == "A":
                            ra, de = pos_ra, pos_de
                        else:
                            ra, de = relay_apparent(ev, z, v_mjd)
                        hit = None
                        for d in dets:
                            if point_in_quad(ra, de, d["corners"]):
                                hit = d
                                break
                        rows_out.append({
                            "channel": ch, "event_id": str(ev["event_id"]),
                            "target_id": tid, "t_ca_mjd": t_ca,
                            "b_min_au": b,
                            "b_min_rsun": b * KM_PER_AU / RSUN_KM,
                            "v_perp_km_s": float(ev["v_perp_km_s"]),
                            "radius_au": r_au,
                            "z_au": z if z else np.nan,
                            "window_days": hi - lo,
                            "visit": visit_id, "band": str(vr[3]),
                            "visit_mjd": v_mjd, "dt_days": dt,
                            "b_at_visit_rsun": b_at * KM_PER_AU / RSUN_KM,
                            "pred_ra": ra, "pred_dec": de,
                            "on_detector": int(hit is not None),
                            "detector": (hit["detector"] if hit else -1),
                            "maglim": (hit["magLim"] if hit and
                                       hit["magLim"] is not None
                                       else np.nan),
                            "seeing": (hit["seeing"] if hit and
                                       hit["seeing"] is not None
                                       else np.nan),
                        })
                        if hit:
                            print(f"  ON-DET [{ch} {r_au}] {tid} "
                                  f"visit {visit_id} {vr[3]} z={z} "
                                  f"det {hit['detector']} "
                                  f"magLim {hit['magLim']}")

    tab = Table(rows=rows_out) if rows_out else Table()
    tab.write(OUT / "coverage_v1_events.ecsv", format="ascii.ecsv",
              overwrite=True)

    # summary: unit = (target, channel-rung); searchable if any
    # (visit, z) row passes the footprint gate with magLim
    summary = {"era_mjd": list(era), "n_visits_total": n_visits,
               "z_grid_au": list(Z_GRID_AU),
               "ladder": {k: list(v) for k, v in LADDER.items()},
               "discovery_deg": DISCOVERY_DEG,
               "a_1au_geometry_only": {
                   "in_era_events": len(events["A_1AU"]),
                   "targets": len(set(str(e["target_id"])
                                      for e in events["A_1AU"]))},
               "channels": {}}
    for ch in ("A", "B"):
        summary["channels"][ch] = {}
        for r_au in LADDER[ch]:
            sel = [r for r in rows_out
                   if r["channel"] == ch and r["radius_au"] == r_au]
            events_here = sorted(set(r["event_id"] for r in sel))
            n_era = sum(1 for e in events[ch]
                        if float(e["b_min_au"]) < r_au)
            units = {}
            for r in sel:
                if r["on_detector"] and not np.isnan(r["maglim"]):
                    k = f"{r['target_id']} {r['band']}"
                    units.setdefault(k, []).append(
                        {"visit": r["visit"], "z_au": r["z_au"],
                         "detector": r["detector"], "maglim": r["maglim"],
                         "dt_days": round(r["dt_days"], 3),
                         "b_at_visit_rsun": round(r["b_at_visit_rsun"],
                                                  3)})
            summary["channels"][ch][str(r_au)] = {
                "in_era_events": n_era,
                "events_with_candidate_visits": len(events_here),
                "searchable_units": units,
            }
    (OUT / "coverage_v1_summary.json").write_text(
        json.dumps(summary, indent=2, default=float) + "\n")
    print(json.dumps(summary["channels"], indent=2, default=float))
    print(f"a_1au geometry-only: {summary['a_1au_geometry_only']}")


if __name__ == "__main__":
    main()
