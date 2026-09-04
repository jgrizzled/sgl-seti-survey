"""Radio crossings coverage intersection, extension v2 (plan §5.8 item 7).

Repeats the VLASS geometry-only construction (report/radio_crossings.md,
decision unchanged: no radio signal search in this repo) on two archives
that revisit fields with per-observation dates:

* CASDA (CSIRO ASKAP Science Data Archive) ObsCore via TAP: every
  Stokes-I continuum restored image (``dataproduct_subtype =
  'cont.restored.t0'``, ``filename LIKE 'image.i.%'``) of any ASKAP
  collection whose ``s_region`` contains the target-channel position —
  RACS (low/mid/high and the low re-epochs), VAST (pilot 2019–2021 and
  the full survey 2022–), EMU and any other project. Each row carries
  the scheduling block (SBID) and its start/stop MJD; SBIDs whose
  cubes lack ``t_min`` (the VAST pilot) are dated from
  ``casda.observation``. One observation per (SBID, field).
* LoTSS DR3 (ASTRON VO TAP, ``lotss_dr3.pointings``): 2,551 HBA
  pointings, each with the mid-MJD of every 8-h observation that fed it
  (``dateallobs``).

Geometry (declared):
  * position = per-event star (channel A) or antipode/relay (channel B)
    ICRS position from crossings/universal_v1; the TAP cone is issued at
    the per-target-channel median position.
  * ASKAP: ``in_image`` = the archive's own s_region containment;
    ``in_footprint`` = separation from the field centre <= 2.25 deg +
    HPBW/2, HPBW = 1.09 lambda / 12 m (closepack36, 0.9 deg pitch), i.e.
    ~3.1 deg at 888 MHz, ~2.8 deg at 1.37 GHz, ~2.7 deg at 1.66 GHz.
  * LoTSS: ``in_image`` = separation <= 2.61 deg (30 % of a 3.96 deg
    FWHM HBA primary beam, the DR mosaic trim); ``in_beam`` = <= 1.98 deg
    (half power).
Time (declared): an observation covers a window if [t_start, t_stop]
  overlaps [t_ca - hd, t_ca + hd] with hd from the flat-chord rung
  radius (no scheduling tolerance — these are actual on-sky intervals);
  LoTSS observations are taken as 8 h centred on the listed mid-MJD.
Era: events with t_ca in 2014-01-01 .. 2028-01-01 (LoTSS from 2014,
  ASKAP from 2019); per-arm archive date ranges are recorded so windows
  beyond the last observation can be separated from uncovered ones.
Rungs: the frozen programme ladder (A 0.1 / 1.0 AU; B 1.2 / 2.5 Rsun,
  0.1 AU). Metadata only — no pixel is touched, no statistic formed.
Raw TAP responses are snapshotted under runs/radio-crossings/v2/.
"""

from __future__ import annotations

import json
import re
import sys
import time as _time
from pathlib import Path

import numpy as np
import requests
from astropy.table import Table
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from sglsurvey.snapshots import SnapshotStore  # noqa: E402

CASDA_TAP = "https://casda.csiro.au/casda_vo_tools/tap/sync"
ASTRON_TAP = "https://vo.astron.nl/__system__/tap/run/tap/sync"
KM_PER_AU = 1.495978707e8
ERA = (56658.0, 61771.0)             # 2014-01-01 .. 2028-01-01
LADDER = {"A": (0.1, 1.0), "B": (0.0056, 0.0116, 0.1)}
ASKAP_DISH_M = 12.0
ASKAP_FOOTPRINT_CORE_DEG = 2.25      # closepack36, 0.9 deg pitch
LOTSS_FWHM_DEG = 3.96
LOTSS_IN_BEAM_DEG = LOTSS_FWHM_DEG / 2.0                         # 1.98
LOTSS_IN_IMAGE_DEG = LOTSS_IN_BEAM_DEG * np.sqrt(np.log(1 / 0.3)
                                                 / np.log(2))    # 2.61
LOTSS_OBS_HALF_DAYS = 4.0 / 24.0
LOTSS_DEC_MIN = -10.0
QUALITY_RANK = {"GOOD": 0, "UNCERTAIN": 1, "NOT_VALIDATED": 2,
                "REJECTED": 3, "BAD": 4}
VALIDATED = ("GOOD", "UNCERTAIN", "NOT_VALIDATED")
OUT = REPO / "surveys" / "radio-crossings" / "results"
RUN = REPO / "runs" / "radio-crossings" / "v2"


# ----------------------------------------------------------------- events
def load_events():
    t = Table.read(REPO / "crossings" / "universal_v1" / "events.ecsv")
    t["mjd"] = Time(list(t["t_ca_utc"]), format="isot", scale="utc").mjd
    side = np.asarray(t["axis_distance_au"])
    era = (t["mjd"] >= ERA[0]) & (t["mjd"] <= ERA[1])
    A = t[(t["link_direction"] == "inbound") & (side > 0) & era]
    B = t[(t["link_direction"] == "outbound") & (side < 0) & era]
    out = {}
    for ch, tab in (("A", A), ("B", B)):
        out[ch] = tab[np.asarray(tab["b_min_au"]) < max(LADDER[ch])]
    return out


def pos_of(ch, ev):
    if ch == "A":
        return float(ev["star_icrs_ra_deg"]), float(ev["star_icrs_dec_deg"])
    return float(ev["relay_icrs_ra_deg"]), float(ev["relay_icrs_dec_deg"])


def half_days(b, r, vp):
    return float(np.sqrt(r * r - b * b) * KM_PER_AU / vp / 86400.0)


def sep_deg(ra1, de1, ra2, de2):
    c1, c2 = np.radians(de1), np.radians(de2)
    dra = np.radians(ra1 - ra2)
    x = (np.sin(c1) * np.sin(c2)
         + np.cos(c1) * np.cos(c2) * np.cos(dra))
    return float(np.degrees(np.arccos(np.clip(x, -1, 1))))


def windows_of(ch, ev):
    b, vp = float(ev["b_min_au"]), float(ev["v_perp_km_s"])
    for r in LADDER[ch]:
        if b >= r:
            continue
        hd = half_days(b, r, vp)
        yield r, hd, float(ev["mjd"]) - hd, float(ev["mjd"]) + hd


def target_channels(events):
    for ch, tab in events.items():
        for tid in sorted(set(str(x) for x in tab["target_id"])):
            sub = tab[np.asarray([str(x) == tid for x in tab["target_id"]])]
            ra0 = float(np.median([pos_of(ch, e)[0] for e in sub]))
            de0 = float(np.median([pos_of(ch, e)[1] for e in sub]))
            yield ch, tid, sub, ra0, de0


# -------------------------------------------------------------------- TAP
def tap_csv(session, store, url, query, maxrec=200000, label=""):
    params = {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv",
              "QUERY": query, "MAXREC": str(maxrec)}
    for attempt in range(4):
        try:
            r = session.get(url, params=params, timeout=600)
            r.raise_for_status()
            text = r.text
            break
        except Exception as exc:
            print(f"  [{label}] attempt {attempt + 1}: {exc}", flush=True)
            _time.sleep(15 * (attempt + 1))
    else:
        return None
    tab = Table.read(text, format="ascii.csv") if text.count("\n") > 1 \
        else Table(names=text.strip().split(","))
    store.store(service_url=url, query=query, request_utc=Time.now().isot,
                response_bytes=r.content, row_count=len(tab),
                http_status=r.status_code)
    return tab


def sval(row, key):
    v = row[key]
    if v is None or (hasattr(v, "mask") and np.ma.is_masked(v)):
        return None
    try:
        if isinstance(v, (float, np.floating)) and np.isnan(v):
            return None
    except Exception:
        pass
    return v


# ------------------------------------------------------------- ASKAP arm
def askap_arm(events, session, store):
    cols = ("obs_collection, obs_id, target_name, filename, s_ra, s_dec, "
            "s_fov, t_min, t_max, t_exptime, em_min, em_max, quality_level, "
            "calib_level")
    rows, per_tc, sbid_dates = [], [], {}
    for ch, tid, sub, ra0, de0 in target_channels(events):
        q = (f"SELECT {cols} FROM ivoa.obscore WHERE dataproduct_type='cube'"
             f" AND dataproduct_subtype='cont.restored.t0'"
             f" AND filename LIKE 'image.i.%'"
             f" AND 1=CONTAINS(POINT('ICRS', {ra0:.5f}, {de0:.5f}), s_region)")
        tab = tap_csv(session, store, CASDA_TAP, q, label=f"casda {ch}/{tid}")
        if tab is None:
            per_tc.append({"channel": ch, "target_id": tid,
                           "status": "query_failed"})
            continue
        # one observation per (SBID, field centre); several cubes per
        # SBID exist (pipeline v1/v2 re-runs, lowres/raw variants) —
        # keep the best-validated one
        obs = {}
        for r in tab:
            sbid = str(r["obs_id"])
            m = re.search(r"([A-Za-z]+_[0-9]{4}[+-][0-9]{2}[A-Z]?|[A-Za-z]+_\d+)",
                          str(r["filename"]))
            field = sval(r, "target_name") or (m.group(1) if m else "?")
            key = (sbid, round(float(r["s_ra"]), 2), round(float(r["s_dec"]), 2))
            qual = str(sval(r, "quality_level") or "")
            if key in obs and QUALITY_RANK.get(qual, 9) >= QUALITY_RANK.get(
                    obs[key]["quality"], 9):
                continue
            obs[key] = {
                "collection": str(r["obs_collection"]), "sbid": sbid,
                "field": str(field), "field_ra": float(r["s_ra"]),
                "field_dec": float(r["s_dec"]),
                "t_min": sval(r, "t_min"), "t_max": sval(r, "t_max"),
                "em_min_m": sval(r, "em_min"),
                "quality": str(sval(r, "quality_level") or ""),
                "filename": str(r["filename"])}
        # date SBIDs whose cubes carry no t_min (VAST pilot)
        missing = sorted({o["sbid"] for o in obs.values()
                          if o["t_min"] is None and o["sbid"].startswith("ASKAP-")})
        for sb in missing:
            n = sb.split("-")[1]
            if n not in sbid_dates:
                t2 = tap_csv(session, store, CASDA_TAP,
                             f"SELECT sbid, obs_start_mjd, obs_end_mjd, obs_program "
                             f"FROM casda.observation WHERE sbid={n}",
                             label=f"casda obs {n}")
                sbid_dates[n] = (None if t2 is None or len(t2) == 0 else
                                 (float(t2[0]["obs_start_mjd"]),
                                  float(t2[0]["obs_end_mjd"])))
            if sbid_dates[n] is not None:
                for o in obs.values():
                    if o["sbid"] == sb and o["t_min"] is None:
                        o["t_min"], o["t_max"] = sbid_dates[n]
                        o["dated_from"] = "casda.observation"
        obs_list = list(obs.values())
        n_dated = sum(1 for o in obs_list if o["t_min"] is not None)
        n_in = 0
        for ev in sub:
            ra_e, de_e = pos_of(ch, ev)
            for r_au, hd, lo, hi in windows_of(ch, ev):
                for o in obs_list:
                    if o["t_min"] is None:
                        continue
                    t0, t1 = float(o["t_min"]), float(o["t_max"])
                    if t1 < lo or t0 > hi:
                        continue
                    s = sep_deg(o["field_ra"], o["field_dec"], ra_e, de_e)
                    lam = o["em_min_m"]
                    freq = 299.792458 / float(lam) if lam else None   # MHz
                    hpbw = np.degrees(1.09 * float(lam) / ASKAP_DISH_M) \
                        if lam else 1.76
                    foot = ASKAP_FOOTPRINT_CORE_DEG + hpbw / 2
                    rows.append({
                        "channel": ch, "target_id": tid,
                        "event_id": str(ev["event_id"]), "radius_au": r_au,
                        "b_min_au": float(ev["b_min_au"]),
                        "t_ca_mjd": float(ev["mjd"]), "window_days": 2 * hd,
                        "collection": o["collection"], "sbid": o["sbid"],
                        "field": o["field"], "obs_start_mjd": t0,
                        "obs_end_mjd": t1,
                        "offset_days": round(0.5 * (t0 + t1) - float(ev["mjd"]), 3),
                        "freq_mhz": round(freq, 1) if freq else 0.0,
                        "sep_deg": round(s, 3),
                        "footprint_half_deg": round(foot, 2),
                        "in_footprint": bool(s <= foot),
                        "quality": o["quality"],
                        "dated_from": o.get("dated_from", "obscore")})
                    n_in += 1
        colls = {}
        for o in obs_list:
            colls[o["collection"]] = colls.get(o["collection"], 0) + 1
        per_tc.append({"channel": ch, "target_id": tid, "status": "ok",
                       "ra": ra0, "dec": de0, "n_events": len(sub),
                       "n_cubes": len(tab), "n_obs": len(obs_list),
                       "n_obs_dated": n_dated, "collections": colls,
                       "n_inwindow_rows": n_in})
        print(f"[askap {ch}] {tid}: {len(tab)} cubes -> {len(obs_list)} obs "
              f"({n_dated} dated), {n_in} in-window rows", flush=True)
    if rows:
        Table(rows=rows).write(OUT / "askap_inwindow_v2.ecsv",
                               format="ascii.ecsv", overwrite=True)
    (OUT / "askap_targets_v2.json").write_text(
        json.dumps(per_tc, indent=1) + "\n")
    return rows, per_tc


# ------------------------------------------------------------- LoTSS arm
def parse_mjd_list(s):
    return [float(x) for x in re.findall(r"[-+]?\d+(?:\.\d+)?", str(s))]


def lotss_arm(events, session, store):
    tab = tap_csv(session, store, ASTRON_TAP,
                  "SELECT pointing_id, centeralpha, centerdelta, dateobs, "
                  "datefirstobs, datelastobs, dateallobs, lofar_obsids "
                  "FROM lotss_dr3.pointings", maxrec=10000, label="lotss")
    if tab is None:
        raise SystemExit("LoTSS pointing table unavailable")
    pts = []
    for r in tab:
        mjds = parse_mjd_list(r["dateallobs"]) or [float(r["dateobs"])]
        pts.append((str(r["pointing_id"]), float(r["centeralpha"]),
                    float(r["centerdelta"]), mjds))
    all_mjd = [m for p in pts for m in p[3]]
    print(f"LoTSS DR3 pointings: {len(pts)}, observations: {len(all_mjd)}, "
          f"MJD {min(all_mjd):.1f}..{max(all_mjd):.1f}")
    rows, per_tc = [], []
    for ch, tid, sub, ra0, de0 in target_channels(events):
        if de0 < LOTSS_DEC_MIN:
            continue
        near = [p for p in pts
                if sep_deg(p[1], p[2], ra0, de0) <= LOTSS_IN_IMAGE_DEG + 0.5]
        n_in = 0
        for ev in sub:
            ra_e, de_e = pos_of(ch, ev)
            for r_au, hd, lo, hi in windows_of(ch, ev):
                for pid, pra, pde, mjds in near:
                    s = sep_deg(pra, pde, ra_e, de_e)
                    if s > LOTSS_IN_IMAGE_DEG:
                        continue
                    for m in mjds:
                        t0, t1 = m - LOTSS_OBS_HALF_DAYS, m + LOTSS_OBS_HALF_DAYS
                        if t1 < lo or t0 > hi:
                            continue
                        rows.append({
                            "channel": ch, "target_id": tid,
                            "event_id": str(ev["event_id"]),
                            "radius_au": r_au, "b_min_au": float(ev["b_min_au"]),
                            "t_ca_mjd": float(ev["mjd"]),
                            "window_days": 2 * hd, "pointing": pid,
                            "obs_mid_mjd": m,
                            "offset_days": round(m - float(ev["mjd"]), 3),
                            "sep_deg": round(s, 3),
                            "in_beam": bool(s <= LOTSS_IN_BEAM_DEG)})
                        n_in += 1
        n_obs = sum(len(p[3]) for p in near
                    if sep_deg(p[1], p[2], ra0, de0) <= LOTSS_IN_IMAGE_DEG)
        per_tc.append({"channel": ch, "target_id": tid, "ra": ra0, "dec": de0,
                       "n_events": len(sub),
                       "n_pointings_in_image": sum(
                           1 for p in near
                           if sep_deg(p[1], p[2], ra0, de0) <= LOTSS_IN_IMAGE_DEG),
                       "n_obs_in_image": n_obs, "n_inwindow_rows": n_in})
        print(f"[lotss {ch}] {tid}: {n_obs} obs over the position, "
              f"{n_in} in-window rows", flush=True)
    if rows:
        Table(rows=rows).write(OUT / "lotss_inwindow_v2.ecsv",
                               format="ascii.ecsv", overwrite=True)
    (OUT / "lotss_targets_v2.json").write_text(
        json.dumps(per_tc, indent=1) + "\n")
    return rows, per_tc, (min(all_mjd), max(all_mjd))


# ------------------------------------------------------------------ main
def main():
    OUT.mkdir(parents=True, exist_ok=True)
    RUN.mkdir(parents=True, exist_ok=True)
    events = load_events()
    session = requests.Session()
    store = SnapshotStore(RUN)
    askap_rows, askap_tc = ([], [])
    if "--skip-askap" not in sys.argv:
        askap_rows, askap_tc = askap_arm(events, session, store)
    lotss_rows, lotss_tc, lotss_range = ([], [], None)
    if "--skip-lotss" not in sys.argv:
        lotss_rows, lotss_tc, lotss_range = lotss_arm(events, session, store)

    def cnt(rows, ch, r, flag=None, coll=None, validated=False):
        sel = [x for x in rows if x["channel"] == ch and x["radius_au"] == r
               and (flag is None or x.get(flag, True))
               and (coll is None or coll in x.get("collection", ""))
               and (not validated or x.get("quality") in VALIDATED)]
        return {"rows": len(sel),
                "events": len({x["event_id"] for x in sel}),
                "targets": len({x["target_id"] for x in sel})}

    summary = {"era_mjd": ERA, "askap": {}, "lotss": {},
               "lotss_obs_mjd_range": lotss_range,
               "n_events": {ch: len(t) for ch, t in events.items()}}
    for ch in LADDER:
        for r in LADDER[ch]:
            summary["askap"][f"{ch}_{r}"] = {
                "in_image": cnt(askap_rows, ch, r),
                "in_footprint": cnt(askap_rows, ch, r, "in_footprint"),
                "in_footprint_validated": cnt(askap_rows, ch, r,
                                              "in_footprint", validated=True),
                "in_footprint_racs": cnt(askap_rows, ch, r, "in_footprint",
                                         "Rapid"),
                "in_footprint_vast": cnt(askap_rows, ch, r, "in_footprint",
                                         "VAST")}
            summary["lotss"][f"{ch}_{r}"] = {
                "in_image": cnt(lotss_rows, ch, r),
                "in_beam": cnt(lotss_rows, ch, r, "in_beam")}
    (OUT / "coverage_v2_summary.json").write_text(
        json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
