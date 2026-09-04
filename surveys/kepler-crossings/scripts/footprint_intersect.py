"""Kepler prime + K2 campaign footprint intersect (plan §5.8 item 8; 2026-09-04).

Geometry-only: does any beam-crossing channel position (A = on-star,
B = antipode/relay) fall on Kepler science silicon during a Kepler
quarter (2009-05 .. 2013-05) or a K2 campaign (2014-03 .. 2018-09)?

Observer: the Kepler spacecraft (``crossings/kepler_v1``, JPL Horizons
-227). Kepler flew an Earth-trailing heliocentric orbit and was
0.05-1.0 AU from Earth over its science life, so the Earth-center
universal list the plan row names is invalid for it at every rung;
the Earth-center list is loaded only for the matched-event shift
diagnostic (like the TESS validation).

Footprints: ``K2fov`` (Kepler/K2 GO office, silicon-level channel
polygons incl. the module-3/7 and module-4 failures) for campaigns
0-19 and the prime field (campaign 1000, roll 20 + 90 x season,
envelope over the four seasons). Campaign date ranges are the *actual*
data ranges from MAST CAOM (``dbo.obspointing`` grouped by
``sequence_number``), snapshotted under runs/kepler-crossings/recon.

Archive check (independent of K2fov): for every target-channel unit
with an in-era event at <= 0.1 AU, a MAST box count of Kepler/K2
products within +-8 deg of the unit position, grouped by campaign.

No signal statistic is formed; nothing here is a search.
"""
from __future__ import annotations

import json
import re
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import requests
from astropy.coordinates import get_body_barycentric
from astropy.table import Table
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
from sglsurvey.snapshots import SnapshotStore  # noqa: E402

RUN = REPO / "runs" / "kepler-crossings" / "recon"
OUT = Path(__file__).resolve().parents[1] / "results" / "footprint_intersect_v1.json"
EVENTS = REPO / "crossings" / "kepler_v1" / "events.ecsv"
EVENTS_EARTH = REPO / "crossings" / "universal_v1" / "events.ecsv"
SC_TABLE = REPO / "crossings" / "observers" / "kepler_sc_ephemeris.npz"
MAST_TAP = "https://mast.stsci.edu/vo-tap/api/v0.1/caom/async"   # sync hits a 60 s gateway timeout

KM_PER_AU = 1.495978707e8
RSUN_AU = 0.00465047
RUNGS = [("B", "1.2Rsun", 1.2 * RSUN_AU), ("B", "2.5Rsun", 2.5 * RSUN_AU),
         ("B", "0.1AU", 0.1), ("A", "0.1AU", 0.1)]
PAD_PIX = 12              # K2fov default silicon-edge padding
DATE_PAD_D = 2.0          # slack around the MAST campaign data ranges
PLAN_CLIP_D = 15.0        # MAST ranges are clipped to the K2fov planned dates +- this
                          # (a C12-era product is tagged sequence 14 in CAOM: its raw
                          # range would start 168 d early)
BOX_HALF_DEG = 8.0        # Kepler FOV half-diagonal ~7.6 deg
PRIME_RA, PRIME_DEC = 290.66666667, 44.5


def _mast_async(adql, timeout):
    """UWS job: submit, poll, fetch the VOTable result."""
    r = requests.post(MAST_TAP, data={"REQUEST": "doQuery", "LANG": "ADQL", "PHASE": "RUN",
                                      "QUERY": adql}, timeout=60, allow_redirects=False)
    r.raise_for_status()
    job = r.headers["location"]
    t0 = time.time()
    while time.time() - t0 < timeout:
        phase = requests.get(job + "/phase", timeout=60).text.strip()
        if phase == "COMPLETED":
            break
        if phase in ("ERROR", "ABORTED"):
            raise RuntimeError(f"job {phase}: " + requests.get(job + "/error", timeout=60).text[:800])
        time.sleep(10)
    else:
        raise TimeoutError(f"job {job} still {phase} after {timeout}s")
    res = requests.get(job + "/results/result", timeout=300)
    res.raise_for_status()
    if "QUERY_STATUS\" value=\"ERROR\"" in res.text:
        raise RuntimeError(res.text[:800])
    return res


def mast(store, adql, timeout=900):
    request_utc = datetime.now(timezone.utc).isoformat()
    for attempt in range(4):
        try:
            r = _mast_async(adql, timeout)
            break
        except Exception as exc:
            if attempt == 3:
                raise
            print(f"  retry after: {str(exc)[:200]}", file=sys.stderr)
            time.sleep(20)
    rows = [[c for c in re.findall(r"<TD>([^<]*)</TD>", tr)]
            for tr in re.findall(r"<TR>.*?</TR>", r.text, flags=re.S)]
    store.store(service_url=MAST_TAP, query=adql, request_utc=request_utc,
                response_bytes=r.content, row_count=len(rows), http_status=r.status_code)
    return rows


def campaign_ranges(store):
    """Actual per-campaign / prime data ranges (MJD) from MAST CAOM."""
    rows = mast(store, "SELECT obs_collection, sequence_number, dataproduct_type, "
                       "MIN(t_min) AS tmin, MAX(t_max) AS tmax, COUNT(*) AS n "
                       "FROM dbo.obspointing WHERE obs_collection IN ('K2','Kepler') "
                       "GROUP BY obs_collection, sequence_number, dataproduct_type "
                       "ORDER BY obs_collection, sequence_number, dataproduct_type")
    out = {}
    for coll, seq, dpt, tmin, tmax, n in rows:
        key = "prime" if coll == "Kepler" else f"c{int(seq)}"
        out.setdefault(key, {})[dpt] = {"t_min_mjd": float(tmin), "t_max_mjd": float(tmax), "n": int(n)}
    return out


def load_fovs():
    from K2fov import fields, fov
    fovs = {}
    for c in range(20):                       # C20 never flew (retired 2018-10-30)
        info = fields.getFieldInfo(c)
        fovs[f"c{c}"] = {"fov": [fields.getKeplerFov(c)], "ra": info["ra"], "dec": info["dec"],
                         "roll": info["roll"], "planned": (info["start"], info["stop"])}
    prime = []
    for season in range(4):
        prime.append(fov.KeplerFov(PRIME_RA, PRIME_DEC,
                                   fov.getFovAngleFromSpacecraftRoll(20.0 + 90.0 * season),
                                   brokenChannels=[]))
    fovs["prime"] = {"fov": prime, "ra": PRIME_RA, "dec": PRIME_DEC, "roll": "20+90*season",
                     "planned": ("2009-05-13", "2013-05-14")}
    return fovs


def on_silicon(entry, ra, dec):
    hits = []
    for i, f in enumerate(entry["fov"]):
        try:
            hit = bool(f.isOnSilicon(ra, dec, padding_pix=PAD_PIX))
        except Exception:
            hit = False
        hits.append(hit)
    return hits


def angsep(ra1, dec1, ra2, dec2):
    r1, d1, r2, d2 = map(np.radians, (ra1, dec1, ra2, dec2))
    c = np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2)
    return float(np.degrees(np.arccos(np.clip(c, -1, 1))))


def sc_antisun(sc, mjd):
    """Anti-sun unit vector (ICRS) and heliocentric distance from the spacecraft."""
    xyz = np.array([np.interp(mjd, sc["mjd_utc"], sc["xyz_au"][:, k]) for k in range(3)])
    sun = get_body_barycentric("sun", Time(mjd, format="mjd", scale="utc")).xyz.to_value("AU")
    v = xyz - sun
    r = float(np.linalg.norm(v))
    return v / r, r


def radec_to_vec(ra, dec):
    r, d = np.radians(ra), np.radians(dec)
    return np.array([np.cos(d) * np.cos(r), np.cos(d) * np.sin(r), np.sin(d)])


def vec_to_radec(v):
    return float(np.degrees(np.arctan2(v[1], v[0])) % 360), float(np.degrees(np.arcsin(v[2])))


def window_half_d(ev, r_au):
    b = float(ev["b_min_au"])
    return np.sqrt(max(r_au * r_au - b * b, 0.0)) * KM_PER_AU / float(ev["v_perp_km_s"]) / 86400.0


def earth_shift(ev, sc_mjd):
    """Matched Earth-center events (same target, direction, nearest t_ca)."""
    te = Table.read(EVENTS_EARTH)
    te_mjd = Time(list(te["t_ca_utc"]), format="isot", scale="utc").mjd
    idx = defaultdict(list)
    for i, (tid, ld) in enumerate(zip(te["target_id"], te["link_direction"])):
        idx[(str(tid), str(ld))].append(i)
    dt, db, unmatched = [], [], 0
    for i, e in enumerate(ev):
        cand = idx[(str(e["target_id"]), str(e["link_direction"]))]
        if not cand:
            unmatched += 1
            continue
        j = min(cand, key=lambda k: abs(te_mjd[k] - sc_mjd[i]))
        d = te_mjd[j] - sc_mjd[i]
        if abs(d) > 120:
            unmatched += 1
            continue
        dt.append(d)
        db.append(float(te["b_min_au"][j]) - float(e["b_min_au"]))
    dt, db = np.array(dt), np.array(db)
    return {"n_matched": int(len(dt)), "n_unmatched_within_120d": unmatched,
            "abs_dt_ca_days": {"median": float(np.median(np.abs(dt))), "p90": float(np.percentile(np.abs(dt), 90)),
                               "max": float(np.abs(dt).max())},
            "abs_db_au": {"median": float(np.median(np.abs(db))), "max": float(np.abs(db).max())},
            "abs_db_rsun_max": float(np.abs(db).max() / RSUN_AU)}


def main():
    RUN.mkdir(parents=True, exist_ok=True)
    store = SnapshotStore(RUN)
    sc = np.load(SC_TABLE)
    fovs = load_fovs()
    ranges = campaign_ranges(store)
    for k in fovs:
        ts = ranges.get(k, {}).get("timeseries")
        if ts is None:
            raise SystemExit(f"no MAST timeseries range for {k}")
        p0, p1 = (Time(d, scale="utc").mjd for d in fovs[k]["planned"])
        lo, hi = max(ts["t_min_mjd"], p0 - PLAN_CLIP_D), min(ts["t_max_mjd"], p1 + PLAN_CLIP_D)
        if (lo, hi) != (ts["t_min_mjd"], ts["t_max_mjd"]):
            print(f"  {k}: MAST range {ts['t_min_mjd']:.1f}-{ts['t_max_mjd']:.1f} clipped to planned+-{PLAN_CLIP_D:.0f} d -> {lo:.1f}-{hi:.1f}")
        fovs[k]["mjd"] = (lo - DATE_PAD_D, hi + DATE_PAD_D)
        fovs[k]["mast_raw_mjd"] = [ts["t_min_mjd"], ts["t_max_mjd"]]
        fovs[k]["n_timeseries"] = ts["n"]
        fovs[k]["ffi"] = ranges[k].get("image") or ranges[k].get("cube")
    print("campaign data ranges from MAST:",
          {k: (round(v["mjd"][0], 1), round(v["mjd"][1], 1)) for k, v in fovs.items()})

    # ---- campaign boresight vs the spacecraft's anti-sun point ------
    campaign_geom = {}
    for k, v in fovs.items():
        lo, hi = v["mjd"]
        samples = np.linspace(lo + DATE_PAD_D, hi - DATE_PAD_D, 25)
        el, anti = [], []
        for m in samples:
            a, _ = sc_antisun(sc, m)
            ang = angsep(*vec_to_radec(a), v["ra"], v["dec"])
            anti.append(ang)
            el.append(180.0 - ang)      # sun elongation of the boresight
        campaign_geom[k] = {"ra": v["ra"], "dec": v["dec"], "roll": v["roll"],
                            "planned": v["planned"], "data_mjd": [round(lo + DATE_PAD_D, 3), round(hi - DATE_PAD_D, 3)],
                            "mast_raw_mjd": v["mast_raw_mjd"],
                            "n_timeseries": v["n_timeseries"], "ffi": v["ffi"],
                            "boresight_sun_elongation_deg": {"start": round(el[0], 1), "mid": round(el[12], 1),
                                                             "end": round(el[-1], 1)},
                            "boresight_to_antisun_min_deg": round(min(anti), 1)}
    print("boresight-to-anti-sun min (deg):",
          {k: v["boresight_to_antisun_min_deg"] for k, v in campaign_geom.items()})

    # ---- events -------------------------------------------------------
    t = Table.read(EVENTS)
    mjd_all = Time(list(t["t_ca_utc"]), format="isot", scale="utc").mjd
    side = np.asarray(t["axis_distance_au"])
    is_a_all = np.asarray(t["link_direction"]) == "inbound"
    not_sunward = np.where(is_a_all, side > 0, side < 0)
    ev = t[not_sunward]
    mjd = mjd_all[not_sunward]
    is_a = is_a_all[not_sunward]
    b_au = np.asarray(ev["b_min_au"])
    ra_pos = np.where(is_a, ev["star_icrs_ra_deg"], ev["relay_icrs_ra_deg"])
    dec_pos = np.where(is_a, ev["star_icrs_dec_deg"], ev["relay_icrs_dec_deg"])
    print(f"{len(t)} kepler_v1 events; {len(ev)} not sunward")

    scope = {}
    for ch, rung, lim in RUNGS:
        m = (is_a if ch == "A" else ~is_a) & (b_au <= lim)
        per = defaultdict(int)
        for tid in ev["target_id"][m]:
            per[str(tid)] += 1
        scope[f"{ch}_{rung}"] = {"events": int(m.sum()), "targets": len(per), "per_target": dict(sorted(per.items()))}
    print("era scope:", {k: (v["events"], v["targets"]) for k, v in scope.items()})

    # per-event: fields active at t_ca (+- the widest window), on-silicon, diagnostics
    per_event, hits_any_time, nearest = [], [], []
    n_active = 0
    for i in range(len(ev)):
        e = ev[i]
        half = window_half_d(e, 0.1) if b_au[i] <= 0.1 else 0.5
        active = [k for k, v in fovs.items() if v["mjd"][0] <= mjd[i] + half and v["mjd"][1] >= mjd[i] - half]
        anti, r_helio = sc_antisun(sc, float(mjd[i]))
        off = angsep(*vec_to_radec(anti), float(ra_pos[i]), float(dec_pos[i]))
        row = {"event_id": str(e["event_id"]), "target_id": str(e["target_id"]),
               "channel": "A" if is_a[i] else "B", "t_ca_utc": str(e["t_ca_utc"]),
               "b_au": float(b_au[i]), "b_rsun": float(e["b_min_solar_radii"]),
               "ra": float(ra_pos[i]), "dec": float(dec_pos[i]),
               "antisun_offset_deg": round(off, 2), "r_helio_au": round(r_helio, 4),
               "fields_active": active, "on_silicon": {}, "boresight_sep_deg": {}}
        if active:
            n_active += 1
        for k in active:
            row["on_silicon"][k] = on_silicon(fovs[k], float(ra_pos[i]), float(dec_pos[i]))
            row["boresight_sep_deg"][k] = round(angsep(float(ra_pos[i]), float(dec_pos[i]), fovs[k]["ra"], fovs[k]["dec"]), 2)
        # sky-only (time ignored): does any field ever contain the position?
        sky = {k: on_silicon(fovs[k], float(ra_pos[i]), float(dec_pos[i])) for k in fovs}
        row["on_silicon_any_time"] = sorted(k for k, v in sky.items() if any(v))
        seps = {k: angsep(float(ra_pos[i]), float(dec_pos[i]), fovs[k]["ra"], fovs[k]["dec"]) for k in fovs}
        kmin = min(seps, key=seps.get)
        row["nearest_field_any_time"] = {"field": kmin, "boresight_sep_deg": round(seps[kmin], 2)}
        per_event.append(row)
        if row["on_silicon_any_time"]:
            hits_any_time.append(row)
        nearest.append(seps[kmin])

    hits = [r for r in per_event if any(any(v) for v in r["on_silicon"].values())]
    min_active_sep = min((s for r in per_event for s in r["boresight_sep_deg"].values()), default=None)
    print(f"events with an active field at t_ca: {n_active}; on-silicon during the field: {len(hits)}; "
          f"on silicon at some other time: {len(hits_any_time)}; min boresight sep while active: {min_active_sep}")

    # ---- Earth-center shift diagnostic --------------------------------
    shift = earth_shift(ev, mjd)
    print("Earth-center matched shift:", shift)

    # ---- archive box counts per <=0.1 AU unit --------------------------
    units = {}
    for i in np.where(b_au <= 0.1)[0]:
        units.setdefault((str(ev["target_id"][i]), "A" if is_a[i] else "B"), []).append(int(i))
    archive = {}
    for (tid, ch), idx in sorted(units.items()):
        ra = float(np.mean(ra_pos[idx])); dec = float(np.mean(dec_pos[idx]))
        dra = BOX_HALF_DEG / max(np.cos(np.radians(dec)), 0.2)
        ra_lo, ra_hi = ra - dra, ra + dra
        if ra_lo < 0 or ra_hi > 360:
            cond = (f"(s_ra > {(ra_lo) % 360:.4f} OR s_ra < {(ra_hi) % 360:.4f})")
        else:
            cond = f"s_ra BETWEEN {ra_lo:.4f} AND {ra_hi:.4f}"
        rows = mast(store, "SELECT obs_collection, sequence_number, COUNT(*) AS n, MIN(t_min) AS tmin, MAX(t_max) AS tmax "
                           f"FROM dbo.obspointing WHERE obs_collection IN ('K2','Kepler') AND {cond} "
                           f"AND s_dec BETWEEN {dec - BOX_HALF_DEG:.4f} AND {dec + BOX_HALF_DEG:.4f} "
                           "GROUP BY obs_collection, sequence_number")
        spread = max(angsep(ra, dec, float(ra_pos[j]), float(dec_pos[j])) for j in idx)
        archive[f"{tid}:{ch}"] = {"ra": ra, "dec": dec, "n_events": len(idx), "position_spread_deg": round(spread, 3),
                                  "box_half_deg": BOX_HALF_DEG,
                                  "products_by_campaign": [{"collection": r[0], "sequence": int(r[1]), "n": int(r[2]),
                                                            "t_min_mjd": float(r[3]), "t_max_mjd": float(r[4])} for r in rows],
                                  "event_mjd": [round(float(mjd[j]), 2) for j in idx]}
        print(f"  {tid}:{ch} ({ra:.2f},{dec:+.2f}) events={len(idx)} products-in-box={[(r[0], r[1], r[2]) for r in rows]}")

    out = {"generated_utc": datetime.now(timezone.utc).isoformat(), "events_input": str(EVENTS.relative_to(REPO)),
           "observer": "kepler-spacecraft (Horizons -227)", "k2fov_padding_pix": PAD_PIX, "date_pad_days": DATE_PAD_D,
           "n_events_total": len(t), "n_events_not_sunward": len(ev), "era_scope": scope,
           "campaigns": campaign_geom, "n_events_with_active_field": n_active,
           "on_silicon_during_field": hits, "on_silicon_any_time": hits_any_time,
           "min_boresight_sep_while_active_deg": min_active_sep,
           "nearest_field_any_time_min_deg": round(float(min(nearest)), 2),
           "antisun_offset_deg_by_rung": {f"{ch}_{rung}": {"max": round(max((r["antisun_offset_deg"] for r in per_event
                                                                           if r["channel"] == ch and r["b_au"] <= lim), default=0.0), 2)}
                                          for ch, rung, lim in RUNGS},
           "earth_center_shift": shift, "archive_box_counts": archive, "events": per_event}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1))
    print("wrote", OUT.relative_to(REPO))


if __name__ == "__main__":
    main()
