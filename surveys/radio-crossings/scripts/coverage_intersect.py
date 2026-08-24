"""Radio crossings coverage intersection (plan §11.2 step 7).

Metadata-only: no signal statistic is formed anywhere. Two archives,
two geometric arms each:

* Breakthrough Listen Open Data API (per-file records with target,
  ra/decl, mjd, telescope, center_freq): one position cone per
  registry target at the star (channel A, on-star/uplink) and one at
  the antipode (channel B, downlink pre-lens). A pointing "covers" a
  position if their separation is under half the primary beam at that
  file's center frequency (HPBW ~ 1.02 lambda/D; D: GBT 100 m,
  Parkes 64 m, MeerKAT 13.5 m dish FoV) and its MJD lies inside a
  crossing window of the given rung. Files are deduplicated to
  observations in 5-minute bins per (target, telescope).
* VLASS (dyn_summary tile table, 4 epochs 2018-2026): the antipode
  (and star) position's tile per event, epoch observing date vs the
  window with a +/-3 d tile-scheduling tolerance (declared).

Windows from crossings/universal_v1 (Earth-center; topocentric fine),
era 2016-01-01 .. 2028-01-01 for BL, VLASS per-tile dates. Rungs: the
frozen programme ladder (A 0.1/1.0 AU; B 1.2/2.5 Rsun + 0.1 AU). Raw
API responses snapshotted under runs/radio-crossings/.
"""

from __future__ import annotations

import json
import sys
import time as _time
from collections import defaultdict
from pathlib import Path

import numpy as np
import requests
from astropy.table import Table
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from sglsurvey.snapshots import SnapshotStore

API = "http://seti.berkeley.edu/opendata/api/query-files"
VLASS_URL = "https://archive-new.nrao.edu/vlass/VLASS_dyn_summary.php"
KM_PER_AU = 1.495978707e8
BL_ERA = (57388.0, 61771.0)          # 2016-01-01 .. 2028-01-01
LADDER = {"A": (0.1, 1.0), "B": (0.0056, 0.0116, 0.1)}
CONE_DEG = 0.5
DISH_M = {"GBT": 100.0, "Parkes": 64.0, "MeerKAT": 13.5}
OBS_BIN_DAYS = 0.0035                # ~5 min
VLASS_TOL_DAYS = 3.0
OUT = REPO / "surveys" / "radio-crossings" / "results"
RUN = REPO / "runs" / "radio-crossings"


def load_events():
    t = Table.read(REPO / "crossings" / "universal_v1" / "events.ecsv")
    t["mjd"] = Time(list(t["t_ca_utc"]), format="isot", scale="utc").mjd
    side = np.asarray(t["axis_distance_au"])
    A = t[(t["link_direction"] == "inbound") & (side > 0)]
    B = t[(t["link_direction"] == "outbound") & (side < 0)]
    return {"A": A, "B": B}


def pos_of(ch, ev):
    if ch == "A":
        return float(ev["star_icrs_ra_deg"]), float(ev["star_icrs_dec_deg"])
    return float(ev["relay_icrs_ra_deg"]), float(ev["relay_icrs_dec_deg"])


def half_days(b, r, vp):
    return float(np.sqrt(r * r - b * b) * KM_PER_AU / vp / 86400.0)


def half_beam_deg(telescope, center_freq_mhz):
    d = DISH_M.get(telescope)
    if d is None or not center_freq_mhz or center_freq_mhz <= 0:
        return 0.15          # conservative default, noted per row
    lam = 299.792458 / float(center_freq_mhz)      # m
    return 0.5 * np.degrees(1.02 * lam / d)


def sep_deg(ra1, de1, ra2, de2):
    c1, c2 = np.radians(de1), np.radians(de2)
    dra = np.radians(ra1 - ra2)
    x = (np.sin(c1) * np.sin(c2)
         + np.cos(c1) * np.cos(c2) * np.cos(dra))
    return float(np.degrees(np.arccos(np.clip(x, -1, 1))))


def bl_cone(session, store, label, ra, dec):
    params = {"target": "", "pos-ra": f"{ra:.5f}",
              "pos-dec": f"{dec:.5f}", "pos-rad": f"{CONE_DEG}",
              "limit": "100000"}
    for attempt in range(3):
        try:
            r = session.get(API, params=params, timeout=300)
            r.raise_for_status()
            doc = r.json()
            break
        except Exception as exc:
            print(f"  [{label}] attempt {attempt+1}: {exc}", flush=True)
            _time.sleep(20)
    else:
        return None
    store.store(service_url=API,
                query="&".join(f"{k}={v}" for k, v in params.items()),
                request_utc=Time.now().isot,
                response_bytes=r.content,
                row_count=len(doc.get("data", [])),
                http_status=r.status_code)
    return doc.get("data", [])


def dedup_obs(files):
    seen = {}
    for f in files:
        try:
            key = (f["target"], f["telescope"],
                   round(float(f["mjd"]) / OBS_BIN_DAYS))
        except Exception:
            continue
        if key not in seen:
            seen[key] = f
    return list(seen.values())


def bl_arm():
    events = load_events()
    session = requests.Session()
    store = SnapshotStore(RUN)
    rows, target_summary = [], []
    for ch, tab in events.items():
        era = (tab["mjd"] >= BL_ERA[0]) & (tab["mjd"] <= BL_ERA[1])
        tab = tab[era]
        tab = tab[np.asarray(tab["b_min_au"]) < max(LADDER[ch])]
        for tid in sorted(set(str(x) for x in tab["target_id"])):
            sub = tab[np.asarray([str(x) == tid
                                  for x in tab["target_id"]])]
            ras = [pos_of(ch, e)[0] for e in sub]
            des = [pos_of(ch, e)[1] for e in sub]
            ra0, de0 = float(np.median(ras)), float(np.median(des))
            files = bl_cone(session, store, f"{ch}/{tid}", ra0, de0)
            if files is None:
                target_summary.append({"channel": ch, "target_id": tid,
                                       "status": "query_failed"})
                continue
            obs = dedup_obs(files)
            n_in = 0
            for ev in sub:
                b, vp = float(ev["b_min_au"]), float(ev["v_perp_km_s"])
                ra_e, de_e = pos_of(ch, ev)
                for r in LADDER[ch]:
                    if b >= r:
                        continue
                    hd = half_days(b, r, vp)
                    lo, hi = float(ev["mjd"]) - hd, float(ev["mjd"]) + hd
                    for o in obs:
                        m = float(o["mjd"])
                        if not (lo <= m <= hi):
                            continue
                        s = sep_deg(float(o["ra"]), float(o["decl"]),
                                    ra_e, de_e)
                        hb = half_beam_deg(o["telescope"],
                                           o.get("center_freq"))
                        if s > CONE_DEG:
                            continue
                        rows.append({
                            "channel": ch, "target_id": tid,
                            "event_id": str(ev["event_id"]),
                            "radius_au": r, "t_ca_mjd": float(ev["mjd"]),
                            "window_days": 2 * hd, "obs_mjd": m,
                            "bl_target": o["target"],
                            "telescope": o["telescope"],
                            "center_freq_mhz": float(
                                o.get("center_freq") or 0.0),
                            "sep_deg": round(s, 4),
                            "half_beam_deg": round(hb, 4),
                            "in_beam": bool(s <= hb),
                        })
                        n_in += 1
            target_summary.append({
                "channel": ch, "target_id": tid, "status": "ok",
                "n_files": len(files), "n_obs": len(obs),
                "n_inwindow_rows": n_in})
            print(f"[{ch}] {tid}: {len(files)} files -> {len(obs)} obs, "
                  f"{n_in} in-window rows", flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    if rows:
        Table(rows=rows).write(OUT / "bl_inwindow_v1.ecsv",
                               format="ascii.ecsv", overwrite=True)
    (OUT / "bl_targets_v1.json").write_text(
        json.dumps(target_summary, indent=1) + "\n")
    return rows, target_summary


def vlass_arm():
    events = load_events()
    resp = requests.get(VLASS_URL, timeout=120)
    resp.raise_for_status()
    RUN.mkdir(parents=True, exist_ok=True)
    (RUN / "vlass_dyn_summary.txt").write_text(resp.text)
    tiles = []
    for line in resp.text.splitlines():
        p = line.split()
        if len(p) >= 7 and p[0].startswith("T") and "t" in p[0]:
            try:
                dec0, dec1 = float(p[1]), float(p[2])
                ra0, ra1 = float(p[3]) * 15, float(p[4]) * 15
                epoch, date = p[5], p[6]
                mjd = Time(date, format="iso", scale="utc").mjd \
                    if date[:2] == "20" else None
            except Exception:
                continue
            if mjd is not None:
                tiles.append((dec0, dec1, ra0, ra1, epoch, float(mjd)))
    print(f"VLASS tiles with dates: {len(tiles)}")

    rows = []
    for ch, tab in events.items():
        tab = tab[np.asarray(tab["b_min_au"]) < max(LADDER[ch])]
        for ev in tab:
            ra_e, de_e = pos_of(ch, ev)
            if de_e < -40.0:
                continue
            b, vp = float(ev["b_min_au"]), float(ev["v_perp_km_s"])
            hits = [t for t in tiles
                    if t[0] <= de_e < t[1] and t[2] <= ra_e < t[3]]
            for r in LADDER[ch]:
                if b >= r:
                    continue
                hd = half_days(b, r, vp)
                lo = float(ev["mjd"]) - hd - VLASS_TOL_DAYS
                hi = float(ev["mjd"]) + hd + VLASS_TOL_DAYS
                for dec0, dec1, ra0, ra1, epoch, mjd in hits:
                    if lo <= mjd <= hi:
                        rows.append({
                            "channel": ch,
                            "target_id": str(ev["target_id"]),
                            "event_id": str(ev["event_id"]),
                            "radius_au": r,
                            "t_ca_mjd": float(ev["mjd"]),
                            "window_days": 2 * hd,
                            "vlass_epoch": epoch,
                            "tile_obs_mjd": mjd,
                            "offset_days": round(mjd - float(ev["mjd"]),
                                                 1)})
    if rows:
        Table(rows=rows).write(OUT / "vlass_inwindow_v1.ecsv",
                               format="ascii.ecsv", overwrite=True)
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    bl_rows, bl_targets = ([], [])
    if "--skip-bl" not in sys.argv:
        bl_rows, bl_targets = bl_arm()
    vl_rows = vlass_arm() if "--skip-vlass" not in sys.argv else []

    def cnt(rows, ch, r, beam_only):
        sel = [x for x in rows if x["channel"] == ch
               and x["radius_au"] == r
               and (x.get("in_beam", True) or not beam_only)]
        return {"rows": len(sel),
                "events": len(set(x["event_id"] for x in sel)),
                "targets": len(set(x["target_id"] for x in sel))}

    summary = {"bl": {}, "vlass": {}}
    for ch in LADDER:
        for r in LADDER[ch]:
            summary["bl"][f"{ch}_{r}"] = {
                "in_beam": cnt(bl_rows, ch, r, True),
                "in_cone_0p5deg": cnt(bl_rows, ch, r, False)}
            summary["vlass"][f"{ch}_{r}"] = cnt(vl_rows, ch, r, False)
    (OUT / "coverage_v1_summary.json").write_text(
        json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
