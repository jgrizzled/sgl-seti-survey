"""High-energy archive reachability recon (plan §5.26; 2026-09-07).

Chandra / XMM-Newton / eROSITA / Swift / Fermi-LAT for (i) opportunistic
high-energy coincidence with the universal crossing windows (Pipeline B,
channels A on-star and B antipode, rungs 1.2 Rsun / 2.5 Rsun / 0.1 AU,
flat-chord windows as in every crossings survey) and (ii) catalogue-level
persistent-source screens on the anti-star corridors (Pipeline A step 1).

Stages (each writes its own JSON under results/ and snapshots every
response verbatim under runs/highenergy-crossings/recon/):

  masters   HEASARC Xamin TAP cones on chanmaster / xmmmaster / swiftmastr
            at the era-mean channel positions of the 7 deep-family targets,
            then a local flat-chord window intersection per rung.
  slew      XSA TAP xsa.v_slew_exposure footprint-polygon CONTAINS at the
            channel positions (XMM slew survey, 2001-), window intersect.
  bat       swiftmastr per-event window queries (all rungs share the
            0.1 AU window) for pointings whose BAT survey/event exposure
            is non-zero within 40 deg of the channel position — the
            wide-field coded-mask coverage of the windows.
  lat       Fermi-LAT weekly spacecraft files for every grazing-family
            window (<= 2.5 Rsun, both channels) and for the 0.1 AU windows
            of one target: in-FoV livetime (theta <= 60 deg, zenith <= 100
            deg, DATA_QUAL > 0, LAT_CONFIG == 1) inside each rung window.
  catalogs  Cones on the 88 anti-star corridor positions and the 7 A-side
            star positions in 5XMM-DR15 (+stack), CSC 2.1, 2SXPS/LSXPS,
            eRASS1 main/hard (HEASARC) + eRODat DR1/DR2 cone + upper-limit
            APIs, BAT 157-month, 4FGL-DR4.
  routes    One data-product route per archive (listing / HEAD only).

Coverage counting only — no event is touched, no statistic is formed.
"""
from __future__ import annotations

import csv
import io
import json
import re
import sys
import time
import urllib.parse
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import requests
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
from sglsurvey.snapshots import SnapshotStore  # noqa: E402

RUN = REPO / "runs" / "highenergy-crossings" / "recon"
OUT_DIR = Path(__file__).resolve().parents[1] / "results"
EVENTS = REPO / "crossings" / "universal_v1" / "events.ecsv"

HEASARC_TAP = "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync"
CSC_TAP = "https://cda.cfa.harvard.edu/csc21tap/sync"
XSA_TAP = "https://nxsa.esac.esa.int/tap-server/tap/sync"
ERODAT = "https://erosita.mpe.mpg.de/erodat"
FERMI_FTP = "https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/weekly"

KM_PER_AU = 1.495978707e8
RSUN_AU = 0.00465047
RUNGS = [("1.2Rsun", 1.2 * RSUN_AU), ("2.5Rsun", 2.5 * RSUN_AU), ("0.1AU", 0.1)]
DEEP = ["gj-1276", "gj-908", "ross-128", "ross-154", "teegarden",
        "van-maanen", "wolf-359"]

# mission eras (UTC) used to bound the event population per archive
ERAS = {
    "chandra": ("1999-07-23", "2026-12-31"),
    "xmm": ("2000-01-01", "2026-12-31"),
    "swift": ("2004-11-20", "2026-12-31"),
    "fermi": ("2008-08-04", "2026-12-31"),
    "erosita_dr1": ("2019-12-12", "2020-06-11"),
    "erosita_dr2": ("2019-12-12", "2021-06-16"),
}
# discovery cone (arcmin) per pointed instrument: FoV half-size + margin
DISC_ARCMIN = {"chanmaster": 25.0, "xmmmaster": 20.0, "swiftmastr": 15.0}
BAT_DISC_DEG = 40.0
LAT_THETA_MAX = 60.0
LAT_ZENITH_MAX = 100.0
CORRIDOR_ARCMIN = 7.0     # 550 AU parallax = 375 arcsec; 7' contains every corridor
STAR_ARCMIN = 2.0
FGL_ARCMIN = 30.0


def now_utc():
    return datetime.now(timezone.utc).isoformat()


def get(store, url, params=None, label="", timeout=300, tries=3, sleep=8,
        method="GET", json_body=None, allow_redirects=True):
    last = None
    for k in range(tries):
        try:
            if method == "POST":
                r = requests.post(url, params=params, json=json_body,
                                  timeout=timeout, allow_redirects=allow_redirects)
            else:
                r = requests.get(url, params=params, timeout=timeout,
                                 allow_redirects=allow_redirects)
            q = urllib.parse.urlencode(params or {})
            if json_body is not None:
                q += " JSON:" + json.dumps(json_body)[:4000]
            store.store(service_url=url, query=q, request_utc=now_utc(),
                        response_bytes=r.content, row_count=None,
                        http_status=r.status_code)
            if r.status_code >= 500:
                raise RuntimeError(f"HTTP {r.status_code}")
            return r.status_code, r
        except Exception as exc:
            last = exc
            print(f"  [{label}] attempt {k+1}: {exc}", file=sys.stderr, flush=True)
            time.sleep(sleep)
    return None, last


def tap_csv(store, base, adql, label, timeout=300):
    status, r = get(store, base, {"REQUEST": "doQuery", "LANG": "ADQL",
                                  "FORMAT": "csv", "QUERY": adql}, label, timeout)
    if status is None or status >= 400:
        return None, f"HTTP {status} {str(getattr(r, 'text', r))[:300]}"
    text = r.text
    if text.lstrip().startswith("<"):
        m = re.search(r'QUERY_STATUS"? value="?ERROR"?>([^<]*)', text)
        return None, (m.group(1) if m else text[:300])
    rows = list(csv.DictReader(io.StringIO(text)))
    return rows, None


def heasarc(store, adql, label):
    """HEASARC Xamin TAP: CSV is not supported by the sync endpoint, use
    FORMAT=text (pipe-separated, trailing row/column count lines)."""
    status, r = get(store, HEASARC_TAP,
                    {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "text",
                     "QUERY": adql}, label)
    if status is None or status >= 400:
        return None, f"HTTP {status}"
    text = r.text
    if text.lstrip().startswith("<"):
        m = re.search(r'QUERY_STATUS"? value="?ERROR"?>([^<]*)', text)
        return None, (m.group(1) if m else text[:300])
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        return [], None
    hdr = [h.strip() for h in lines[0].split("|")]
    rows = []
    for ln in lines[1:]:
        if ln.startswith("Number of "):
            continue
        parts = [p.strip() for p in ln.split("|")]
        if len(parts) != len(hdr):
            continue
        rows.append(dict(zip(hdr, parts)))
    return rows, None


def ffloat(x, default=None):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def sep_deg(ra1, de1, ra2, de2):
    r1, d1, r2, d2 = map(np.radians, (ra1, de1, ra2, de2))
    c = np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(r1 - r2)
    return float(np.degrees(np.arccos(np.clip(c, -1.0, 1.0))))


# ---------------------------------------------------------------- events
def load_events():
    t = Table.read(EVENTS)
    t["mjd"] = Time(list(t["t_ca_utc"]), format="isot", scale="utc").mjd
    side = np.asarray(t["axis_distance_au"])
    ld = np.asarray(t["link_direction"])
    ok = np.asarray(t["validity"]) == "valid"
    A = t[(ld == "inbound") & (side > 0) & ok]
    B = t[(ld == "outbound") & (side < 0) & ok]
    A = A[np.asarray(A["b_min_au"]) < 0.1]
    B = B[np.asarray(B["b_min_au"]) < 0.1]
    return {"A": A, "B": B}


def source_radec(ch, ev):
    if ch == "A":
        return float(ev["star_icrs_ra_deg"]), float(ev["star_icrs_dec_deg"])
    return float(ev["relay_icrs_ra_deg"]), float(ev["relay_icrs_dec_deg"])


def windows_for(ev):
    b = float(ev["b_min_au"]); vp = float(ev["v_perp_km_s"]); tc = float(ev["mjd"])
    out = {}
    for name, lim in RUNGS:
        if b >= lim:
            out[name] = None
        else:
            hw = np.sqrt(lim * lim - b * b) * KM_PER_AU / vp / 86400.0
            out[name] = (tc - hw, tc + hw, hw)
    return out


def positions(events):
    """Era-mean channel positions per (channel, target) plus per-target
    anti-star corridor positions for all 88 endpoints."""
    pos = {}
    for ch, tab in events.items():
        for tid in sorted(set(map(str, tab["target_id"]))):
            s = tab[tab["target_id"] == tid]
            ra = float(np.median([source_radec(ch, e)[0] for e in s]))
            de = float(np.median([source_radec(ch, e)[1] for e in s]))
            pos[(ch, tid)] = (ra, de, len(s))
    return pos


def all_endpoints():
    t = Table.read(EVENTS)
    out = {}
    for tid in sorted(set(map(str, t["target_id"]))):
        s = t[t["target_id"] == tid]
        ra = float(np.median(s["star_icrs_ra_deg"]))
        de = float(np.median(s["star_icrs_dec_deg"]))
        out[tid] = {"star": (ra, de), "anti": ((ra + 180.0) % 360.0, -de)}
    return out


def in_era(mjd, era):
    lo, hi = [Time(x, format="iso", scale="utc").mjd for x in era]
    return lo <= mjd <= hi


# ---------------------------------------------------------------- masters
MASTER_COLS = {
    "chanmaster": 'obsid, status, name, ra, dec, "time", detector, grating, exposure, type, public_date, data_mode',
    "xmmmaster": 'obsid, status, name, ra, dec, "time", end_time, duration, pn_time, mos1_time, mos2_time, pn_mode, public_date',
    "swiftmastr": 'obsid, name, ra, dec, start_time, stop_time, xrt_exposure, xrt_expo_pc, xrt_expo_wt, uvot_exposure, bat_exposure, bat_expo_sv, bat_expo_ev, archive_date',
}


def stage_masters(store, events, pos):
    out = {"stage": "masters", "utc": now_utc(), "cones": [], "units": []}
    for (ch, tid), (ra, de, n) in sorted(pos.items()):
        for table, cols in MASTER_COLS.items():
            r_deg = DISC_ARCMIN[table] / 60.0
            adql = (f"SELECT {cols} FROM {table} WHERE "
                    f"CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra:.6f}, {de:.6f}, {r_deg:.5f}))=1")
            rows, err = heasarc(store, adql, f"{table} {ch} {tid}")
            rec = {"channel": ch, "target": tid, "table": table, "ra": ra,
                   "dec": de, "radius_arcmin": DISC_ARCMIN[table],
                   "n_rows": None if rows is None else len(rows), "error": err,
                   "rows": rows or []}
            out["cones"].append(rec)
            print(f"[masters] {table:10s} {ch} {tid:11s} -> "
                  f"{'ERR ' + str(err) if rows is None else len(rows)}", flush=True)
            if not rows:
                continue
            # local window intersection, all events of this (ch, target)
            tab = events[ch]
            sub = tab[tab["target_id"] == tid]
            for row in rows:
                if table == "swiftmastr":
                    t0 = ffloat(row.get("start_time")); t1 = ffloat(row.get("stop_time"))
                    expo = ffloat(row.get("xrt_exposure"), 0.0)
                elif table == "xmmmaster":
                    t0 = ffloat(row.get("time")); t1 = ffloat(row.get("end_time"))
                    expo = ffloat(row.get("duration"), 0.0)
                    if t1 is None and t0 is not None:
                        t1 = t0 + (expo or 0.0) / 86400.0
                else:
                    t0 = ffloat(row.get("time"))
                    expo = ffloat(row.get("exposure"), 0.0)
                    t1 = None if t0 is None else t0 + (expo or 0.0) / 86400.0
                if t0 is None:
                    continue
                pra, pde = ffloat(row.get("ra")), ffloat(row.get("dec"))
                for ev in sub:
                    w = windows_for(ev)
                    for rung, win in w.items():
                        if win is None:
                            continue
                        w0, w1 = win[0], win[1]
                        if t1 >= w0 and t0 <= w1:
                            era, ede = source_radec(ch, ev)
                            out["units"].append({
                                "channel": ch, "target": tid, "table": table,
                                "rung": rung, "event_id": str(ev["event_id"]),
                                "t_ca_utc": str(ev["t_ca_utc"]),
                                "b_min_rsun": float(ev["b_min_au"]) / RSUN_AU,
                                "window_utc": [Time(w0, format="mjd").isot,
                                               Time(w1, format="mjd").isot],
                                "obsid": row.get("obsid"), "name": row.get("name"),
                                "obs_start_utc": Time(t0, format="mjd").isot,
                                "obs_stop_utc": Time(t1, format="mjd").isot,
                                "offset_days": float(t0 - float(ev["mjd"])),
                                "overlap_s": float(max(0.0, min(t1, w1) - max(t0, w0)) * 86400.0),
                                "exposure_s": expo,
                                "pointing_offset_arcmin": (
                                    None if pra is None else sep_deg(pra, pde, era, ede) * 60.0),
                                "row": row})
    n_units = len(out["units"])
    print(f"[masters] in-window pointed units: {n_units}", flush=True)
    (OUT_DIR / "recon_masters_v0.json").write_text(json.dumps(out, indent=1))
    return out



# ---------------------------------------------------------------- XMM slew
def stage_slew(store, events, pos):
    """XMM-Newton slew survey: every slew exposure whose footprint polygon
    contains the era-mean channel position (XSA TAP, pgsphere CONTAINS),
    intersected with the windows. Slew subimages are ~10 s per position at
    ~90 deg/h; the window test uses the exposure start/end (an upper bound
    on the in-window time; the exact pass time needs the exposure image)."""
    out = {"stage": "slew", "utc": now_utc(), "positions": [], "units": []}
    for (ch, tid), (ra, de, n) in sorted(pos.items()):
        adql = ("SELECT slew_observation_id, slew_exposure_id, instrument, filter, start_utc, end_utc, ra, dec "
                "FROM xsa.v_slew_exposure WHERE "
                f"1=CONTAINS(POINT('ICRS',{ra:.6f},{de:.6f}), slew_exposure_fov_spoly)")
        rows, err = tap_csv(store, XSA_TAP, adql, f"slew {ch} {tid}")
        rec = {"channel": ch, "target": tid, "ra": ra, "dec": de,
               "n_rows": None if rows is None else len(rows), "error": err, "rows": rows or []}
        out["positions"].append(rec)
        sub = events[ch][events[ch]["target_id"] == tid]
        nearest = None
        for row in rows or []:
            try:
                t0 = Time(row["start_utc"].strip(), format="isot", scale="utc").mjd
                t1 = Time(row["end_utc"].strip(), format="isot", scale="utc").mjd
            except Exception:
                continue
            for ev in sub:
                dt = 0.5 * (t0 + t1) - float(ev["mjd"])
                if nearest is None or abs(dt) < abs(nearest[0]):
                    nearest = (dt, row["slew_observation_id"], row["start_utc"])
                for rung, win in windows_for(ev).items():
                    if win is None:
                        continue
                    if t1 >= win[0] and t0 <= win[1]:
                        out["units"].append({"channel": ch, "target": tid, "rung": rung,
                                             "event_id": str(ev["event_id"]), "t_ca_utc": str(ev["t_ca_utc"]),
                                             "b_min_rsun": float(ev["b_min_au"]) / RSUN_AU,
                                             "offset_days": float(0.5 * (t0 + t1) - float(ev["mjd"])), "row": row})
        rec["nearest_days"] = None if nearest is None else nearest[0]
        rec["nearest"] = None if nearest is None else nearest[1:]
        print(f"[slew] {ch} {tid:11s} -> {rec['n_rows']} exposures, nearest {rec['nearest_days']}", flush=True)
    print(f"[slew] in-window slew units: {len(out['units'])}", flush=True)
    (OUT_DIR / "recon_slew_v0.json").write_text(json.dumps(out, indent=1))
    return out

# ---------------------------------------------------------------- BAT
def stage_bat(store, events):
    out = {"stage": "bat", "utc": now_utc(), "disc_deg": BAT_DISC_DEG,
           "events": []}
    era = ERAS["swift"]
    n_q = 0
    for ch, tab in events.items():
        for ev in tab:
            if not in_era(float(ev["mjd"]), era):
                continue
            w = windows_for(ev)
            w0, w1 = w["0.1AU"][0], w["0.1AU"][1]
            ra, de = source_radec(ch, ev)
            adql = ("SELECT obsid, name, ra, dec, start_time, stop_time, bat_exposure, "
                    "bat_expo_sv, bat_expo_ev, xrt_exposure FROM swiftmastr WHERE "
                    f"start_time <= {w1:.5f} AND stop_time >= {w0:.5f} AND "
                    f"(bat_expo_sv > 0 OR bat_expo_ev > 0) AND "
                    f"CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra:.5f}, {de:.5f}, {BAT_DISC_DEG}))=1")
            rows, err = heasarc(store, adql, f"bat {ch} {ev['target_id']} {ev['t_ca_utc'][:10]}")
            n_q += 1
            rec = {"channel": ch, "target": str(ev["target_id"]),
                   "event_id": str(ev["event_id"]), "t_ca_utc": str(ev["t_ca_utc"]),
                   "b_min_rsun": float(ev["b_min_au"]) / RSUN_AU, "ra": ra, "dec": de,
                   "error": err, "n_rows": None if rows is None else len(rows),
                   "rungs": {}}
            for rung, win in w.items():
                if win is None:
                    continue
                r0, r1 = win[0], win[1]
                secs = 0.0; secs30 = 0.0; secs20 = 0.0; n = 0; best = None
                for row in rows or []:
                    t0 = ffloat(row.get("start_time")); t1 = ffloat(row.get("stop_time"))
                    if t0 is None or t1 is None:
                        continue
                    ov = max(0.0, min(t1, r1) - max(t0, r0)) * 86400.0
                    if ov <= 0:
                        continue
                    pra, pde = ffloat(row.get("ra")), ffloat(row.get("dec"))
                    off = sep_deg(pra, pde, ra, de)
                    # scale the BAT exposure by the fraction of the pointing inside the window
                    bexp = (ffloat(row.get("bat_expo_sv"), 0.0) or 0.0) + (ffloat(row.get("bat_expo_ev"), 0.0) or 0.0)
                    span = max((t1 - t0) * 86400.0, 1.0)
                    e = bexp * min(1.0, ov / span)
                    n += 1; secs += e
                    if off <= 30.0:
                        secs30 += e
                    if off <= 20.0:
                        secs20 += e
                    if best is None or off < best[0]:
                        best = (off, row.get("obsid"), Time(t0, format="mjd").isot, e)
                rec["rungs"][rung] = {"window_utc": [Time(r0, format="mjd").isot, Time(r1, format="mjd").isot],
                                      "half_width_d": win[2], "n_pointings": n,
                                      "bat_s_40deg": secs, "bat_s_30deg": secs30,
                                      "bat_s_20deg": secs20,
                                      "best": None if best is None else
                                      {"offset_deg": best[0], "obsid": best[1], "start": best[2], "bat_s": best[3]}}
            out["events"].append(rec)
            g = rec["rungs"].get("2.5Rsun")
            print(f"[bat] {ch} {ev['target_id']:11s} {ev['t_ca_utc'][:10]} b={rec['b_min_rsun']:.2f} "
                  f"rows={rec['n_rows']} 0.1AU s30={rec['rungs']['0.1AU']['bat_s_30deg']:.0f}"
                  + (f" 2.5Rs s30={g['bat_s_30deg']:.0f}" if g else ""), flush=True)
    out["n_queries"] = n_q
    (OUT_DIR / "recon_bat_v0.json").write_text(json.dumps(out, indent=1))
    return out


# ---------------------------------------------------------------- LAT
def fermi_week_table(store):
    rows, err = heasarc(store, 'SELECT week_number, "time", end_time, met_start, met_end FROM fermilweek', "fermilweek")
    if rows is None:
        raise RuntimeError(err)
    weeks = []
    for r in rows:
        weeks.append((int(r["week_number"]), ffloat(r["time"]), ffloat(r["end_time"])))
    return sorted(weeks)


def fermi_sc_file(store, week):
    name = f"lat_spacecraft_weekly_w{week:03d}_p310_v001.fits"
    path = RUN / "fermi_spacecraft" / name
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        url = f"{FERMI_FTP}/spacecraft/{name}"
        for k in range(3):
            try:
                r = requests.get(url, timeout=600)
                if r.status_code == 200:
                    path.write_bytes(r.content)
                    store.store(service_url=url, query="", request_utc=now_utc(),
                                response_bytes=r.content[:0] + name.encode(), row_count=None,
                                http_status=200)
                    break
                print(f"  [fermi sc w{week}] HTTP {r.status_code}", file=sys.stderr, flush=True)
            except Exception as exc:
                print(f"  [fermi sc w{week}] {exc}", file=sys.stderr, flush=True)
            time.sleep(8)
    if not path.exists():
        return None
    return path


MET_EPOCH_MJD = 51910.0 + 7.428703703703703e-4   # 2001-01-01T00:00:00 UTC in MJD(UTC) ~ 51910.00074


def met_to_mjd(met):
    # Fermi MET counts SI seconds from 2001-01-01 00:00:00 UTC (TT-UTC offset ignored; leap seconds
    # since 2001 total 5 s, below the 30-s sampling of the spacecraft file)
    return 51910.0 + met / 86400.0


def lat_livetime(sc_tables, ra, de, w0, w1):
    """In-FoV LAT livetime (s) between MJD w0 and w1 for a position."""
    tot = 0.0; n = 0
    for tab in sc_tables:
        t0 = met_to_mjd(tab["START"]); t1 = met_to_mjd(tab["STOP"])
        m = (t1 > w0) & (t0 < w1)
        if not m.any():
            continue
        t0 = t0[m]; t1 = t1[m]
        frac = np.clip((np.minimum(t1, w1) - np.maximum(t0, w0)) / np.maximum(t1 - t0, 1e-9), 0, 1)
        theta = np.array([sep_deg(a, b, ra, de) for a, b in zip(tab["RA_SCZ"][m], tab["DEC_SCZ"][m])])
        zen = np.array([sep_deg(a, b, ra, de) for a, b in zip(tab["RA_ZENITH"][m], tab["DEC_ZENITH"][m])])
        good = (theta <= LAT_THETA_MAX) & (zen <= LAT_ZENITH_MAX) & (tab["DATA_QUAL"][m] > 0) & (tab["LAT_CONFIG"][m] == 1)
        lt = tab["LIVETIME"][m] * frac
        tot += float(lt[good].sum()); n += int(good.sum())
    return tot, n


def stage_lat(store, events, all_01au=False):
    """all_01au=False: every grazing-family window + the van-maanen 0.1 AU
    windows (recon v0). all_01au=True: the 0.1 AU windows of every event
    of all seven targets (the open-item extension), written to
    recon_lat01_v0.json."""
    out = {"stage": "lat01" if all_01au else "lat", "utc": now_utc(), "theta_max_deg": LAT_THETA_MAX,
           "zenith_max_deg": LAT_ZENITH_MAX, "events": []}
    weeks = fermi_week_table(store)
    out["n_weeks"] = len(weeks)
    out["week_span_utc"] = [Time(weeks[0][1], format="mjd").isot, Time(weeks[-1][2], format="mjd").isot]
    era = ERAS["fermi"]
    cache = {}

    def tables_for(w0, w1):
        need = [wk for wk, a, b in weeks if b >= w0 and a <= w1]
        tabs = []
        for wk in need:
            if wk not in cache:
                p = fermi_sc_file(store, wk)
                if p is None:
                    cache[wk] = None
                else:
                    with fits.open(p) as h:
                        d = h["SC_DATA"].data
                        cache[wk] = {k: np.array(d[k]) for k in
                                     ("START", "STOP", "RA_SCZ", "DEC_SCZ", "RA_ZENITH", "DEC_ZENITH",
                                      "LIVETIME", "DATA_QUAL", "LAT_CONFIG")}
            if cache[wk] is not None:
                tabs.append(cache[wk])
        return need, tabs

    for ch, tab in events.items():
        for ev in tab:
            if not in_era(float(ev["mjd"]), era):
                continue
            b_rs = float(ev["b_min_au"]) / RSUN_AU
            grazing = b_rs < 2.5
            if all_01au:
                grazing = False
            # every grazing-family event; 0.1 AU windows only for van-maanen (sample)
            elif not grazing and str(ev["target_id"]) != "van-maanen":
                continue
            w = windows_for(ev)
            ra, de = source_radec(ch, ev)
            rec = {"channel": ch, "target": str(ev["target_id"]), "event_id": str(ev["event_id"]),
                   "t_ca_utc": str(ev["t_ca_utc"]), "b_min_rsun": b_rs, "ra": ra, "dec": de, "rungs": {}}
            for rung, win in w.items():
                if win is None:
                    continue
                if not grazing and rung != "0.1AU":
                    continue
                need, tabs = tables_for(win[0], win[1])
                lt, n = lat_livetime(tabs, ra, de, win[0], win[1])
                rec["rungs"][rung] = {"window_utc": [Time(win[0], format="mjd").isot, Time(win[1], format="mjd").isot],
                                      "half_width_d": win[2], "weeks": need, "n_weeks_found": len(tabs),
                                      "livetime_s": lt, "n_intervals": n,
                                      "window_s": 2 * win[2] * 86400.0,
                                      "duty": lt / (2 * win[2] * 86400.0)}
            out["events"].append(rec)
            s = " ".join(f"{k}:{v['livetime_s']:.0f}s({v['duty']:.2f})" for k, v in rec["rungs"].items())
            print(f"[lat] {ch} {ev['target_id']:11s} {ev['t_ca_utc'][:10]} b={b_rs:.2f} {s}", flush=True)
    (OUT_DIR / ("recon_lat01_v0.json" if all_01au else "recon_lat_v0.json")).write_text(json.dumps(out, indent=1))
    return out


# ---------------------------------------------------------------- catalogs
HEASARC_CATS = {
    "xmmssc": ("5XMM-DR15 detections", "srcid, detid, obsid, ra, dec, error_radius, ep_8_flux, sc_ep_8_flux, sc_var_flag, \"time\", ep_8_det_ml"),
    "xmmstack": ("5XMM-DR15 stacked", "srcid, ra, dec, radec_err, ep_flux, ep_det_ml, n_obs, n_contrib, var_prob"),
    "swift2sxps": ("2SXPS", "name, ra, dec, err90, exposure, rate0, flux_pow, detflag, fieldflag"),
    "swiftlsxps": ("LSXPS", "name, ra, dec, err90, exposure, rate0, flux_pow, detflag, fieldflag, end_time"),
    "erass1main": ("eRASS1 main", "name, ra, dec, radec_err, det_like_0, ml_flux_1, ml_rate_1, ext"),
    "erass1hard": ("eRASS1 hard", "name, ra, dec, radec_err, det_like_3, ml_flux_3"),
    "swbat157m": ("BAT 157-month", "name, ra, dec, snr, flux, counterpart_name"),
    "fermilpsc": ("4FGL-DR4", "name, ra, dec, semi_major_axis_95, significance, flux_1000, assoc_name_1, class_1"),
}



def csc_box(store, ra, de, r_am, label):
    r_deg = r_am / 60.0
    cosd = max(np.cos(np.radians(de)), 0.05)
    adql = ("SELECT name, ra, dec, err_ellipse_r0, flux_aper_b, significance, var_inter_prob_b "
            "FROM csc21.master_source WHERE "
            f"ra BETWEEN {ra - r_deg / cosd:.6f} AND {ra + r_deg / cosd:.6f} AND "
            f"dec BETWEEN {de - r_deg:.6f} AND {de + r_deg:.6f}")
    status, r = get(store, CSC_TAP, {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": adql}, label)
    if status is None or status >= 400:
        return None, f"HTTP {status}"
    try:
        from astropy.io.votable import parse_single_table
        tab = parse_single_table(io.BytesIO(r.content)).to_table()
    except Exception as exc:
        return None, f"votable: {exc}"
    rows = []
    for row in tab:
        d = {k: (row[k].item() if hasattr(row[k], "item") else row[k]) for k in tab.colnames}
        d["name"] = str(d["name"])
        d["_sep_arcmin"] = sep_deg(float(d["ra"]), float(d["dec"]), ra, de) * 60.0
        if d["_sep_arcmin"] <= r_am:
            rows.append(d)
    return rows, None

def cat_columns(store, table):
    rows, err = heasarc(store, f"SELECT column_name FROM tap_schema.columns WHERE table_name='{table}'", f"cols {table}")
    return set(r["column_name"].strip('"') for r in rows) if rows else set()


def stage_catalogs(store, pos):
    out = {"stage": "catalogs", "utc": now_utc(), "columns": {}, "cones": [], "erosita_ul": None,
           "erosita_scs": [], "csc": []}
    ends = all_endpoints()
    # pin the column lists (docs vs live)
    for table, (label, cols) in HEASARC_CATS.items():
        have = cat_columns(store, table)
        want = [c.strip().strip('"') for c in cols.split(",")]
        use = [c for c in want if c in have]
        out["columns"][table] = {"requested": want, "available": use,
                                 "missing": [c for c in want if c not in have]}
        HEASARC_CATS[table] = (label, ", ".join(('"time"' if c == "time" else c) for c in use))
    targets = []
    for tid, d in ends.items():
        targets.append(("corridor", tid, d["anti"][0], d["anti"][1]))
    for tid in DEEP:
        targets.append(("star", tid, ends[tid]["star"][0], ends[tid]["star"][1]))
    for kind, tid, ra, de in targets:
        for table, (label, cols) in HEASARC_CATS.items():
            r_am = FGL_ARCMIN if table == "fermilpsc" else (CORRIDOR_ARCMIN if kind == "corridor" else STAR_ARCMIN)
            if table == "swbat157m":
                r_am = 12.0
            adql = (f"SELECT {cols} FROM {table} WHERE CONTAINS(POINT('ICRS', ra, dec), "
                    f"CIRCLE('ICRS', {ra:.6f}, {de:.6f}, {r_am / 60.0:.6f}))=1")
            rows, err = heasarc(store, adql, f"{table} {kind} {tid}")
            for r in rows or []:
                r["_sep_arcmin"] = sep_deg(ffloat(r["ra"]), ffloat(r["dec"]), ra, de) * 60.0
            out["cones"].append({"kind": kind, "target": tid, "table": table, "ra": ra, "dec": de,
                                 "radius_arcmin": r_am, "n_rows": None if rows is None else len(rows),
                                 "error": err, "rows": rows or []})
        # CSC 2.1 TAP: the service rejects ADQL geometry (CIRCLE/CONTAINS
        # "not found") and ignores FORMAT=csv -> RA/Dec box + VOTable parse +
        # local separation cut
        r_am = CORRIDOR_ARCMIN if kind == "corridor" else STAR_ARCMIN
        rows, err = csc_box(store, ra, de, r_am, f"csc {kind} {tid}")
        out["csc"].append({"kind": kind, "target": tid, "ra": ra, "dec": de, "radius_arcmin": r_am,
                           "n_rows": None if rows is None else len(rows), "error": err, "rows": rows or []})
        # eRODat cone (DR2 main = eRASS:3; DR1 main = eRASS1)
        for cat in ("DR1_Main", "DR2_Main"):
            status, r = get(store, f"{ERODAT}/catalogue/SCS",
                            {"CAT": cat, "RA": f"{ra:.6f}", "DEC": f"{de:.6f}", "SR": f"{r_am / 60.0:.6f}", "VERB": 1},
                            f"erodat {cat} {kind} {tid}")
            n = None; err = None
            if status == 200:
                txt = r.text
                n = txt.count("<TR>")
                m = re.search(r'name="Error" value="([^"]*)"', txt)
                err = m.group(1) if m else None
            out["erosita_scs"].append({"kind": kind, "target": tid, "cat": cat, "ra": ra, "dec": de,
                                       "radius_arcmin": r_am, "http": status, "n_rows": n, "error": err})
        nz = {t: c["n_rows"] for t, c in [(c["table"], c) for c in out["cones"] if c["target"] == tid and c["kind"] == kind]}
        print(f"[cat] {kind:8s} {tid:11s} " + " ".join(f"{k}={v}" for k, v in nz.items())
              + f" csc={out['csc'][-1]['n_rows']} eDR1={out['erosita_scs'][-2]['n_rows']} eDR2={out['erosita_scs'][-1]['n_rows']}", flush=True)
    # eROSITA upper limits at every corridor + star position, DR1 eRASS1 + DR2 eRASS:3
    q = []
    for kind, tid, ra, de in targets:
        for drs, band in (("DR1_eRASS1", "024"), ("DR1_eRASS1", "021"), ("DR2_eRASSc3", "024")):
            q.append({"ra": ra, "dec": de, "band": band, "dr_survey": drs})
    status, r = get(store, f"{ERODAT}/upperlimit/service_multi", label="erodat UL", method="POST", json_body=q)
    if status == 200:
        try:
            res = r.json()
            lim = res.get("limits", [])
            out["erosita_ul"] = {"http": status, "error": res.get("error"), "n": len(lim),
                                 "rows": [dict(t={"kind": k, "target": tid}, **l) for (k, tid, _, _), l in
                                          zip([t for t in targets for _ in range(3)], lim)]}
        except Exception as exc:
            out["erosita_ul"] = {"http": status, "error": str(exc), "text": r.text[:500]}
    else:
        out["erosita_ul"] = {"http": status, "error": str(getattr(r, "text", r))[:500]}
    (OUT_DIR / "recon_catalogs_v0.json").write_text(json.dumps(out, indent=1))
    return out


# ---------------------------------------------------------------- routes
def stage_routes(store, masters):
    out = {"stage": "routes", "utc": now_utc(), "probes": []}

    def probe(label, url, method="GET", expect=None, timeout=120):
        try:
            r = requests.request(method, url, timeout=timeout, allow_redirects=True, stream=True)
            head = b""
            try:
                for chunk in r.iter_content(4096):
                    head += chunk
                    if len(head) >= 4096:
                        break
            finally:
                r.close()
            store.store(service_url=url, query=method, request_utc=now_utc(), response_bytes=head,
                        row_count=None, http_status=r.status_code)
            rec = {"label": label, "url": url, "http": r.status_code,
                   "content_type": r.headers.get("Content-Type"),
                   "content_length": r.headers.get("Content-Length"),
                   "final_url": r.url, "head": head[:300].decode("utf-8", "replace")}
        except Exception as exc:
            rec = {"label": label, "url": url, "error": str(exc)}
        out["probes"].append(rec)
        print(f"[route] {label}: {rec.get('http', rec.get('error'))} {rec.get('content_length', '')}", flush=True)
        return rec

    # Chandra: pick an in-window obsid if any, else the nearest chanmaster row to ross-128
    chan = [u for u in masters["units"] if u["table"] == "chanmaster"]
    cone = [c for c in masters["cones"] if c["table"] == "chanmaster" and c["n_rows"]]
    obsid = None
    if chan:
        obsid = chan[0]["obsid"]
    elif cone:
        obsid = cone[0]["rows"][0]["obsid"]
    if obsid:
        # cda.harvard.edu/cdaftp 404s; the products are served from the cxc host
        probe("chandra cdaftp primary listing", f"https://cxc.cfa.harvard.edu/cdaftp/byobsid/{str(obsid)[-1]}/{obsid}/primary/")
        probe("chandra cxctap ObsCore", f"https://cda.harvard.edu/cxctap/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&QUERY="
              + urllib.parse.quote(f"SELECT TOP 5 obs_id, t_min, t_max, t_exptime, s_ra, s_dec, access_url FROM ivoa.obscore WHERE obs_id='{obsid}'"))
    # XMM: AIO servlet for one obsid
    xmm = [u for u in masters["units"] if u["table"] == "xmmmaster"]
    conex = [c for c in masters["cones"] if c["table"] == "xmmmaster" and c["n_rows"]]
    xobs = xmm[0]["obsid"] if xmm else (conex[0]["rows"][0]["obsid"] if conex else None)
    if xobs:
        probe("xmm nxsa AIO PPS listing", f"https://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno={xobs}&level=PPS&name=PIEVLI", method="HEAD")
        probe("xmm XSA TAP public obs", f"{XSA_TAP}?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&QUERY="
              + urllib.parse.quote(f"SELECT TOP 5 observation_id, start_utc, end_utc, ra, dec, duration FROM v_public_observations WHERE observation_id='{xobs}'"))
    # Swift: FTP obs tree for one obsid
    sw = [u for u in masters["units"] if u["table"] == "swiftmastr"]
    cones = [c for c in masters["cones"] if c["table"] == "swiftmastr" and c["n_rows"]]
    if sw or cones:
        row = sw[0]["row"] if sw else cones[0]["rows"][0]
        sobs = row["obsid"]; t0 = Time(ffloat(row["start_time"]), format="mjd").datetime
        probe("swift heasarc obs tree", f"https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/{t0:%Y_%m}/{sobs}/")
        probe("swift heasarc xrt event dir", f"https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/{t0:%Y_%m}/{sobs}/xrt/event/")
        probe("swift heasarc bat survey dir", f"https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/{t0:%Y_%m}/{sobs}/bat/survey/")
    probe("swift UKSSDC swifttools API root", "https://www.swift.ac.uk/API/")
    probe("swift LSXPS upper-limit server page", "https://www.swift.ac.uk/LSXPS/ulserv.php")
    # Fermi
    probe("fermi weekly photon HEAD", f"{FERMI_FTP}/photon/lat_photon_weekly_w700_p305_v001.fits", method="HEAD")
    probe("fermi LAT data server form", "https://fermi.gsfc.nasa.gov/cgi-bin/ssc/LAT/LATDataQuery.cgi")
    # eROSITA DR1 skytile + download tree
    r = probe("erosita skytile api (van-maanen antipode)", f"{ERODAT}/skyview/skytile_search_api/?RA=192.3056&DEC=-5.3572&RAD=0")
    tile = None
    try:
        j = requests.get(r["url"], timeout=60).json()
        tile = j["tiles"][0]
        out["erosita_tile_example"] = tile
    except Exception as exc:
        out["erosita_tile_example"] = str(exc)
    if isinstance(tile, dict):
        n = str(tile.get("srvmap"))          # tile RRRDDD -> /DDD/RRR/ (eRODat basket convention)
        probe("erosita DR1 download tree (tile)", f"{ERODAT}/data/download/{n[3:]}/{n[:3]}/")
        probe("erosita DR1 download EXP_010 (tile)", f"{ERODAT}/data/download/{n[3:]}/{n[:3]}/EXP_010/")
        probe("erosita DR1 event list HEAD (tile)", f"{ERODAT}/data/download/{n[3:]}/{n[:3]}/EXP_010/em01_{n}_020_EventList_c010.fits.gz", method="HEAD")
    (OUT_DIR / "recon_routes_v0.json").write_text(json.dumps(out, indent=1))
    return out


def main(stages):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    store = SnapshotStore(RUN)
    events = load_events()
    pos = positions(events)
    for ch in events:
        print(f"[events] {ch}: {len(events[ch])} events b<0.1 AU, "
              f"{len(set(map(str, events[ch]['target_id'])))} targets", flush=True)
    masters = None
    if "masters" in stages:
        masters = stage_masters(store, events, pos)
    if "slew" in stages:
        stage_slew(store, events, pos)
    if "bat" in stages:
        stage_bat(store, events)
    if "lat" in stages:
        stage_lat(store, events)
    if "lat01" in stages:
        stage_lat(store, events, all_01au=True)
    if "catalogs" in stages:
        stage_catalogs(store, pos)
    if "routes" in stages:
        if masters is None:
            masters = json.loads((OUT_DIR / "recon_masters_v0.json").read_text())
        stage_routes(store, masters)


if __name__ == "__main__":
    main(set(sys.argv[1:]) or {"masters", "slew", "bat", "lat", "catalogs", "routes"})
