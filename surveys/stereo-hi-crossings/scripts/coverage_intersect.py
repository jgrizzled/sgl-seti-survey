"""STEREO-A HI-1 coverage intersect (pre-freeze stage).

Inventory of the HI-1A synoptic archive against the 0.1 AU sunward
events of `crossings/stereoa_v1`:

1. day lists: RAL level-1 tree (CGI, the search substrate) and the
   NASA SSC level-2 mirror (plain HTTP; identical file stems);
2. per-day frame lists from the SSC mirror (snapshotted HTML), giving
   every frame's nominal start time from its file stem;
3. per-day pointing: the FITS header of the first frame of each day
   via an HTTP Range request on the SSC mirror (CRVAL1/2 = HPLN/HPLT
   of the CCD centre, DATE-AVG, EXPTIME, N_IMAGES, NMISSING) — the
   roll/pointing history over the mission;
4. per event: frames whose start time falls in the in-beam arc, with
   the source's (HPLN, HPLT) tested against that day's footprint
   (centre from the header, half-widths from the recon frame), and
   frames in the star-fixed baseline transit (HPLN in [-12, -6] deg);
5. L1/L2 parity spot-check on a seeded sample of days (RAL CGI).

Outputs: results/coverage_v1.json; snapshots under
runs/stereo-hi-crossings/coverage/.
"""

from __future__ import annotations

import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.time import Time

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hi_geometry as G

REPO = Path(__file__).resolve().parents[3]
EVENTS = REPO / "crossings" / "stereoa_v1" / "events.ecsv"
CENSUS = REPO / "surveys" / "stereo-hi-crossings" / "results" / "stereoa_census_v1.json"
OUT = REPO / "surveys" / "stereo-hi-crossings" / "results" / "coverage_v1.json"
RUNS = REPO / "runs" / "stereo-hi-crossings" / "coverage"
LIST = RUNS / "listings"
HEAD = RUNS / "day_headers"

SSC = "https://stereo-ssc.nascom.nasa.gov/data/ins_data/secchi_hi/L2/a/img/hi_1/"
RAL_CGI = "https://www.stereo.rl.ac.uk/cgi-bin/data.py"
RAL_L1 = "lz/L1/a/img/hi_1"
HALF_W_HPLN = (G.HI1_HPLN[1] - G.HI1_HPLN[0]) / 2   # from the recon frame
HALF_W_HPLT = (G.HI1_HPLT[1] - G.HI1_HPLT[0]) / 2
ARC_ABS_HPLN = (4.05, 6.0)        # in-beam + in-FOV arc, either side of the Sun
BASE_ABS_HPLN = (6.0, 12.0)       # star-fixed baseline transit, same side
STEM = re.compile(r'href="((\d{8})_(\d{6})_24h1A_br01\.fts)"')


def http_get(url: str, headers: dict | None = None, timeout: int = 60,
             tries: int = 4) -> bytes | None:
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers=headers or {})
            return urllib.request.urlopen(req, timeout=timeout).read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            time.sleep(2 + 2 * k)
        except Exception:
            time.sleep(2 + 2 * k)
    return None


def ral_post(target: str, timeout: int = 120) -> bytes | None:
    data = urllib.parse.urlencode({"source": "SECCHI", "target": target}).encode()
    for k in range(3):
        try:
            return urllib.request.urlopen(urllib.request.Request(RAL_CGI, data=data),
                                          timeout=timeout).read()
        except Exception:
            time.sleep(3 + 3 * k)
    return None


def day_lists() -> tuple[list[str], list[str]]:
    LIST.mkdir(parents=True, exist_ok=True)
    p = LIST / "ssc_days.html"
    if not p.exists():
        p.write_bytes(http_get(SSC, timeout=120) or b"")
    ssc = sorted(set(re.findall(r'href="(\d{8})/"', p.read_text())))
    q = LIST / "ral_l1_days.html"
    if not q.exists():
        q.write_bytes(ral_post(RAL_L1) or b"")
    ral = sorted(set(re.findall(rf'{RAL_L1}/(\d{{8}})"', q.read_text())))
    return ssc, ral


def day_frames(day: str) -> list[tuple[str, float]]:
    """[(stem, mjd_start)] for one day from the SSC listing snapshot."""
    p = LIST / f"ssc_{day}.html"
    if not p.exists():
        raw = http_get(SSC + day + "/")
        p.write_bytes(raw if raw is not None else b"<!--404-->")
    out = []
    for fname, d, t in STEM.findall(p.read_text(errors="replace")):
        iso = f"{d[:4]}-{d[4:6]}-{d[6:]}T{t[:2]}:{t[2:4]}:{t[4:]}"
        out.append((fname, Time(iso, scale="utc").mjd))
    return out


def day_header(day: str, frames: list[tuple[str, float]]) -> dict | None:
    """Header keywords of the first frame of the day (Range request)."""
    HEAD.mkdir(parents=True, exist_ok=True)
    p = HEAD / f"{day}.json"
    if p.exists():
        return json.loads(p.read_text())
    if not frames:
        return None
    raw = http_get(SSC + day + "/" + frames[0][0], headers={"Range": "bytes=0-17279"})
    if raw is None:
        return None
    txt = raw[:17280].decode("latin-1")
    cards = {}
    for i in range(0, len(txt) - 79, 80):
        c = txt[i:i + 80]
        if c.startswith("END "):
            break
        if "=" in c[:10]:
            key = c[:8].strip()
            val = c[10:].split("/")[0].strip().strip("'").strip()
            cards[key] = val
    want = ["DATE-AVG", "EXPTIME", "N_IMAGES", "NMISSING", "CRVAL1", "CRVAL2",
            "CRVAL1A", "CRVAL2A", "CDELT1", "PC1_1A", "PC1_2A", "SUMMED", "NAXIS1",
            "NAXIS2", "ATT_FILE", "VERSION", "IMGSEQ", "BIASMEAN", "RAVG"]
    rec = {k: cards.get(k) for k in want}
    rec["frame"] = frames[0][0]
    p.write_text(json.dumps(rec) + "\n")
    return rec


def parity_check(days: list[str], n: int = 40, seed: int = 20260904) -> dict:
    rng = np.random.default_rng(seed)
    sample = sorted(rng.choice(days, size=min(n, len(days)), replace=False))
    rows = []
    for d in sample:
        raw = ral_post(f"{RAL_L1}/{d}")
        l1 = set(re.findall(r'(\d{8}_\d{6})_14h1A\.fts', (raw or b"").decode(errors="replace")))
        l2 = {s[:15] for s, _ in day_frames(d)}
        rows.append({"day": d, "n_l1": len(l1), "n_l2": len(l2),
                     "l1_not_l2": len(l1 - l2), "l2_not_l1": len(l2 - l1)})
    return {"sampled_days": len(rows), "rows": rows,
            "total_l1": sum(r["n_l1"] for r in rows), "total_l2": sum(r["n_l2"] for r in rows),
            "l1_only": sum(r["l1_not_l2"] for r in rows), "l2_only": sum(r["l2_not_l1"] for r in rows)}


def main() -> None:
    ssc_days, ral_days = day_lists()
    print(f"SSC L2 days: {len(ssc_days)} ({ssc_days[0]}..{ssc_days[-1]}); "
          f"RAL L1 days: {len(ral_days)} ({ral_days[0]}..{ral_days[-1]})")
    with ThreadPoolExecutor(6) as ex:
        frames_by_day = dict(zip(ssc_days, ex.map(day_frames, ssc_days)))
    print(f"frames listed: {sum(len(v) for v in frames_by_day.values())}")
    with ThreadPoolExecutor(6) as ex:
        hdr_by_day = dict(zip(ssc_days, ex.map(lambda d: day_header(d, frames_by_day[d]), ssc_days)))
    n_hdr = sum(1 for v in hdr_by_day.values() if v)
    print(f"day headers: {n_hdr}")

    # per-day pointing (HPLN/HPLT of the CCD centre)
    pointing = {}
    for d, h in hdr_by_day.items():
        if h and h.get("CRVAL1") not in (None, ""):
            try:
                pointing[d] = (float(h["CRVAL1"]), float(h["CRVAL2"]))
            except ValueError:
                pass
    cx = np.array([v[0] for v in pointing.values()])
    cy = np.array([v[1] for v in pointing.values()])
    nominal = (np.abs(np.abs(cx) - 14.0) < 0.5) & (np.abs(cy) < 2.5)
    off_days = sorted(d for d, v in pointing.items()
                      if not (abs(abs(v[0]) - 14.0) < 0.5 and abs(v[1]) < 2.5))
    east_days = int((cx < 0).sum())
    # pointing transitions (sign of the centre HPLN)
    trans, prev = [], None
    for d in sorted(pointing):
        side = "east" if pointing[d][0] < 0 else "west"
        if side != prev:
            trans.append({"day": d, "side": side, "centre": pointing[d]})
            prev = side

    # all frame times, day-indexed
    all_frames = sorted((m, d, s) for d, v in frames_by_day.items() for s, m in v)
    fm = np.array([x[0] for x in all_frames])
    fday = [x[1] for x in all_frames]

    cen = json.loads(CENSUS.read_text())
    t = Table.read(EVENTS)
    by_id = {str(r["event_id"]): r for r in t}
    rows = []
    for r in cen["visibility_0.1AU"]["rows"]:
        ev = by_id[r["event_id"]]
        e = {"star_icrs_ra_deg": float(ev["star_icrs_ra_deg"]),
             "star_icrs_dec_deg": float(ev["star_icrs_dec_deg"])}
        s_hat = G.source_direction(e, r["channel"])
        mjd_ca = r["mjd_ca"]
        lo, hi = mjd_ca - 13.0, mjd_ca + 13.0
        i0, i1 = np.searchsorted(fm, lo), np.searchsorted(fm, hi)
        n_arc = n_base = n_arc_nominal = 0
        arc_days, sides = set(), set()
        if i1 > i0:
            mj = fm[i0:i1]
            hpln, hplt = G.helioprojective(s_hat, mj)
            b_e = G.impact_parameter_au(e, mj)
            for k in range(i1 - i0):
                d = fday[i0 + k]
                c = pointing.get(d)
                if c is None:
                    continue
                inside = (abs(hpln[k] - c[0]) <= HALF_W_HPLN - G.HI1_EDGE_MARGIN_DEG
                          and abs(hplt[k] - c[1]) <= HALF_W_HPLT - G.HI1_EDGE_MARGIN_DEG)
                if not inside:
                    continue
                if b_e[k] <= 0.1 and ARC_ABS_HPLN[0] <= abs(hpln[k]) <= ARC_ABS_HPLN[1]:
                    n_arc += 1
                    arc_days.add(d)
                    sides.add("pre_tca" if mj[k] < mjd_ca else "post_tca")
                    if abs(abs(c[0]) - 14.0) < 0.5 and abs(c[1]) < 2.5:
                        n_arc_nominal += 1
                elif BASE_ABS_HPLN[0] < abs(hpln[k]) <= BASE_ABS_HPLN[1]:
                    n_base += 1
        rows.append({**{k: r[k] for k in ("event_id", "target_id", "channel", "t_ca_utc",
                                          "b_min_rsun", "visible_hours")},
                     "n_frames_arc": n_arc, "n_frames_arc_nominal_pointing": n_arc_nominal,
                     "n_frames_baseline": n_base, "arc_days": sorted(arc_days),
                     "side": "/".join(sorted(sides)) if sides else None,
                     "covered": n_arc >= 10})

    agg = {}
    for r in rows:
        a = agg.setdefault(r["channel"], {}).setdefault(r["target_id"], {
            "events": 0, "covered": 0, "arc_frames": [], "base_frames": []})
        a["events"] += 1
        a["covered"] += r["covered"]
        a["arc_frames"].append(r["n_frames_arc"])
        a["base_frames"].append(r["n_frames_baseline"])
    for ch in agg:
        for tid, a in agg[ch].items():
            a["median_arc_frames"] = float(np.median(a.pop("arc_frames")))
            a["median_baseline_frames"] = float(np.median(a.pop("base_frames")))

    # era gaps: runs of >= 3 consecutive calendar days with no frames
    days_with = sorted(d for d, v in frames_by_day.items() if v)
    mjd_days = np.array([Time(f"{d[:4]}-{d[4:6]}-{d[6:]}").mjd for d in days_with])
    gaps = []
    for a, b in zip(mjd_days[:-1], mjd_days[1:]):
        if b - a >= 4:
            gaps.append({"from": Time(a + 1, format="mjd").iso[:10],
                         "to": Time(b - 1, format="mjd").iso[:10], "days": int(b - a - 1)})
    parity = parity_check(ral_days)

    per_year = {}
    for m, d, _ in all_frames:
        per_year[d[:4]] = per_year.get(d[:4], 0) + 1
    out = {
        "stage": "coverage_v1",
        "events_input": str(EVENTS.relative_to(REPO)),
        "archive": {"ssc_l2_days": len(ssc_days), "ssc_first_day": ssc_days[0],
                    "ssc_last_day": ssc_days[-1], "ral_l1_days": len(ral_days),
                    "ral_first_day": ral_days[0], "ral_last_day": ral_days[-1],
                    "n_frames_listed": int(len(all_frames)), "frames_per_year": per_year,
                    "day_headers_read": n_hdr},
        "pointing": {"n_days": len(pointing), "nominal_days": int(nominal.sum()),
                     "off_nominal_days": len(off_days), "east_days": east_days,
                     "west_days": int(len(pointing) - east_days),
                     "abs_hpln_centre_median": float(np.median(np.abs(cx))),
                     "hplt_centre_median": float(np.median(cy)),
                     "hplt_centre_range": [float(cy.min()), float(cy.max())],
                     "transitions": trans,
                     "off_nominal_day_list": off_days[:400]},
        "gaps_ge_3d": gaps,
        "l1_l2_parity": parity,
        "arc_definition": {"abs_hpln_deg": ARC_ABS_HPLN, "baseline_abs_hpln_deg": BASE_ABS_HPLN,
                           "edge_margin_deg": G.HI1_EDGE_MARGIN_DEG,
                           "footprint_half_widths_deg": [HALF_W_HPLN, HALF_W_HPLT]},
        "aggregates": agg, "rows": rows,
    }
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: out[k] for k in ("archive", "pointing", "gaps_ge_3d")}, indent=1)[:4000])
    print(json.dumps(parity, indent=1)[:600])
    for ch in agg:
        for tid, a in agg[ch].items():
            print(ch, tid, a)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
