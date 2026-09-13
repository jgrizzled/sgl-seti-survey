"""Find a multi-epoch numbered-asteroid positive control inside the
PTF corridor exposures themselves (the ZTF/DECam (60000) control has
only 10 PTF frames, 4 surviving the ZP gate — PTF's coverage of any
one asteroid's track is too sparse). SkyBoT cone searches at every
usable exposure of the near-ecliptic pilot corridors, numbered
asteroids with V in [vmin, vmax] and ephemeris error <= 1", counted
per object across frames. Writes runs/ptf/control_v1/recon.json with
the ranked candidates and the chosen one's frame list.

Usage: uv run python surveys/ptf/scripts/asteroid_control_recon.py
           [--corridors vanmaanen ross128] [--vmin 18.5] [--vmax 20.5]
"""

from __future__ import annotations

import argparse
import json
import sys
import threading
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import requests
from astropy.coordinates import Angle

from sglsurvey.adapters.irsa_ptf import PtfNominalFootprint
from sglsurvey.records import Observation, read_records
from sglsurvey.snapshots import SnapshotStore

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ptf_corridors import CORRIDOR_OF  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
RUNS = REPO / "runs" / "ptf"
WORK = RUNS / "control_v1"
SKYBOT = "https://vo.imcce.fr/webservices/skybot/skybotconesearch_query.php"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--corridors", nargs="*", default=["vanmaanen", "ross128"])
    ap.add_argument("--vmin", type=float, default=18.5)
    ap.add_argument("--vmax", type=float, default=20.5)
    a = ap.parse_args()
    obs = {}
    for r in read_records(RUNS / "coarse_v1" / "records" / "observation.jsonl"):
        r["corners_icrs_deg"] = tuple(map(tuple, r["corners_icrs_deg"]))
        obs[r["observation_id"]] = Observation(**r)
    use = set()
    for r in read_records(RUNS / "precise_v1" / "records" / "intersection_evaluation.jsonl"):
        if r["usable"] in ("usable", "partial") and CORRIDOR_OF[r["endpoint_id"]] in a.corridors:
            use.add(r["observation_id"])
    frames = sorted(use, key=lambda o: obs[o].t_mid_mjd_utc)
    print(f"{len(frames)} usable exposures in {a.corridors}", flush=True)
    WORK.mkdir(parents=True, exist_ok=True)
    store = SnapshotStore(WORK)
    lock = threading.Lock()

    def one(oid):
        o = obs[oid]
        w = o.wcs
        params = {"EPOCH": f"{o.t_mid_mjd_utc + 2400000.5:.6f}", "RA": f"{w['crval1']:.5f}",
                  "DEC": f"{w['crval2']:.5f}", "SR": "0.7", "-mime": "text", "-loc": "675",
                  "-filter": "120"}
        try:
            resp = requests.get(SKYBOT, params=params, timeout=120)
        except Exception as exc:
            return oid, []
        with lock:
            store.store(service_url=SKYBOT, query="&".join(f"{k}={v}" for k, v in params.items()),
                        request_utc=datetime.now(timezone.utc).isoformat(),
                        response_bytes=resp.content, row_count=None, http_status=resp.status_code)
        rows = []
        for line in resp.text.splitlines():
            if line.startswith(("#", "-")) or "|" not in line:
                continue
            f = [x.strip() for x in line.split("|")]
            try:
                num = int(f[0]); ra = Angle(f[2] + " hours").degree
                dec = Angle(f[3] + " degrees").degree; mv = float(f[5]); err = float(f[6])
            except (ValueError, IndexError):
                continue
            if err > 1.0 or not (a.vmin <= mv <= a.vmax):
                continue
            if PtfNominalFootprint(o, pad_arcsec=-60.0).contains(ra, dec):
                rows.append({"num": num, "name": f[1], "ra": ra, "dec": dec, "mv": mv})
        return oid, rows

    seen = defaultdict(list)
    with ThreadPoolExecutor(6) as ex:
        for k, (oid, rows) in enumerate(ex.map(one, frames)):
            for r in rows:
                seen[r["num"]].append({"observation_id": oid, "band": obs[oid].band,
                                       "mjd": obs[oid].t_mid_mjd_utc, "name": r["name"],
                                       "ra": r["ra"], "dec": r["dec"], "mv": r["mv"]})
            if (k + 1) % 50 == 0:
                print(f"  {k + 1}/{len(frames)}", flush=True)
    ranked = sorted(seen.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    out = {"corridors": a.corridors, "v_range": [a.vmin, a.vmax], "n_frames_queried": len(frames),
           "candidates": [{"num": n, "name": h[0]["name"], "n_frames": len(h),
                           "n_R": sum(1 for x in h if x["band"] == "R"),
                           "n_g": sum(1 for x in h if x["band"] == "g"),
                           "n_nights": len({round(x["mjd"]) for x in h}),
                           "mv": [round(x["mv"], 2) for x in h]} for n, h in ranked[:15]],
           "frames": {str(n): h for n, h in ranked[:5]}}
    (WORK / (f"recon_v{a.vmin}-{a.vmax}.json" if a.vmin > 20.5 else "recon.json")).write_text(json.dumps(out, indent=1))
    for c in out["candidates"][:10]:
        print(c)


if __name__ == "__main__":
    main()
