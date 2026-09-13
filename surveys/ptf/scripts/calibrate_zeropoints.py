"""Per-frame star-calibrated zero points through the IDENTICAL matched
filter (PS1 lesson 1; hypotheses §4): for every fetched stage-2 cutout,
build the flux map and measure the corridor's PS1 DR2 calibrators
(ptf_calib rules); zp_star = median(m_pred + 2.5 log10 F_mf) under the
>= 5-star / <= 0.2-mag scatter gate. Frames failing the gate get
zp_star = null and are unusable (never header-MAGZPT fallback).

Writes runs/ptf/zeropoints.jsonl (one row per observation id).

Usage: uv run python surveys/ptf/scripts/calibrate_zeropoints.py
       [--corridors c1 c2 ...] [--limit N] [--workers N]
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time as _time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from sglsurvey.adapters.irsa_ptf import MASK_FATAL_TEMPLATE
from sglsurvey.photometry import build_flux_map_ptf
from sglsurvey.records import read_records

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ptf_calib import corridor_calibrators, star_zeropoint  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
RUNS = REPO / "runs" / "ptf"
CUT = RUNS / "products" / "cut"
MSK = RUNS / "products" / "msk"
OUT = RUNS / "zeropoints.jsonl"


def one(args):
    m, band, mjd = args
    rec = {"observation_id": m["observation_id"], "band": band,
           "zp_star": None, "n_stars": 0, "zp_mad": None,
           "n_in_frame": 0, "magzpt_header": None}
    try:
        fm = build_flux_map_ptf(CUT / m["files"]["sci"], MSK / m["msk"],
                                MASK_FATAL_TEMPLATE, band, mjd)
    except Exception as exc:
        rec["error"] = str(exc)[:120]
        return rec
    rec["magzpt_header"] = fm.magzp
    rec["fwhm_arcsec"] = round(float(fm.fwhm_arcsec), 3)
    rec["bg_sigma"] = round(float(fm.bg_sigma), 3)
    cors = m["corridor"] if isinstance(m["corridor"], list) else [m["corridor"]]
    stars = np.vstack([corridor_calibrators(c, band) for c in cors]) \
        if cors else np.empty((0, 3))
    if len(stars):
        stars = np.unique(stars, axis=0)
    zp, n, mad, ninb = star_zeropoint(fm, stars)
    rec.update(zp_star=(round(zp, 4) if zp is not None else None),
               n_stars=n, zp_mad=(round(mad, 4) if mad is not None else None),
               n_in_frame=ninb)
    return rec


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--corridors", nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--workers", type=int, default=6)
    args = ap.parse_args()
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    only = set(args.corridors) if args.corridors else None
    obs = {r["observation_id"]: r for r in read_records(
        RUNS / "coarse_v1" / "records" / "observation.jsonl")}
    rows = []
    for line in open(CUT / "manifest.jsonl"):
        m = json.loads(line)
        if "sci" in m.get("files", {}) and m.get("msk"):
            rows.append(m)
    if only:
        rows = [m for m in rows if set(m["corridor"]) & only]
    done = set()
    if OUT.exists():
        done = {json.loads(l)["observation_id"] for l in open(OUT)}
    todo = [m for m in rows if m["observation_id"] not in done]
    if args.limit:
        todo = todo[:args.limit]
    print(f"{len(rows)} exposures, {len(todo)} to calibrate", flush=True)
    jobs = [(m, obs[m["observation_id"]]["band"],
             obs[m["observation_id"]]["t_mid_mjd_utc"]) for m in todo]
    t0 = _time.monotonic()
    with OUT.open("a") as fh, ProcessPoolExecutor(args.workers) as ex:
        for i, rec in enumerate(ex.map(one, jobs, chunksize=4)):
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            if (i + 1) % 50 == 0:
                print(f"  {i + 1}/{len(todo)} "
                      f"({(i + 1) / (_time.monotonic() - t0):.2f}/s)",
                      flush=True)
    recs = [json.loads(l) for l in open(OUT)]
    n_ok = sum(1 for r in recs if r.get("zp_star") is not None)
    print(f"done: {n_ok} calibrated of {len(recs)}")


if __name__ == "__main__":
    main()
