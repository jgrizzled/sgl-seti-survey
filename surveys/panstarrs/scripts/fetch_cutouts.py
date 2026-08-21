"""Fetch image + weight fitscut cutouts for every usable/partial
precise-pass warp (PS1 port of surveys/ztf/scripts/fetch_cutouts.py).
Cutout centre/size come from the precise pass's cutout index (per-warp
size from the locus extent); the full skycell mask from the precise pass
is aligned to the cutout by CRPIX in the photometry builder. Writes
runs/panstarrs/products/cut/manifest.jsonl.

Usage: uv run python surveys/panstarrs/scripts/fetch_cutouts.py [--workers N]
"""

from __future__ import annotations

import argparse
import json
import sys
import time as _time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

from sglsurvey.adapters.base import CutoutSpec
from sglsurvey.adapters.mast_ps1 import Ps1WarpAdapter
from sglsurvey.records import Observation, read_records

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ps1_corridors import CORRIDOR_OF  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
COARSE_DIR = REPO / "runs" / "panstarrs" / "coarse_v1"
PRECISE_DIR = REPO / "runs" / "panstarrs" / "precise_v1"
CUT_DIR = REPO / "runs" / "panstarrs" / "products" / "cut"
KINDS = ("img", "wt")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--corridors", nargs="*", default=None,
                    help="restrict to these corridors (batching)")
    args = ap.parse_args()
    only = set(args.corridors) if args.corridors else None

    obs_by_id = {}
    for r in read_records(COARSE_DIR / "records" / "observation.jsonl"):
        r["corners_icrs_deg"] = tuple(map(tuple, r["corners_icrs_deg"]))
        obs_by_id[r["observation_id"]] = Observation(**r)
    usable = {}
    for r in read_records(PRECISE_DIR / "records"
                          / "intersection_evaluation.jsonl"):
        if r["stage"] == "precise" and r["usable"] in ("usable", "partial"):
            usable.setdefault(r["observation_id"], set()).add(
                CORRIDOR_OF[r["endpoint_id"]])
    cut_index = {}
    with (PRECISE_DIR / "records" / "cutout_index.jsonl").open() as fh:
        for line in fh:
            rec = json.loads(line)
            cut_index[rec["observation_id"]] = rec

    CUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest_path = CUT_DIR / "manifest.jsonl"
    done = set()
    if manifest_path.exists():
        with manifest_path.open() as fh:
            done = {json.loads(l)["observation_id"] for l in fh}
    if only is not None:
        usable = {o: cs for o, cs in usable.items() if cs & only}
    todo = [o for o in usable if o not in done and o in cut_index]
    px = [cut_index[o]["size_pix"] for o in todo]
    est_mb = sum(2 * 4 * p * p for p in px) / 1e6
    print(f"{len(usable)} usable warps, {len(todo)} to fetch "
          f"(~{est_mb:.0f} MB, sizes {min(px) if px else 0}-{max(px) if px else 0} px)",
          flush=True)

    def work(oid):
        obs = obs_by_id[oid]
        ci = cut_index[oid]
        cutout = CutoutSpec(ra_deg=ci["center_ra_deg"],
                            dec_deg=ci["center_dec_deg"],
                            size_pix=ci["size_pix"])
        adapter = Ps1WarpAdapter(session=requests.Session())
        files, missing = {}, []
        for kind in KINDS:
            for attempt in (1, 2, 3):
                try:
                    ps = adapter.fetch(obs, [kind], CUT_DIR, cutout=cutout)
                    files[kind] = ps.products[0].path.name
                    break
                except FileNotFoundError:
                    missing.append(kind)
                    break
                except Exception:
                    if attempt == 3:
                        missing.append(kind)
                    _time.sleep(2.0 * attempt)
        return {"observation_id": oid, "corridor": sorted(usable[oid]),
                "files": files, "missing": missing,
                "msk": Path(ci["path"]).name if ci["path"] else None,
                "header": ci.get("header"),
                "center_ra_deg": cutout.ra_deg,
                "center_dec_deg": cutout.dec_deg,
                "size_pix": cutout.size_pix}

    t0 = _time.monotonic()
    n = 0
    with ThreadPoolExecutor(max_workers=args.workers) as ex, \
            manifest_path.open("a") as fh:
        futs = {ex.submit(work, o): o for o in todo}
        for fut in as_completed(futs):
            rec = fut.result()
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            n += 1
            if n % 50 == 0:
                print(f"  {n}/{len(todo)} ({n / (_time.monotonic() - t0):.2f}/s)",
                      flush=True)
    print("done")


if __name__ == "__main__":
    main()
