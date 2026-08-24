"""Fetch single-CCD image + weight-map HDUs for every usable/partial
precise-pass exposure (DECam port of surveys/panstarrs/scripts/
fetch_cutouts.py).

CCD selection: the retrieval unit is the CCD (``?hdus=``,
EXTNAME-asserted), and the CCDs to fetch are those the *usable locus
crossing* actually touches — recomputed here from the precise-pass
hits, located on the exposure's dqmask. Selecting by the locus-arc
centre instead loses exposures whose centre sits in a chip gap (51 of
213 pilot exposures) and clips multi-CCD arcs. The precise pass's full
dqmask supplies the per-exposure HDU map and is re-fetched only if
purged. Writes runs/decam/products/cut/manifest.jsonl
(one row per exposure; files keyed "<kind>:<EXTNAME>").

Usage: uv run python surveys/decam/scripts/fetch_cutouts.py [--workers N]
"""

from __future__ import annotations

import argparse
import json
import sys
import time as _time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import warnings

import numpy as np
import requests
from astropy.time import Time

# astropy's warning logger mutates shared state and crashes under
# concurrent FITS opens in ThreadPoolExecutor workers; suppressing the
# warnings avoids that code path entirely (decam recon addendum).
warnings.filterwarnings("ignore")
try:
    from astropy.logger import log as _astropy_log
    _astropy_log.disabled = True
except Exception:
    pass

from sglseti import Role, adaptive_locus, load_target_registry

from sglsurvey.adapters.noirlab_decam import (DecamExactFootprint,
                                              DecamInstcalAdapter)
from sglsurvey.geometry import GeometryContext
from sglsurvey.records import Observation, read_records

sys.path.insert(0, str(Path(__file__).resolve().parent))
from decam_corridors import CORRIDOR_OF  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
COARSE_DIR = REPO / "runs" / "decam" / "coarse_v1"
PRECISE_DIR = REPO / "runs" / "decam" / "precise_v1"
DQ_DIR = REPO / "runs" / "decam" / "products" / "dqmask"
CUT_DIR = REPO / "runs" / "decam" / "products" / "cut"
REGISTRY_PATH = REPO / "registries" / "pilot_wise_2026.yaml"
KINDS = ("image", "wtmap")
LOCUS_TOLERANCE_ARCSEC = 1.0
ROLE_BY_NAME = {r.value: r for r in Role}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--corridors", nargs="*", default=None,
                    help="restrict to these corridors (batching)")
    args = ap.parse_args()
    only = set(args.corridors) if args.corridors else None

    registry = load_target_registry(REGISTRY_PATH)
    ctx = GeometryContext.decam_default()

    obs_by_id = {}
    for r in read_records(COARSE_DIR / "records" / "observation.jsonl"):
        r["corners_icrs_deg"] = tuple(map(tuple, r["corners_icrs_deg"]))
        obs_by_id[r["observation_id"]] = Observation(**r)
    usable_pairs = defaultdict(set)   # obs_id -> {(endpoint, role)}
    for r in read_records(PRECISE_DIR / "records"
                          / "intersection_evaluation.jsonl"):
        if r["stage"] == "precise" and r["usable"] in ("usable", "partial"):
            usable_pairs[r["observation_id"]].add(
                (r["endpoint_id"], r["role"]))
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
            for line in fh:
                m = json.loads(line)
                if m.get("schema") != "v1":
                    continue  # v0 rows (centre-CCD selection) redone
                dq_ok = (not m.get("dqmask")) or (
                    DQ_DIR / m["dqmask"]).exists()
                if dq_ok and m["files"] and not m["missing"] and all(
                        (CUT_DIR / f).exists() for f in m["files"].values()):
                    done.add(m["observation_id"])
    if only is not None:
        usable_pairs = {
            o: ps for o, ps in usable_pairs.items()
            if {CORRIDOR_OF[e] for e, _ in ps} & only}
    todo = [o for o in usable_pairs if o not in done and o in cut_index]
    print(f"{len(usable_pairs)} usable exposures, {len(todo)} to fetch",
          flush=True)

    def work(oid):
        obs = obs_by_id[oid]
        ci = cut_index[oid]
        adapter = DecamInstcalAdapter(session=requests.Session())
        files, missing, ccds = {}, [], []
        # dqmask (re-fetch if purged) + covered-locus CCD selection
        try:
            dq = adapter.fetch(obs, ["dqmask"], DQ_DIR)
        except Exception as exc:
            return {"schema": "v1", "observation_id": oid,
                    "corridor": sorted({CORRIDOR_OF[e]
                                        for e, _ in usable_pairs[oid]}),
                    "files": {}, "missing": ["dqmask"], "ccds": [],
                    "error": str(exc)[:200], "dqmask": None,
                    "header": ci.get("header")}
        t_mid = Time(obs.t_mid_mjd_utc, format="mjd")
        with DecamExactFootprint(dq.products[0].path) as fp:
            hit_ccds, onccd_ccds = set(), set()
            for endpoint_id, role_name in sorted(usable_pairs[oid]):
                al = adaptive_locus(
                    target=registry[endpoint_id],
                    role=ROLE_BY_NAME[role_name],
                    observation_time=t_mid, observer=ctx.observer,
                    relay_range=ctx.relay_range,
                    tolerance_arcsec=LOCUS_TOLERANCE_ARCSEC,
                    ephemeris=ctx.ephemeris, model=ctx.model)
                pts = np.array([[p.icrs_ra_deg, p.icrs_dec_deg]
                                for p in al.points])
                for name, _x, _y, us in fp.locate(pts):
                    if name is not None:
                        onccd_ccds.add(name)
                        if us:
                            hit_ccds.add(name)
        # marginal partial cells (usable_fraction 0) graze only masked
        # pixels at the 1" locus sampling — fall back to any on-CCD
        # crossing so stage 2 sees whatever pixels exist
        ccds = sorted(hit_ccds or onccd_ccds)
        if not ccds:
            missing.append("no-usable-ccd")
        for name in ccds:
            for kind in KINDS:
                for attempt in (1, 2, 3):
                    try:
                        ps = adapter.fetch(obs, [kind], CUT_DIR,
                                           extname=name,
                                           dqmask_dir=DQ_DIR)
                        files[f"{kind}:{name}"] = \
                            ps.products[0].path.name
                        break
                    except FileNotFoundError:
                        missing.append(f"{kind}:{name}")
                        break
                    except Exception:
                        if attempt == 3:
                            missing.append(f"{kind}:{name}")
                        _time.sleep(2.0 * attempt)
        return {"schema": "v1", "observation_id": oid,
                "corridor": sorted({CORRIDOR_OF[e]
                                    for e, _ in usable_pairs[oid]}),
                "files": files, "missing": missing, "ccds": ccds,
                "dqmask": dq.products[0].path.name,
                "header": ci.get("header"),
                "center_ra_deg": ci.get("center_ra_deg"),
                "center_dec_deg": ci.get("center_dec_deg"),
                "size_pix": ci.get("size_pix")}

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
            if n % 25 == 0:
                print(f"  {n}/{len(todo)} "
                      f"({n / (_time.monotonic() - t0):.2f}/s)", flush=True)
    print("done")


if __name__ == "__main__":
    main()
