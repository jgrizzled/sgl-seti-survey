"""Fetch -int and -unc cutouts covering the corridor arc for every
usable precise-pass frame-band (stage-2 input, plan §4.6).

Cutout geometry per frame: center and size are computed from the union
of the member endpoint x role locus arcs at the frame epoch, with
margin for the PSF, the frozen residual-motion bound (1 "/yr over the
mission), and the search padding. One int+unc cutout pair per unique
frame-band serves all member endpoint-roles of its corridor.

Idempotent: existing files are kept; a manifest records geometry.

Usage: uv run python surveys/wise/scripts/fetch_cutouts.py
"""

from __future__ import annotations

import json
import time as _time
from collections import defaultdict
from pathlib import Path

import numpy as np
import requests
from astropy.time import Time

from sglseti import Role, load_target_registry

from sglsurvey.geometry import GeometryContext, locus_radec
from sglsurvey.records import Observation, read_records

REPO = Path(__file__).resolve().parents[3]
COARSE_DIR = REPO / "runs" / "wise" / "coarse_v1"
PRECISE_DIR = REPO / "runs" / "wise" / "precise_v1"
CUT_DIR = REPO / "runs" / "wise" / "products" / "cut"
REGISTRY_PATH = REPO / "registries" / "pilot_wise_2026.yaml"

CORRIDOR_OF = {
    "barnard-star": "barnard", "ross-154": "ross154",
    "lalande-21185": "lalande", "alpha-cen-a": "alphacen",
    "alpha-cen-b": "alphacen", "sirius-a": "sirius",
    "sirius-b": "sirius",
}
MEMBERS = defaultdict(list)
for e, c in CORRIDOR_OF.items():
    MEMBERS[c].append(e)

MARGIN_ARCSEC = 45.0   # PSF (12") + |mu|*14yr (14") + padding (10") + buffer
PIX_ARCSEC = 2.75
MAX_SIZE_PIX = 240
LOCUS_BIN_DAYS = 0.5


def main() -> None:
    registry = load_target_registry(REGISTRY_PATH)
    ctx = GeometryContext.wise_coarse_default()
    session = requests.Session()
    CUT_DIR.mkdir(parents=True, exist_ok=True)

    obs_by_id = {}
    for r in read_records(COARSE_DIR / "records" / "observation.jsonl"):
        r["corners_icrs_deg"] = tuple(map(tuple, r["corners_icrs_deg"]))
        obs_by_id[r["observation_id"]] = Observation(**r)

    # Frame-bands with any usable endpoint-role, plus which corridor.
    frame_corridor = {}
    for r in read_records(PRECISE_DIR / "records"
                          / "intersection_evaluation.jsonl"):
        if r["stage"] == "precise" and r["usable"] in ("usable", "partial"):
            frame_corridor[r["observation_id"]] = \
                CORRIDOR_OF[r["endpoint_id"]]

    manifest_path = CUT_DIR / "manifest.jsonl"
    done = set()
    if manifest_path.exists():
        done = {json.loads(l)["observation_id"]
                for l in open(manifest_path)}
    todo = [oid for oid in sorted(
        frame_corridor, key=lambda o: obs_by_id[o].t_mid_mjd_utc)
        if oid not in done]
    print(f"{len(frame_corridor)} usable frame-bands, {len(todo)} to fetch",
          flush=True)

    locus_cache: dict = {}

    def arc_points(corridor, mjd):
        key = (corridor, round(mjd / LOCUS_BIN_DAYS))
        if key not in locus_cache:
            t = Time(key[1] * LOCUS_BIN_DAYS, format="mjd")
            pts = []
            for e in MEMBERS[corridor]:
                for role in (Role.RX, Role.TX):
                    p, _ = locus_radec(ctx, registry[e], role, t,
                                       tolerance_arcsec=2.0)
                    pts.append(p)
            locus_cache[key] = np.concatenate(pts)
        return locus_cache[key]

    t0 = _time.monotonic()
    n_err = 0
    with open(manifest_path, "a") as mf:
        for i, oid in enumerate(todo):
            obs = obs_by_id[oid]
            pts = arc_points(frame_corridor[oid], obs.t_mid_mjd_utc)
            cra = float(np.mean(pts[:, 0]))
            cdec = float(np.mean(pts[:, 1]))
            cosd = np.cos(np.deg2rad(cdec))
            half_deg = float(np.max(np.hypot(
                (pts[:, 0] - cra) * cosd, pts[:, 1] - cdec)))
            size_pix = min(MAX_SIZE_PIX, int(np.ceil(
                2 * (half_deg * 3600 + MARGIN_ARCSEC) / PIX_ARCSEC)))
            row = {"observation_id": oid,
                   "corridor": frame_corridor[oid],
                   "center_radec": [cra, cdec], "size_pix": size_pix,
                   "files": {}}
            ok = True
            for kind in ("int", "unc"):
                url = obs.products[kind]["url"]
                name = f"cut-{url.rsplit('/', 1)[-1]}"
                path = CUT_DIR / name
                if not path.exists():
                    for attempt in (1, 2, 3):
                        try:
                            resp = session.get(
                                url, params={
                                    "center": f"{cra},{cdec}",
                                    "size": f"{size_pix}pix"},
                                timeout=120)
                            resp.raise_for_status()
                            path.write_bytes(resp.content)
                            break
                        except Exception as exc:
                            if attempt == 3:
                                print(f"  FAIL {oid} {kind}: {exc}",
                                      flush=True)
                                ok = False
                            else:
                                _time.sleep(2.0 * attempt)
                if ok:
                    row["files"][kind] = name
            if ok:
                mf.write(json.dumps(row) + "\n")
                mf.flush()
            else:
                n_err += 1
            if (i + 1) % 250 == 0:
                rate = (i + 1) / (_time.monotonic() - t0)
                print(f"  {i + 1}/{len(todo)} ({rate:.1f}/s, "
                      f"{n_err} errors)", flush=True)
    print(f"done: {len(todo)} attempted, {n_err} errors", flush=True)


if __name__ == "__main__":
    main()
