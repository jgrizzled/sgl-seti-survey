"""PTF v2 positive control (hypotheses §6): recover a known numbered
main-belt asteroid through the full chain — IBE discovery at the
Horizons position, sci + dmask cutouts, the PTF matched filter, the
in-frame PS1 DR2 star calibration (ptf_calib rules) and the v2 ring-48
rule (sglsurvey.control.rescore). Predicted magnitudes = Horizons V +
solar colours (g = V + 0.25; R_AB = V - 0.36 + 0.21 = V - 0.15).

Two passes: --stage fetch discovers frames and writes the frozen
exposure list configs/asteroid_control_v1.json plus the cutouts under
runs/ptf/control_v1/; --stage score re-scores that list with the
profile's rule (frozen clip and no-clip). Output:
runs/ptf/v2/control/<asteroid>{,_noclip}/summary.json.

Usage: uv run python surveys/ptf/scripts/asteroid_control_v2.py
           --stage fetch [--asteroid 60000] [--vmax 20.6] [--max-nights 80]
       uv run python surveys/ptf/scripts/asteroid_control_v2.py --stage score
"""

from __future__ import annotations

import argparse
import json
import sys
import time as _time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from profile import PROFILE as P  # noqa: E402
from ptf_calib import calibrators_from_rows, star_zeropoint  # noqa: E402

from sglsurvey.adapters.base import ConeRegion, CutoutSpec, MjdRange  # noqa: E402
from sglsurvey.adapters.irsa_ptf import (MASK_FATAL_TEMPLATE,  # noqa: E402
                                         PtfLevel1Adapter)
from sglsurvey.adapters.mast_ps1 import MEAN_COLS, catalog_cone  # noqa: E402
from sglsurvey.photometry import build_flux_map_ptf  # noqa: E402
from sglsurvey.records import Observation, read_records  # noqa: E402
from sglsurvey.snapshots import SnapshotStore  # noqa: E402
from sglsurvey.control import rescore  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
CFG = REPO / "surveys" / "ptf" / "configs" / "asteroid_control_v1.json"
WORK = REPO / "runs" / "ptf" / "control_v1"
CUT_DIR = WORK / "cut"
HORIZONS = "https://ssd.jpl.nasa.gov/api/horizons.api"
JD_MJD = 2400000.5
CUTOUT_PIX = 384     # calibrator-density floor measured at the crossings dev
CAL_CONE_DEG = 0.12
COLOUR = {"g": 0.25, "R": -0.15}


def horizons(session, command, **kw):
    params = {"format": "text", "COMMAND": f"'{command};'", "OBJ_DATA": "NO",
              "MAKE_EPHEM": "YES", "EPHEM_TYPE": "OBSERVER", "CENTER": "'675'",
              "QUANTITIES": "'1,9,23'", "ANG_FORMAT": "DEG", "CSV_FORMAT": "YES",
              "EXTRA_PREC": "YES"}
    params.update(kw)
    r = session.get(HORIZONS, params=params, timeout=120)
    r.raise_for_status()
    body = r.text.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines()
    rows = []
    for line in body:
        c = [x.strip() for x in line.split(",")]
        rows.append({"date": c[0], "ra": float(c[3]), "dec": float(c[4]),
                     "v": float(c[5]) if c[5] not in ("n.a.", "") else np.nan,
                     "elong": float(c[7]) if len(c) > 7 and c[7] not in ("n.a.", "") else np.nan})
    return rows, r.content, r.url


class _LockedStore:
    """SnapshotStore wrapper serialising the record append across threads."""

    def __init__(self, store, lock):
        self._s, self._l = store, lock

    def store(self, **kw):
        with self._l:
            return self._s.store(**kw)


def fetch(args):
    WORK.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    store = SnapshotStore(WORK)
    adapter = PtfLevel1Adapter(session=session)
    if args.recon:
        # frames from asteroid_control_recon.py (SkyBoT over the pilot
        # corridor exposures): the coarse Observation records by id
        rec = json.loads((WORK / args.recon_file).read_text())["frames"][str(args.asteroid)]
        by_id = {}
        for r in read_records(REPO / "runs" / "ptf" / "coarse_v1" / "records" / "observation.jsonl"):
            r["corners_icrs_deg"] = tuple(map(tuple, r["corners_icrs_deg"]))
            by_id[r["observation_id"]] = r
        frames = [Observation(**by_id[x["observation_id"]]) for x in rec]
        print(f"{len(frames)} recon frames for ({args.asteroid})", flush=True)
        _fetch_frames(args, session, store, adapter, frames)
        return
    daily, raw, url = horizons(session, args.asteroid, START_TIME="'2009-03-01'",
                               STOP_TIME="'2015-02-01'", STEP_SIZE="'1d'")
    store.store(service_url=url, query="daily", request_utc=datetime.now(timezone.utc).isoformat(),
                response_bytes=raw, row_count=len(daily), http_status=200)
    nights = [d for d in daily if np.isfinite(d["v"]) and d["v"] <= args.vmax
              and d["elong"] >= 100.0 and d["dec"] > -28.0]
    print(f"{len(daily)} daily epochs, {len(nights)} nights with V<={args.vmax}, elong>=100", flush=True)
    if len(nights) > args.max_nights:
        idx = np.linspace(0, len(nights) - 1, args.max_nights).astype(int)
        nights = [nights[i] for i in idx]
    from concurrent.futures import ThreadPoolExecutor
    import threading
    lock = threading.Lock()

    def one(d):
        t = datetime.strptime(d["date"][:11], "%Y-%b-%d")
        mjd0 = (t - datetime(1858, 11, 17)).days
        ad = PtfLevel1Adapter(session=requests.Session())
        try:
            got = list(ad.discover(ConeRegion(d["ra"], d["dec"], 0.02),
                                   MjdRange(mjd0, mjd0 + 1.0), _LockedStore(store, lock)))
        except Exception as exc:
            print(f"  discover fail {d['date']}: {exc}", flush=True)
            got = []
        return got

    frames = []
    with ThreadPoolExecutor(8) as ex:
        for got in ex.map(one, nights):
            frames.extend(got)
    print(f"{len(frames)} PTF frames on those nights", flush=True)
    if not frames:
        return
    _fetch_frames(args, session, store, adapter, frames)


def _fetch_frames(args, session, store, adapter, frames):
    frames.sort(key=lambda o: o.t_mid_mjd_utc)
    if len(frames) > args.max_frames:
        idx = np.linspace(0, len(frames) - 1, args.max_frames).astype(int)
        frames = [frames[i] for i in idx]
    eph = {}
    for i in range(0, len(frames), 80):
        batch = frames[i:i + 80]
        tlist = " ".join(f"{o.t_mid_mjd_utc + JD_MJD:.6f}" for o in batch)
        rows, raw, url = horizons(session, args.asteroid, TLIST=f"'{tlist}'", TLIST_TYPE="JD")
        store.store(service_url=url, query=f"tlist batch {i}",
                    request_utc=datetime.now(timezone.utc).isoformat(),
                    response_bytes=raw, row_count=len(rows), http_status=200)
        for o, r in zip(batch, rows):
            eph[o.observation_id] = r
    # PS1 calibrator cones, one per 0.1-deg cell along the track
    cal_cells = {}
    exposures = []
    for k, o in enumerate(frames):
        e = eph[o.observation_id]
        cut = CutoutSpec(e["ra"], e["dec"], CUTOUT_PIX)
        try:
            sci = adapter.fetch(o, ["sci"], CUT_DIR, cutout=cut)
            msk = adapter.fetch(o, ["msk"], CUT_DIR, cutout=cut)
        except Exception as exc:
            print(f"  fetch fail {o.observation_id}: {exc}", flush=True)
            continue
        cell = (round(e["ra"] / 0.1), round(e["dec"] / 0.1))
        if cell not in cal_cells:
            rows, snaps = catalog_cone("mean", e["ra"], e["dec"], CAL_CONE_DEG, MEAN_COLS, store)
            name = f"cal_{cell[0]}_{cell[1]}.json"
            (WORK / name).write_text(json.dumps({"center": [e["ra"], e["dec"]],
                                                 "radius_deg": CAL_CONE_DEG,
                                                 "snapshots": snaps, "rows": rows}))
            cal_cells[cell] = name
        exposures.append({
            "observation_id": o.observation_id, "native_key": o.native_key, "band": o.band,
            "mjd_mid": o.t_mid_mjd_utc, "exptime_s": o.exptime_s, "products": o.products,
            "quality_flags": o.quality_flags, "ra_pred": e["ra"], "dec_pred": e["dec"],
            "v_pred": e["v"], "sci": sci.products[0].path.name, "msk": msk.products[0].path.name,
            "calibrators": cal_cells[cell]})
        if (k + 1) % 20 == 0:
            print(f"  {k + 1}/{len(frames)}", flush=True)
    CFG.parent.mkdir(parents=True, exist_ok=True)
    cfg_path = CFG if args.asteroid == "798452" or not args.recon else CFG.with_name(f"asteroid_control_{args.asteroid}.json")
    cfg_path.write_text(json.dumps({
        "asteroid": args.asteroid, "site": "675 (Palomar)", "vmax": args.vmax,
        "cutout_pix": CUTOUT_PIX, "colours": COLOUR, "n_exposures": len(exposures),
        "selected_utc": datetime.now(timezone.utc).isoformat(), "exposures": exposures}, indent=1))
    print(f"frozen exposure list: {len(exposures)} exposures -> {cfg_path}")


def score(args):
    cfg_path = CFG if args.asteroid == "798452" or not args.recon else CFG.with_name(f"asteroid_control_{args.asteroid}.json")
    cfg = json.loads(cfg_path.read_text())
    frames, maps, positions, vpred = [], [], [], {b: [] for b in P.bands}
    zps = []
    cal_cache = {}
    for e in cfg["exposures"]:
        band = e["band"]
        try:
            fm = build_flux_map_ptf(CUT_DIR / e["sci"], CUT_DIR / e["msk"], MASK_FATAL_TEMPLATE,
                                    band, e["mjd_mid"])
        except Exception as exc:
            print(f"  fluxmap fail {e['observation_id']}: {exc}", flush=True)
            continue
        key = (e["calibrators"], band)
        if key not in cal_cache:
            rows = json.loads((WORK / e["calibrators"]).read_text())["rows"]
            cal_cache[key] = np.array(calibrators_from_rows(rows, band)).reshape(-1, 3)
        stars = cal_cache[key].copy()
        if band == "R" and len(stars):
            stars[:, 2] += 0.21   # Cousins R Vega -> AB, as in the survey profile
        zp, n, mad, _ = star_zeropoint(fm, stars)
        if zp is None:
            continue
        fm.magzp = zp
        f0, v0, g0 = fm.sample(e["ra_pred"], e["dec_pred"])
        if not np.isfinite(f0[0]):
            continue
        zps.append({"observation_id": e["observation_id"], "band": band, "zp_star": zp,
                    "n_stars": n, "zp_mad": mad,
                    "single_snr": float(f0[0] / np.sqrt(v0[0])) if v0[0] > 0 else None})
        frames.append({"band": band, "mjd": e["mjd_mid"], "v_pred": e["v_pred"]})
        maps.append(fm)
        positions.append((e["ra_pred"], e["dec_pred"]))
        vpred[band].append(e["v_pred"])
    pred = {b: float(-2.5 * np.log10(np.mean(10 ** (-0.4 * np.array(v)))) + COLOUR[b])
            for b, v in vpred.items() if v}
    extra = {"source": "configs/asteroid_control_v1.json (IBE + Horizons site 675, snapshotted)",
             "colours": COLOUR, "n_calibrated": len(frames),
             "zp_scatter_median": float(np.median([z["zp_mad"] for z in zps])) if zps else None,
             "single_snr_median": float(np.nanmedian([z["single_snr"] for z in zps if z["single_snr"] is not None])) if zps else None,
             "single_frames_above_5sigma": int(sum(1 for z in zps if (z["single_snr"] or 0) > 5))}
    (WORK / "frames_scored.jsonl").write_text("\n".join(json.dumps(z) for z in zps) + "\n")
    label = cfg["asteroid"]
    rescore(P, frames, maps, positions, P.run_dir / "control" / label, label, predicted_mag=pred, extra=extra)
    rescore(P, frames, maps, positions, P.run_dir / "control" / f"{label}_noclip",
            f"{label} (no single-epoch clip)", predicted_mag=pred,
            extra={**extra, "note": "the frozen clip removes frames where a bright mover is individually detected"},
            clip_sigma=None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["fetch", "score"], required=True)
    ap.add_argument("--asteroid", default="60000")
    ap.add_argument("--vmax", type=float, default=20.6)
    ap.add_argument("--max-frames", type=int, default=120)
    ap.add_argument("--max-nights", type=int, default=250)
    ap.add_argument("--recon-file", default="recon.json")
    ap.add_argument("--recon", action="store_true",
                    help="take frames from runs/ptf/control_v1/recon.json")
    args = ap.parse_args()
    fetch(args) if args.stage == "fetch" else score(args)


if __name__ == "__main__":
    main()
