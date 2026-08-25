"""Stack-regime positive-control fetch for ZTF (notes/learnings.md §10
item 1): cutouts and ring-grid samples of a known numbered asteroid on
nights where its predicted V is FAINTER than the single-exposure limit,
so no frame detects it individually (nothing for the frozen |S_e| <= 5
clip to remove) and only the trajectory-weighted stack can — the regime
the SGL search actually operates in. Counterpart of the catalogue-layer
(60000) control (runs/ztf/v2/control/), which the clip demotes by
construction.

Fetch + sample only; the v2-rule scoring (per archive and joint) is
`surveys/joint/scripts/asteroid_control_stack.py`. Default target
(220000): V in [22.0, 22.6] on ~128 ZTF-era nights (Horizons probe
2026-08-25), ~1-1.5 mag below the zr single-frame 5-sigma limit.

Outputs runs/ztf/v2/control_stack/{cut/, frames.jsonl,
samples_ring.npz, snapshots}: per ok frame the (49, 5, 5) flux/var/
good-frac samples — the 48-offset ring plus the real track, each on the
±2" residual grid — on the ZP 25 scale, the frame weight for the
per-frame cap, and the single-frame S/N at the ephemeris position.

Usage: uv run python surveys/ztf/scripts/asteroid_stack_fetch.py
           [--asteroid 220000] [--vmin 22.0] [--vmax 22.6] [--max-nights 90]
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
from asteroid_fetch import horizons, JD_MJD, CUTOUT_PIX  # noqa: E402

from sglsurvey.adapters.base import ConeRegion, CutoutSpec, MjdRange  # noqa: E402
from sglsurvey.adapters.irsa_ztf import ZtfExactFootprint, ZtfSciAdapter  # noqa: E402
from sglsurvey.photometry import build_flux_map_ztf  # noqa: E402
from sglsurvey.snapshots import SnapshotStore  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
RUN_DIR = REPO / "runs" / "ztf" / "v2" / "control_stack"
CUT_DIR = RUN_DIR / "cut"
RESID_GRID = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])


def sample_ring(fm, ra, dec, scale, offsets):
    """(len(offsets), 5, 5) f/v/g at ra/dec + offsets + residual grid."""
    nr, nt = len(RESID_GRID), len(offsets)
    cosd = np.cos(np.deg2rad(dec))
    off = np.empty((nt, nr, nr, 2))
    off[..., 0] = RESID_GRID[None, :, None] + offsets[:, 0][:, None, None]
    off[..., 1] = RESID_GRID[None, None, :] + offsets[:, 1][:, None, None]
    ra_q = ra + off[..., 0] / 3600.0 / cosd
    dec_q = dec + off[..., 1] / 3600.0
    f, v, g = fm.sample(ra_q.ravel(), dec_q.ravel())
    return (f.reshape(nt, nr, nr) * scale, np.minimum(v.reshape(nt, nr, nr) * scale * scale, 1e30),
            g.reshape(nt, nr, nr))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--asteroid", default="220000")
    ap.add_argument("--vmin", type=float, default=22.0)
    ap.add_argument("--vmax", type=float, default=22.6)
    ap.add_argument("--max-nights", type=int, default=90)
    a = ap.parse_args()
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    store = SnapshotStore(RUN_DIR)
    adapter = ZtfSciAdapter(session=session)
    offsets = np.vstack([[0.0, 0.0], np.array(P.ring)])

    daily, raw, url = horizons(session, a.asteroid, START_TIME="'2018-03-01'",
                               STOP_TIME="'2026-08-01'", STEP_SIZE="'1d'")
    store.store(service_url=url, query="daily", request_utc=datetime.now(timezone.utc).isoformat(),
                response_bytes=raw, row_count=len(daily), http_status=200)
    nights = [d for d in daily if np.isfinite(d["v"]) and a.vmin <= d["v"] <= a.vmax
              and d["elong"] >= 95.0 and d["dec"] > -28.0]
    print(f"{len(daily)} daily epochs, {len(nights)} nights with {a.vmin}<=V<={a.vmax}, elong>=95", flush=True)
    if len(nights) > a.max_nights:
        idx = np.linspace(0, len(nights) - 1, a.max_nights).astype(int)
        nights = [nights[i] for i in idx]

    frames = []
    for d in nights:
        t = datetime.strptime(d["date"][:11], "%Y-%b-%d")
        mjd0 = (t - datetime(1858, 11, 17)).days
        try:
            for obs in adapter.discover(ConeRegion(d["ra"], d["dec"], 0.02), MjdRange(mjd0, mjd0 + 1.0), store):
                if obs.quality_flags.get("bad_quality"):
                    continue
                frames.append(obs)
        except Exception as exc:
            print(f"  discover fail {d['date']}: {exc}", flush=True)
    print(f"{len(frames)} public good-quality frames on those nights", flush=True)
    if not frames:
        return

    frames.sort(key=lambda o: o.t_mid_mjd_utc)
    eph = {}
    for i in range(0, len(frames), 80):
        batch = frames[i:i + 80]
        tlist = " ".join(f"{o.t_mid_mjd_utc + JD_MJD:.6f}" for o in batch)
        rows, raw, url = horizons(session, a.asteroid, TLIST=f"'{tlist}'", TLIST_TYPE="JD")
        store.store(service_url=url, query=f"tlist batch {i}", request_utc=datetime.now(timezone.utc).isoformat(),
                    response_bytes=raw, row_count=len(rows), http_status=200)
        for o, r in zip(batch, rows):
            eph[o.observation_id] = r

    recs, F, V, G, mjds, vpred, bands, wf, snr = [], [], [], [], [], [], [], [], []
    t0 = _time.monotonic()
    for k, o in enumerate(frames):
        e = eph[o.observation_id]
        cut = CutoutSpec(e["ra"], e["dec"], CUTOUT_PIX)
        try:
            sci = adapter.fetch(o, ["sci"], CUT_DIR, cutout=cut)
            msk = adapter.fetch(o, ["msk"], CUT_DIR, cutout=cut)
        except Exception as exc:
            recs.append({"observation_id": o.observation_id, "status": f"fetch fail: {exc}"})
            continue
        try:
            diff = adapter.fetch(o, ["diff"], CUT_DIR, cutout=cut)
            diff_path = diff.products[0].path
        except Exception:
            diff_path = None
        try:
            fm = build_flux_map_ztf(sci.products[0].path, diff_path, msk.products[0].path,
                                    ZtfExactFootprint.FATAL_MASK, o.band, o.t_mid_mjd_utc)
        except Exception as exc:
            recs.append({"observation_id": o.observation_id, "status": f"fluxmap fail: {exc}"})
            continue
        if fm.magzp is None:
            continue
        f0, v0, g0 = fm.sample(e["ra"], e["dec"])
        if not np.isfinite(f0[0]):
            recs.append({"observation_id": o.observation_id, "status": "position outside usable cutout"})
            continue
        scale = 10.0 ** ((P.zp_ref - fm.magzp) / 2.5)
        sel = np.isfinite(fm.var) & (fm.good_frac >= P.min_good_frac) & (fm.var > 0)
        if sel.sum() < 100:
            recs.append({"observation_id": o.observation_id, "status": "too few usable pixels"})
            continue
        Fk, Vk, Gk = sample_ring(fm, e["ra"], e["dec"], scale, offsets)
        F.append(Fk.astype(np.float32)); V.append(Vk.astype(np.float32)); G.append(Gk.astype(np.float32))
        mjds.append(o.t_mid_mjd_utc); vpred.append(e["v"]); bands.append(o.band)
        wf.append(1.0 / (np.median(fm.var[sel]) * scale * scale))
        s1 = float(f0[0] / np.sqrt(v0[0]))
        snr.append(s1)
        recs.append({"observation_id": o.observation_id, "status": "ok", "band": o.band,
                     "mjd": o.t_mid_mjd_utc, "ra": e["ra"], "dec": e["dec"], "v_pred": e["v"],
                     "diff": diff_path is not None, "single_snr": s1, "magzp": float(fm.magzp),
                     "seeing": float(o.quality_flags.get("seeing") or 0.0), "native_key": o.native_key})
        if (k + 1) % 20 == 0:
            print(f"  {k + 1}/{len(frames)} ({(k + 1) / (_time.monotonic() - t0):.2f}/s)", flush=True)
    with (RUN_DIR / "frames.jsonl").open("w") as fh:
        for r in recs:
            fh.write(json.dumps(r) + "\n")
    np.savez_compressed(RUN_DIR / "samples_ring.npz", f=np.array(F), v=np.array(V), g=np.array(G),
                        mjd=np.array(mjds), v_pred=np.array(vpred), band=np.array(bands),
                        frame_w=np.array(wf), single_snr=np.array(snr),
                        resid_grid=RESID_GRID, offsets=offsets,
                        meta=json.dumps({"survey": "ztf", "asteroid": a.asteroid, "zp_ref": P.zp_ref,
                                         "vmin": a.vmin, "vmax": a.vmax}))
    s = np.array(snr)
    print(f"done: {len(mjds)} ok frames; single-frame S/N median {np.median(s):.2f}, "
          f"max {s.max():.2f}, >5 sigma: {(s > 5).sum()}", flush=True)


if __name__ == "__main__":
    main()
