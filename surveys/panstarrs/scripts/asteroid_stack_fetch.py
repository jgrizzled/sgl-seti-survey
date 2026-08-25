"""Stack-regime positive-control fetch for PS1 (notes/learnings.md §10
item 1): warp cutouts and ring-grid samples of a known numbered
asteroid on nights where its predicted V is fainter than the
single-warp limit — nothing for the frozen |S_e| <= 5 clip to remove,
only the stack can recover it. Port of the retired v1 control fetch
(git history: surveys/panstarrs/scripts/asteroid_control.py) with the
faint-window selection and the 48-ring sample grid.

Flux scale = header FPA.ZP + the per-filter median (ZP_star - FPA.ZP)
offset from the v1 star calibration (runs/panstarrs/calib_v1/
zeropoints.jsonl) — the per-warp star calibration used by the survey
tensors needs a corridor star catalogue that does not exist at the
asteroid's positions; the transfer offset is accurate to ~0.05 mag and
declared in the summary.

Fetch + sample only; scoring is
`surveys/joint/scripts/asteroid_control_stack.py`. Default target
(220000): V in [22.0, 22.6] on ~67 PS1-era nights (Horizons probe
2026-08-25).

Outputs runs/panstarrs/v2/control_stack/{cut/, frames.jsonl,
samples_ring.npz, snapshots}.

Usage: uv run python surveys/panstarrs/scripts/asteroid_stack_fetch.py
           [--asteroid 220000] [--vmin 22.0] [--vmax 22.6] [--max-nights 120]
"""

from __future__ import annotations

import argparse
import json
import sys
import time as _time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from profile import PROFILE as P  # noqa: E402

from sglsurvey.adapters.base import CutoutSpec  # noqa: E402
from sglsurvey.adapters.mast_ps1 import Ps1ExactFootprint, Ps1WarpAdapter  # noqa: E402
from sglsurvey.photometry import build_flux_map_ps1  # noqa: E402
from sglsurvey.snapshots import SnapshotStore  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
RUN_DIR = REPO / "runs" / "panstarrs" / "v2" / "control_stack"
CUT_DIR = RUN_DIR / "cut"
MSK_DIR = REPO / "runs" / "panstarrs" / "products" / "msk"
ZP_LOG = REPO / "runs" / "panstarrs" / "calib_v1" / "zeropoints.jsonl"
HORIZONS = "https://ssd.jpl.nasa.gov/api/horizons.api"
JD_MJD = 2400000.5
CUTOUT_PIX = 240  # 60" at 0.25"/pix
RESID_GRID = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])


def horizons(session, command, **kw):
    params = {"format": "text", "COMMAND": f"'{command};'", "OBJ_DATA": "NO", "MAKE_EPHEM": "YES",
              "EPHEM_TYPE": "OBSERVER", "CENTER": "'F51'", "QUANTITIES": "'1,9,23'",
              "ANG_FORMAT": "DEG", "CSV_FORMAT": "YES", "EXTRA_PREC": "YES"}
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


def sample_ring(fm, ra, dec, scale, offsets):
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
    ap.add_argument("--max-nights", type=int, default=120)
    a = ap.parse_args()
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    store = SnapshotStore(RUN_DIR)
    adapter = Ps1WarpAdapter(session=session)
    offsets = np.vstack([[0.0, 0.0], np.array(P.ring)])

    acc = defaultdict(list)
    for line in ZP_LOG.read_text().splitlines():
        r = json.loads(line)
        if r["n_cal"] >= 5 and r["zp_hdr"] is not None:
            acc[r["band"]].append(r["zp_star"] - r["zp_hdr"])
    zp_off = {b: float(np.median(v)) for b, v in acc.items()}
    print(f"ZP_star - FPA.ZP offsets by filter: {zp_off}", flush=True)

    daily, raw, url = horizons(session, a.asteroid, START_TIME="'2009-05-01'",
                               STOP_TIME="'2014-04-30'", STEP_SIZE="'1d'")
    store.store(service_url=url, query="daily", request_utc=datetime.now(timezone.utc).isoformat(),
                response_bytes=raw, row_count=len(daily), http_status=200)
    nights = [d for d in daily if np.isfinite(d["v"]) and a.vmin <= d["v"] <= a.vmax
              and d["elong"] >= 90.0 and d["dec"] > -30.0]
    print(f"{len(daily)} daily epochs, {len(nights)} nights with {a.vmin}<=V<={a.vmax}, elong>=90", flush=True)
    if len(nights) > a.max_nights:
        idx = np.linspace(0, len(nights) - 1, a.max_nights).astype(int)
        nights = [nights[i] for i in idx]

    frames, seen = [], set()
    for d in nights:
        t = datetime.strptime(d["date"][:11], "%Y-%b-%d")
        mjd0 = (t - datetime(1858, 11, 17)).days
        try:
            for row in adapter.list_warps(d["ra"], d["dec"], store):
                mjd = float(row["mjd"])
                if not (mjd0 - 0.5 <= mjd < mjd0 + 1.5):
                    continue
                if row["filename"] in seen:
                    continue
                seen.add(row["filename"])
                frames.append(adapter._observation(row))
        except Exception as exc:
            print(f"  listing fail {d['date']}: {exc}", flush=True)
    print(f"{len(frames)} warps on those nights", flush=True)
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
            img = adapter.fetch(o, ["img"], CUT_DIR, cutout=cut)
            msk = adapter.fetch(o, ["msk"], MSK_DIR)
        except Exception as exc:
            recs.append({"observation_id": o.observation_id, "status": f"fetch fail: {exc}"})
            continue
        try:
            wt = adapter.fetch(o, ["wt"], CUT_DIR, cutout=cut)
            wt_path = wt.products[0].path
        except Exception:
            wt_path = None
        try:
            fm = build_flux_map_ps1(img.products[0].path, wt_path, msk.products[0].path,
                                    Ps1ExactFootprint.FATAL_MASK, o.band, o.t_mid_mjd_utc)
        except Exception as exc:
            recs.append({"observation_id": o.observation_id, "status": f"fluxmap fail: {exc}"})
            continue
        if fm.magzp is None:
            continue
        fm.magzp = fm.magzp + zp_off.get(o.band, 0.0)
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
                     "single_snr": s1, "magzp": float(fm.magzp), "native_key": o.native_key})
        if (k + 1) % 20 == 0:
            print(f"  {k + 1}/{len(frames)} ({(k + 1) / (_time.monotonic() - t0):.2f}/s)", flush=True)
    with (RUN_DIR / "frames.jsonl").open("w") as fh:
        for r in recs:
            fh.write(json.dumps(r) + "\n")
    np.savez_compressed(RUN_DIR / "samples_ring.npz", f=np.array(F), v=np.array(V), g=np.array(G),
                        mjd=np.array(mjds), v_pred=np.array(vpred), band=np.array(bands),
                        frame_w=np.array(wf), single_snr=np.array(snr),
                        resid_grid=RESID_GRID, offsets=offsets,
                        meta=json.dumps({"survey": "ps1", "asteroid": a.asteroid, "zp_ref": P.zp_ref,
                                         "vmin": a.vmin, "vmax": a.vmax, "zp_offsets": zp_off}))
    s = np.array(snr)
    if len(s):
        print(f"done: {len(mjds)} ok frames; single-frame S/N median {np.median(s):.2f}, "
              f"max {s.max():.2f}, >5 sigma: {(s > 5).sum()}", flush=True)


if __name__ == "__main__":
    main()
