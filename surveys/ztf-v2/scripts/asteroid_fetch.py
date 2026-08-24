"""Positive control for the ZTF pilot (hypotheses §9): recover a known
numbered main-belt asteroid through the same cutout -> matched-filter
flux map -> trajectory-sample -> weighted-stack machinery used for the
SGL tracks, with the JPL Horizons ephemeris playing the role of the
trajectory model.

Default target: (60000) 1999 RH27-class object with V ~ 19.1-20.9 at
Palomar during 2018-2026 — at or below the ZTF single-exposure limit,
so single frames are marginal and the stack must do the work.

Outputs runs/ztf/control_v1/{frames.jsonl, samples.npz, summary.json}.

Usage: uv run python surveys/ztf/scripts/asteroid_control.py
           [--asteroid 60000] [--max-nights 80] [--vmax 20.2]
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import time as _time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import requests

from sglsurvey.adapters.base import CutoutSpec, MjdRange
from sglsurvey.adapters.irsa_ztf import ZtfExactFootprint, ZtfSciAdapter
from sglsurvey.adapters.base import ConeRegion
from sglsurvey.photometry import build_flux_map_ztf
from sglsurvey.snapshots import SnapshotStore

REPO = Path(__file__).resolve().parents[3]
RUN_DIR = REPO / "runs" / "ztf-v2" / "control_fetch"
CUT_DIR = RUN_DIR / "cut"
HORIZONS = "https://ssd.jpl.nasa.gov/api/horizons.api"
JD_MJD = 2400000.5
CUTOUT_PIX = 120
ZP_REF = 25.0
# 5x5 grid of +-2" residual offsets mimics the mu grid of the SGL stack.
RESID_GRID = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
OFFSETS = [(0.0, 0.0), (20.0, 0.0), (-20.0, 0.0), (30.0, 0.0),
           (-30.0, 0.0), (40.0, 0.0), (-40.0, 0.0), (0.0, 25.0),
           (0.0, -25.0)]


def horizons(session, command, **kw):
    params = {"format": "text", "COMMAND": f"'{command};'",
              "OBJ_DATA": "NO", "MAKE_EPHEM": "YES",
              "EPHEM_TYPE": "OBSERVER", "CENTER": "'I41'",
              "QUANTITIES": "'1,9,23'", "ANG_FORMAT": "DEG",
              "CSV_FORMAT": "YES", "EXTRA_PREC": "YES"}
    params.update(kw)
    r = session.get(HORIZONS, params=params, timeout=120)
    r.raise_for_status()
    txt = r.text
    body = txt.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines()
    rows = []
    for line in body:
        c = [x.strip() for x in line.split(",")]
        rows.append({"date": c[0], "ra": float(c[3]), "dec": float(c[4]),
                     "v": float(c[5]) if c[5] not in ("n.a.", "") else np.nan,
                     "elong": float(c[7]) if len(c) > 7 and c[7] not in ("n.a.", "") else np.nan})
    return rows, r.content, r.url


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--asteroid", default="60000")
    ap.add_argument("--max-nights", type=int, default=80)
    ap.add_argument("--vmax", type=float, default=20.2)
    args = ap.parse_args()
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    store = SnapshotStore(RUN_DIR)
    adapter = ZtfSciAdapter(session=session)

    # 1. daily ephemeris -> candidate nights
    daily, raw, url = horizons(session, args.asteroid,
                               START_TIME="'2018-03-01'",
                               STOP_TIME="'2026-08-01'", STEP_SIZE="'1d'")
    store.store(service_url=url, query="daily", request_utc=datetime.now(
        timezone.utc).isoformat(), response_bytes=raw, row_count=len(daily),
        http_status=200)
    nights = [d for d in daily if np.isfinite(d["v"]) and d["v"] <= args.vmax
              and d["elong"] >= 100.0 and d["dec"] > -28.0]
    print(f"{len(daily)} daily epochs, {len(nights)} nights with "
          f"V<={args.vmax}, elong>=100", flush=True)
    if len(nights) > args.max_nights:
        idx = np.linspace(0, len(nights) - 1, args.max_nights).astype(int)
        nights = [nights[i] for i in idx]

    # 2. IBE frames per night (public, good quality), containing the
    #    predicted position within a small box
    frames = []
    for d in nights:
        t = datetime.strptime(d["date"][:11], "%Y-%b-%d")
        mjd0 = (t - datetime(1858, 11, 17)).days  # 0h UT
        region = ConeRegion(d["ra"], d["dec"], 0.02)
        rng = MjdRange(mjd0, mjd0 + 1.0)
        try:
            for obs in adapter.discover(region, rng, store):
                if obs.quality_flags.get("bad_quality"):
                    continue
                frames.append(obs)
        except Exception as exc:
            print(f"  discover fail {d['date']}: {exc}", flush=True)
    print(f"{len(frames)} public good-quality frames on those nights",
          flush=True)
    if not frames:
        return

    # 3. exact ephemeris at every frame mid-time (batched TLIST, JD)
    frames.sort(key=lambda o: o.t_mid_mjd_utc)
    eph = {}
    for i in range(0, len(frames), 80):
        batch = frames[i:i + 80]
        tlist = " ".join(f"{o.t_mid_mjd_utc + JD_MJD:.6f}" for o in batch)
        rows, raw, url = horizons(session, args.asteroid,
                                  TLIST=f"'{tlist}'", TLIST_TYPE="JD")
        store.store(service_url=url, query=f"tlist batch {i}",
                    request_utc=datetime.now(timezone.utc).isoformat(),
                    response_bytes=raw, row_count=len(rows), http_status=200)
        for o, r in zip(batch, rows):
            eph[o.observation_id] = r

    # 4. cutouts + flux maps + samples
    nr, nt = len(RESID_GRID), len(OFFSETS)
    recs, F, V, G, mjds, vpred, bands, seeing = [], [], [], [], [], [], [], []
    t0 = _time.monotonic()
    for k, o in enumerate(frames):
        e = eph[o.observation_id]
        cut = CutoutSpec(e["ra"], e["dec"], CUTOUT_PIX)
        try:
            sci = adapter.fetch(o, ["sci"], CUT_DIR, cutout=cut)
            msk = adapter.fetch(o, ["msk"], CUT_DIR, cutout=cut)
        except Exception as exc:
            recs.append({"observation_id": o.observation_id, "status":
                         f"fetch fail: {exc}"})
            continue
        try:
            diff = adapter.fetch(o, ["diff"], CUT_DIR, cutout=cut)
            diff_path = diff.products[0].path
        except Exception:
            diff_path = None
        try:
            fm = build_flux_map_ztf(sci.products[0].path, diff_path,
                                    msk.products[0].path,
                                    ZtfExactFootprint.FATAL_MASK, o.band,
                                    o.t_mid_mjd_utc)
        except Exception as exc:
            recs.append({"observation_id": o.observation_id,
                         "status": f"fluxmap fail: {exc}"})
            continue
        if fm.magzp is None:
            continue
        # containment check: predicted position inside the cutout
        f0, v0, g0 = fm.sample(e["ra"], e["dec"])
        if not np.isfinite(f0[0]):
            recs.append({"observation_id": o.observation_id,
                         "status": "position outside usable cutout"})
            continue
        scale = 10.0 ** ((ZP_REF - fm.magzp) / 2.5)
        cosd = np.cos(np.deg2rad(e["dec"]))
        Fk = np.empty((nt, nr, nr), np.float32)
        Vk = np.empty_like(Fk)
        Gk = np.empty_like(Fk)
        for ti, (dra, ddec) in enumerate(OFFSETS):
            ra_q = (e["ra"] + (RESID_GRID[:, None] + dra) / 3600.0 / cosd
                    + 0 * RESID_GRID[None, :])
            dec_q = (e["dec"] + (RESID_GRID[None, :] + ddec) / 3600.0
                     + 0 * RESID_GRID[:, None])
            f, v, g = fm.sample(ra_q.ravel(), dec_q.ravel())
            Fk[ti] = f.reshape(nr, nr) * scale
            Vk[ti] = v.reshape(nr, nr) * scale * scale
            Gk[ti] = g.reshape(nr, nr)
        F.append(Fk); V.append(Vk); G.append(Gk)
        mjds.append(o.t_mid_mjd_utc); vpred.append(e["v"]); bands.append(o.band)
        seeing.append(o.quality_flags.get("seeing") or 0.0)
        recs.append({"observation_id": o.observation_id, "status": "ok",
                     "band": o.band, "mjd": o.t_mid_mjd_utc,
                     "v_pred": e["v"], "diff": diff_path is not None,
                     "single_snr": float(f0[0] / np.sqrt(v0[0])),
                     "seeing": float(o.quality_flags.get("seeing") or 0.0),
                     "single_mag": (float(fm.magzp - 2.5 * np.log10(f0[0]))
                                    if f0[0] > 0 else None),
                     "magzp": float(fm.magzp),
                     "native_key": o.native_key})
        if (k + 1) % 20 == 0:
            print(f"  {k + 1}/{len(frames)} ({(k + 1) / (_time.monotonic() - t0):.2f}/s)",
                  flush=True)
    with (RUN_DIR / "frames.jsonl").open("w") as fh:
        for r in recs:
            fh.write(json.dumps(r) + "\n")
    F, V, G = np.array(F), np.array(V), np.array(G)
    mjds, vpred, seeing = np.array(mjds), np.array(vpred), np.array(seeing)
    bands = np.array(bands)
    np.savez_compressed(RUN_DIR / "samples.npz", f=F, v=V, g=G, mjd=mjds,
                        v_pred=vpred, band=bands, seeing=seeing,
                        resid_grid=RESID_GRID, offsets=np.array(OFFSETS))

    # 5. stack exactly as the calibration does (weights, cap 20x median)
    summary = {"asteroid": args.asteroid, "n_frames_ok": int(len(mjds))}
    for band in sorted(set(bands.tolist())):
        sel = bands == band
        f, v, g = F[sel], V[sel], G[sel]
        valid = np.isfinite(f) & np.isfinite(v) & (v > 0) & (g >= 0.7)
        w = np.where(valid, 1.0 / np.where(v > 0, v, 1.0), 0.0)
        wmed = np.nanmedian(np.where(w > 0, w, np.nan), axis=0)
        w = np.minimum(w, 20.0 * np.nan_to_num(wmed, nan=np.inf))
        A = (np.where(valid, f, 0.0) * w).sum(axis=0)
        B = w.sum(axis=0)
        S = A / np.sqrt(np.where(B > 0, B, np.inf))
        T = float(np.nanmax(S[1:]))
        i_real = np.unravel_index(np.nanargmax(S[0]), S[0].shape)
        flux = A[0][i_real] / B[0][i_real]
        mag = ZP_REF - 2.5 * np.log10(flux) if flux > 0 else None
        # expected stack mag: flux-weighted mean of predicted V (AB r~V+0.2 ignored)
        fpred = 10 ** (-0.4 * vpred[sel])
        single = np.array([r["single_snr"] for r in recs
                           if r.get("status") == "ok" and r["band"] == band])
        summary[band] = {
            "n_epochs": int(sel.sum()),
            "stack_S_real": float(S[0][i_real]),
            "stack_S_center": float(S[0][2, 2]),
            "threshold_8_controls": T,
            "control_maxima": np.round(np.nanmax(S[1:], axis=(1, 2)), 2).tolist(),
            "peak_offset_arcsec": [float(RESID_GRID[i_real[0]]),
                                   float(RESID_GRID[i_real[1]])],
            "recovered_mag_zp_scale": mag,
            "predicted_V_mean_flux": float(-2.5 * np.log10(fpred.mean())),
            "predicted_V_median": float(np.median(vpred[sel])),
            "single_frame_snr_median": float(np.median(single)),
            "single_frames_above_5sigma": int((single > 5).sum()),
            "detected": bool(S[0][i_real] > T),
        }
        print(f"[{band}] N={sel.sum()} S_real={S[0][i_real]:.1f} "
              f"(center {S[0][2,2]:.1f}) T={T:.2f} peak offset "
              f"{summary[band]['peak_offset_arcsec']}\" recovered mag "
              f"{mag if mag is None else round(mag, 2)} vs predicted V "
              f"{summary[band]['predicted_V_mean_flux']:.2f}; single-frame "
              f"median S/N {np.median(single):.1f}, "
              f"{(single > 5).sum()} frames >5 sigma", flush=True)
    (RUN_DIR / "summary.json").write_text(json.dumps(summary, indent=2,
                                                     default=float))


if __name__ == "__main__":
    main()
