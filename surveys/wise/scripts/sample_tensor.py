"""Stage-7 sampling pass: store the full per-epoch trajectory sample
tensor for every endpoint x role, for the real trajectory plus 8
offset-control trajectories (plan §4.7 matched controls).

Per pair, writes runs/wise/calib_v1/tensors/<endpoint>__<role>.npz:
  mjd (E,), band_idx (E,) [1..4], phase (E,) [0/1],
  f (9, E, 64, 5, 5) float32   matched-filter flux at common ZP 20
  v (9, E, 64, 5, 5) float32   variance at common ZP
  g (9, E, 64, 5, 5) float16   PSF good-pixel weight fraction
Trajectory 0 is the real (z, mu) grid; 1..8 are the declared offsets.

All downstream calibration (nulls, thresholds, injections, constraint
curves) runs from these arrays without touching FITS again.

Usage: uv run python surveys/wise/scripts/sample_tensor.py [--only-missing]
  --only-missing  build tensors only for endpoint x role pairs without
                  an existing .npz (incremental batches); existing
                  tensors are left untouched.
"""

from __future__ import annotations

import json
import time as _time
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.time import Time

from sglseti import Role, adaptive_locus, load_target_registry

from sglsurvey.adapters.irsa_wise import WiseExactFootprint
from sglsurvey.geometry import GeometryContext
from sglsurvey.photometry import build_flux_map
from sglsurvey.records import read_records

REPO = Path(__file__).resolve().parents[3]
COARSE_DIR = REPO / "runs" / "wise" / "coarse_v1"
PRECISE_DIR = REPO / "runs" / "wise" / "precise_v1"
CUT_DIR = REPO / "runs" / "wise" / "products" / "cut"
MSK_DIR = REPO / "runs" / "wise" / "products" / "msk"
OUT_DIR = REPO / "runs" / "wise" / "calib_v1" / "tensors"
REGISTRY_PATH = REPO / "registries" / "pilot_wise_2026.yaml"

# Uniform in q = 1/z: the locus position scales as parallax ~ 1/z, so
# this gives constant ~5.6" on-sky spacing between adjacent nodes —
# just under the W1 PSF FWHM — where a log grid leaves ~17" gaps at
# low z and oversamples high z.
Z_GRID = 1.0 / np.linspace(1.0 / 550.0, 1.0 / 10000.0, 64)
MU_GRID = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
T0_MJD = 57800.0
ZP_REF = 20.0
LOCUS_BIN_DAYS = 0.5
# (dRA, dDec) arcsec control-trajectory offsets; 0 = real.
OFFSETS = [(0.0, 0.0), (20.0, 0.0), (-20.0, 0.0), (30.0, 0.0),
           (-30.0, 0.0), (40.0, 0.0), (-40.0, 0.0), (0.0, 25.0),
           (0.0, -25.0)]

from wise_corridors import CORRIDOR_OF

BAND_IDX = {"W1": 1, "W2": 2, "W3": 3, "W4": 4}


def main() -> None:
    import sys
    only_missing = "--only-missing" in sys.argv[1:]
    registry = load_target_registry(REGISTRY_PATH)
    ctx = GeometryContext.wise_coarse_default()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    obs_by_id = {r["observation_id"]: r for r in read_records(
        COARSE_DIR / "records" / "observation.jsonl")}
    usable = defaultdict(list)
    for r in read_records(PRECISE_DIR / "records"
                          / "intersection_evaluation.jsonl"):
        if r["stage"] == "precise" and r["usable"] in ("usable", "partial"):
            usable[(r["endpoint_id"], r["role"])].append(
                r["observation_id"])
    if only_missing:
        usable = {k: v for k, v in usable.items()
                  if not (OUT_DIR / f"{k[0]}__{k[1]}.npz").exists()}
        print(f"--only-missing: {len(usable)} endpoint-role pairs to build",
              flush=True)
    manifest = {json.loads(l)["observation_id"]: json.loads(l)
                for l in open(CUT_DIR / "manifest.jsonl")}

    corridor_frames = defaultdict(set)
    for (e, role), oids in usable.items():
        corridor_frames[CORRIDOR_OF[e]].update(oids)
    phase_ref = {}
    for c, oids in corridor_frames.items():
        doys = np.array([obs_by_id[o]["t_mid_mjd_utc"] % 365.25
                         for o in sorted(oids)])
        phase_ref[c] = float(np.median(doys[:50]))

    def phase_of(corridor, mjd):
        d = (mjd % 365.25) - phase_ref[corridor]
        d = (d + 182.625) % 365.25 - 182.625
        return 0 if abs(d) < 91.3 else 1

    zcache: dict = {}

    def points_at_zgrid(endpoint, role, mjd):
        key = (endpoint, role, round(mjd / LOCUS_BIN_DAYS))
        if key not in zcache:
            al = adaptive_locus(
                target=registry[endpoint], role=Role(role),
                observation_time=Time(key[2] * LOCUS_BIN_DAYS,
                                      format="mjd"),
                observer=ctx.observer, relay_range=ctx.relay_range,
                tolerance_arcsec=2.0, ephemeris=ctx.ephemeris,
                model=ctx.model)
            zs = np.array([p.z_au for p in al.points])
            ra = np.array([p.icrs_ra_deg for p in al.points])
            dec = np.array([p.icrs_dec_deg for p in al.points])
            q = 1.0 / zs
            o = np.argsort(q)
            qg = 1.0 / Z_GRID
            zcache[key] = np.stack([np.interp(qg, q[o], ra[o]),
                                    np.interp(qg, q[o], dec[o])], axis=1)
        return zcache[key]

    pairs_of_frame = defaultdict(list)
    for p, oids in usable.items():
        for o in oids:
            pairs_of_frame[o].append(p)

    rows = defaultdict(list)   # pair -> [(mjd, band, phase, f, v, g)]
    frames = sorted(pairs_of_frame,
                    key=lambda o: obs_by_id[o]["t_mid_mjd_utc"])
    nz, nm, nt = len(Z_GRID), len(MU_GRID), len(OFFSETS)
    t0 = _time.monotonic()
    n = 0
    for oid in frames:
        row = manifest.get(oid)
        obs = obs_by_id[oid]
        if row is None or "int" not in row.get("files", {}):
            continue
        msk_name = obs["products"]["msk"]["url"].rsplit("/", 1)[-1]
        try:
            fm = build_flux_map(
                CUT_DIR / row["files"]["int"],
                CUT_DIR / row["files"]["unc"],
                MSK_DIR / msk_name, obs["band"], obs["t_mid_mjd_utc"],
                WiseExactFootprint.FATAL_MASK,
                magzp=obs["quality_flags"].get("magzp"))
        except Exception as exc:
            print(f"  fluxmap FAIL {oid}: {exc}", flush=True)
            continue
        n += 1
        mjd = obs["t_mid_mjd_utc"]
        zp = fm.magzp if fm.magzp else ZP_REF
        scale = 10.0 ** ((ZP_REF - zp) / 2.5)
        dt_yr = (mjd - T0_MJD) / 365.25
        ph = phase_of(row["corridor"], mjd)
        for pair in pairs_of_frame[oid]:
            endpoint, role = pair
            base = points_at_zgrid(endpoint, role, mjd)
            cosd = np.cos(np.deg2rad(base[:, 1]))
            dmu = MU_GRID * dt_yr / 3600.0
            ra_all = np.broadcast_to(
                base[:, 0][:, None, None]
                + dmu[None, :, None] / cosd[:, None, None],
                (nz, nm, nm))
            dec_all = np.broadcast_to(
                base[:, 1][:, None, None] + dmu[None, None, :],
                (nz, nm, nm))
            F = np.empty((nt, nz, nm, nm), dtype=np.float32)
            V = np.empty_like(F)
            G = np.empty((nt, nz, nm, nm), dtype=np.float16)
            for ti, (dra, ddec) in enumerate(OFFSETS):
                ra_q = ra_all + dra / 3600.0 / cosd[:, None, None]
                dec_q = dec_all + ddec / 3600.0
                f, v, g = fm.sample(ra_q.ravel(), dec_q.ravel())
                F[ti] = (f.reshape(nz, nm, nm) * scale).astype(np.float32)
                V[ti] = (v.reshape(nz, nm, nm) * scale * scale
                         ).astype(np.float32)
                G[ti] = g.reshape(nz, nm, nm).astype(np.float16)
            rows[pair].append((mjd, BAND_IDX[obs["band"]], ph, F, V, G))
        if n % 250 == 0:
            print(f"  {n}/{len(frames)} ({n / (_time.monotonic() - t0):.1f}/s)",
                  flush=True)

    for (endpoint, role), rs in rows.items():
        rs.sort(key=lambda r: r[0])
        np.savez_compressed(
            OUT_DIR / f"{endpoint}__{role}.npz",
            z_grid=Z_GRID, mu_grid=MU_GRID, t0_mjd=T0_MJD,
            offsets=np.array(OFFSETS),
            mjd=np.array([r[0] for r in rs]),
            band_idx=np.array([r[1] for r in rs], dtype=np.uint8),
            phase=np.array([r[2] for r in rs], dtype=np.uint8),
            f=np.stack([r[3] for r in rs], axis=1),
            v=np.stack([r[4] for r in rs], axis=1),
            g=np.stack([r[5] for r in rs], axis=1))
        print(f"wrote {endpoint}__{role}: {len(rs)} epochs", flush=True)


if __name__ == "__main__":
    main()
