"""Stage-7 sampling pass for the ZTF pilot: per-epoch trajectory sample
tensors for every endpoint x role, real trajectory plus 8 offset
controls (same design as surveys/wise/scripts/sample_tensor.py).

Writes runs/ztf/calib_v1/tensors/<endpoint>__<role>.npz:
  mjd (E,), band_idx (E,) [1 g, 2 r, 3 i], phase (E,) [0/1],
  seeing (E,), f/v (9, E, NZ, 5, 5) float32 at common ZP, g float16,
  dfrac (9, E, NZ, 5, 5) float16 = matched-filter weight fraction drawn
  from difference-image pixels (hybrid search image, see photometry).

ZTF-specific choices: NZ = 192 nodes uniform in 1/z (1.85" spacing,
under the ~2" seeing); T0 = MJD 59800 (baseline midpoint); fluxes from
the difference image (science-image fallback flagged per epoch).

Processes one corridor at a time (memory), writing each endpoint x role
tensor as soon as its corridor is done.

Usage: uv run python surveys/ztf/scripts/sample_tensor.py [--only-missing]
"""

from __future__ import annotations

import json
import sys
import time as _time
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.time import Time

from sglseti import Role, adaptive_locus, load_target_registry

from sglsurvey.adapters.irsa_ztf import ZtfExactFootprint
from sglsurvey.geometry import GeometryContext
from sglsurvey.photometry import build_flux_map_ztf
from sglsurvey.records import read_records

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ztf_corridors import CORRIDOR_OF  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
COARSE_DIR = REPO / "runs" / "ztf" / "coarse_v1"
PRECISE_DIR = REPO / "runs" / "ztf" / "precise_v1"
CUT_DIR = REPO / "runs" / "ztf" / "products" / "cut"
MSK_DIR = REPO / "runs" / "ztf" / "products" / "msk"
OUT_DIR = REPO / "runs" / "ztf" / "calib_v1" / "tensors"
REGISTRY_PATH = REPO / "registries" / "pilot_wise_2026.yaml"

NZ = 192
Z_GRID = 1.0 / np.linspace(1.0 / 550.0, 1.0 / 10000.0, NZ)
MU_GRID = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
T0_MJD = 59800.0
ZP_REF = 25.0
LOCUS_BIN_DAYS = 0.05
OFFSETS = [(0.0, 0.0), (20.0, 0.0), (-20.0, 0.0), (30.0, 0.0),
           (-30.0, 0.0), (40.0, 0.0), (-40.0, 0.0), (0.0, 25.0),
           (0.0, -25.0)]
BAND_IDX = {"zg": 1, "zr": 2, "zi": 3}


def main() -> None:
    only_missing = "--only-missing" in sys.argv[1:]
    registry = load_target_registry(REGISTRY_PATH)
    ctx = GeometryContext.ztf_default()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    obs_by_id = {r["observation_id"]: r for r in read_records(
        COARSE_DIR / "records" / "observation.jsonl")}
    usable = defaultdict(list)
    for r in read_records(PRECISE_DIR / "records"
                          / "intersection_evaluation.jsonl"):
        if r["stage"] == "precise" and r["usable"] in ("usable", "partial"):
            usable[(r["endpoint_id"], r["role"])].append(r["observation_id"])
    manifest = {}
    with (CUT_DIR / "manifest.jsonl").open() as fh:
        for line in fh:
            rec = json.loads(line)
            manifest[rec["observation_id"]] = rec

    corridor_frames = defaultdict(set)
    for (e, role), oids in usable.items():
        corridor_frames[CORRIDOR_OF[e]].update(oids)
    phase_ref = {}
    for c, oids in corridor_frames.items():
        doys = np.array(sorted(obs_by_id[o]["t_mid_mjd_utc"] % 365.25
                               for o in oids))
        # reference = circular median of observing day-of-year
        ang = doys / 365.25 * 2 * np.pi
        phase_ref[c] = float((np.arctan2(np.sin(ang).mean(),
                                         np.cos(ang).mean())
                              % (2 * np.pi)) / (2 * np.pi) * 365.25)

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
                observation_time=Time(key[2] * LOCUS_BIN_DAYS, format="mjd"),
                observer=ctx.observer, relay_range=ctx.relay_range,
                tolerance_arcsec=0.5, ephemeris=ctx.ephemeris,
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
        if only_missing and (OUT_DIR / f"{p[0]}__{p[1]}.npz").exists():
            continue
        for o in oids:
            pairs_of_frame[o].append(p)
    frames_by_corridor = defaultdict(list)
    for o, pairs in pairs_of_frame.items():
        frames_by_corridor[CORRIDOR_OF[pairs[0][0]]].append(o)
    nz, nm, nt = len(Z_GRID), len(MU_GRID), len(OFFSETS)
    t0 = _time.monotonic()
    n = n_diff = 0
    for corridor, frames in sorted(frames_by_corridor.items()):
        frames.sort(key=lambda o: obs_by_id[o]["t_mid_mjd_utc"])
        rows = defaultdict(list)
        zcache.clear()
        for oid in frames:
            row = manifest.get(oid)
            obs = obs_by_id[oid]
            if row is None or "sci" not in row.get("files", {}) or not row["msk"]:
                continue
            diff = (CUT_DIR / row["files"]["diff"]
                    if "diff" in row["files"] else None)
            try:
                fm = build_flux_map_ztf(
                    CUT_DIR / row["files"]["sci"], diff, MSK_DIR / row["msk"],
                    ZtfExactFootprint.FATAL_MASK, obs["band"],
                    obs["t_mid_mjd_utc"])
            except Exception as exc:
                print(f"  fluxmap FAIL {oid}: {exc}", flush=True)
                continue
            if fm.magzp is None:
                continue
            n += 1
            n_diff += diff is not None
            mjd = obs["t_mid_mjd_utc"]
            scale = 10.0 ** ((ZP_REF - fm.magzp) / 2.5)
            dt_yr = (mjd - T0_MJD) / 365.25
            for pair in pairs_of_frame[oid]:
                endpoint, role = pair
                ph = phase_of(CORRIDOR_OF[endpoint], mjd)
                base = points_at_zgrid(endpoint, role, mjd)
                cosd = np.cos(np.deg2rad(base[:, 1]))
                dmu = MU_GRID * dt_yr / 3600.0
                ra_all = np.broadcast_to(
                    base[:, 0][:, None, None]
                    + dmu[None, :, None] / cosd[:, None, None], (nz, nm, nm))
                dec_all = np.broadcast_to(
                    base[:, 1][:, None, None] + dmu[None, None, :], (nz, nm, nm))
                F = np.empty((nt, nz, nm, nm), dtype=np.float32)
                V = np.empty_like(F)
                G = np.empty((nt, nz, nm, nm), dtype=np.float16)
                D = np.empty((nt, nz, nm, nm), dtype=np.float16)
                for ti, (dra, ddec) in enumerate(OFFSETS):
                    ra_q = ra_all + dra / 3600.0 / cosd[:, None, None]
                    dec_q = dec_all + ddec / 3600.0
                    f, v, g = fm.sample(ra_q.ravel(), dec_q.ravel())
                    F[ti] = (f.reshape(nz, nm, nm) * scale).astype(np.float32)
                    V[ti] = (v.reshape(nz, nm, nm) * scale * scale
                             ).astype(np.float32)
                    G[ti] = g.reshape(nz, nm, nm).astype(np.float16)
                    D[ti] = fm.sample_aux(ra_q.ravel(), dec_q.ravel()
                                          ).reshape(nz, nm, nm).astype(np.float16)
                rows[pair].append((mjd, BAND_IDX[obs["band"]], ph,
                                   obs["quality_flags"].get("seeing") or 0.0,
                                   diff is not None, F, V, G, D))
            if n % 200 == 0:
                print(f"  {n} maps ({n / (_time.monotonic() - t0):.1f}/s)",
                      flush=True)
        for (endpoint, role), rs in rows.items():
            rs.sort(key=lambda r: r[0])
            np.savez_compressed(
                OUT_DIR / f"{endpoint}__{role}.npz",
                z_grid=Z_GRID, mu_grid=MU_GRID, t0_mjd=T0_MJD, zp_ref=ZP_REF,
                offsets=np.array(OFFSETS),
                mjd=np.array([r[0] for r in rs]),
                band_idx=np.array([r[1] for r in rs], dtype=np.uint8),
                phase=np.array([r[2] for r in rs], dtype=np.uint8),
                seeing=np.array([r[3] for r in rs], dtype=np.float32),
                is_diff=np.array([r[4] for r in rs], dtype=np.uint8),
                f=np.stack([r[5] for r in rs], axis=1),
                v=np.stack([r[6] for r in rs], axis=1),
                g=np.stack([r[7] for r in rs], axis=1),
                dfrac=np.stack([r[8] for r in rs], axis=1))
            print(f"wrote {endpoint}__{role}: {len(rs)} epochs "
                  f"[{corridor}]", flush=True)
        del rows
    print(f"{n} flux maps built ({n_diff} from difference images)")


if __name__ == "__main__":
    main()
