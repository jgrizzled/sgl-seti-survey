"""Trajectory-aware forced-photometry stack (plan §4.6).

For every endpoint x role, sample the matched-filter flux maps of all
usable frames at the positions predicted by a (z, mu_resid) grid —
64 log-spaced relay distances x a 5x5 residual-proper-motion grid
inside the frozen |mu| <= 1 "/yr bound — and accumulate inverse-
variance-weighted stacks per band. This is the joint relay-distance +
bounded-residual-motion fit, sensitive below single-exposure depth.

Built-in checks:
  - parallax-phase-split stacks (the static-star veto validated in
    screen_v1: a real relay contributes at BOTH alternating phases);
  - an offset-trajectory negative control (+35" RA) run through the
    identical machinery.

Fluxes are converted to a common zero point (magzp -> 20.0) before
stacking. Emits an AnalysisRun record and a stacks .npz; constraints
wait for injection calibration (plan §4.7).

Usage: uv run python surveys/wise/scripts/forced_stack.py
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import time as _time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy.time import Time

from sglseti import (Role, adaptive_locus, canonical_json,
                     load_target_registry, stable_hash, stable_id)

from sglsurvey.adapters.irsa_wise import WiseExactFootprint
from sglsurvey.geometry import GeometryContext
from sglsurvey.photometry import build_flux_map
from sglsurvey.records import AnalysisRun, append_records, read_records

REPO = Path(__file__).resolve().parents[3]
COARSE_DIR = REPO / "runs" / "wise" / "coarse_v1"
PRECISE_DIR = REPO / "runs" / "wise" / "precise_v1"
CUT_DIR = REPO / "runs" / "wise" / "products" / "cut"
MSK_DIR = REPO / "runs" / "wise" / "products" / "msk"
RUN_DIR = REPO / "runs" / "wise" / "stack_v1"
REGISTRY_PATH = REPO / "registries" / "pilot_wise_2026.yaml"
HYPOTHESES_PATH = REPO / "surveys" / "wise" / "hypotheses.md"
HYPOTHESIS_VERSION = "wise-hypotheses-v1.0"

Z_GRID = np.geomspace(550.0, 10000.0, 64)
MU_GRID = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])   # arcsec/yr, each axis
T0_MJD = 57800.0                                   # mu reference epoch
CONTROL_OFFSET_ARCSEC = 35.0
ZP_REF = 20.0
LOCUS_BIN_DAYS = 0.5
MIN_GOOD_FRAC = 0.7

CORRIDOR_OF = {
    "barnard-star": "barnard", "ross-154": "ross154",
    "lalande-21185": "lalande", "alpha-cen-a": "alphacen",
    "alpha-cen-b": "alphacen", "sirius-a": "sirius",
    "sirius-b": "sirius",
}


def main() -> None:
    registry = load_target_registry(REGISTRY_PATH)
    ctx = GeometryContext.wise_coarse_default()
    hyp_hash = "sha256:" + hashlib.sha256(
        HYPOTHESES_PATH.read_bytes()).hexdigest()

    obs_by_id = {r["observation_id"]: r for r in read_records(
        COARSE_DIR / "records" / "observation.jsonl")}
    usable = defaultdict(list)    # (endpoint, role) -> [obs_id]
    ixn_ids = []
    for r in read_records(PRECISE_DIR / "records"
                          / "intersection_evaluation.jsonl"):
        if r["stage"] == "precise" and r["usable"] in ("usable", "partial"):
            usable[(r["endpoint_id"], r["role"])].append(
                r["observation_id"])
            ixn_ids.append(r["intersection_id"])

    manifest = {json.loads(l)["observation_id"]: json.loads(l)
                for l in open(CUT_DIR / "manifest.jsonl")}

    # Phase assignment per corridor from usable epoch day-of-year.
    corridor_frames = defaultdict(set)
    for (e, role), oids in usable.items():
        corridor_frames[CORRIDOR_OF[e]].update(oids)
    phase_ref = {}
    for c, oids in corridor_frames.items():
        doys = np.array([obs_by_id[o]["t_mid_mjd_utc"] % 365.25
                         for o in oids])
        phase_ref[c] = float(np.median(doys[:50]))

    def phase_of(corridor, mjd):
        d = (mjd % 365.25) - phase_ref[corridor]
        d = (d + 182.625) % 365.25 - 182.625
        return 0 if abs(d) < 91.3 else 1

    # Locus interpolated at Z_GRID per (endpoint, role, halfday).
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

    nz, nm = len(Z_GRID), len(MU_GRID)
    pairs = sorted(usable)
    bands = ["W1", "W2", "W3", "W4"]

    def acc():
        return {b: {"A": np.zeros((nz, nm, nm)),
                    "B": np.zeros((nz, nm, nm)),
                    "Aph": np.zeros((2, nz, nm, nm)),
                    "Bph": np.zeros((2, nz, nm, nm)),
                    "n": 0} for b in bands}

    stacks = {p: acc() for p in pairs}
    controls = {p: acc() for p in pairs}
    epoch_rows = defaultdict(list)   # (pair) -> per-epoch mu=0 slices

    frames = sorted(corridor_frames["barnard"] | corridor_frames["ross154"]
                    | corridor_frames["lalande"]
                    | corridor_frames["alphacen"]
                    | corridor_frames["sirius"])
    pairs_of_frame = defaultdict(list)
    for p, oids in usable.items():
        for o in oids:
            pairs_of_frame[o].append(p)

    t0 = _time.monotonic()
    n_built = 0
    n_skip = 0
    for oid in frames:
        row = manifest.get(oid)
        obs = obs_by_id[oid]
        if row is None or "int" not in row.get("files", {}):
            n_skip += 1
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
            n_skip += 1
            continue
        n_built += 1
        band = obs["band"]
        mjd = obs["t_mid_mjd_utc"]
        zp = fm.magzp if fm.magzp else ZP_REF
        scale = 10.0 ** ((ZP_REF - zp) / 2.5)
        dt_yr = (mjd - T0_MJD) / 365.25
        ph = phase_of(row["corridor"], mjd)

        for pair in pairs_of_frame[oid]:
            endpoint, role = pair
            base = points_at_zgrid(endpoint, role, mjd)   # (nz, 2)
            cosd = np.cos(np.deg2rad(base[:, 1]))
            # trajectory positions: base + mu*dt, all (nz, nm, nm)
            dmu = MU_GRID * dt_yr / 3600.0   # deg offsets
            ra_all = (base[:, 0][:, None, None]
                      + dmu[None, :, None] / cosd[:, None, None])
            dec_all = (base[:, 1][:, None, None]
                       + dmu[None, None, :])
            for target, offset in ((stacks, 0.0),
                                   (controls, CONTROL_OFFSET_ARCSEC)):
                ra_q = np.broadcast_to(
                    ra_all + offset / 3600.0 / cosd[:, None, None],
                    (nz, nm, nm))
                dec_q = np.broadcast_to(dec_all, (nz, nm, nm))
                f, v, g = fm.sample(ra_q.ravel(), dec_q.ravel())
                f = f.reshape(nz, nm, nm) * scale
                v = v.reshape(nz, nm, nm) * scale * scale
                g = g.reshape(nz, nm, nm)
                ok = np.isfinite(f) & np.isfinite(v) & (v > 0) \
                    & (g >= MIN_GOOD_FRAC)
                w = np.where(ok, 1.0 / np.where(v > 0, v, 1.0), 0.0)
                a = target[pair][band]
                a["A"] += np.where(ok, f, 0.0) * w
                a["B"] += w
                a["Aph"][ph] += np.where(ok, f, 0.0) * w
                a["Bph"][ph] += w
                a["n"] += 1
                if offset == 0.0:
                    i0 = nm // 2
                    epoch_rows[pair].append(
                        (mjd, band, ph,
                         f[:, i0, i0], v[:, i0, i0], g[:, i0, i0]))
        if n_built % 250 == 0:
            rate = n_built / (_time.monotonic() - t0)
            print(f"  {n_built}/{len(frames)} flux maps ({rate:.1f}/s)",
                  flush=True)

    print(f"built {n_built} flux maps, skipped {n_skip}", flush=True)

    # ---- results ----------------------------------------------------
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    out = {}
    summary_lines = []
    for pair in pairs:
        e, role = pair
        for b in bands:
            a = stacks[pair][b]
            c = controls[pair][b]
            if a["B"].max() <= 0:
                continue
            with np.errstate(divide="ignore", invalid="ignore"):
                S = a["A"] / np.sqrt(np.where(a["B"] > 0, a["B"], np.inf))
                Sc = c["A"] / np.sqrt(np.where(c["B"] > 0, c["B"],
                                               np.inf))
                Sph = a["Aph"] / np.sqrt(np.where(a["Bph"] > 0, a["Bph"],
                                                  np.inf))
            i = np.unravel_index(np.nanargmax(S), S.shape)
            key = f"{e}/{role}/{b}"
            out[key] = {"S": S, "S_control": Sc, "S_phase": Sph,
                        "A": a["A"], "B": a["B"]}
            summary_lines.append(
                f"{e:15s} {role} {b}: max S={S[i]:6.2f} at "
                f"z={Z_GRID[i[0]]:7.1f} AU mu=({MU_GRID[i[1]]:+.1f},"
                f"{MU_GRID[i[2]]:+.1f}) | phase S=("
                f"{Sph[0][i]:5.2f},{Sph[1][i]:5.2f}) | "
                f"control max={np.nanmax(Sc):5.2f}")

    npz_path = RUN_DIR / "stacks.npz"
    np.savez_compressed(
        npz_path, z_grid=Z_GRID, mu_grid=MU_GRID, t0_mjd=T0_MJD,
        **{k.replace("/", "__") + "__" + q: v
           for k, d in out.items() for q, v in d.items()})
    ep_path = RUN_DIR / "epoch_fluxes.npz"
    np.savez_compressed(ep_path, **{
        f"{e}__{r}": np.array(
            [(m, {"W1": 1, "W2": 2, "W3": 3, "W4": 4}[b], ph)
             for m, b, ph, *_ in rows]) for (e, r), rows
        in epoch_rows.items()},
        **{f"{e}__{r}__flux": np.array([f for *_, f, v, g in rows])
           for (e, r), rows in epoch_rows.items()},
        **{f"{e}__{r}__var": np.array([v for *_, f, v, g in rows])
           for (e, r), rows in epoch_rows.items()})

    config = {
        "pipeline_id": "wise-forced-stack",
        "z_grid": [float(Z_GRID[0]), float(Z_GRID[-1]), nz],
        "mu_grid_arcsec_yr": MU_GRID.tolist(), "t0_mjd": T0_MJD,
        "zp_ref": ZP_REF, "min_good_frac": MIN_GOOD_FRAC,
        "control_offset_arcsec": CONTROL_OFFSET_ARCSEC,
        "psf_model": "gaussian-fwhm-per-band",
        "background_model": "sigma-clipped-cutout-median",
        "locus_tolerance_arcsec": 2.0, "locus_bin_days": LOCUS_BIN_DAYS,
    }
    try:
        commit = subprocess.run(
            ["git", "-C", str(REPO.parent / "sglseti"), "rev-parse",
             "--short", "HEAD"], capture_output=True, text=True,
            timeout=10).stdout.strip()
    except Exception:
        commit = "unknown"
    run = AnalysisRun(
        analysis_run_id=stable_id("run", {
            "config": json.loads(canonical_json(config)),
            "observation_set": stable_hash(sorted(
                {o for p in pairs for o in usable[p]})),
            "intersection_set": stable_hash(sorted(ixn_ids)),
            "registry": registry.source_hash,
            "hypothesis": HYPOTHESIS_VERSION}),
        pipeline_id="wise-forced-stack", pipeline_version="0.1.0",
        config=config,
        observation_set_hash=stable_hash(sorted(
            {o for p in pairs for o in usable[p]})),
        intersection_set_hash=stable_hash(sorted(ixn_ids)),
        registry_source_hash=registry.source_hash,
        hypothesis_version=HYPOTHESIS_VERSION,
        environment={"sglseti_commit": commit,
                     "hypothesis_hash": hyp_hash},
        started_utc=datetime.now(timezone.utc).isoformat(),
        finished_utc=datetime.now(timezone.utc).isoformat(),
        output_files={"stacks": str(npz_path.name),
                      "epoch_fluxes": str(ep_path.name)},
    )
    append_records(RUN_DIR / "records" / "analysis_run.jsonl", [run])

    print(f"\nanalysis_run: {run.analysis_run_id}")
    print("\n=== stack summary (max significance per pair x band) ===")
    for line in summary_lines:
        print(line)


if __name__ == "__main__":
    main()
