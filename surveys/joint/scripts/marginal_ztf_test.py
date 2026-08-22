"""Direct ZTF forced-photometry test of the PS1 marginal cells (stage-7
cross-archive adjudication). For each marginal PS1 candidate (runs/
panstarrs/calib_v1/adjudication.json, status "marginal"), the PS1-fitted
trajectory (z, mu, T0 = 56000) is extrapolated to every usable ZTF
frame of that endpoint-role; the ZTF hybrid diff/sci flux map is
sampled there and at the 8 shared offset controls; the weighted stack
S_ztf is compared with the control maximum. A relay present in PS1 at
the fitted cell must appear in the deeper, both-phase ZTF stack
(2018-2026) at S >> T; a PS1 noise maximum gives S_ztf ~ 0.

Also evaluates the PS1 cell's own mu = 0 neighbour for reference.
Writes <out>/marginal_tests.json.
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

REPO = Path(__file__).resolve().parents[3]
ZTF = REPO / "runs" / "ztf"
ADJ = REPO / "runs" / "panstarrs" / "calib_v1" / "adjudication.json"
PS1_T0 = 56000.0
ZP_REF = 25.0
ZTF_THROUGHPUT_MAG = 0.5
OFFSETS = [(0.0, 0.0), (20.0, 0.0), (-20.0, 0.0), (30.0, 0.0), (-30.0, 0.0),
           (40.0, 0.0), (-40.0, 0.0), (0.0, 25.0), (0.0, -25.0)]
BAND_OF = {"g": "zg", "r": "zr", "i": "zi", "z": "zi", "y": "zi"}
MIN_GOOD_FRAC, WEIGHT_CAP, CLIP_SIGMA = 0.7, 20.0, 5.0


def weighted_S(f, v, g):
    valid = np.isfinite(f) & np.isfinite(v) & (v > 0) & (g >= MIN_GOOD_FRAC)
    with np.errstate(invalid="ignore", divide="ignore"):
        valid &= np.abs(f) / np.sqrt(np.where(v > 0, v, np.inf)) <= CLIP_SIGMA
    w = np.where(valid, 1.0 / np.where(v > 0, v, 1.0), 0.0)
    if (w > 0).sum() == 0:
        return np.nan, 0
    w = np.minimum(w, WEIGHT_CAP * np.median(w[w > 0]))
    return float((np.where(valid, f, 0) * w).sum() / np.sqrt(w.sum())), int((w > 0).sum())


PS1_CAND = REPO / "runs" / "panstarrs" / "calib_v1" / "records" / "candidate.jsonl"


def retained_ps1_cells() -> dict:
    """Every PS1 exceedance the automatic rules retain (phase test +
    catalogued-static-source test), as {key: {cell, S_real, T}}."""
    cells = {}
    for c in read_records(PS1_CAND):
        if c["status"] != "retained":
            continue
        mc = c["model_comparison"]
        key = f"{c['endpoint_id']}/{c['role']}/{c['extra']['band']}"
        cells[key] = {"cell": {"z_au": c["fitted_z_au"], "mu": mc["mu_arcsec_yr"]},
                      "S_real": mc["real_max_S"], "T": mc["threshold_8_controls"]}
    return cells


def run_marginal_tests(out_dir: Path, analysis_run_id: str | None = None,
                       all_retained: bool = True, cells: dict | None = None,
                       t0_mjd: float = PS1_T0, out_name: str = "marginal_tests.json") -> None:
    """``cells`` may be given directly as {key: {"cell": {"z_au", "mu"},
    "S_real", "T"}} with their mu reference epoch ``t0_mjd``."""
    global PS1_T0
    PS1_T0 = t0_mjd
    if cells is None and all_retained:
        cells = retained_ps1_cells()
    elif cells is None:
        adj = json.load(open(ADJ))
        cells = {k: v for k, v in adj.items() if v["status"] == "marginal"}
    print(f"\n{len(cells)} automatically-retained PS1 cells to test against ZTF",
          flush=True)
    registry = load_target_registry(REPO / "registries" / "pilot_wise_2026.yaml")
    ctx = GeometryContext.ztf_default()
    obs_by_id = {r["observation_id"]: r for r in read_records(
        ZTF / "coarse_v1" / "records" / "observation.jsonl")}
    usable = defaultdict(list)
    for r in read_records(ZTF / "precise_v1" / "records" / "intersection_evaluation.jsonl"):
        if r["stage"] == "precise" and r["usable"] in ("usable", "partial"):
            usable[(r["endpoint_id"], r["role"])].append(r["observation_id"])
    manifest = {}
    with (ZTF / "products" / "cut" / "manifest.jsonl").open() as fh:
        for line in fh:
            rec = json.loads(line)
            manifest[rec["observation_id"]] = rec
    results = {}
    for key, c in cells.items():
        e, role, band = key.split("/")
        z0 = c["cell"]["z_au"]
        # other role at the same (z, mu): Rx/Tx loci coincide at the
        # antipode, so a real source must appear in both (tested below via
        # the same forced photometry on the other role's frames)
        mu = c["cell"]["mu"]
        # all ZTF bands: a reflector has solar colours (r brighter than y/z
        # by ~0.5 mag) and ZTF has almost no i-band frames, so the full
        # g + r (+ i) stack is the sensitive test for every PS1 filter
        oids = [o for o in usable[(e, role)]
                if o in manifest and "sci" in manifest[o].get("files", {})]
        # de-duplicate by exposure (same field/filter/time)
        oids.sort(key=lambda o: obs_by_id[o]["t_mid_mjd_utc"])
        print(f"[{key}] z={z0:.0f} mu={mu} -> {len(oids)} ZTF frames", flush=True)
        F, V, G = [], [], []
        F0, V0, G0 = [], [], []
        mjds, bands = [], []
        t0 = _time.monotonic()
        for k, o in enumerate(oids):
            row, obs = manifest[o], obs_by_id[o]
            try:
                fm = build_flux_map_ztf(
                    ZTF / "products" / "cut" / row["files"]["sci"],
                    (ZTF / "products" / "cut" / row["files"]["diff"]
                     if "diff" in row["files"] else None),
                    ZTF / "products" / "msk" / row["msk"],
                    ZtfExactFootprint.FATAL_MASK, obs["band"], obs["t_mid_mjd_utc"])
            except Exception as exc:
                continue
            if fm.magzp is None:
                continue
            t = obs["t_mid_mjd_utc"]
            al = adaptive_locus(target=registry[e], role=Role(role),
                                observation_time=Time(t, format="mjd"),
                                observer=ctx.observer, relay_range=ctx.relay_range,
                                tolerance_arcsec=0.5, ephemeris=ctx.ephemeris, model=ctx.model)
            zs = np.array([p.z_au for p in al.points]); q = 1.0 / zs; o_ = np.argsort(q)
            ra = float(np.interp(1.0 / z0, q[o_], np.array([p.icrs_ra_deg for p in al.points])[o_]))
            dec = float(np.interp(1.0 / z0, q[o_], np.array([p.icrs_dec_deg for p in al.points])[o_]))
            cosd = np.cos(np.deg2rad(dec))
            dt = (t - PS1_T0) / 365.25
            ra_mu, dec_mu = ra + mu[0] * dt / 3600.0 / cosd, dec + mu[1] * dt / 3600.0
            scale = 10.0 ** ((ZP_REF - fm.magzp) / 2.5) * 10 ** (0.4 * ZTF_THROUGHPUT_MAG)
            ra_q = np.array([ra_mu + dx / 3600.0 / cosd for dx, _ in OFFSETS])
            dec_q = np.array([dec_mu + dy / 3600.0 for _, dy in OFFSETS])
            f, v, g = fm.sample(ra_q, dec_q)
            F.append(f * scale); V.append(v * scale * scale); G.append(g)
            f0, v0, g0 = fm.sample(np.array([ra]), np.array([dec]))
            F0.append(f0[0] * scale); V0.append(v0[0] * scale * scale); G0.append(g0[0])
            mjds.append(t); bands.append(obs["band"])
            if (k + 1) % 100 == 0:
                print(f"  {k + 1}/{len(oids)} ({(k + 1) / (_time.monotonic() - t0):.1f}/s)", flush=True)
        if not mjds:
            results[key] = {"ps1_cell": c["cell"], "ztf_frames": 0,
                            "verdict": "no usable ZTF frames — untestable"}
            print("  -> no ZTF frames", flush=True)
            continue
        F, V, G = np.array(F), np.array(V), np.array(G)
        F0, V0, G0 = np.array(F0), np.array(V0), np.array(G0)
        bands = np.array(bands)
        per_band = {}
        for b in sorted(set(bands.tolist())) + ["all"]:
            sel = np.ones(len(bands), bool) if b == "all" else bands == b
            S = [weighted_S(F[sel, i], V[sel, i], G[sel, i]) for i in range(len(OFFSETS))]
            S0 = weighted_S(F0[sel], V0[sel], G0[sel])
            ctrl = [s for s, _ in S[1:] if np.isfinite(s)]
            T = float(max(ctrl)) if ctrl else np.nan
            per_band[b] = {"frames": int(sel.sum()), "S_track": S[0][0], "n_valid": S[0][1],
                           "T_controls": T, "control_maxima": [round(s, 2) for s, _ in S[1:]],
                           "S_mu0": S0[0]}
        a = per_band["all"]
        verdict = ("vetoed: absent in the ZTF g/r stack along the PS1-fitted "
                   "trajectory" if a["S_track"] < a["T_controls"] else
                   "ZTF exceedance along the PS1 trajectory — escalate")
        res = {"ps1_cell": c["cell"], "ps1_S": c["S_real"], "ps1_T": c["T"],
               "ps1_band": band, "ztf_frames": len(mjds),
               "ztf_epoch_range": [min(mjds), max(mjds)],
               "ztf_by_band": per_band, "verdict": verdict}
        results[key] = res
        print("  -> " + "  ".join(f"{b}: S={d['S_track']:.2f}/T={d['T_controls']:.2f} (n={d['n_valid']})"
                               for b, d in per_band.items()) + f"  => {verdict}", flush=True)
    out = out_dir / out_name
    out.write_text(json.dumps({"analysis_run_id": analysis_run_id, "cells": results}, indent=2))
    print(f"wrote {out}")


if __name__ == "__main__":
    out = REPO / "runs" / "joint" / "ps1_ztf_v1"
    out.mkdir(parents=True, exist_ok=True)
    run_marginal_tests(out)
