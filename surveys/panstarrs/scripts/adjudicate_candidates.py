"""Stage-7 adjudication of Candidate records that survive the automatic
phase test in injection_calibrate.py (PS1 scale-up).

For every retained candidate, four independent tests at the fitted
(z, mu) cell, all computable from the sample tensors and the DR2 mean
table already snapshotted by catalog_screen.py:

1. split-half: weighted-stack S from the early half and late half of
   the epochs (a persistent source is significant in both; a single
   outlying epoch or a static star sampled at one season is not);
2. static-star: DR2 mean-table object within STAR_RADIUS of the
   trajectory position at the median major-phase epoch (a static star
   at the phase-0 position reproduces a high phase-0 S and a chance
   phase-1 S);
3. other-bands: S at the same (z, mu) cell in the other filters of the
   same endpoint-role, against each filter's own control threshold (a
   reflector or self-luminous source is achromatic to first order);
4. grid-edge: fitted z within 2 nodes of 550/10,000 AU or |mu| at the
   grid limit (noise maxima pile up at search-box edges).

Verdicts: vetoed if the split-half test fails (either half < 2 sigma)
or a catalogued star sits on the major-phase position; otherwise
"marginal" with the supporting evidence listed. Writes
runs/panstarrs/calib_v1/adjudication.json and appends superseding
Candidate records (same candidate_id, status updated) to
runs/panstarrs/calib_v1/records/candidate_adjudicated.jsonl.

Usage: uv run python surveys/panstarrs/scripts/adjudicate_candidates.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.time import Time

from sglseti import Role, adaptive_locus, load_target_registry

from sglsurvey.geometry import GeometryContext
from sglsurvey.records import Candidate, append_records, read_records

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ps1_corridors import CORRIDOR_OF  # noqa: E402
import importlib.util as _ilu  # the WISE scripts dir also has a sample_tensor

_spec = _ilu.spec_from_file_location(
    "ps1_sample_tensor", Path(__file__).resolve().parent / "sample_tensor.py")
_st = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(_st)
load_mean_stars, BAND_IDX = _st.load_mean_stars, _st.BAND_IDX

REPO = Path(__file__).resolve().parents[3]
CAL = REPO / "runs" / "panstarrs" / "calib_v1"
STAR_RADIUS_ARCSEC = 2.0
MIN_GOOD_FRAC, WEIGHT_CAP, CLIP_SIGMA, MIN_EPOCHS = 0.7, 20.0, 5.0, 5
BAND_NAME = {v: k for k, v in BAND_IDX.items()}


def stack_S(f, v, g, sel, cell):
    """Weighted-stack significance at one cell over epochs ``sel`` with
    the calibration's weighting (cap, clip)."""
    fz = f[sel][:, cell[0], cell[1], cell[2]]
    vz = v[sel][:, cell[0], cell[1], cell[2]]
    gz = g[sel][:, cell[0], cell[1], cell[2]]
    ok = np.isfinite(fz) & np.isfinite(vz) & (vz > 0) & (gz >= MIN_GOOD_FRAC)
    with np.errstate(invalid="ignore", divide="ignore"):
        ok &= np.abs(fz) / np.sqrt(np.where(vz > 0, vz, np.inf)) <= CLIP_SIGMA
    w = np.where(ok, 1.0 / np.where(vz > 0, vz, 1.0), 0.0)
    if (w > 0).sum() == 0:
        return np.nan, 0
    w = np.minimum(w, WEIGHT_CAP * np.median(w[w > 0]))
    B = w.sum()
    return float((np.where(ok, fz, 0) * w).sum() / np.sqrt(B)), int((w > 0).sum())


def main() -> None:
    registry = load_target_registry(REPO / "registries" / "pilot_wise_2026.yaml")
    ctx = GeometryContext.ps1_default()
    tr = json.load(open(CAL / "threshold_report.json"))
    cands = [c for c in read_records(CAL / "records" / "candidate.jsonl")
             if c["status"] == "retained"]
    print(f"{len(cands)} retained candidates to adjudicate")
    report, superseded = {}, []
    stars_cache = {}
    for c in cands:
        e, role, band = c["endpoint_id"], c["role"], c["extra"]["band"]
        key = f"{e}/{role}/{band}"
        r = tr[key]
        d = np.load(CAL / "tensors" / f"{e}__{role}.npz")
        z_grid, mu_grid = d["z_grid"], d["mu_grid"]
        zi = int(np.argmin(np.abs(z_grid - r["real_max_z"])))
        mi = int(np.argmin(np.abs(mu_grid - r["real_max_mu"][0])))
        mj = int(np.argmin(np.abs(mu_grid - r["real_max_mu"][1])))
        cell = (zi, mi, mj)
        b = BAND_IDX[band]
        eb = d["band_idx"] == b
        f, v, g = d["f"][0], d["v"][0], np.asarray(d["g"][0], dtype=np.float32)
        mjd, phase = d["mjd"], d["phase"]
        # 1. split-half
        idx = np.flatnonzero(eb)
        order = idx[np.argsort(mjd[idx])]
        half = len(order) // 2
        sel_a = np.zeros(len(mjd), bool); sel_a[order[:half]] = True
        sel_b = np.zeros(len(mjd), bool); sel_b[order[half:]] = True
        S_a, n_a = stack_S(f, v, g, sel_a, cell)
        S_b, n_b = stack_S(f, v, g, sel_b, cell)
        # 2. static star at the major-phase position
        maj = 0 if r["phase_n"]["0"] >= r["phase_n"]["1"] else 1
        t_maj = float(np.median(mjd[eb & (phase == maj)]))
        al = adaptive_locus(target=registry[e], role=Role(role),
                            observation_time=Time(t_maj, format="mjd"),
                            observer=ctx.observer, relay_range=ctx.relay_range,
                            tolerance_arcsec=0.5, ephemeris=ctx.ephemeris,
                            model=ctx.model)
        zs = np.array([p.z_au for p in al.points])
        ra = np.array([p.icrs_ra_deg for p in al.points])
        dec = np.array([p.icrs_dec_deg for p in al.points])
        q = 1.0 / zs; o = np.argsort(q)
        ra0 = float(np.interp(1.0 / z_grid[zi], q[o], ra[o]))
        dec0 = float(np.interp(1.0 / z_grid[zi], q[o], dec[o]))
        dt = (t_maj - float(d["t0_mjd"])) / 365.25
        cosd = np.cos(np.deg2rad(dec0))
        ra0 += mu_grid[mi] * dt / 3600.0 / cosd
        dec0 += mu_grid[mj] * dt / 3600.0
        corr = CORRIDOR_OF[e]
        if corr not in stars_cache:
            stars_cache[corr] = load_mean_stars(corr)
        st = stars_cache[corr]
        sep = np.hypot((st["ra"] - ra0) * cosd, st["dec"] - dec0) * 3600.0
        near = np.flatnonzero(sep < STAR_RADIUS_ARCSEC)
        star = None
        if len(near):
            k = near[np.argmin(sep[near])]
            star = {"sep_arcsec": round(float(sep[k]), 2),
                    "ndet": int(st["ndet"][k]),
                    "mags": {bb: round(float(st[f"{bb}_mag"][k]), 2)
                             for bb in BAND_IDX if st[f"{bb}_mag"][k] > 0}}
        # 3. other bands at the same cell
        others = {}
        for bb in sorted(set(d["band_idx"].tolist())):
            if bb == b:
                continue
            S_o, n_o = stack_S(f, v, g, d["band_idx"] == bb, cell)
            T_o = tr.get(f"{e}/{role}/{BAND_NAME[bb]}", {}).get("T")
            others[BAND_NAME[bb]] = {"S": None if np.isnan(S_o) else round(S_o, 2),
                                     "n": n_o, "T": T_o}
        # 4. grid edge
        edge = (zi <= 1 or zi >= len(z_grid) - 2
                or abs(mu_grid[mi]) == mu_grid.max()
                or abs(mu_grid[mj]) == mu_grid.max())
        # verdict
        reasons = []
        if not (S_a >= 2.0 and S_b >= 2.0):
            reasons.append(f"split-half fails (early S={S_a:.1f}/n={n_a}, "
                           f"late S={S_b:.1f}/n={n_b}): not persistent")
        if star is not None:
            reasons.append(f"catalogued DR2 star {star['sep_arcsec']}\" from the "
                           f"major-phase track position (static source)")
        if edge:
            reasons.append("peak at the z/mu search-box edge (noise-maximum "
                           "signature; not a veto on its own)")
        status = "vetoed" if (not (S_a >= 2.0 and S_b >= 2.0) or star) else "marginal"
        report[key] = {"cell": {"z_au": float(z_grid[zi]),
                                "mu": [float(mu_grid[mi]), float(mu_grid[mj])]},
                       "S_real": r["real_max_S"], "T": r["T"],
                       "phase_S": r["phase_S"], "phase_n": r["phase_n"],
                       "split_half": {"early": [round(S_a, 2), n_a],
                                      "late": [round(S_b, 2), n_b]},
                       "track_position_major_phase": [ra0, dec0, t_maj],
                       "nearest_star": star, "other_bands": others,
                       "grid_edge": bool(edge), "status": status,
                       "reasons": reasons}
        superseded.append(Candidate(**{**c, "status": status,
                                       "veto_reason": ("; ".join(reasons)
                                                       if status == "vetoed"
                                                       else None),
                                       "extra": {**c["extra"],
                                                 "adjudication": "stage-7 split-half / static-star / other-band / grid-edge",
                                                 "supersedes_status": "retained"}}))
        print(f"{key:28s} S={r['real_max_S']:.2f} T={r['T']:.2f} z={z_grid[zi]:.0f} "
              f"mu={mu_grid[mi]:+.1f},{mu_grid[mj]:+.1f} | halves {S_a:.1f}/{S_b:.1f} "
              f"| star {star['sep_arcsec'] if star else '-'} | other "
              + " ".join(f"{k}:{v['S']}" for k, v in others.items())
              + f" | edge={edge} -> {status.upper()}")
    (CAL / "adjudication.json").write_text(json.dumps(report, indent=2))
    append_records(CAL / "records" / "candidate_adjudicated.jsonl", superseded)


if __name__ == "__main__":
    main()
