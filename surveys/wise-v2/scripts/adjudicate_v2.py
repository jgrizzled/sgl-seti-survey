"""Step F / G / H: candidates under the frozen rule, calibrated vetoes,
annotations, W3/W4 confirmation procedure (hypotheses v2.0 §3.6–3.7).

Reads the set's null-ensemble analysis (candidates = cells with
R~ >= R~_FWER and a stable null) and, for each candidate, from the
tensor and the corridor's CatWISE2020 catalogue:

  rejection tests (the only grounds for rejection, each calibrated):
    a. flux-consistent catalogued static source / halo;
    b. held-out-epoch prediction test (early refit, 2022-24 forced
       photometry) with the null false-pass and injection true-pass
       rates quoted from the ledger;
    c. W3/W4 confirmation procedure — W1/W2 forced photometry at the
       same node on the cryo frames against the 300 K W2/W3 flux-ratio
       prediction, post-cryo W1/W2 persistence, catalogue counterpart.
  annotations (never grounds for rejection): parallax-phase split and
    ratio, single-epoch share, companion-band significance, cryo visit
    count, nearest catalogue source, grid-edge fit.

Every candidate becomes a Candidate record with status
``rejected:<test>`` or ``retained-ambiguous``, its global p-value, its
rank statement and the selection-function loss of the test that fired
(from the set's completeness ledger). Nothing discretionary.

Usage: uv run python surveys/wise-v2/scripts/adjudicate_v2.py --set dev|confirmatory
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v2common as C  # noqa: E402

from sglseti import stable_id  # noqa: E402
from sglsurvey import inject, nulls  # noqa: E402
from sglsurvey.photometry import _gaussian_kernel  # noqa: E402
from sglsurvey.records import Candidate, append_records  # noqa: E402
from sglsurvey.vetting import (flux_consistency, holdout_prediction_test,  # noqa: E402
                               load_catwise_vizier, parallax_phase_test,
                               radial_response_table)

MASKS = ("primary", "strict", "loose")
F = {"S_max": 0, "iz": 1, "imu": 2, "S_phase0": 3, "S_phase1": 4, "n_epochs": 5, "top_share": 6}
H = {"S_early": 0, "f_early": 1, "S_late": 2, "f_late": 3, "n_late": 4, "iz_early": 5, "imu_early": 6}
CRYO_END_MJD = 55600.0
_PRF: dict = {}


def response_table(band, pix):
    key = (band, round(pix, 3))
    if key not in _PRF:
        prf = inject.PRFGrid.load(band, C.PRF_DIR / band.lower())
        fwhm_pix = C.PSF_FWHM[band] / pix
        k = _gaussian_kernel(fwhm_pix, int(np.ceil(2.5 * fwhm_pix)))
        r, v = radial_response_table(prf.mean_template(), k)
        _PRF[key] = (r * pix, v)
    return _PRF[key]


def bb300_ratio(band_num: str, band_den: str) -> float:
    """F_nu ratio of a 300 K blackbody between two WISE bands (at lambda_iso)."""
    def bnu(lam_um):
        x = 1.4388e4 / (lam_um * 300.0)
        return (1.0 / lam_um) ** 3 / np.expm1(x)
    return bnu(inject.LAMBDA_ISO_UM[band_num]) / bnu(inject.LAMBDA_ISO_UM[band_den])


def node_stack(d, b, node, epoch_sel):
    """(S, f, sigma_f, n) of the real trajectory at ``node`` over the
    epochs of band ``b`` selected by ``epoch_sel`` (boolean over all E)."""
    eb = (d["band_idx"] == b) & epoch_sel
    if eb.sum() < 1:
        return np.nan, np.nan, np.nan, 0
    f = d["f"][0, eb][:, node[0], node[1], node[2]][:, None]
    v = d["v"][0, eb][:, node[0], node[1], node[2]][:, None]
    g = d["g"][0, eb][:, node[0], node[1], node[2]][:, None]
    S, A, B, n, _ = nulls.stack_S(f, v, g, return_parts=True, min_epochs=1)
    if B[0] <= 0:
        return np.nan, np.nan, np.nan, int(n[0])
    return float(S[0]), float(A[0] / B[0]), float(1 / np.sqrt(B[0])), int(n[0])


def adjudicate(key: str, cinfo: dict, mask: str, catalog, completeness: dict, null_cal: dict) -> Candidate:
    endpoint, role, band = key.split("/")
    b = C.BAND_IDX[band]
    mi = MASKS.index(mask)
    # cross-track variants: adjudicate the offset with the largest S_max
    best_path, best_S = C.TENSOR_DIR / f"{endpoint}__{role}.npz", -np.inf
    for tp in [best_path] + sorted(C.TENSOR_DIR.glob(f"{endpoint}__{role}__xt*.npz")):
        with np.load(tp) as dv:
            s = float(np.nan_to_num(dv["summary"][mi, 0, b - 1, 0], nan=-np.inf))
        if s > best_S:
            best_path, best_S = tp, s
    d = np.load(best_path)
    eb = d["band_idx"] == b
    ok_all = d["masks"][:, mi]
    ok = ok_all[eb]
    f0, v0, g0 = d["f"][0, eb], d["v"][0, eb], d["g"][0, eb]
    mjd, phase = d["mjd"][eb], d["phase"][eb]
    S, A, B, n, w = nulls.stack_S(f0, v0, g0, ok, return_parts=True)
    s_max, node = nulls.grid_max(S)
    iz, im1, im2 = node
    z = float(C.Z_GRID[iz]); mu = (float(C.MU_GRID[im1]), float(C.MU_GRID[im2]))
    we = w[:, iz, im1, im2]
    fe = np.nan_to_num(f0[:, iz, im1, im2].astype(np.float64))
    contrib = fe * we
    S_ph = {}
    for p in (0, 1):
        sel = ok & (phase == p)
        Bq = we[sel].sum()
        S_ph[p] = float(contrib[sel].sum() / np.sqrt(Bq)) if Bq > 0 else np.nan
    annotations = {
        "parallax_phase": parallax_phase_test(S_ph, {p: int((ok & (phase == p)).sum()) for p in (0, 1)}),
        "top_epoch_share": float(contrib.max() / A[node]) if A[node] > 0 else None,
        "n_epochs": int(n[node]), "cryo_visit_cell": bool(band in ("W3", "W4")),
        "grid_edge": {"z": bool(iz in (0, len(C.Z_GRID) - 1)),
                      "mu": bool(im1 in (0, 4) or im2 in (0, 4))},
        "fitted": {"z_au": z, "mu_arcsec_yr": mu, "S_max": float(s_max),
                   "xt_arcsec": float(d["xt_arcsec"]) if "xt_arcsec" in d.files else 0.0},
    }
    # companion-band significance at the same node (annotation)
    other = {"W1": 2, "W2": 1, "W3": 1, "W4": 1}[band]
    S_o, f_o, sf_o, n_o = node_stack(d, other, node, ok_all)
    annotations["companion_band"] = {"band": C.BAND_NAME[other], "S": S_o, "n": n_o,
                                     "significance_ratio": (float(s_max / S_o) if np.isfinite(S_o) and S_o > 0 else None)}
    # --- calibrated tests -----------------------------------------------------
    rejection = None
    tests = {}
    # a. flux-consistent static source / halo
    centre = tuple(d["centre"]); cosd = np.cos(np.deg2rad(centre[1]))
    dt_yr = (mjd - C.T0_MJD) / 365.25
    d0 = d["d0"][eb]
    tr_ra = centre[0] + (d0[:, iz, 0] + mu[0] * dt_yr) / 3600.0 / cosd
    tr_dec = centre[1] + (d0[:, iz, 1] + mu[1] * dt_yr) / 3600.0
    if catalog is not None and len(catalog.ra):
        pix = float(np.median(d["pix_arcsec"][eb]))
        rt = response_table(band, pix)
        fc = flux_consistency(catalog, tr_ra, tr_dec, we * ok, phase, C.ZP_REF, rt[0], rt[1],
                              float(s_max), S_by_phase=S_ph, search_arcsec=C.STATIC_SEARCH_ARCSEC,
                              factor=C.STATIC_FACTOR, min_ndet=C.CATWISE_MIN_NDET)
        tests["static_flux_consistent"] = fc
        annotations["nearest_catalogue_source"] = fc["sources"][0] if fc["sources"] else None
        if fc["consistent"]:
            rejection = "static_flux_consistent"
    else:
        tests["static_flux_consistent"] = {"available": False}
    # b. held-out-epoch prediction (early refit; from the build summaries)
    h = d["holdout"][mi, 0, b - 1]
    hres = holdout_prediction_test(h[H["S_early"]], h[H["f_early"]], h[H["S_late"]], h[H["f_late"]],
                                   int(h[H["n_late"]]) if np.isfinite(h[H["n_late"]]) else 0,
                                   C.HOLDOUT_MIN_S_LATE, C.HOLDOUT_MIN_FLUX_RATIO)
    hres["calibration"] = {"null_false_pass_rate": null_cal.get("false_pass_rate"),
                           "null_false_pass_ci68": null_cal.get("false_pass_ci68")}
    # v2.1: an annotation with its measured pass rates, never a rejection
    annotations["holdout_prediction"] = hres
    # c. W3/W4 confirmation procedure
    if band in ("W3", "W4"):
        cryo = d["mjd"] < CRYO_END_MJD
        f_w34 = A[node] / B[node]
        m_w34 = C.ZP_REF - 2.5 * np.log10(f_w34) if f_w34 > 0 else np.nan
        conf = {"f_node": float(f_w34), "mag_node": float(m_w34)}
        thermal_rejected = None
        for ob in ("W2", "W1"):
            S_c, f_c, sf_c, n_c = node_stack(d, C.BAND_IDX[ob], node, ok_all & cryo)
            if np.isfinite(m_w34):
                fnu = inject.fnu_from_vega_mag(band, m_w34, "blackbody_300K") * bb300_ratio(ob, band)
                m_pred = inject.vega_mag_from_fnu(ob, fnu, "blackbody_300K")
                f_pred = 10 ** (0.4 * (C.ZP_REF - m_pred))
            else:
                m_pred = f_pred = np.nan
            conf[f"{ob}_cryo"] = {"S": S_c, "f": f_c, "sigma_f": sf_c, "n": n_c,
                                  "f_pred_300K": float(f_pred), "mag_pred_300K": float(m_pred),
                                  "absent_at_3sigma": bool(np.isfinite(f_c) and np.isfinite(f_pred)
                                                           and f_c + 3 * sf_c < f_pred)}
            if ob == "W2" and conf[f"{ob}_cryo"]["absent_at_3sigma"]:
                thermal_rejected = True
        S_p, f_p, sf_p, n_p = node_stack(d, 2, node, ok_all & ~cryo)
        conf["W2_postcryo"] = {"S": S_p, "f": f_p, "sigma_f": sf_p, "n": n_p}
        conf["thermal_interpretation_rejected"] = thermal_rejected
        tests["w34_confirmation"] = conf
        if rejection is None and thermal_rejected:
            rejection = "w34_confirmation_no_w2"
    d.close()
    loss = None
    comp = completeness.get(key)
    if comp and rejection:
        tot = comp.get("n_threshold") or 0
        loss = {"injections_threshold_recovered": tot,
                "lost_to_this_test": comp.get("vetoes_fired", {}).get(rejection, 0),
                "fraction": (comp.get("vetoes_fired", {}).get(rejection, 0) / tot) if tot else None}
    status = f"rejected:{rejection}" if rejection else "retained-ambiguous"
    return Candidate(
        candidate_id=stable_id("cnd", {"cell": key, "hypothesis": C.HYPOTHESIS_VERSION,
                                       "freeze": C.freeze_hash(), "mask": mask}),
        analysis_run_id=cinfo.get("analysis_run_id", "see-null-ensemble"),
        endpoint_id=endpoint, role=role, observation_ids=(),
        fitted_z_au=z, fitted_residual_motion={"mu_ra": mu[0], "mu_dec": mu[1]},
        model_comparison={"S_max": float(s_max), "T": cinfo["T"], "R": cinfo["R"], "q95": cinfo["q95"],
                          "R_norm": cinfo["R_norm"], "tests": tests},
        status=status, veto_reason=rejection, rejection_test=rejection,
        global_p_value=cinfo["global_p"], rank_statement=cinfo["rank"],
        annotations=annotations,
        extra={"band": band, "set": cinfo.get("set"), "mask": mask, "selection_function_loss": loss})


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=["dev", "confirmatory"], required=True)
    ap.add_argument("--mask", default="primary")
    a = ap.parse_args()
    ne_all = json.loads((C.NULL_DIR / f"{a.set}_null_ensemble.json").read_text())
    ne = ne_all[a.mask]
    comp_path = C.RUN_DIR / "completeness" / f"{a.set}_completeness.json"
    completeness = json.loads(comp_path.read_text())["cells"] if comp_path.exists() else {}
    cands = []
    cats = {}
    for key in ne["candidates"]:
        cinfo = dict(ne["cells"][key]); cinfo["set"] = a.set
        endpoint = key.split("/")[0]
        d = np.load(C.TENSOR_DIR / f"{endpoint}__{key.split('/')[1]}.npz")
        corridor = str(d["corridor"]); centre = tuple(d["centre"]); d.close()
        if corridor not in cats:
            try:
                cats[corridor] = load_catwise_vizier(centre[0], centre[1], 0.1, C.RUN_DIR / "catwise")
            except Exception as exc:
                print(f"  catalogue unavailable for {corridor}: {exc}")
                cats[corridor] = None
        c = adjudicate(key, cinfo, a.mask, cats[corridor], completeness, ne["holdout_null_calibration"])
        cands.append(c)
        print(f"{key:26s} R~={cinfo['R_norm']:.2f} p_global={cinfo['global_p']:.3f} -> {c.status}"
              f"  phase={c.annotations['parallax_phase']['signature']} "
              f"S_ph={[round(x, 1) for x in c.annotations['parallax_phase']['S_by_phase'].values()]}", flush=True)
    out = C.RUN_DIR / "records" / f"candidate_{a.set}.jsonl"
    if out.exists():
        out.unlink()
    append_records(out, cands)
    summary = {"set": a.set, "mask": a.mask, "R_fwer_norm": ne["fwer"]["R_fwer"],
               "n_candidates": len(cands),
               "status_counts": {s: sum(1 for c in cands if c.status == s) for s in set(c.status for c in cands)},
               "retained_ambiguous": [c.endpoint_id + "/" + c.role + "/" + c.extra["band"]
                                      for c in cands if c.status == "retained-ambiguous"],
               "written": str(out), "utc": datetime.now(timezone.utc).isoformat()}
    (C.RUN_DIR / "records" / f"adjudication_{a.set}.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
