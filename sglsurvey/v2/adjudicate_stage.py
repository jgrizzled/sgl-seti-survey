"""v2 engine: candidates under the frozen rule, calibrated vetoes and
annotations (profile-parameterised port of
surveys/wise-v2/scripts/adjudicate_v2.py). The W3/W4 confirmation
procedure is WISE-specific and stays in the WISE script; a profile may
declare a cross-band consistency hook later."""

from __future__ import annotations

import json
from datetime import datetime, timezone

import numpy as np

from sglseti import stable_id

from sglsurvey import nulls
from sglsurvey.photometry import _gaussian_kernel
from sglsurvey.records import Candidate, append_records
from sglsurvey.v2.inject_stage import oversampled_template
from sglsurvey.vetting import (flux_consistency, holdout_prediction_test, parallax_phase_test,
                               radial_response_table)

MASKS = ("primary", "strict", "loose")
F = {"S_max": 0, "iz": 1, "imu": 2, "S_phase0": 3, "S_phase1": 4, "n_epochs": 5, "top_share": 6}
H = {"S_early": 0, "f_early": 1, "S_late": 2, "f_late": 3, "n_late": 4, "iz_early": 5, "imu_early": 6}
_PRF: dict = {}



def response_table(P, band, fm_like, fwhm_arcsec, pix):
    key = (band, round(fwhm_arcsec, 2), round(pix, 3))
    if key not in _PRF:
        psf = P.make_psf(band, fm_like)
        fwhm_pix = fwhm_arcsec / pix
        k = _gaussian_kernel(fwhm_pix, int(np.ceil(2.5 * fwhm_pix)))
        r, v = radial_response_table(oversampled_template(psf), k)
        _PRF[key] = (r * pix, v)
    return _PRF[key]


class _FmLike:
    def __init__(self, fwhm_arcsec):
        self.fwhm_arcsec = fwhm_arcsec


def node_stack(P, d, b, node, epoch_sel):
    eb = (d["band_idx"] == b) & epoch_sel
    if eb.sum() < 1:
        return np.nan, np.nan, np.nan, 0
    f = d["f"][0, eb][:, node[0], node[1], node[2]][:, None]
    v = d["v"][0, eb][:, node[0], node[1], node[2]][:, None]
    g = d["g"][0, eb][:, node[0], node[1], node[2]][:, None]
    S, A, B, n, _ = nulls.stack_S(f, v, g, return_parts=True, min_epochs=1, frame_cap=d["frame_cap"][eb],
                                  clip_sigma=P.clip_sigma)
    if B[0] <= 0:
        return np.nan, np.nan, np.nan, int(n[0])
    return float(S[0]), float(A[0] / B[0]), float(1 / np.sqrt(B[0])), int(n[0])


def adjudicate(P, key: str, cinfo: dict, mask: str, catalog, completeness: dict, null_cal: dict) -> Candidate:
    endpoint, role, band = key.split("/")
    b = P.band_idx[band]
    mi = MASKS.index(mask)
    # cross-track variants: adjudicate the offset with the largest S_max
    best_path, best_S = P.tensor_dir / f"{endpoint}__{role}.npz", -np.inf
    for tp in [best_path] + sorted(P.tensor_dir.glob(f"{endpoint}__{role}__xt*.npz")):
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
    fcap = d["frame_cap"][eb]
    S, A, B, n, w = nulls.stack_S(f0, v0, g0, ok, return_parts=True, frame_cap=fcap, clip_sigma=P.clip_sigma)
    s_max, node = nulls.grid_max(S)
    iz, im1, im2 = node
    z = float(P.z_grid[iz]); mu = (float(P.mu_grid[im1]), float(P.mu_grid[im2]))
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
        "n_epochs": int(n[node]),
        "grid_edge": {"z": bool(iz in (0, len(P.z_grid) - 1)),
                      "mu": bool(im1 in (0, len(P.mu_grid) - 1) or im2 in (0, len(P.mu_grid) - 1))},
        "fitted": {"z_au": z, "mu_arcsec_yr": mu, "S_max": float(s_max),
                   "xt_arcsec": float(d["xt_arcsec"]) if "xt_arcsec" in d.files else 0.0},
    }
    # other-band significance at the same node (annotation)
    annotations["other_bands"] = {}
    for ob, obi in P.band_idx.items():
        if obi == b:
            continue
        S_o, f_o, sf_o, n_o = node_stack(P, d, obi, node, ok_all)
        annotations["other_bands"][ob] = {"S": S_o, "n": n_o}
    # --- calibrated tests -----------------------------------------------------
    rejection = None
    tests = {}
    # a. flux-consistent static source / halo
    centre = tuple(d["centre"]); cosd = np.cos(np.deg2rad(centre[1]))
    dt_yr = (mjd - P.t0_mjd) / 365.25
    d0 = d["d0"][eb]
    tr_ra = centre[0] + (d0[:, iz, 0] + mu[0] * dt_yr) / 3600.0 / cosd
    tr_dec = centre[1] + (d0[:, iz, 1] + mu[1] * dt_yr) / 3600.0
    catalog = catalog.get(band) if isinstance(catalog, dict) else catalog
    if catalog is not None and len(catalog.ra):
        pix = float(np.nanmedian(d["pix_arcsec"][eb]))
        rt = response_table(P, band, _FmLike(float(np.nanmedian(d["fwhm_arcsec"][eb]))),
                            float(np.nanmedian(d["fwhm_arcsec"][eb])), pix)
        esc = P.static_epoch_scale(d, np.where(eb)[0], node) if P.static_epoch_scale else None
        fc = flux_consistency(catalog, tr_ra, tr_dec, we * ok, phase, P.zp_ref, rt[0], rt[1],
                              float(s_max), S_by_phase=S_ph, search_arcsec=P.static_search_arcsec,
                              factor=P.static_factor, min_ndet=P.catalog_min_ndet, epoch_scale=esc)
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
                                   P.holdout_min_s_late, P.holdout_min_flux_ratio)
    hres["calibration"] = {"null_false_pass_rate": null_cal.get("false_pass_rate"),
                           "null_false_pass_ci68": null_cal.get("false_pass_ci68")}
    # v2.1: an annotation with its measured pass rates, never a rejection
    annotations["holdout_prediction"] = hres
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
        candidate_id=stable_id("cnd", {"cell": key, "hypothesis": P.hypothesis_version,
                                       "freeze": P.freeze_hash(), "mask": mask}),
        analysis_run_id=cinfo.get("analysis_run_id", "see-null-ensemble"),
        endpoint_id=endpoint, role=role, observation_ids=(),
        fitted_z_au=z, fitted_residual_motion={"mu_ra": mu[0], "mu_dec": mu[1]},
        model_comparison={"S_max": float(s_max), "T": cinfo["T"], "R": cinfo["R"], "q95": cinfo["q95"],
                          "R_norm": cinfo["R_norm"], "tests": tests},
        status=status, veto_reason=rejection, rejection_test=rejection,
        global_p_value=cinfo["global_p"], rank_statement=cinfo["rank"],
        annotations=annotations,
        extra={"band": band, "set": cinfo.get("set"), "mask": mask, "selection_function_loss": loss})


def run(P, set_name: str, mask: str = "primary") -> None:
    ne_all = json.loads((P.null_dir / f"{set_name}_null_ensemble.json").read_text())
    ne = ne_all[mask]
    comp_path = P.run_dir / "completeness" / f"{set_name}_completeness.json"
    completeness = json.loads(comp_path.read_text())["cells"] if comp_path.exists() else {}
    cands, cats = [], {}
    for key in ne["candidates"]:
        cinfo = dict(ne["cells"][key]); cinfo["set"] = set_name
        endpoint, role = key.split("/")[:2]
        with np.load(P.tensor_dir / f"{endpoint}__{role}.npz") as d:
            corridor = str(d["corridor"]); centre = tuple(d["centre"])
        if corridor not in cats:
            cats[corridor] = P.catalog(corridor, centre) if P.catalog else None
        c = adjudicate(P, key, cinfo, mask, cats[corridor], completeness, ne["holdout_null_calibration"])
        cands.append(c)
        print(f"{key:26s} R~={cinfo['R_norm']:.2f} p_global={cinfo['global_p']:.3f} -> {c.status}  "
              f"phase={c.annotations['parallax_phase']['signature']}", flush=True)
    out = P.run_dir / "records" / f"candidate_{set_name}.jsonl"
    if out.exists():
        out.unlink()
    append_records(out, cands)
    summary = {"survey": P.name, "set": set_name, "mask": mask, "R_fwer_norm": ne["fwer"]["R_fwer"],
               "n_candidates": len(cands),
               "status_counts": {s: sum(1 for c in cands if c.status == s) for s in set(c.status for c in cands)},
               "retained_ambiguous": [f"{c.endpoint_id}/{c.role}/{c.extra['band']}" for c in cands if c.status == "retained-ambiguous"],
               "written": str(out), "utc": datetime.now(timezone.utc).isoformat()}
    (P.run_dir / "records" / f"adjudication_{set_name}.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))
