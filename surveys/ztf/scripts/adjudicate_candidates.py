"""Stage-7 adjudication of retained Candidate records (plan §3.4 items
5-6: competing-model tests and held-out-epoch prediction).

For every Candidate with status "retained" in calib_v1 we recompute the
stack at its peak cell (same weights, cap, clip as the calibration) and
apply, in order:

  A. few-epoch dominance: the top-3 epochs carry > 50 % of the weighted
     sum A (the v0.2 weight cap bounds weights, not flux outliers below
     the 5-sigma clip);
  B. static-position test: epochs with S_e > 1.5 cluster on the sky
     (rms < 3") while the predicted track spans > 10" — a fixed source
     crossed by the parallax ellipse, not a track follower;
  C. held-out consistency: odd/even epoch split and first/second-half
     time split must each keep S >= 2 in both halves (a persistent relay
     predicts every subset);
  D. difference-image regime check: if the cell is science-image
     dominated (diff weight < 0.3) a static star is the default
     explanation unless B passes.

  E. survey-wide trials: the look-elsewhere null (look_elsewhere.json,
     leave-one-out control S/T ratios over every pair-band search) gives
     a global per-draw p for the candidate's S/T; with N_trials pair-band
     searches the expected number of chance exceedances at that strength
     is N_trials x p. Vetoed when that expectation is >= 1.

Anything surviving A-E is left "retained" with the diagnostics attached
for manual inspection (thumbnails are the next step). Writes updated
Candidate records to calib_v1/records/candidate_adjudicated.jsonl and a
human-readable adjudication.md.

Usage: uv run python surveys/ztf/scripts/adjudicate_candidates.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from astropy.time import Time

from sglseti import Role, adaptive_locus, load_target_registry

from sglsurvey.geometry import GeometryContext
from sglsurvey.records import Candidate, append_records

REPO = Path(__file__).resolve().parents[3]
CAL = REPO / "runs" / "ztf" / "calib_v1"
MIN_GOOD_FRAC, WEIGHT_CAP, CLIP = 0.7, 20.0, 5.0
BAND_IDX = {"zg": 1, "zr": 2, "zi": 3}


def stack(fb, vb, gb, sel=None):
    valid = np.isfinite(fb) & np.isfinite(vb) & (vb > 0) & (gb >= MIN_GOOD_FRAC)
    with np.errstate(all="ignore"):
        valid &= np.abs(fb) / np.sqrt(np.where(vb > 0, vb, np.inf)) <= CLIP
    if sel is not None:
        valid &= sel
    w = np.where(valid, 1.0 / np.where(vb > 0, vb, 1.0), 0.0)
    if (w > 0).any():
        w = np.minimum(w, WEIGHT_CAP * np.median(w[w > 0]))
    A, B = (np.where(valid, fb, 0.0) * w).sum(), w.sum()
    S = A / np.sqrt(B) if B > 0 else np.nan
    return S, A, B, w, valid


def main() -> None:
    registry = load_target_registry(REPO / "registries" / "pilot_wise_2026.yaml")
    ctx = GeometryContext.ztf_default()
    cands = [json.loads(l) for l in (CAL / "records" / "candidate.jsonl").open()]
    retained = [c for c in cands if c["status"] == "retained"]
    print(f"{len(cands)} candidates, {len(retained)} retained")
    le = json.loads((CAL / "look_elsewhere.json").read_text())
    n_trials = le["summary"]["n_pair_bands"]
    out, lines = [], ["# Stage-7 adjudication of retained candidates", "",
                      f"Look-elsewhere null: {n_trials} pair-band searches; real "
                      f"exceedances {le['summary']['real_exceedances']} vs "
                      f"{le['summary']['expected_exceedances_from_controls']} expected "
                      f"from leave-one-out controls.", ""]
    for c in retained:
        e, r, band = c["endpoint_id"], c["role"], c["extra"]["band"]
        d = np.load(CAL / "tensors" / f"{e}__{r}.npz")
        z_grid, mu_grid, mjd = d["z_grid"], d["mu_grid"], d["mjd"]
        eb = d["band_idx"] == BAND_IDX[band]
        zi = int(np.argmin(np.abs(z_grid - c["fitted_z_au"])))
        mi = int(np.argmin(np.abs(mu_grid - c["model_comparison"]["mu_arcsec_yr"][0])))
        mj = int(np.argmin(np.abs(mu_grid - c["model_comparison"]["mu_arcsec_yr"][1])))
        fb = d["f"][0, eb, zi, mi, mj].astype(float)
        vb = d["v"][0, eb, zi, mi, mj].astype(float)
        gb = d["g"][0, eb, zi, mi, mj].astype(float)
        db = np.nan_to_num(d["dfrac"][0, eb, zi, mi, mj].astype(float))
        mb = mjd[eb]
        S, A, B, w, valid = stack(fb, vb, gb)
        contrib = np.where(valid, fb * w, 0.0)
        top3 = float(np.sort(contrib)[::-1][:3].sum() / A) if A > 0 else np.nan
        s_e = np.where(valid, fb * np.sqrt(w), 0.0)
        diff_w = float((db * w).sum() / B) if B > 0 else np.nan
        # predicted positions of this cell per epoch
        pos = []
        for m in mb:
            al = adaptive_locus(target=registry[e], role=Role(r),
                                observation_time=Time(float(m), format="mjd"),
                                observer=ctx.observer, relay_range=ctx.relay_range,
                                tolerance_arcsec=1.0, ephemeris=ctx.ephemeris,
                                model=ctx.model)
            zs = np.array([p.z_au for p in al.points])
            k = int(np.argmin(np.abs(zs - z_grid[zi])))
            pos.append((al.points[k].icrs_ra_deg, al.points[k].icrs_dec_deg))
        pos = np.array(pos)
        cosd = np.cos(np.deg2rad(pos[:, 1].mean()))
        sig = valid & (s_e > 1.5)
        def rms(p):
            return float(np.hypot((p[:, 0] - p[:, 0].mean()) * cosd,
                                  p[:, 1] - p[:, 1].mean()).std() * 3600) if len(p) > 1 else np.nan
        sky_rms, track_rms = rms(pos[sig]), rms(pos[valid])
        order = np.argsort(mb)
        odd = np.zeros(len(mb), bool); odd[order[::2]] = True
        half = np.zeros(len(mb), bool); half[order[:len(order) // 2]] = True
        S_odd, S_even = stack(fb, vb, gb, odd)[0], stack(fb, vb, gb, ~odd)[0]
        S_h1, S_h2 = stack(fb, vb, gb, half)[0], stack(fb, vb, gb, ~half)[0]
        diag = {"S_recomputed": round(float(S), 2), "n_valid": int(valid.sum()),
                "top3_fraction_of_A": round(top3, 2), "n_sig_epochs": int(sig.sum()),
                "sig_epoch_sky_rms_arcsec": None if np.isnan(sky_rms) else round(sky_rms, 2),
                "track_rms_arcsec": None if np.isnan(track_rms) else round(track_rms, 2),
                "S_odd_even": [round(float(S_odd), 2), round(float(S_even), 2)],
                "S_first_second_half": [round(float(S_h1), 2), round(float(S_h2), 2)],
                "diff_image_weight": round(diff_w, 2)}
        reasons = []
        if top3 > 0.5:
            reasons.append(f"A: top-3 epochs carry {top3:.0%} of the stack")
        if sig.sum() >= 2 and not np.isnan(sky_rms) and sky_rms < 3.0 and track_rms > 10.0:
            reasons.append(f"B: significant epochs cluster at one sky position (rms {sky_rms:.1f}\" vs track {track_rms:.0f}\")")
        if min(S_odd, S_even) < 2.0:
            reasons.append(f"C: odd/even split fails ({S_odd:.1f}, {S_even:.1f})")
        if min(S_h1, S_h2) < 2.0:
            reasons.append(f"C: first/second-half split fails ({S_h1:.1f}, {S_h2:.1f})")
        gp = le["summary"]["per_candidate_global_p"].get(f"{e}/{r}/{band}")
        if gp is not None:
            diag["global_p_per_draw"] = gp
            diag["expected_chance_exceedances_at_this_strength"] = round(gp * n_trials, 1)
            if gp * n_trials >= 1.0:
                reasons.append(f"E: not significant after survey-wide trials (global p={gp:.3f}, "
                               f"expected {gp * n_trials:.1f} chance exceedances this strong among {n_trials} searches)")
        if diff_w < 0.3 and not reasons:
            reasons.append(f"D: science-image regime (diff weight {diff_w:.2f}) with no track-following evidence beyond the stack")
        status = "vetoed" if reasons else "retained"
        out.append(Candidate(**{**c, "status": status,
                                "veto_reason": "; ".join(reasons) if reasons else None,
                                "extra": {**c["extra"], "adjudication": diag,
                                          "adjudication_stage": "stage-7 A-E"}}))
        lines.append(f"## {e} / {r} / {band} — **{status}**")
        lines.append(f"z = {c['fitted_z_au']:.0f} AU, µ = {c['model_comparison']['mu_arcsec_yr']}, "
                     f"S = {c['model_comparison']['real_max_S']:.2f} vs T = {c['model_comparison']['threshold_8_controls']:.2f}")
        lines.append("```\n" + json.dumps(diag, indent=1) + "\n```")
        if reasons:
            lines.append("Veto: " + "; ".join(reasons))
        lines.append("")
        print(f"{e:16s} {r} {band}: {status:8s} {'; '.join(reasons)}")
    append_records(CAL / "records" / "candidate_adjudicated.jsonl", out)
    (CAL / "adjudication.md").write_text("\n".join(lines))
    n_ret = sum(o.status == "retained" for o in out)
    print(f"\n{len(out)} adjudicated: {len(out) - n_ret} vetoed, {n_ret} retained")


if __name__ == "__main__":
    main()
