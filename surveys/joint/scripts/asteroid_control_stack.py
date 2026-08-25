"""Stack-regime asteroid positive control — scoring (notes/learnings.md
§10 item 1; joint hypotheses v3.0 dev-stage validation).

Scores the ring-grid samples fetched by
surveys/{ztf,panstarrs}/scripts/asteroid_stack_fetch.py with the frozen
v2 rule — per-frame cap (20x band-median frame weight), single-epoch
clip |S_e| <= 5 ACTIVE, 48-offset ring null, R~ = R/q95 — for each
archive alone and for the joint g/r/i family combination
S_joint = (A_ztf + A_ps1)/sqrt(B_ztf + B_ps1), i.e. the exact statistic
of surveys/joint/joint.py with the Horizons ephemeris playing the role
of the trajectory model.

The point of the control: every frame's single-frame S/N is below the
clip (verified and reported), so the clip removes nothing and only the
trajectory-weighted stack can recover the source — the regime the SGL
search operates in. Contrast the catalogue-layer (60000) control
(runs/ztf/v2/control/), which the clip demotes by construction.

Outputs: runs/ztf/v2/control/<n>_stack/summary.json,
runs/panstarrs/v2/control/<n>_stack/summary.json,
runs/joint/v3/control/<n>_stack/summary.json.

Usage: uv run python surveys/joint/scripts/asteroid_control_stack.py
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "surveys" / "ztf"))
import importlib.util  # noqa: E402


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


ZTF = _load("ztf_profile_c", REPO / "surveys" / "ztf" / "profile.py").PROFILE
PS1 = _load("ps1_profile_c", REPO / "surveys" / "panstarrs" / "profile.py").PROFILE

from sglsurvey import nulls  # noqa: E402

SAMPLES = {"ztf": REPO / "runs" / "ztf" / "v2" / "control_stack" / "samples_ring.npz",
           "ps1": REPO / "runs" / "panstarrs" / "v2" / "control_stack" / "samples_ring.npz"}
OUT = {"ztf": REPO / "runs" / "ztf" / "v2" / "control",
       "ps1": REPO / "runs" / "panstarrs" / "v2" / "control",
       "joint": REPO / "runs" / "joint" / "v3" / "control"}
#: solar colours vs V (v1/v2 control convention)
COLOUR = {"g": 0.25, "r": -0.19, "i": -0.35, "z": -0.44, "y": -0.47,
          "zg": 0.25, "zr": -0.19, "zi": -0.35}
FAMILIES = {"g": {"ps1": "g", "ztf": "zg"}, "r": {"ps1": "r", "ztf": "zr"}, "i": {"ps1": "i", "ztf": "zi"}}
CLIP = 5.0
RESID_CENTRE = (2, 2)


def band_parts(d, band, P):
    """Per-trajectory (S, A, B, n) with the archive's frozen cap and
    clip, plus per-band bookkeeping, for one archive's samples. Bands
    with fewer than P.min_epochs frames still return their stack sums
    (for the joint family combination, where the OTHER archive supplies
    the epochs) but are flagged unusable for standalone scoring."""
    sel = np.asarray(d["band"]) == band
    if sel.sum() < 1:
        return None
    f, v, g = d["f"][sel], d["v"][sel], d["g"][sel]      # (E, 49, 5, 5)
    wf = d["frame_w"][sel]
    cap = P.weight_cap * float(np.median(wf))
    fc = np.full(int(sel.sum()), cap, np.float32)
    nt = f.shape[1]
    S = np.empty((nt, 5, 5)); A = np.empty_like(S); B = np.empty_like(S); N = np.empty_like(S)
    for t in range(nt):
        St, At, Bt, nt_, _ = nulls.stack_S(f[:, t], v[:, t], g[:, t], return_parts=True,
                                           frame_cap=fc, clip_sigma=CLIP, min_epochs=P.min_epochs)
        S[t], A[t], B[t], N[t] = St, At, Bt, nt_
    vp = np.asarray(d["v_pred"])[sel]
    snr = np.asarray(d["single_snr"])[sel]
    pred = float(-2.5 * np.log10(np.mean(10 ** (-0.4 * vp))) + COLOUR.get(band, 0.0))
    return {"S": S, "A": A, "B": B, "n_grid": N, "n": int(sel.sum()), "pred": pred,
            "scoreable": bool(sel.sum() >= P.min_epochs),
            "snr_median": float(np.median(snr)), "snr_max": float(snr.max()),
            "n_above_clip": int((snr > CLIP).sum())}


def score(S, A, B, designated, norm_quantile, zp_ref):
    mx = np.array([nulls.grid_max(S[t])[0] for t in range(S.shape[0])])
    R, T = nulls.exceedance_ratios(mx, designated)
    s_max, node = nulls.grid_max(S[0])
    ring_R = R[1:]
    q95 = float(np.nanquantile(ring_R, norm_quantile))
    flux = A[0][node] / B[0][node] if B[0][node] > 0 else np.nan
    mag = zp_ref - 2.5 * np.log10(flux) if flux > 0 else None
    return {"S_max_real": float(s_max), "S_centre": float(S[0][RESID_CENTRE]),
            "T_designated": float(T), "R_real": float(R[0]),
            "ring_R_median": float(np.nanmedian(ring_R)), "ring_R_q95": q95,
            "R_norm_real": float(R[0] / q95) if q95 > 0 else None,
            "rank_statement": nulls.rank_statement(R[0], ring_R),
            "peak_offset_arcsec": [float(x - 2.0) for x in node],
            "recovered_stack_mag": mag, "detected": bool(R[0] > 1.0)}


def main():
    data, meta = {}, {}
    for k, p in SAMPLES.items():
        if not p.exists():
            raise SystemExit(f"missing {p} — run the {k} asteroid_stack_fetch first")
        data[k] = np.load(p)
        meta[k] = json.loads(str(data[k]["meta"]))
        meta[k]["samples_sha256"] = "sha256:" + hashlib.sha256(p.read_bytes()).hexdigest()
    assert meta["ztf"]["asteroid"] == meta["ps1"]["asteroid"]
    label = meta["ztf"]["asteroid"]
    designated_z = np.array([1 + i for i in ZTF.designated])
    parts = {"ztf": {}, "ps1": {}}

    for k, P in (("ztf", ZTF), ("ps1", PS1)):
        summary = {"asteroid": label, "kind": "stack-regime positive control",
                   "rule": f"v2 ring-48 / per-frame cap / clip |S_e| <= {CLIP} ACTIVE",
                   "source": meta[k], "colours_vs_V": COLOUR}
        for band in sorted(set(np.asarray(data[k]["band"]).tolist())):
            bp = band_parts(data[k], band, P)
            if bp is None:
                continue
            parts[k][band] = bp
            if not bp["scoreable"]:
                summary[band] = {"n_epochs": bp["n"], "status": "insufficient_epochs_for_standalone_score",
                                 "single_frame_snr_median": bp["snr_median"],
                                 "single_frame_snr_max": bp["snr_max"], "frames_above_clip": bp["n_above_clip"],
                                 "note": "stack sums still enter the joint family combination"}
                print(f"[{k} {band}] N={bp['n']} < min_epochs — joint-only", flush=True)
                continue
            sc = score(bp["S"], bp["A"], bp["B"], designated_z, P.norm_quantile, P.zp_ref)
            summary[band] = {**sc, "n_epochs": bp["n"], "predicted_mag": bp["pred"],
                             "throughput_dmag": (sc["recovered_stack_mag"] - bp["pred"])
                             if sc["recovered_stack_mag"] is not None else None,
                             "single_frame_snr_median": bp["snr_median"],
                             "single_frame_snr_max": bp["snr_max"],
                             "frames_above_clip": bp["n_above_clip"],
                             "stack_regime": bp["n_above_clip"] == 0}
            print(f"[{k} {band}] N={bp['n']} S={sc['S_max_real']:.1f} R={sc['R_real']:.2f} "
                  f"R~={sc['R_norm_real']:.2f} mag {sc['recovered_stack_mag'] and round(sc['recovered_stack_mag'], 2)} "
                  f"vs pred {bp['pred']:.2f}; single S/N med {bp['snr_median']:.2f} max {bp['snr_max']:.2f} "
                  f"(>clip: {bp['n_above_clip']})", flush=True)
        out = OUT[k] / f"{label}_stack"
        out.mkdir(parents=True, exist_ok=True)
        (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float))

    # joint family combination: S = (A_ztf + A_ps1)/sqrt(B_ztf + B_ps1),
    # each archive's sums under its own frozen cap and clip
    jsummary = {"asteroid": label, "kind": "stack-regime positive control (joint g/r/i families)",
                "rule": "joint S = (A_ztf + A_ps1)/sqrt(B_ztf + B_ps1); per-archive caps and clip",
                "source": {k: meta[k] for k in meta}}
    for fam, bmap in FAMILIES.items():
        bz, bp_ = parts["ztf"].get(bmap["ztf"]), parts["ps1"].get(bmap["ps1"])
        if bz is None and bp_ is None:
            continue
        use = [x for x in (bz, bp_) if x is not None]
        A = sum(x["A"] for x in use); B = sum(x["B"] for x in use)
        N = sum(x["n_grid"] for x in use)
        with np.errstate(all="ignore"):
            S = A / np.sqrt(np.where(B > 0, B, np.inf))
        S = np.where(N >= 5, S, np.nan)
        sc = score(S, A, B, designated_z, ZTF.norm_quantile, 25.0)
        fz = sum(10 ** (-0.4 * x["pred"]) for x in use) / len(use)
        jsummary[fam] = {**sc, "archives": {"ztf": bz is not None and bz["n"], "ps1": bp_ is not None and bp_["n"]},
                         "predicted_mag_mean_flux": float(-2.5 * np.log10(fz)),
                         "throughput_dmag": (sc["recovered_stack_mag"] - float(-2.5 * np.log10(fz)))
                         if sc["recovered_stack_mag"] is not None else None}
        print(f"[joint {fam}] S={sc['S_max_real']:.1f} R={sc['R_real']:.2f} R~={sc['R_norm_real']:.2f} "
              f"mag {sc['recovered_stack_mag'] and round(sc['recovered_stack_mag'], 2)}", flush=True)
    out = OUT["joint"] / f"{label}_stack"
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.json").write_text(json.dumps(jsummary, indent=2, default=float))


if __name__ == "__main__":
    main()
