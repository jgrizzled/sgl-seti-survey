"""Threshold freeze v1 (hypotheses.md D4/D6): the analytic Poisson threshold
T = -log10(alpha_trial) for 150 trials at FWER 0.05, validated per trial on
the pseudo-window ensembles (off-window photons only): the dispersion gate
(> 10 % of pseudo-windows with p < 0.05 -> constraint_only) and the
ensemble exceedance count at T. Real in-window statistics are computed
in memory by the shared routine but are NOT written here — only the
ensemble diagnostics are (the confirmatory stage re-computes everything).
Writes results/thresholds_v1.json and ../thresholds.md."""
from __future__ import annotations

import json
import hashlib

import numpy as np
from astropy.time import Time

import hlib as H
import recon_scan as R


def main():
    events = R.load_events()
    units = H.all_units(events)
    out = {"utc": Time.now().isot, "T": H.FREEZE["T"], "alpha_trial": H.FREEZE["alpha_trial"],
           "n_trials": H.FREEZE["n_trials"], "seed": H.FREEZE["seed"], "trials": [], "units": []}
    for u in units:
        ph = H.photons(u["key"], u["ra"], u["dec"])
        ws = H.unit_windows(u["channel"], u["target"], u["rung"], events)
        excl = [w["event_id"] for w in ws if H.is_forced_dev(u["channel"], u["target"], w["t_ca_utc"])]
        a = H.unit_analysis(u["channel"], u["target"], u["rung"], ph, u["ra"], u["dec"], ws, exclude_event_ids=excl)
        d = a["diag"]
        rec = {"unit": f"{u['channel']} {u['target']} {u['rung']}", "dev": u["dev"], "n_searchable": a["n_searchable"],
               "n_photons_gated": int(ph["n_gated"][0])}
        for lane in ("L", "H"):
            ens = np.array(d["stack_ens_" + lane]) if d["stack_ens_" + lane] else np.array([])
            frac = d["disp_frac_" + lane]
            gate = (frac is not None) and (not np.isnan(frac)) and frac <= H.FREEZE["dispersion_gate"][1]
            # amendment v1.2: per-trial ensemble exceedance gate — the pseudo-window
            # exceedances at T must not exceed max(3, 3 x expected)
            n_pw = d["n_pw_" + lane]; exc = d["pw_exceed_T_" + lane]
            exc_gate = exc <= max(3.0, 3.0 * n_pw * H.FREEZE["alpha_trial"])
            gate = gate and exc_gate
            rec[lane] = {"n_pw": d["n_pw_" + lane], "disp_frac": frac, "disp_ok": bool(gate), "exceedance_gate_ok": bool(exc_gate),
                         "pw_exceed_T": d["pw_exceed_T_" + lane],
                         "stack_ens_n": int(len(ens)), "stack_ens_max": float(ens.max()) if len(ens) else None,
                         "stack_ens_p95": float(np.percentile(ens, 95)) if len(ens) else None,
                         "stack_ens_exceed_T": int((ens > H.FREEZE["T"]).sum()) if len(ens) else 0}
            for stat in ("S_stack_" + lane, "S_event_" + lane):
                out["trials"].append({"unit": rec["unit"], "statistic": stat, "dev": u["dev"], "T": H.FREEZE["T"],
                                      "status": ("calibrated" if gate else "constraint_only")})
        out["trials"].append({"unit": rec["unit"], "statistic": "S_burst", "dev": u["dev"], "T": H.FREEZE["T"],
                              "status": "calibrated" if (rec["L"]["disp_ok"] and rec["H"]["disp_ok"]) else "constraint_only"})
        # per-window rates for the record (off-window quantities)
        rec["windows"] = [{"t_ca_utc": w["t_ca_utc"], "searchable": w["searchable"], "n_pw": w["n_pw"],
                           "rate_L": w.get("rate_L"), "rate_H": w.get("rate_H"), "rate_U": w.get("rate_U"),
                           "lam_L": w.get("lam_L"), "lam_H": w.get("lam_H"), "expo_L": w["expo_L"], "expo_H": w["expo_H"]}
                          for w in a["windows"]]
        out["units"].append(rec)
        print(f"[thr] {rec['unit']:28s} L disp {rec['L']['disp_frac']:.3f} ({rec['L']['n_pw']} pw, exc {rec['L']['pw_exceed_T']}, stack-ens max {rec['L']['stack_ens_max']}) "
              f"H disp {rec['H']['disp_frac']:.3f} (exc {rec['H']['pw_exceed_T']}, max {rec['H']['stack_ens_max']})", flush=True)
    conf = [t for t in out["trials"] if not t["dev"]]
    out["n_confirmatory_trials"] = len(conf)
    out["n_constraint_only"] = sum(1 for t in conf if t["status"] == "constraint_only")
    tot_pw = sum(r["L"]["n_pw"] + r["H"]["n_pw"] for r in out["units"])
    tot_exc = sum(r["L"]["pw_exceed_T"] + r["H"]["pw_exceed_T"] for r in out["units"])
    out["ensemble_exceedance"] = {"n_members": tot_pw, "n_exceed_T": tot_exc,
                                  "expected": tot_pw * H.FREEZE["alpha_trial"],
                                  "ok": tot_exc <= max(2.0 * tot_pw * H.FREEZE["alpha_trial"], 3)}
    txt = json.dumps(H.clean(out), indent=1)
    (H.RES / "thresholds_v1.json").write_text(txt)
    (H.RES / "threshold_freeze_v1_sha.txt").write_text(hashlib.sha256(txt.encode()).hexdigest() + "\n")
    print(f"confirmatory trials {len(conf)} (constraint_only {out['n_constraint_only']}); ensemble members {tot_pw}, "
          f"exceed T {tot_exc} vs expected {out['ensemble_exceedance']['expected']:.2f}")


if __name__ == "__main__":
    main()
