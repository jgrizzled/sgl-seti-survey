"""Blind confirmatory run v1 (hypotheses.md D4/D6, amendment v1.1): the 29
confirmatory units x 5 statistics at the frozen T, with the frozen veto
ladder applied to every exceedance (4FGL association, strict re-run with
doubled Moon/Sun exclusion + zenith 90 + theta 50, attitude/data-quality
census, recurrence). Writes results/confirmatory_v1.json / .md and the
per-unit analyses under results/units_v1/."""
from __future__ import annotations

import json
from collections import Counter

import numpy as np
from astropy.time import Time

import hlib as H
import recon_scan as R

FGL_VAR_THRESHOLD = 18.48


def fgl_sources(store_dir):
    """4FGL-DR4 sources within 3 deg of every channel position (HEASARC TAP)."""
    cache = H.RUN / "fgl_dr4_near_positions.json"
    if cache.exists():
        return json.loads(cache.read_text())
    from sglsurvey.snapshots import SnapshotStore
    store = SnapshotStore(H.RUN / "catalogs")
    pos = R.positions(R.load_events())
    out = {}
    for (ch, tid), (ra, de, n) in sorted(pos.items()):
        adql = ("SELECT name, ra, dec, semi_major_axis_95, flux_1_100_gev, variability_index, frac_variability, assoc_name "
                f"FROM fermilpsc WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra:.5f}, {de:.5f}, 3.0))=1")
        rows, err = R.heasarc(store, adql, f"4fgl {ch} {tid}")
        for r in rows or []:
            r["_sep"] = R.sep_deg(R.ffloat(r["ra"]), R.ffloat(r["dec"]), ra, de)
        out[f"{ch}-{tid}"] = rows or []
    cache.write_text(json.dumps(out, indent=1))
    return out


def adjudicate(u, a, stat, ph, ws, fgl):
    """Frozen veto ladder for one exceedance; returns a disposition dict."""
    T = H.FREEZE["T"]
    lane = stat[-1] if stat.endswith(("_L", "_H")) else "U"
    win_utc = a["diag"].get(f"{stat}_window") if stat.startswith("S_event") else (
        a["diag"]["S_burst_rec"]["window"] if stat == "S_burst" and a["diag"]["S_burst_rec"] else None)
    rec = {"unit": f"{u['channel']} {u['target']} {u['rung']}", "statistic": stat, "value": a["stats"][stat], "window": win_utc, "steps": []}
    # (1) 4FGL association: sources inside the aperture; variable ones
    r_ap = H.FREEZE["lanes"][lane]["r_ap"]
    near = [s for s in fgl.get(u["key"], []) if s["_sep"] <= r_ap]
    var = [s for s in near if (R.ffloat(s.get("variability_index")) or 0) > FGL_VAR_THRESHOLD]
    rec["steps"].append({"step": "4FGL", "in_aperture": [(s["name"], round(s["_sep"], 2), s.get("variability_index")) for s in near],
                         "variable": [s["name"] for s in var]})
    # excess-photon centroid for event/burst statistics
    if win_utc:
        w = [x for x in a["windows"] if x["t_ca_utc"] == win_utc]
        if w:
            w = w[0]
            n, t = H.count(ph, w["mjd0"], w["mjd1"], lane)
            lo = np.searchsorted(ph["mjd"], w["mjd0"]); hi = np.searchsorted(ph["mjd"], w["mjd1"])
            e = ph["ENERGY"][lo:hi]; s = ph["sep"][lo:hi]
            L = H.FREEZE["lanes"][lane]
            m = (e >= L["emin"]) & (e < L["emax"]) & (s <= L["r_ap"])
            if m.sum():
                cra = float(np.mean(ph["RA"][lo:hi][m])); cde = float(np.mean(ph["DEC"][lo:hi][m]))
                d_fgl = [(x["name"], round(R.sep_deg(R.ffloat(x["ra"]), R.ffloat(x["dec"]), cra, cde), 2)) for x in fgl.get(u["key"], [])]
                d_fgl = sorted(d_fgl, key=lambda x: x[1])[:3]
                rec["steps"].append({"step": "centroid", "n": int(m.sum()), "ra": cra, "dec": cde, "nearest_4FGL": d_fgl,
                                     "photon_mjd": [float(x) for x in ph["mjd"][lo:hi][m]], "photon_E": [float(x) for x in e[m]]})
    # (2) strict re-run
    a2 = H.unit_analysis(u["channel"], u["target"], u["rung"], ph, u["ra"], u["dec"], ws, strict=True,
                         exclude_event_ids=[w["event_id"] for w in ws if H.is_forced_dev(u["channel"], u["target"], w["t_ca_utc"])])
    rec["steps"].append({"step": "strict", "value": a2["stats"].get(stat), "survives": bool(a2["stats"].get(stat, 0) > T)})
    # (4) attitude / data quality inside the window
    if win_utc and w:
        iv = H.intervals(u["ra"], u["dec"], w["mjd0"], w["mjd1"])
        rec["steps"].append({"step": "attitude", "modes": Counter(iv["mode"][iv["gate"]].tolist()) if iv is not None else {},
                             "rock_range": [float(iv["rock"][iv["gate"]].min()), float(iv["rock"][iv["gate"]].max())] if iv is not None and iv["gate"].any() else None,
                             "theta_min": float(iv["theta"][iv["gate"]].min()) if iv is not None and iv["gate"].any() else None})
    # (5) recurrence: other windows' per-window S in this lane
    if stat.startswith("S_event"):
        others = sorted([(H.nlog10_sf(x["n_" + lane], x["lam_" + lane]), x["t_ca_utc"][:10]) for x in a["windows"]
                         if x["searchable"] and x["t_ca_utc"] != win_utc], reverse=True)[:3]
        rec["steps"].append({"step": "recurrence", "next_best_windows": others})
    # disposition
    if var or (rec["steps"][1]["step"] == "centroid" and rec["steps"][1]["nearest_4FGL"] and rec["steps"][1]["nearest_4FGL"][0][1] <= (0.5 if lane == "H" else 1.5)):
        rec["disposition"] = "vetoed_known_source"
    elif not rec["steps"][-2 if stat.startswith("S_event") else -1]["survives"] if any(s["step"] == "strict" for s in rec["steps"]) else False:
        rec["disposition"] = "vetoed_strict_gate"
    else:
        rec["disposition"] = "retained_ambiguous"
    return rec


def main():
    events = R.load_events()
    T = H.FREEZE["T"]
    thr = json.loads((H.RES / "thresholds_v1.json").read_text())
    status = {(t["unit"], t["statistic"]): t["status"] for t in thr["trials"]}
    fgl = fgl_sources(None)
    (H.RES / "units_v1").mkdir(exist_ok=True)
    out = {"utc": Time.now().isot, "T": T, "n_trials": 0, "units": [], "exceedances": [], "expected_exceedances": None}
    md = ["# Blind confirmatory run v1 — 2026-09-07", "", f"T = {T:.3f}; 29 units × 5 statistics.", "",
          "| Unit | n_win | S_stack_L | S_event_L | S_stack_H | S_event_H | S_burst | n_L / λ_L | n_H / λ_H |", "|---|---|---|---|---|---|---|---|---|"]
    n_trials = 0; n_exc = 0
    for u in H.all_units(events):
        if u["dev"]:
            continue
        ph = H.photons(u["key"], u["ra"], u["dec"])
        ws = H.unit_windows(u["channel"], u["target"], u["rung"], events)
        excl = [w["event_id"] for w in ws if H.is_forced_dev(u["channel"], u["target"], w["t_ca_utc"])]
        a = H.unit_analysis(u["channel"], u["target"], u["rung"], ph, u["ra"], u["dec"], ws, exclude_event_ids=excl)
        name = f"{u['channel']} {u['target']} {u['rung']}"
        rec = {"unit": name, "n_searchable": a["n_searchable"], "stats": a["stats"], "diag": {k: v for k, v in a["diag"].items() if not k.startswith("stack_ens")},
               "trials": {}}
        for stat, val in a["stats"].items():
            st = status.get((name, stat), "calibrated")
            n_trials += 1
            exc = (val is not None) and (not np.isnan(val)) and val > T
            rec["trials"][stat] = {"value": val, "status": st, "exceeds": bool(exc)}
            if exc:
                n_exc += 1
                adj = adjudicate(u, a, stat, ph, ws, fgl)
                adj["trial_status"] = st
                out["exceedances"].append(adj)
                print(f"[conf] EXCEEDANCE {name} {stat} = {val:.2f} -> {adj['disposition']}", flush=True)
        out["units"].append(rec)
        (H.RES / "units_v1" / f"{u['channel']}_{u['target']}_{u['rung']}.json").write_text(json.dumps(H.clean(a), indent=1))
        d = a["diag"]
        md.append(f"| {name} | {a['n_searchable']} | {a['stats']['S_stack_L']:.2f} | {a['stats']['S_event_L']:.2f} | {a['stats']['S_stack_H']:.2f} | "
                  f"{a['stats']['S_event_H']:.2f} | {a['stats']['S_burst']:.2f} | {d['n_tot_L']} / {d['lam_tot_L']:.1f} | {d['n_tot_H']} / {d['lam_tot_H']:.2f} |")
        print(f"[conf] {name:28s} " + " ".join(f"{k}={v:.2f}" for k, v in a["stats"].items()), flush=True)
    out["n_trials"] = n_trials
    out["expected_exceedances"] = n_trials * H.FREEZE["alpha_trial"]
    out["n_exceedances"] = n_exc
    disp = Counter(e["disposition"] for e in out["exceedances"])
    md += ["", f"Trials {n_trials}; exceedances {n_exc} vs {out['expected_exceedances']:.3f} expected (FWER 0.05). Dispositions: {dict(disp)}", ""]
    for e in out["exceedances"]:
        md.append(f"- **{e['unit']} {e['statistic']} = {e['value']:.2f}** (window {e['window']}): {e['disposition']}; steps: {json.dumps(H.clean(e['steps']))[:600]}")
    (H.RES / "confirmatory_v1.json").write_text(json.dumps(H.clean(out), indent=1))
    (H.RES / "confirmatory_v1.md").write_text("\n".join(md))
    print(f"trials {n_trials} exceedances {n_exc} expected {out['expected_exceedances']:.3f} {dict(disp)}")


if __name__ == "__main__":
    main()
