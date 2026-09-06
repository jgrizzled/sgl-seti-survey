"""Post-blind adjudication analysis (freeze §6 ladder, chord-shape and
background tests): for every exceeding unit of a search, recompute the
per-event source excess under (a) the frozen same-frame-median
differential and (b) a curvature-corrected differential — a quadratic
in orbit-latitude offset fitted to the 8 control patches on each frame
and evaluated at the source — with split halves. A unit whose raw z is
positive in every event, correlated with |orbit latitude| of the source,
and removed by the curvature term is adjudicated 'ridge-curvature
systematic' (the HI-1 finding H1 in WISPR's orbit-latitude coordinate).
This is a DIAGNOSTIC on the blind result, not a re-search: the frozen
statistics and thresholds stand.

    python adjudicate_ridge.py confirmatory v1
"""
from __future__ import annotations
import json, sys
from collections import defaultdict
from pathlib import Path
import numpy as np
from astropy.time import Time

sys.path.insert(0, str(Path(__file__).resolve().parent))
import psp_geometry as G
import series as S
import wispr_lib as W


def main(role: str, version: str = "v1"):
    res = json.loads((S.SURV / "results" / f"{role}_search_{version}.json").read_text())
    idx = json.loads((S.SER / f"index_{role}_v1.json").read_text())
    patches, evmeta = idx["patches"], idx["events"]
    ccoef = res["colour_system"]["coeff"]
    ra_, rb_ = res["template_response"]["a"], res["template_response"]["b"]
    kq = res["bright_star_class"]["per_q_group"]
    qgroup = {}
    for grp in S.CFG["era"]["q_groups_rsun"]:
        a_, b_ = grp.split("-")
        for n_ in range(int(a_[1:]), int(b_[1:]) + 1):
            qgroup[f"E{n_:02d}"] = grp
    exc_units = [u for u, r in res["units"].items() if r["status"] == "searched"
                 and any(r[s]["exceeds"] for s in ("S_stack", "S_event", "S_pulse"))]
    # gather epochs for those units only (last record per frame)
    last = {}
    with open(S.SER / f"measurements_{role}_v1.jsonl") as fh:
        for line in fh:
            try:
                rec = json.loads(line)
            except Exception:
                continue
            last[rec["frame"]] = rec
    gp = S.CFG["gates"]["epoch"]
    epochs = defaultdict(dict)
    for rec in last.values():
        if "unusable" in rec or "error" in rec:
            continue
        for m in rec["records"]:
            if m["u"] in exc_units and m["k"] == "arc" and m["valid"] >= gp["aperture_finite_frac_min"] and np.isfinite(m["flux"] or np.nan):
                epochs[(m["u"], m["e"], round(m["mjd"], 5))][m["p"]] = m["flux"]
    # vetoes as in reduce (proximity)
    mjds = np.array(sorted({k[2] for k in epochs}))
    pdirs = G.planet_directions(mjds) if len(mjds) else {}
    lim_prox = 2 * np.sin(np.radians(gp["body_proximity_px"] * W.PIX_ARCSEC / 3600) / 2)
    bright = {b: np.radians(r) for b, r in gp["bright_body_proximity_deg"].items()}
    excl = {e: [Time(x["utc"]).mjd for x in xs] for e, xs in S.CFG["recon_exclusions"].items()}
    out = {}
    for u in exc_units:
        r = res["units"][u]
        rows = []
        for key, pt in patches.items():
            if not key.startswith(u + "|"):
                continue
            eid = key.split("|")[1]
            ev = r["events"].get(eid, {})
            if ev.get("status") != "included":
                continue
            offs = np.array(pt["ladder"], float)
            T = {pi: ra_ + rb_ * W.template_flux(pt["stars"][pi], ccoef) for pi in range(9)}
            if u in kq:
                T[0] *= kq[u].get(qgroup.get(evmeta[key]["encounter"]), 1.0)
            pv = np.array([G.radec_to_vec(a, b) for a, b in pt["radec"]])
            Draw, Dq, mj = [], [], []
            for (uu, ee, m), by_p in epochs.items():
                if uu != u or ee != eid or 0 not in by_p:
                    continue
                if eid in excl and any(abs(m - x) < 0.003 for x in excl[eid]):
                    continue
                i = min(int(np.searchsorted(mjds, m)), len(mjds) - 1)
                bad = False
                for body, dirs in pdirs.items():
                    if np.any(np.linalg.norm(pv - dirs[i], axis=1) < lim_prox):
                        bad = True
                        break
                    if body in bright and np.any(np.arccos(np.clip(pv @ dirs[i], -1, 1)) < bright[body]):
                        bad = True
                        break
                if bad:
                    continue
                cs = [pi for pi in range(1, 9) if pi in by_p]
                if len(cs) < 5:
                    continue
                E = {pi: by_p[pi] - T[pi] for pi in [0] + cs}
                ec = np.array([E[c] for c in cs])
                oc = np.array([offs[c - 1] for c in cs])
                Draw.append(E[0] - float(np.median(ec)))
                A = np.c_[np.ones(len(cs)), oc, oc ** 2]
                coef, *_ = np.linalg.lstsq(A, ec, rcond=None)
                Dq.append(E[0] - float(coef[0]))
                mj.append(m)
            if len(Draw) < 10:
                continue
            Draw, Dq = np.array(Draw), np.array(Dq)
            sd_r, sd_q = S._robust_sd(Draw), S._robust_sd(Dq)
            n = len(Draw)
            h = n // 2
            rows.append({"event": eid, "encounter": evmeta[key]["encounter"], "n": n,
                         "bet_src": round(pt["bet0"], 2), "lam_src": round(pt["lam0"], 2),
                         "z_raw": round(float(np.mean(Draw) / (sd_r / np.sqrt(n))), 2),
                         "z_quad": round(float(np.mean(Dq) / (sd_q / np.sqrt(n))), 2) if sd_q > 0 else None,
                         "split_raw": [round(float(np.mean(Draw[:h]) / (sd_r / np.sqrt(h))), 1), round(float(np.mean(Draw[h:]) / (sd_r / np.sqrt(n - h))), 1)],
                         "mean_raw": round(float(np.mean(Draw)), 3), "mean_quad": round(float(np.mean(Dq)), 3),
                         "sd_raw": round(sd_r, 3)})
        zr = np.array([x["z_raw"] for x in rows]); zq = np.array([x["z_quad"] for x in rows if x["z_quad"] is not None])
        out[u] = {"exceeds": {s: r[s]["exceeds"] for s in ("S_stack", "S_event", "S_pulse")},
                  "S": {s: (r[s]["S"], r[s]["T"]) for s in ("S_stack", "S_event", "S_pulse")},
                  "n_events": len(rows), "z_raw_median": round(float(np.median(zr)), 2) if len(zr) else None,
                  "z_quad_median": round(float(np.median(zq)), 2) if len(zq) else None,
                  "S_stack_raw": round(float(np.sum(zr) / np.sqrt(len(zr))), 2) if len(zr) else None,
                  "S_stack_quad": round(float(np.sum(zq) / np.sqrt(len(zq))), 2) if len(zq) else None,
                  "frac_events_positive_raw": round(float(np.mean(zr > 0)), 2) if len(zr) else None,
                  "bet_src_median": round(float(np.median([x["bet_src"] for x in rows])), 2) if rows else None,
                  "events": rows}
        print(u, {k: v for k, v in out[u].items() if k != "events"}, flush=True)
    (S.SURV / "results" / f"adjudication_ridge_{role}_{version}.json").write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main(*sys.argv[1:3])
