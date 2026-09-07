"""Completeness v1 (hypotheses.md D9): photon-count injections on the real
per-window quantities of every confirmatory unit.

  S_stack : persistent flux F (ph cm^-2 s^-1 in the lane) in every
            searchable window -> extra counts ~ Poisson(F * expo_k[cm^2 s] * c_lane);
            F90 = flux at which 90 % of 200 realizations give S_stack > T.
  S_event : the same in the best-exposure window only (and the median window).
  S_burst : a 1-ks pulse of fluence Phi (ph cm^-2, union band) at a random gated
            time of the best-exposure window -> photons added to the real list,
            burst statistic re-run; Phi90 likewise.

Flux -> power through the rung cone with the frozen convention
P = F_E * pi * b_rung^2 (F_E = F * <E>, Gamma = 2), and the Gamma +/- 0.3 and
+/-10 % exposure systematics as multiplicative factors on <E> and F.
Writes results/completeness_v1.json / .md."""
from __future__ import annotations

import json

import numpy as np
from astropy.time import Time
from scipy.stats import poisson

import hlib as H

N_REAL = 200
F_GRID = np.geomspace(1e-9, 1e-5, 41)     # ph cm^-2 s^-1
PHI_GRID = np.geomspace(1e-6, 1e0, 61)    # ph cm^-2 (1-ks pulse)


def f90(grid, frac, level=0.9):
    frac = np.asarray(frac)
    i = np.argmax(frac >= level)
    if frac[i] < level:
        return None
    if i == 0:
        return float(grid[0])
    # log interpolation between the bracketing grid points
    x0, x1 = np.log10(grid[i - 1]), np.log10(grid[i]); y0, y1 = frac[i - 1], frac[i]
    return float(10 ** (x0 + (level - y0) * (x1 - x0) / max(y1 - y0, 1e-9)))


def main():
    T = H.FREEZE["T"]
    rng = np.random.default_rng(H.FREEZE["seed"])
    cont = {lane: H.lane_containment(lane, 0.8) for lane in ("L", "H", "U")}
    out = {"utc": Time.now().isot, "T": T, "n_real": N_REAL, "containment": cont, "units": []}
    md = ["# Completeness v1 — 2026-09-07", "", f"T = {T:.3f}; {N_REAL} realizations per grid point; containment L {cont['L']:.2f}, H {cont['H']:.2f}, U {cont['U']:.2f}.", "",
          "| Unit | F90 stack L | F90 event L (best / median win) | F90 stack H | F90 event H (best) | Φ90 / Φ50 burst (1 ks, ph cm⁻²) | P90 stack L (W) | P90 event L (W) |",
          "|---|---|---|---|---|---|---|---|"]
    for f in sorted((H.RES / "units_v1").glob("*.json")):
        a = json.loads(f.read_text())
        name = f"{a['channel']} {a['target']} {a['rung']}"
        sw = [w for w in a["windows"] if w["searchable"] and w["n_pw"] >= 4]
        rec = {"unit": name, "rung": a["rung"], "n_searchable": len(sw)}
        if not sw:
            out["units"].append(rec); continue
        for lane in ("L", "H"):
            n = np.array([w["n_" + lane] for w in sw]); lam = np.array([w["lam_" + lane] for w in sw])
            expo_cm2s = np.array([w["expo_" + lane] for w in sw]) * 1e4
            # S_stack
            frac = []
            for F in F_GRID:
                mu = F * expo_cm2s * cont[lane]
                inj = rng.poisson(mu, size=(N_REAL, len(sw)))
                ntot = n.sum() + inj.sum(axis=1)
                s = -np.log10(np.maximum(poisson.sf(ntot - 1, lam.sum()), 1e-300))
                frac.append(float((s > T).mean()))
            rec["F90_stack_" + lane] = f90(F_GRID, frac); rec["frac_stack_" + lane] = frac
            # S_event at best-exposure and median-exposure windows
            order = np.argsort(expo_cm2s)
            for tag, k in (("best", order[-1]), ("median", order[len(order) // 2])):
                frac = []
                for F in F_GRID:
                    mu = F * expo_cm2s[k] * cont[lane]
                    inj = rng.poisson(mu, size=N_REAL)
                    s = -np.log10(np.maximum(poisson.sf(n[k] + inj - 1, lam[k]), 1e-300))
                    frac.append(float((s > T).mean()))
                rec[f"F90_event_{lane}_{tag}"] = f90(F_GRID, frac)
                rec[f"event_{lane}_{tag}_window"] = sw[k]["t_ca_utc"][:10]
                rec[f"event_{lane}_{tag}_expo_cm2s"] = float(expo_cm2s[k])
            rec["expo_stack_" + lane + "_cm2s"] = float(expo_cm2s.sum())
        # S_burst: 1-ks pulse in the best-exposure (lane U) window
        ph = H.photons(f"{a['channel']}-{a['target']}", *[None, None][:0]) if False else None
        pos = json.loads((H.PHOT / "positions.json").read_text())[f"{a['channel']}-{a['target']}"]
        ph = H.photons(f"{a['channel']}-{a['target']}", pos["ra"], pos["dec"])
        w = max(sw, key=lambda x: x["expo_U"])
        iv = H.intervals(pos["ra"], pos["dec"], w["mjd0"], w["mjd1"])
        good = np.where(iv["gate"] & (iv["live"] > 0))[0]
        frac = []
        rate_U = w["rate_U"]
        for Phi in PHI_GRID:
            hits = 0
            for _ in range(max(N_REAL // 4, 40)):
                i = rng.choice(good)
                t0 = H.met_to_mjd(iv["start"][i])
                # pulse spans 1 ks from t0: exposure-weighted A_eff over the covered intervals
                sel = (iv["start"] >= iv["start"][i]) & (iv["start"] < iv["start"][i] + H.FREEZE["boxcar_s"]) & iv["gate"]
                a_eff = (iv["expo_U"][sel].sum() / max(iv["live"][sel].sum(), 1e-9)) * 1e4    # cm^2
                live = iv["live"][sel].sum()
                mu = Phi * a_eff * cont["U"] * (live / H.FREEZE["boxcar_s"])
                k = rng.poisson(mu)
                if k == 0:
                    continue
                # add k photons at random times inside the pulse, inside the union aperture
                tt = t0 + rng.uniform(0, H.FREEZE["boxcar_s"] / 86400.0, size=k)
                ph2 = {"mjd": np.sort(np.r_[ph["mjd"], tt]), }
                idx = np.argsort(np.r_[ph["mjd"], tt])
                ph2["ENERGY"] = np.r_[ph["ENERGY"], np.full(k, 500.0)][idx]
                ph2["sep"] = np.r_[ph["sep"], np.full(k, 0.3)][idx]
                s, _ = H.burst_stat(ph2, pos["ra"], pos["dec"], w["mjd0"], w["mjd1"], rate_U, iv=iv)
                hits += s > T
            frac.append(hits / max(N_REAL // 4, 40))
        rec["Phi90_burst"] = f90(PHI_GRID, frac); rec["Phi50_burst"] = f90(PHI_GRID, frac, 0.5); rec["burst_window"] = w["t_ca_utc"][:10]
        rec["burst_frac_max"] = float(max(frac))
        # power (W) through the rung cone
        for key in ("F90_stack_L", "F90_event_L_best", "F90_stack_H", "F90_event_H_best"):
            lane = "L" if "_L" in key else "H"
            rec["P90_" + key[4:]] = H.power_w(rec[key], lane, a["rung"]) if rec.get(key) else None
        out["units"].append(rec)
        fmt = lambda x: "—" if x is None else f"{x:.1e}"
        md.append(f"| {name} | {fmt(rec['F90_stack_L'])} | {fmt(rec['F90_event_L_best'])} / {fmt(rec['F90_event_L_median'])} | {fmt(rec['F90_stack_H'])} | "
                  f"{fmt(rec['F90_event_H_best'])} | {fmt(rec['Phi90_burst'])} / {fmt(rec['Phi50_burst'])} | {fmt(rec.get('P90_stack_L'))} | {fmt(rec.get('P90_event_L_best'))} |")
        print(f"[compl] {name:28s} F90 stack L {fmt(rec['F90_stack_L'])} event L {fmt(rec['F90_event_L_best'])} stack H {fmt(rec['F90_stack_H'])} burst {fmt(rec['Phi90_burst'])}", flush=True)
    # systematics
    out["systematics"] = {"gamma_pm0.3_meanE_factor_L": [float(np.log(10) / (1 / 100 - 1 / 1000) * 0 + 1)],
                          "note": "see md"}
    g = H.FREEZE["gamma"]
    def meanE(gam, emin, emax):
        if abs(gam - 2) < 1e-9:
            return np.log(emax / emin) / (1 / emin - 1 / emax)
        return ((emax ** (2 - gam) - emin ** (2 - gam)) / (2 - gam)) / ((emax ** (1 - gam) - emin ** (1 - gam)) / (1 - gam))
    sysd = {}
    for lane in ("L", "H"):
        emin, emax = H.FREEZE["lanes"][lane]["emin"], H.FREEZE["lanes"][lane]["emax"]
        base = meanE(2.0, emin, emax)
        sysd[lane] = {"meanE_MeV": base, "gamma1.7_factor": meanE(1.7, emin, emax) / base, "gamma2.3_factor": meanE(2.3, emin, emax) / base,
                      "exposure_pm": 0.10}
    out["systematics"] = sysd
    md += ["", "Systematics on the power conversion: " + "; ".join(f"lane {k}: ⟨E⟩ ×{v['gamma1.7_factor']:.2f} (Γ=1.7) / ×{v['gamma2.3_factor']:.2f} (Γ=2.3), exposure ±10 %" for k, v in sysd.items())]
    (H.RES / "completeness_v1.json").write_text(json.dumps(H.clean(out), indent=1))
    (H.RES / "completeness_v1.md").write_text("\n".join(md))


if __name__ == "__main__":
    main()
