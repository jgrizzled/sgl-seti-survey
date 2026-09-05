"""Completeness (freeze §6): injections into the real null series against
the fixed confirmatory thresholds, and the power conversion.

Inputs: runs/.../series_<role>_<version>.json (per unit-event: arc
epochs and the same-frame differential D per patch, from reduce) and
results/<role>_search_<version>.json (T per unit and statistic).

Injection model: a source of ZP-normalized flux F (V_eq = ZP_REF -
2.5 log10 F for a B-V 0.65 source) adds F*R to D_src at every arc
epoch (d = 1 chord; R = the stamp response measured in the run, drawn
per event from the measured ratio distribution); for S_pulse, F*R at
one random arc epoch of one random event. 100 draws per flux level;
recovery = fraction of draws with S > max(T, 0). m90 = the faintest
grid flux with recovery >= 0.9 (interpolated in log flux).

Power (declared conversion, LASCO/ATLAS anchors): F_line = 3631e-26 *
10^(-0.4 m) * dnu [W m^-2] with dnu = c * W_eff / lambda_c^2 (HI-1:
lambda_c 680 nm, W_eff 100 nm; a +-0.3 mag band-conversion systematic
is declared). S1: P_cone = F_line * pi * (0.1 AU)^2 (the beam at the
observer, post-Sun, relay at 550 AU). S2: P_tx(10 m) = F_line * pi *
(1.22 lambda_c / 10 m * d_star)^2 (diffraction-limited footprint at
the Sun), d_star from the registry parallax.

Output: results/completeness_<version>.json, results/power_limits_<version>.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hi_lib as H
import series as S

REPO = Path(__file__).resolve().parents[3]
SURV = REPO / "surveys" / "stereo-hi-crossings"
SER = REPO / "runs" / "stereo-hi-crossings" / "series"

FLUX_GRID = np.array([0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0])
N_DRAWS = 100
C_LIGHT, AU_M, PC_M = 2.99792458e8, 1.495978707e11, 3.0857e16
LAMBDA_C, W_EFF = 680e-9, 100e-9
DNU = C_LIGHT * W_EFF / LAMBDA_C ** 2


def f_line(m_veq: float) -> float:
    return 3631e-26 * 10 ** (-0.4 * m_veq) * DNU


def power_W(channel: str, tid: str, m: float, plx_mas: float) -> dict:
    F = f_line(m)
    out = {"F_line_W_m2": F, "P_cone_W": F * np.pi * (0.1 * AU_M) ** 2}
    if channel == "S2":
        d = 1000.0 / plx_mas * PC_M
        out["P_tx_10m_W"] = F * np.pi * (1.22 * LAMBDA_C / 10.0 * d) ** 2
    return out


def stats_from(series: dict, inject: float = 0.0, ratios=None, rng=None,
               pulse: tuple | None = None) -> dict | None:
    """(S_stack, S_event, S_pulse) of the source patch under injection."""
    zs, pulses = [], []
    for i, (key, se) in enumerate(series.items()):
        d = np.array(se["D"]["0"], float)
        if inject:
            r = float(rng.choice(ratios)) if ratios is not None and len(ratios) else 1.0
            d = d + inject * r
        if pulse is not None and pulse[0] == i:
            r = float(rng.choice(ratios)) if ratios is not None and len(ratios) else 1.0
            d = d.copy()
            d[pulse[1]] += pulse[2] * r
        sd = S._robust_sd(d)
        if not sd > 0 or len(d) < S.CFG["gates"]["event"]["min_arc_epochs"]:
            continue
        zs.append(float(np.mean(d) / (sd / np.sqrt(len(d)))))
        pulses.append(float(np.max(d - np.median(d)) / sd))
    if len(zs) < S.CFG["gates"]["unit"]["min_events_for_stack"]:
        return None
    z = np.array(zs)
    return {"S_stack": float(np.sum(z) / np.sqrt(len(z))), "S_event": float(np.max(z)),
            "S_pulse": float(np.max(pulses))}


def m90_from_curve(grid, curve):
    c = np.asarray(curve)
    if not (c >= 0.9).any():
        return None
    i = int(np.argmax(c >= 0.9))
    if i == 0:
        return float(grid[0])
    # log-interpolate between the bracketing grid points
    x0, x1 = np.log10(grid[i - 1]), np.log10(grid[i])
    y0, y1 = c[i - 1], c[i]
    return float(10 ** (x0 + (0.9 - y0) / max(y1 - y0, 1e-9) * (x1 - x0)))


def main(role: str = "confirmatory", version: str = "v1") -> None:
    search = json.loads((SURV / "results" / f"{role}_search_{version}.json").read_text())
    series_all = json.loads((SER / f"series_{role}_{version}.json").read_text())
    rng = np.random.default_rng(S.SEED)
    # stamp-response ratios measured in the run
    ratios = []
    with open(SER / f"measurements_{role}_v1.jsonl") as fh:
        for line in fh:
            rec = json.loads(line)
            ratios.extend(s["ratio"] for s in rec.get("stamps", []))
    ratios = np.array([r for r in ratios if 0.3 < r < 1.5])
    R_med = float(np.median(ratios)) if len(ratios) else 1.0
    from sglseti import load_target_registry
    reg = load_target_registry(REPO / "registries" / "pilot_wise_2026.yaml")

    results, power = {}, {}
    for ukey, unit in search["units"].items():
        if unit.get("status") != "searched":
            results[ukey] = {"status": unit.get("status")}
            continue
        tid, uname = ukey.split("|")
        channel = uname[:2]
        series = {k: v for k, v in series_all.items() if k.startswith(ukey + "|") and "0" in v["D"]}
        # recovery threshold = max(T, observed null S, 0): an injected source
        # must rise above the adjudicated systematic where the null already
        # exceeds T (otherwise m90 is degenerate at the grid floor)
        T = {s: max(unit[s]["T"], unit[s]["S"], 0.0) for s in ("S_stack", "S_event", "S_pulse")}
        rec = {"n_events_series": len(series), "thresholds": T,
               "null_S": {s: unit[s]["S"] for s in ("S_stack", "S_event", "S_pulse")}, "curves": {}}
        curves = {"S_stack": [], "S_event": []}
        for F in FLUX_GRID:
            hits = {"S_stack": 0, "S_event": 0}
            for _ in range(N_DRAWS):
                st = stats_from(series, inject=F, ratios=ratios, rng=rng)
                if st is None:
                    continue
                for s in hits:
                    hits[s] += st[s] > T[s]
            for s in hits:
                curves[s].append(hits[s] / N_DRAWS)
        # pulse: single-epoch injection at a random epoch of a random event
        pgrid = np.concatenate([FLUX_GRID * 3, [40.0, 60.0, 80.0, 120.0, 160.0]])
        pc = []
        keys = list(series)
        for F in pgrid:
            hits = 0
            for _ in range(N_DRAWS):
                i = int(rng.integers(len(keys)))
                j = int(rng.integers(len(series[keys[i]]["D"]["0"])))
                st = stats_from(series, ratios=ratios, rng=rng, pulse=(i, j, F))
                if st and st["S_pulse"] > T["S_pulse"]:
                    hits += 1
            pc.append(hits / N_DRAWS)
        rec["curves"] = {"flux_grid": FLUX_GRID.tolist(), "S_stack": curves["S_stack"],
                         "S_event": curves["S_event"], "pulse_flux_grid": pgrid.tolist(), "S_pulse": pc}
        plx = reg[tid].astrometry.parallax_mas
        power[ukey] = {}
        for s, grid, curve in (("S_stack", FLUX_GRID, curves["S_stack"]),
                               ("S_event", FLUX_GRID, curves["S_event"]),
                               ("S_pulse", pgrid, pc)):
            f90 = m90_from_curve(grid, curve)
            if f90 is not None and unit[s]["S"] > unit[s]["T"] and curve[0] >= 0.9:
                # the null already exceeds T and any injection "recovers":
                # completeness is not measurable for this statistic
                rec[f"m90_note_{s}"] = "degenerate: observed null exceeds T (adjudicated systematic)"
                f90 = None
            rec[f"m90_flux_{s}"] = f90
            if f90:
                m = H.ZP_REF - 2.5 * np.log10(f90)
                rec[f"m90_Veq_{s}"] = round(m, 2)
                power[ukey][s] = {"m90_Veq": round(m, 2), **power_W(channel, tid, m, plx)}
        results[ukey] = rec
        print(ukey, {k: v for k, v in rec.items() if k.startswith("m90_Veq")}, flush=True)

    out = {"stage": f"completeness_{version}", "role": role, "seed": S.SEED,
           "stamp_response": {"R_median": R_med, "n": int(len(ratios)),
                              "mad": S._robust_sd(ratios) if len(ratios) else None},
           "flux_grid": FLUX_GRID.tolist(), "n_draws": N_DRAWS,
           "thresholds_source": f"{role}_search_{version}.json (fixed)",
           "band": {"lambda_c_nm": 680, "w_eff_nm": 100, "dnu_Hz": DNU,
                    "systematic_mag": 0.3, "note": "V_eq of a B-V 0.65 source; 532 nm line out of band"},
           "units": results}
    (SURV / "results" / f"completeness_{version}.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    (SURV / "results" / f"power_limits_{version}.json").write_text(json.dumps(power, indent=1, default=float) + "\n")
    print("wrote completeness/power_limits", version)


if __name__ == "__main__":
    main(*(sys.argv[1:3]))
