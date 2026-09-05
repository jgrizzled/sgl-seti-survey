"""Render results/<family>_search_v1.json (+ completeness_<family>_v1.json if
present) into results/<family>_v1.md: unit table, exceedance ledger,
trial bookkeeping, power limits."""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
fam = sys.argv[1]
res = json.load(open(HERE / "results" / f"{fam}_search_v1.json"))
cpath = HERE / "results" / f"completeness_{fam}_v1.json"
comp = json.load(open(cpath)) if cpath.exists() else None
L = []
P = L.append
P(f"# Spectral-archive family — {fam} results (v1.1, run {res['run_utc'][:16]} UTC; frozen config {res.get('frozen_config')} sha {res['frozen_sha'][:12]})\n")
s = res["summary"]
P(f"Units {s['units']} ({s['searched']} searched), trials {s['trials']}, exceedances {s['exceedances']} vs {s['expected_exceedances']:.2f} expected; dispositions {s['dispositions']}.\n")

P("## Combos (null ensembles)\n")
P("| combo | nulls used | T_line per cell | null S median / p90 | σ/photon | stellar mask | status |")
P("|---|---|---|---|---|---|---|")
for k, c in res["combos"].items():
    if c.get("status") != "ok":
        P(f"| {k} | {c.get('n_ok', '')} | — | — | — | — | {c['status']} |"); continue
    T = ", ".join(f"{cl} {v:.1f}" for cl, v in c["T_line"].items())
    nd = ", ".join(f"{cl} {d['median']:.1f}/{d['p90']:.1f}" for cl, d in c["null_S_distribution"].items())
    sp = ", ".join(f"{cl} {v:.2f}" for cl, v in c["sigma_over_photon"].items())
    P(f"| {k} | {c['n_nulls']} | {T} | {nd} | {sp} | {c.get('stellar_mask_frac', 0):.3f} | ok |")

P("\n## Units\n")
P("| unit | b (R☉) | rungs | usable / n | S_line (S / T per cell) | S_coadd (S / T per cell) | exceed. | expected |")
P("|---|---|---|---|---|---|---|---|")
for u in res["units"]:
    if u.get("status") != "searched":
        P(f"| {u['unit_id']} | {u['b_rsun']:.2f} | {'+'.join(u['rungs'])} | {u.get('n_usable', 0)} / {len(u['spectra'])} | — | — | — | {u.get('status')} |"); continue
    sl = "; ".join(f"{cl} {v['S']:.1f}/{v['T']:.1f}{'**' if v['exceeds'] else ''}" for cl, v in u["S_line"].items() if v)
    sc = "; ".join(f"{cl} {v['S']:.1f}/{v['T']:.1f}{'**' if v['exceeds'] else ''}" for cl, v in u.get("S_coadd", {}).items() if v and v["S"] is not None) or "—"
    P(f"| {u['unit_id']} | {u['b_rsun']:.2f} | {'+'.join(u['rungs'])} | {u['n_usable']} / {len(u['spectra'])} | {sl} | {sc} | {len(u['exceedances'])} | {u['expected_exceedances']:.2f} |")

P("\n## Exceedance ledger (automatic ladder; manual notes in the report)\n")
P("| unit | statistic | cell | S / T | λ_bary (nm) | spectrum (UTC) | FWHM ratio | ladder evidence | disposition |")
P("|---|---|---|---|---|---|---|---|---|")
for u in res["units"]:
    for e in u.get("exceedances", []):
        ev = e["evidence"]; sh = ev.get("shape") or {}
        extra = []
        if ev.get("stellar_line"): extra.append(f"stellar {ev['stellar_line']['line']} dv {ev['stellar_line']['dv_km_s']:.0f}")
        if ev.get("sky_line"): extra.append(f"sky {ev['sky_line']['line']}")
        if ev.get("sigma_ratio"): extra.append(f"σ-band ×{ev['sigma_ratio']:.1f}")
        if ev.get("recurrence_in_unit") is not None: extra.append(f"recurs in {len(ev['recurrence_in_unit'])} other spectra")
        if ev.get("mask_frac"): extra.append(f"mask frac {ev['mask_frac']:.2f}")
        P(f"| {u['unit_id']} | {e['statistic']} | {e['cell']} | {e['S']:.1f} / {e['T']:.1f} | {ev['lambda_bary_vac_nm']:.3f} | {e.get('spectrum','')[:34]} {e.get('utc','')[:16]} | {sh.get('fwhm_ratio', float('nan')):.2f} | {'; '.join(extra)} | {e['disposition']} |")

if comp:
    P("\n## Completeness and power limits (A90 = 90 % recovery amplitude in units of the local pseudo-continuum; P through the rung cone at the cell reference wavelength; range over the cell in brackets)\n")
    P("| unit | cell | T_line | A90 best / median | λ_ref (nm) | F_λ (erg s⁻¹ cm⁻² Å⁻¹) | P per rung (W) |")
    P("|---|---|---|---|---|---|---|")
    for u in comp["units"]:
        for cl, c in u["cells"].items():
            if c.get("A90_best") is None:
                P(f"| {u['unit_id']} | {cl} | {c['T_line']:.1f} | not recovered ≤ 2.0 | {c['lambda_ref_nm']:.0f} | {c['Flambda_ref_erg_s_cm2_A']:.2e} | — |"); continue
            pw = "; ".join(f"{r} {v['at_ref']:.3g} [{v['range_over_cell'][0]:.2g}–{v['range_over_cell'][1]:.2g}]" for r, v in c["power_W"].items())
            P(f"| {u['unit_id']} | {cl} | {c['T_line']:.1f} | {c['A90_best']:.3f} / {c['A90_median']:.3f} | {c['lambda_ref_nm']:.0f} | {c['Flambda_ref_erg_s_cm2_A']:.2e} | {pw} |")
(HERE / "results" / f"{fam}_v1.md").write_text("\n".join(L) + "\n")
print("\n".join(L))
