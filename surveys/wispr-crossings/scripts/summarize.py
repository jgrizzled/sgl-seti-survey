"""Write results/<role>_<version>.md — the per-unit table of a search
stage (S/T per statistic, exceedances, event counts) from the reduce
json, plus the completeness/power table when present."""

from __future__ import annotations

import json
import sys
from pathlib import Path

SURV = Path(__file__).resolve().parents[1]


def main(role: str, version: str = "v1") -> None:
    d = json.loads((SURV / "results" / f"{role}_search_{version}.json").read_text())
    lines = [f"# {role} search {version} — unit table", "",
             f"Construction: {d['construction']}; frames usable {d['frames']['usable']}, "
             f"unusable {d['frames']['unusable']}; median matched stars {d['frames']['median_n_match']}, "
             f"astrometric rms {d['frames']['median_astrom_rms_px']} px, ZP0 {d['frames']['median_zp0']:.3f}, "
             f"ZP MAD {d['frames']['median_zp_mad']:.3f} mag. Epoch vetoes: {d['epoch_vetoes']}. "
             f"Stamp response {d['stamp_response']['median_ratio']:.3f} ± {d['stamp_response']['mad']:.3f} "
             f"(n {d['stamp_response']['n']}); colour coefficient {d['colour_system']['coeff']:.3f} (n {d['colour_system']['n']}).", "",
             f"**Trials {d['trials']['searched']}, exceedances {d['trials']['exceedances']}, "
             f"expected control crossings {d['trials']['expected_control_crossings']}.**", "",
             "| unit | status | n_ev | z med / MAD | S_stack S / T | S_event S / T | S_pulse S / T |",
             "|---|---|---|---|---|---|---|"]
    for u, r in d["units"].items():
        if r["status"] != "searched":
            lines.append(f"| {u} | {r['status']} | {r.get('n_included', '')} | | | | |")
            continue
        p0 = r["per_patch"]["0"]
        cells = []
        for s in ("S_stack", "S_event", "S_pulse"):
            x = r[s]
            cells.append(f"{x['S']:.2f} / {x['T']:.2f}{' **EXC**' if x['exceeds'] else ''}")
        lines.append(f"| {u} | searched | {p0['n_events']} | {p0['z_median']:+.2f} / {p0['z_mad']:.2f} | "
                     + " | ".join(cells) + " |")
    comp = SURV / "results" / f"completeness_{version}.json"
    if comp.exists() and role == "confirmatory":
        c = json.loads(comp.read_text())
        pw = json.loads((SURV / "results" / f"power_limits_{version}.json").read_text())
        lines += ["", "## Completeness (m90, V-equivalent) and power", "",
                  "| unit | S_stack m90 / P | S_event m90 / P | S_pulse m90 / P |", "|---|---|---|---|"]
        for u, r in c["units"].items():
            if "m90_Veq_S_stack" not in r and "m90_Veq_S_event" not in r:
                continue
            cells = []
            for s in ("S_stack", "S_event", "S_pulse"):
                p = pw.get(u, {}).get(s)
                if p:
                    P = p.get("P_tx_10m_W", p["P_cone_W"])
                    cells.append(f"{p['m90_Veq']:.2f} / {P:.2e} W")
                else:
                    cells.append("—")
            lines.append(f"| {u} | " + " | ".join(cells) + " |")
    (SURV / "results" / f"{role}_{version}.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main(*sys.argv[1:3])
