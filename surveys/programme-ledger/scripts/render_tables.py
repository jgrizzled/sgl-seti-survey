"""Render the programme ledger v2 into the Markdown tables used by
report/programme_ledger.md (written to results/programme_ledger_v2_tables.md).

Usage: uv run python surveys/programme-ledger/scripts/render_tables.py
"""
from __future__ import annotations

import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from astropy.table import Table

ROOT = Path(__file__).resolve().parents[3]
RESULTS = ROOT / "surveys/programme-ledger/results"
RUN_DIR = ROOT / "runs/programme-ledger/v2"

STATE_ORDER = ["searched", "retained_ambiguous", "vetoed_known_source", "constraint_only",
               "ledger_only", "structurally_open", "no_survey", "not_constrainable"]
STATE_GLYPH = {"searched": "S", "retained_ambiguous": "S*", "vetoed_known_source": "Sv",
               "constraint_only": "C", "ledger_only": "L", "structurally_open": "open",
               "no_survey": "–", "not_constrainable": "×"}
CROSSING_CELLS = [("A", "0.1 AU"), ("A", "1.0 AU"),
                  ("B", "1.2 Rsun"), ("B", "2.5 Rsun"), ("B", "0.1 AU"),
                  ("S1", "1.2 Rsun"), ("S1", "2.5 Rsun"), ("S1", "0.1 AU"),
                  ("S2", "0.1 AU"), ("S2", "0.90 AU"), ("S2", "0.95 AU"), ("S2", "1.0 AU")]


def _sig(x, n=2):
    """n significant figures, plain decimal notation."""
    if x == 0:
        return "0"
    d = n - int(np.floor(np.log10(abs(x)))) - 1
    return f"{round(x, d):.{max(d, 0)}f}".rstrip("0").rstrip(".") if d > 0 else f"{int(round(x, d))}"


def fmt_power(mw):
    if mw is None or not np.isfinite(mw):
        return ""
    if mw >= 1000:
        return f"{_sig(mw/1000)} GW"
    if mw >= 1:
        return f"{_sig(mw)} MW"
    if mw >= 1e-3:
        return f"{_sig(mw*1e3)} kW"
    return f"{_sig(mw*1e6)} W"


def main() -> int:
    tab = Table.read(RUN_DIR / "programme_ledger_v2.ecsv")
    matrix = json.loads((RESULTS / "programme_ledger_v2_matrix.json").read_text())
    summary = json.loads((RESULTS / "programme_ledger_v2_summary.json").read_text())
    out = []

    # --- 1. survey index -------------------------------------------------
    out.append("### Table 1 — survey index\n")
    out.append("| survey | plan § | pipeline | rows | candidates | states |")
    out.append("|---|---|---|---|---|---|")
    for sid, s in sorted(summary["per_survey"].items(), key=lambda kv: (str(kv[1]["plan_section"]), kv[0])):
        states = ", ".join(f"{k} {v}" for k, v in sorted(s["statuses"].items(), key=lambda kv: STATE_ORDER.index(kv[0]) if kv[0] in STATE_ORDER else 99))
        out.append(f"| `{sid}` | {s['plan_section']} | {s['pipeline'] or '–'} | {s['rows']} | {s['candidates'] if s['candidates'] is not None else '–'} | {states} |")
    out.append("")

    # --- 2. status census -------------------------------------------------
    out.append("### Table 2 — row census by coverage state and channel\n")
    chans = ["corridor", "A", "B", "S1", "S2"]
    out.append("| state | " + " | ".join(chans) + " | total |")
    out.append("|---|" + "---|" * (len(chans) + 1))
    for st in STATE_ORDER:
        counts = [int(np.sum((tab["status"] == st) & (tab["channel"] == c))) for c in chans]
        if sum(counts):
            out.append(f"| {st} | " + " | ".join(str(c) for c in counts) + f" | {sum(counts)} |")
    out.append("")

    # --- 3. per-target crossings matrix (best state + best power) --------
    out.append("### Table 3 — per-target crossings matrix (best state; best transmitter-power limit where one is published)\n")
    out.append("Glyphs: S searched (calibrated depth), C constraint-only, L ledger-only, open = structurally open, – no survey, "
               "× not constrainable; a trailing * marks a retained-ambiguous exceedance on the cell and v an exceedance vetoed as a "
               "known source (S* / Sv alone: the cell has only that disposition). The power is the lowest published transmitter-power "
               "limit among the cell's searched or constraint-only rows, whatever the band or statistic. "
               "Blank = no report speaks to the cell for that target (programme-wide rows are in Table 4).\n")
    hdr = "| target | " + " | ".join(f"{c} {r}" for c, r in CROSSING_CELLS) + " |"
    out.append(hdr)
    out.append("|---|" + "---|" * len(CROSSING_CELLS))
    targets = sorted(t for t in matrix if t != "programme")
    for t in targets:
        cells = []
        for c, r in CROSSING_CELLS:
            m = matrix[t].get(f"{c}|{r}")
            if not m:
                cells.append("")
                continue
            g = STATE_GLYPH[m["best_state"]]
            if m["best_state"] not in ("retained_ambiguous", "vetoed_known_source"):
                if "retained_ambiguous" in m["states"]:
                    g += "*"
                if "vetoed_known_source" in m["states"]:
                    g += "v"
            p = fmt_power(m["best_power_mw"])
            cells.append(f"{g} {p}".strip())
        if any(cells):
            out.append(f"| {t} | " + " | ".join(cells) + " |")
    out.append("")

    # --- 4. programme-wide rows -------------------------------------------
    out.append("### Table 4 — programme-wide statements (rows that speak to the whole target list)\n")
    out.append("| survey | channel | rung | band | state | limit | power | n | source |")
    out.append("|---|---|---|---|---|---|---|---|---|")
    prog = tab[tab["target_id"] == "programme"]
    for r in sorted(prog, key=lambda r: (str(r["plan_section"]), str(r["channel"]), str(r["rung"]))):
        lim = ""
        if np.isfinite(r["limit_value"]):
            lim = f"{r['limit_value']:g} {r['limit_unit']}"
        elif np.isfinite(r["limit_hi"]):
            lim = f"{r['limit_lo']:g}–{r['limit_hi']:g} {r['limit_unit']}"
        pw = fmt_power(r["power_mw"]) if np.isfinite(r["power_mw"]) else (
            f"{fmt_power(r['power_mw_lo'])}–{fmt_power(r['power_mw_hi'])}" if np.isfinite(r["power_mw_hi"]) else "")
        n = r["n_events"] if r["n_events"] >= 0 else ""
        out.append(f"| `{r['survey_id']}` | {r['channel']} | {r['rung']} | {r['band']} | {r['status']} | {lim} | {pw} | {n} | {r['source']} |")
    out.append("")

    # --- 5. corridor (Pipeline A) per-endpoint depths ---------------------
    out.append("### Table 5 — corridor (Pipeline A) depth per endpoint: deepest published persistent-source m90 per survey (band)\n")
    out.append("WISE depths are Vega magnitudes; all others AB. A blank means the endpoint has no searched corridor cell in that survey.\n")
    cor = tab[(tab["channel"] == "corridor") & (tab["status"] == "searched") & (tab["target_id"] != "programme")]
    best = defaultdict(dict)
    for r in cor:
        v = r["limit_value"] if np.isfinite(r["limit_value"]) else None
        if v is None:
            continue
        sid = str(r["survey_id"])
        if sid not in best[r["target_id"]] or v > best[r["target_id"]][sid][0]:
            best[r["target_id"]][sid] = (float(v), str(r["band"]))
    surveys = ["wise_survey", "ztf_survey", "ps1_survey", "joint_ps1_ztf_wise", "spherex_survey", "spherex_joint6", "decam_survey"]
    surveys = [k for k in surveys if any(k in d for d in best.values())]
    out.append("| endpoint | " + " | ".join(f"`{k}`" for k in surveys) + " | n surveys |")
    out.append("|---|" + "---|" * (len(surveys) + 1))
    for t in sorted(best):
        out.append(f"| {t} | " + " | ".join(f"{best[t][k][0]:.1f} ({best[t][k][1]})" if k in best[t] else "" for k in surveys) + f" | {len(best[t])} |")
    out.append("")

    # --- 6. open-cell census ---------------------------------------------
    out.append("### Table 6 — open and unconstrainable cells (crossings channels), by rung\n")
    out.append("| channel | rung | structurally_open rows | targets | not_constrainable rows | no_survey rows |")
    out.append("|---|---|---|---|---|---|")
    for c, r in CROSSING_CELLS:
        sel = (tab["channel"] == c) & (tab["rung"] == r)
        so = tab[sel & (tab["status"] == "structurally_open")]
        nc = int(np.sum(sel & (tab["status"] == "not_constrainable")))
        ns = int(np.sum(sel & (tab["status"] == "no_survey")))
        tg = sorted({str(x) for x in so["target_id"] if x != "programme"})
        if len(so) or nc or ns:
            out.append(f"| {c} | {r} | {len(so)} | {len(tg)} | {nc} | {ns} |")
    out.append("")

    # --- 7. candidates/ambiguous ----------------------------------------
    out.append("### Table 7 — retained-ambiguous and vetoed exceedances\n")
    out.append("| survey | target | channel | rung | band | state | notes |")
    out.append("|---|---|---|---|---|---|---|")
    amb = tab[(tab["status"] == "retained_ambiguous") | (tab["status"] == "vetoed_known_source")]
    for r in sorted(amb, key=lambda r: (str(r["survey_id"]), str(r["target_id"]))):
        out.append(f"| `{r['survey_id']}` | {r['target_id']} | {r['channel']} | {r['rung']} | {r['band']} | {r['status']} | {str(r['notes'])[:160]} |")
    out.append("")

    (RESULTS / "programme_ledger_v2_tables.md").write_text("\n".join(out))
    print(f"wrote {RESULTS / 'programme_ledger_v2_tables.md'} ({len(out)} lines)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
