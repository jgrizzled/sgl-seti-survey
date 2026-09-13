"""Generate the harvest YAML for report/ptf_survey.md (PTF/iPTF corridor
survey, Pipeline A, v2 design) — faithful transcription of the
threshold-completeness Constraint records plus aggregate rows computed
from the null ensembles and records (the same numbers the report
quotes). Same row/format conventions as pipeline_a_spherex_decam.py.

Usage: python pipeline_a_ptf.py [out_dir]
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path("/home/jgreene/code/sgl-seti-survey")
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "runs/programme-ledger/v2/harvest_generated"
FIELDS = ["target_id", "channel", "rung", "band", "substrate", "status", "limit_kind",
          "limit_value", "limit_unit", "limit_range", "power_mw", "power_mw_range",
          "n_events", "n_trials", "epoch_range", "source", "notes"]
REC = "runs/ptf/v2/records/constraint.jsonl"
SUB = "PTF/iPTF level-1 epochal images (Palomar P48, 2009–2015), scie-direct Moffat matched-filter stack, PS1-star flux scale, band {band}"


def load(path):
    return [json.loads(l) for l in open(ROOT / path)]


def rung(z):
    f = lambda x: str(int(x)) if float(x).is_integer() else str(x)
    return f"{f(z[0])}-{f(z[1])} AU"


def row(**kw):
    d = {k: None for k in FIELDS}
    d.update(kw)
    assert set(d) == set(FIELDS), set(d) ^ set(FIELDS)
    return d


def yaml_block_row(d):
    lines = []
    for i, k in enumerate(FIELDS):
        v = d[k]
        vs = ("[" + ", ".join(json.dumps(x, ensure_ascii=False) for x in v) + "]") if isinstance(v, list) \
            else json.dumps(v, ensure_ascii=False)
        lines.append(("  - " if i == 0 else "    ") + f"{k}: {vs}")
    return "\n".join(lines)


def yaml_flow_row(d):
    return "  - " + json.dumps(d, ensure_ascii=False)


def write(path, header, block_rows, flow_rows, comment):
    out = [f"# {c}" for c in comment]
    for k, v in header.items():
        out.append(f"{k}: [" + ", ".join(v) + "]" if isinstance(v, list) else f"{k}: {json.dumps(v, ensure_ascii=False)}")
    out.append("rows:")
    out.append("  # --- aggregate rows (report tables) ---")
    out += [yaml_block_row(r) for r in block_rows]
    out.append("  # --- per-endpoint rows (one per endpoint x role x band x z interval; threshold-completeness records) ---")
    out += [yaml_flow_row(r) for r in flow_rows]
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / path).write_text("\n".join(out) + "\n")
    print(path, "rows:", len(block_rows) + len(flow_rows), "(aggregate", len(block_rows), "+ per-endpoint", len(flow_rows), ")")


recs = load(REC)
thr = [r for r in recs if r["completeness_kind"] == "threshold"]
fc_by = {(r["endpoint_id"], r["role"], r["band"], tuple(r["z_interval_au"]), r["extra"]["set"]): r
         for r in recs if r["completeness_kind"] == "final_candidate"}
flow = []
for r in sorted(thr, key=lambda r: (r["extra"]["set"] != "confirmatory", r["endpoint_id"], r["role"], r["band"], r["z_interval_au"][0])):
    ep, role, band, z, st = r["endpoint_id"], r["role"], r["band"], tuple(r["z_interval_au"]), r["extra"]["set"]
    f = fc_by.get((ep, role, band, z, st))
    if f is None:
        fc_note = "no final-candidate record"
    elif f["kind"] == "recovery_curve":
        fc_note = f"final-candidate persistent m90 {f['flux_limit']['per_model_m90'].get('persistent')} (record {f['constraint_id']})"
    else:
        fc_note = f"final-candidate curve not_constrainable ({f['extra'].get('not_constrainable_reason')}; {f['constraint_id']})"
    base = dict(target_id=ep, channel="corridor", rung=rung(z), band=band, substrate=SUB.format(band=band), n_events=1)
    if r["kind"] == "recovery_curve":
        fl = r["flux_limit"]; pm = fl["per_model_m90"]
        if "persistent" in pm:
            val = pm["persistent"]
            head = f"role={role}; set={st}; threshold-completeness persistent m90 (record {r['constraint_id']})"
        else:
            val = fl["value"]
            head = (f"role={role}; set={st}; NO persistent-model recovery curve — value is the record's worst-model "
                    f"({fl['worst_temporal_model']}) threshold m90 (record {r['constraint_id']})")
        notes = head + "; per-model threshold m90 " + ", ".join(f"{k} {v}" for k, v in pm.items()) + \
            f"; ci_68 {r.get('ci_68')}; n_injections {r['n_injections']}; bright limit {r['extra'].get('bright_limit_mag')}; " + fc_note
        flow.append(row(**base, status="searched", limit_kind="m90", limit_value=val, limit_unit="AB mag",
                        source=f"{REC} {r['constraint_id']}", notes=notes))
    else:
        reason = r["extra"].get("not_constrainable_reason"); bl = r["extra"].get("bright_limit_mag")
        why = "reason null_unstable (void ring null)" if reason == "null_unstable" else \
            (f"reason {reason} (sub-floor epoch counts, single-phase cadence, or bright limits" +
             (f"; single-epoch clip limit {bl:.2f} AB" if bl is not None else "") + ")")
        flow.append(row(**base, status="structurally_open", source=f"{REC} {r['constraint_id']}",
                        notes=f"role={role}; set={st}; record kind not_constrainable, {why}; n_injections {r['n_injections']}; {fc_note}"))

# aggregates from the null ensembles + records
block = []
nulls = {s: json.load(open(ROOT / f"runs/ptf/v2/nulls/{s}_null_ensemble.json"))["primary"] for s in ("confirmatory", "dev")}
for st in ("confirmatory", "dev"):
    n = nulls[st]
    for band in ("g", "R"):
        rr = [r for r in thr if r["extra"]["set"] == st and r["band"] == band and r["kind"] == "recovery_curve"
              and "persistent" in r["flux_limit"]["per_model_m90"]]
        if not rr:
            continue
        vals = [r["flux_limit"]["per_model_m90"]["persistent"] for r in rr]
        fcs = [fc_by[(r["endpoint_id"], r["role"], r["band"], tuple(r["z_interval_au"]), st)] for r in rr]
        fcv = [f["flux_limit"]["per_model_m90"]["persistent"] for f in fcs
               if f["kind"] == "recovery_curve" and "persistent" in f["flux_limit"]["per_model_m90"]]
        others = []
        for m in ("flicker", "visit", "block"):
            v = [r["flux_limit"]["per_model_m90"][m] for r in thr if r["extra"]["set"] == st and r["band"] == band
                 and r["kind"] == "recovery_curve" and m in r["flux_limit"]["per_model_m90"]]
            others.append(f"{m} {np.median(v):.2f} ({len(v)})" if v else f"{m} — (0)")
        label = ("confirmatory set (blind, once)" if st == "confirmatory" else "development set (4 pilot corridors forced + drawn)")
        block.append(row(target_id="programme", channel="corridor", rung="550-10000 AU", band=band, substrate=SUB.format(band=band),
                         status="searched", limit_kind="m90", limit_value=round(float(np.median(vals)), 3), limit_unit="AB mag", n_events=len(vals),
                         source=f"report/ptf_survey.md §Completeness; runs/ptf/v2/nulls/{st}_null_ensemble.json; {REC}",
                         notes=f"{label}: {n['n_endpoints']} endpoints / {n['n_cells']} cells, R̃_FWER {n['fwer']['R_fwer']:.3f}, "
                               f"{n['n_real_R_gt_1']} R>1 vs {n['expected_R_gt_1_from_ring']:.1f} expected, {n['n_null_unstable']} void, "
                               f"{n['n_candidates']} candidates; median persistent m90 over cells × z intervals, threshold completeness "
                               f"{np.median(vals):.2f} AB (n = {len(vals)}) / final-candidate {np.median(fcv):.2f} AB (n = {len(fcv)}); "
                               f"other temporal models threshold median (n): " + ", ".join(others)))
n_nc = sum(1 for r in thr if r["extra"]["set"] == "confirmatory" and r["kind"] != "recovery_curve")
n_thr = sum(1 for r in thr if r["extra"]["set"] == "confirmatory")
block.append(row(target_id="programme", channel="corridor", rung="550-10000 AU", band="any", substrate="PTF g/R confirmatory cells",
                 status="structurally_open", n_events=n_nc, source=f"report/ptf_survey.md §Coverage; {REC}",
                 notes=f"{n_nc} of {n_thr} confirmatory threshold constraints are not_constrainable (sub-floor epoch counts, single-phase "
                       f"cadence, or bright limits); per-endpoint rows carry each cell×interval"))
ov = json.load(open(ROOT / "surveys/ptf/targets/overlay_v1.json"))
below = [r for r in ov["rows"] if not r["searchable"]]
block.append(row(target_id=sorted(e for r in below for e in r["endpoints"]), channel="corridor", rung="550-10000 AU", band="any",
                 substrate="PTF level-1 epochal images", status="ledger_only", n_events=len(below),
                 source="surveys/ptf/targets/overlay_v1.md; report/ptf_survey.md §Coverage",
                 notes="corridors below the frozen epoch floor (best cell < 5 usable exposures) or with no PTF coverage at all — "
                       "coverage-without-statistic, not part of the search family (hypotheses v1.0 §8.2)"))
c = nulls["confirmatory"]; d = nulls["dev"]
header = {
    "survey_id": "ptf_survey", "report": "report/ptf_survey.md", "plan_section": "4.13", "pipeline": "A",
    "archives": ["ptf"], "hypothesis_version": "v1.0", "era": "2009-03 – 2015-01 (closed archive; no epoch hold-out)",
    "run_dir": "runs/ptf/v2", "records_dir": "runs/ptf/v2/records", "candidates": 0,
    "candidate_notes": (f"0 candidates: confirmatory {c['n_cells']} cells (R̃_FWER {c['fwer']['R_fwer']:.3f}; {c['n_real_R_gt_1']} R>1 vs "
                        f"{c['expected_R_gt_1_from_ring']:.1f} expected; {c['n_null_unstable']} void), development {d['n_cells']} cells "
                        f"(R̃_FWER {d['fwer']['R_fwer']:.3f}); positive controls (798452) and (388125) recovered under the frozen rule."),
}
write("ptf_survey.yaml", header, block, flow, [
    "Harvest of report/ptf_survey.md (PTF/iPTF corridor survey, v2 design, blind confirmatory run). Statuses per SCHEMA.md.",
    "Per-endpoint limit_value = THRESHOLD-completeness persistent m90 from runs/ptf/v2/records/constraint.jsonl; the final-candidate",
    "(post-veto) curve is quoted in notes. Cells with no persistent-model curve carry the record's worst-model m90 (flagged in notes).",
    "Record kind not_constrainable is mapped to status structurally_open; below-floor corridors to ledger_only.",
])
