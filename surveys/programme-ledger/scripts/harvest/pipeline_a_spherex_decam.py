"""Generate harvest YAML for the three Pipeline A reports (faithful transcription).

Per-endpoint rows are transcribed one-per-line as JSON flow mappings (valid YAML)
from the threshold-completeness Constraint records; header and aggregate rows are
block style. Every number is copied from a record or a cited results file.
"""
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path("/home/jgreene/code/sgl-seti-survey")
import sys
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "runs/programme-ledger/v2/harvest_generated"

FIELDS = ["target_id", "channel", "rung", "band", "substrate", "status", "limit_kind",
          "limit_value", "limit_unit", "limit_range", "power_mw", "power_mw_range",
          "n_events", "n_trials", "epoch_range", "source", "notes"]


def load(path):
    return [json.loads(l) for l in open(ROOT / path)]


def rung(z):
    f = lambda x: str(int(x)) if float(x).is_integer() else str(x)
    return f"{f(z[0])}-{f(z[1])} AU"


def r3(x):
    return None if x is None else round(float(x), 3)


def row(**kw):
    d = {k: None for k in FIELDS}
    d.update(kw)
    assert set(d) == set(FIELDS), set(d) ^ set(FIELDS)
    return d


def yaml_block_row(d):
    lines = []
    for i, k in enumerate(FIELDS):
        v = d[k]
        if isinstance(v, list):
            vs = "[" + ", ".join(json.dumps(x, ensure_ascii=False) for x in v) + "]"
        else:
            vs = json.dumps(v, ensure_ascii=False)
        lines.append(("  - " if i == 0 else "    ") + f"{k}: {vs}")
    return "\n".join(lines)


def yaml_flow_row(d):
    return "  - " + json.dumps(d, ensure_ascii=False)


def write(path, header, block_rows, flow_rows, comment):
    out = [f"# {c}" for c in comment]
    for k, v in header.items():
        if isinstance(v, list):
            out.append(f"{k}: [" + ", ".join(v) + "]")
        else:
            out.append(f"{k}: {json.dumps(v, ensure_ascii=False)}")
    out.append("rows:")
    out.append("  # --- aggregate rows (report tables) ---")
    out += [yaml_block_row(r) for r in block_rows]
    out.append("  # --- per-endpoint rows (one per endpoint x role x band x z interval; threshold-completeness records) ---")
    out += [yaml_flow_row(r) for r in flow_rows]
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / path).write_text("\n".join(out) + "\n")
    print(path, "rows:", len(block_rows) + len(flow_rows), "(aggregate", len(block_rows), "+ per-endpoint", len(flow_rows), ")")


# --------------------------------------------------------------------------------------
# SPHEREx v2 per-detector survey
# --------------------------------------------------------------------------------------
corr = json.load(open(ROOT / "surveys/spherex/results/template_absorption_corrected.json"))
corr_pd = {(r["cell"], r["set"], tuple(r["z_interval_au"])): r for r in corr["families"]["per-detector"]["rows"]}
corr_j6 = {(r["cell"], r["set"], tuple(r["z_interval_au"])): r for r in corr["families"]["joint6"]["rows"]}
corr_sum_pd = corr["families"]["per-detector"]["summary"]
corr_sum_j6 = corr["families"]["joint6"]["summary"]

CORR_JSON = "surveys/spherex/results/template_absorption_corrected.json"
SPX_REC = "runs/spherex/v2/records/constraint.jsonl"
J6_REC = "runs/spherex/joint6/records/constraint.jsonl"
DEC_REC = "runs/decam/v2/records/constraint.jsonl"


def spherex_rows(records_path, corr_map, corr_family, substrate_fmt):
    recs = [r for r in load(records_path) if r["completeness_kind"] == "threshold"]
    rows = []
    for r in sorted(recs, key=lambda r: (r["extra"]["set"] != "confirmatory", r["endpoint_id"], r["role"], r["band"], r["z_interval_au"][0])):
        ep, role, band, z, st = r["endpoint_id"], r["role"], r["band"], tuple(r["z_interval_au"]), r["extra"]["set"]
        cell = f"{ep}/{role}/{band}"
        base = dict(target_id=ep, channel="corridor", rung=rung(z), band=band,
                    substrate=substrate_fmt.format(band=band), n_events=1)
        if r["kind"] == "recovery_curve":
            c = corr_map[(cell, st, z)]
            fl = r["flux_limit"]
            pm = fl["per_model_m90"]
            n_inj = r["n_injections"]
            common = (f"role={role}; set={st}; corrected for static-template absorption (f_abs {c['f_abs']:.3f}, "
                      f"Δm {c['delta_mag']:.3f} mag; report-level, frozen record unchanged/flagged)")
            if c["m90_persistent_corrected"] is not None:
                notes = (f"{common}; persistent m90 uncorrected {c['m90_persistent']} (record {r['constraint_id']}, "
                         f"threshold = final-candidate curve, no veto); worst temporal model {c['worst_model']} m90 "
                         f"{c['m90']} → corrected {r3(c['m90_corrected'])}; per-model m90 uncorrected "
                         + ", ".join(f"{k} {v}" for k, v in pm.items()) + f"; n_injections {n_inj}")
                val = r3(c["m90_persistent_corrected"])
            else:
                notes = (f"{common}; NO persistent-model recovery curve for this cell×interval — value is the worst-model "
                         f"({c['worst_model']}) m90 corrected; uncorrected {c['m90']} (record {r['constraint_id']}); "
                         f"per-model m90 uncorrected " + ", ".join(f"{k} {v}" for k, v in pm.items()) + f"; n_injections {n_inj}")
                val = r3(c["m90_corrected"])
            rows.append(row(**base, status="searched", limit_kind="m90", limit_value=val, limit_unit="AB mag",
                            source=f"{CORR_JSON} families.{corr_family}.rows cell={cell} set={st} z_interval_au={list(z)}; {records_path} {r['constraint_id']}",
                            notes=notes))
        else:
            reason = r["extra"].get("not_constrainable_reason")
            q95 = r["extra"].get("q95")
            if reason == "null_unstable":
                why = "reason null_unstable (void ring null: heavy-tail / unstable; counted among the report's void cells)"
            elif reason == "insufficient_recovery_fit":
                bl = r["extra"].get("bright_limit_mag")
                why = ("reason insufficient_recovery_fit (too few epochs on, or fit brighter than the single-epoch clip limit"
                       + (f" {bl:.2f} AB" if bl is not None else "") + ")")
            elif reason is None and q95 is None:
                why = "no reason field in the joint record; q95 null (ring null void — report: heavy-tail or q95 ≤ 0 degenerate)"
            else:
                why = f"no reason field in the joint record; q95 {q95:.3f} (ring null valid; recovery fit not obtained)"
            rows.append(row(**base, status="structurally_open",
                            source=f"{records_path} {r['constraint_id']}",
                            notes=f"role={role}; set={st}; record kind not_constrainable, {why}; n_injections {r['n_injections']}"))
    return rows


# ---- spherex_survey.yaml ----
SPX_SUB = "SPHEREx QR2 spectral-image cutouts, v3 static template subtracted, exposure-PSF matched-filter stack, detector {band}"
flow = spherex_rows(SPX_REC, corr_pd, "per-detector", SPX_SUB)

# aggregate rows: corrected persistent medians (report/spherex_survey.md deferred-controls paragraph for confirmatory;
# corrected JSON summary for dev), n = cell x z-interval persistent recovery curves from surveys/spherex/results/report_tables.md
spx_conf = {"D1": (20.08, 19.30, 952, 18.79, "block"), "D2": (20.07, 19.20, 920, 18.73, "block"), "D3": (20.24, 19.43, 952, 18.61, "visit"),
            "D4": (20.36, 19.44, 960, 18.64, "visit"), "D5": (19.58, 18.75, 952, 17.87, "visit"), "D6": (19.05, 18.24, 928, 17.34, "visit")}
spx_dev = {"D1": (20.18, 19.34, 416), "D2": (20.07, 19.29, 408), "D3": (20.34, 19.60, 384),
           "D4": (20.52, 19.68, 392), "D5": (19.58, 18.93, 384), "D6": (19.17, 18.44, 392)}
# report/spherex_survey.md completeness table: worst-of-four = min over the four temporal models per detector
spx_models = {"D1": "flicker 19.26, visit 18.44, block 18.79", "D2": "flicker 19.14, visit 18.38, block 18.73",
              "D3": "flicker 19.38, visit 18.61, block 18.71", "D4": "flicker 19.51, visit 18.64, block 18.97",
              "D5": "flicker 18.69, visit 17.87, block 18.19", "D6": "flicker 18.17, visit 17.34, block 17.34"}
block = []
for b, (unc, cor, n, *_rest) in spx_conf.items():
    block.append(row(target_id="programme", channel="corridor", rung="550-10000 AU", band=b, substrate=SPX_SUB.format(band=b),
                     status="searched", limit_kind="m90", limit_value=cor, limit_unit="AB mag", n_events=n,
                     source="report/spherex_survey.md §Completeness table + §Results 'Deferred controls' paragraph; surveys/spherex/results/report_tables.md",
                     notes=f"confirmatory set (61 endpoints / 732 cells, blind), median persistent flat-Fν m90 over cells × z intervals: v2 frozen {unc} AB, corrected for static-template absorption (median 49 % flux absorbed, Δm 0.70 median) {cor} AB — the corrected value supersedes; other temporal models (uncorrected medians): {spx_models[b]}; n = {n} cell×interval persistent recovery curves"))
for b, (unc, cor, n) in spx_dev.items():
    block.append(row(target_id="programme", channel="corridor", rung="550-10000 AU", band=b, substrate=SPX_SUB.format(band=b),
                     status="searched", limit_kind="m90", limit_value=cor, limit_unit="AB mag", n_events=n,
                     source=f"surveys/spherex/results/report_tables.md §Development set; {CORR_JSON} families.per-detector.summary dev/{b}",
                     notes=f"development set (27 endpoints / 323 cells): median persistent m90 v2 frozen {unc} AB, template-absorption corrected {cor} AB (JSON summary median_m90_persistent_corrected, rounded); n = {n} cell×interval persistent recovery curves"))
block.append(row(target_id="programme", channel="corridor", rung="550-10000 AU", band="any",
                 substrate="SPHEREx QR2 per-detector cells (D1–D6)", status="structurally_open", n_events=18,
                 source="report/spherex_survey.md §Results table; surveys/spherex/results/report_tables.md §Confirmatory set",
                 notes="18 confirmatory cells void (heavy-tail / unstable ring null; 10 heavy-tail) — by band D1 2, D2 6, D3 2, D4 1, D5 2, D6 5; per-endpoint rows below carry them as not_constrainable/null_unstable"))
block.append(row(target_id="programme", channel="corridor", rung="550-10000 AU", band="any",
                 substrate="SPHEREx QR2 per-detector cells (D1–D6)", status="structurally_open", n_events=9,
                 source="report/spherex_survey.md §Results table; surveys/spherex/results/report_tables.md §Development set",
                 notes="9 development cells void (7 heavy-tail) — by band D1 1, D2 2, D3 3, D4 1, D5 2; 620 of 16,896 constraint records are not_constrainable (void null, too few epochs on, or fit brighter than the single-epoch clip limit)"))
header = {
    "survey_id": "spherex_survey", "report": "report/spherex_survey.md", "plan_section": "4.3", "pipeline": "A",
    "archives": ["spherex"], "hypothesis_version": "v2.0 (depths corrected report-level by the v2.1 deferred controls, report/spherex_joint6.md)",
    "era": "QR2 (SPHEREx quick release 2; 2–3 observing seasons)", "run_dir": "runs/spherex/v2", "records_dir": "runs/spherex/v2/records",
    "candidates": 0,
    "candidate_notes": "0 candidates: confirmatory 732 cells (R̃_FWER 2.314, 74 R>1 vs 81.3 expected), development 323 cells (R̃_FWER 2.353); no retained-ambiguous items; no catalogue veto in the v2.0 rule.",
}
write("spherex_survey.yaml", header, block, flow, [
    "Harvest of report/spherex_survey.md (SPHEREx v2 per-detector survey, QR2). Statuses per SCHEMA.md.",
    "Per-endpoint limit_value = template-absorption-CORRECTED persistent m90 from " + CORR_JSON + " (report-level correction;",
    "frozen records in runs/spherex/v2/records/constraint.jsonl unchanged, flagged). Threshold = final-candidate records (no veto).",
    "Record kind not_constrainable is mapped to status structurally_open (no usable depth: void null / too few epochs / bright fit).",
])

# ---- spherex_joint6.yaml ----
J6_SUB = "SPHEREx QR2 six-detector joint stack S_J = ΣA_b/√ΣB_b over D1–D6 (flat-Fν SED) from the v2 per-detector accumulators"
flow = spherex_rows(J6_REC, corr_j6, "joint6", J6_SUB)
block = [
    row(target_id="programme", channel="corridor", rung="550-10000 AU", band="J6", substrate=J6_SUB, status="searched",
        limit_kind="m90", limit_value=20.02, limit_unit="AB mag", n_events=944,
        source="report/spherex_joint6.md §Summary + 'Corrected limits' table; surveys/spherex/results/joint6/report_tables.md",
        notes="confirmatory set (61 endpoints / 122 cells, blind, R̃_FWER 2.249, 11 R>1 vs 12.7 expected): median persistent flat-Fν m90 20.93 AB frozen → 20.02 AB after the template-absorption correction (median Δm 0.69, p90 1.31); worst-of-four temporal model 19.15 → 18.35; other models uncorrected flicker 20.06, visit 19.27, block 19.52; n = 944 cell×interval persistent recovery curves (1,888 of 1,952 confirmatory constraints are recovery curves)"),
    row(target_id="programme", channel="corridor", rung="550-10000 AU", band="J6", substrate=J6_SUB, status="searched",
        limit_kind="m90", limit_value=20.69, limit_unit="AB mag", n_events=386,
        source="report/spherex_joint6.md §Joint cell — results table; surveys/spherex/results/template_absorption_corrected.md §joint6",
        notes="development set (27 endpoints / 54 cells, R̃_FWER 1.589, 0 candidates under all three masks): median persistent m90 21.39 AB frozen → 20.69 AB corrected; worst model 19.46 → 18.77; 772 of 864 dev constraints are recovery curves"),
    row(target_id="programme", channel="corridor", rung="550-10000 AU", band="J6", substrate=J6_SUB, status="structurally_open", n_events=3,
        source="report/spherex_joint6.md §Joint cell — results table",
        notes="3 confirmatory joint cells void (heavy-tail); per-endpoint rows below: gj-667-c/rx, lacaille-8760/rx, ross-248/rx carry q95 null in the records"),
    row(target_id="programme", channel="corridor", rung="550-10000 AU", band="J6", substrate=J6_SUB, status="structurally_open", n_events=4,
        source="report/spherex_joint6.md §Joint cell — results table + 'Dev-driven amendment (a)'",
        notes="4 development joint cells void (2 heavy-tail, 2 degenerate q95 ≤ 0 → void by amendment v2.1(a)); proxima-cen/rx named as degenerate (over-subtracted Galactic-plane field, first-pass q95 = −12.3); records give no reason field"),
]
header = {
    "survey_id": "spherex_joint6", "report": "report/spherex_joint6.md", "plan_section": "4.3", "pipeline": "A",
    "archives": ["spherex"], "hypothesis_version": "v2.1 (joint family amendment incl. dev-driven amendment (a); v2.0 freeze unchanged)",
    "era": "QR2 (SPHEREx quick release 2)", "run_dir": "runs/spherex/joint6", "records_dir": "runs/spherex/joint6/records",
    "candidates": 0,
    "candidate_notes": "0 candidates in the joint family (confirmatory 122 cells, dev 54); top confirmatory cell sigma-dra/rx R̃ 1.34 (global p 0.56); 82-eri/rx flat-Fν χ² 15.2/5 is annotation only; survey-wide FWER across per-detector + joint families ≤ 0.10.",
}
write("spherex_joint6.yaml", header, block, flow, [
    "Harvest of report/spherex_joint6.md (SPHEREx deferred controls: six-detector joint cell + template-absorption control).",
    "Per-endpoint limit_value = template-absorption-CORRECTED persistent m90 (joint6 family) from " + CORR_JSON + ";",
    "frozen records in runs/spherex/joint6/records/constraint.jsonl unchanged, flagged. Threshold = final-candidate records (no veto).",
    "Record kind not_constrainable is mapped to status structurally_open. Joint records carry no not_constrainable_reason field.",
])

# --------------------------------------------------------------------------------------
# DECam v2
# --------------------------------------------------------------------------------------
DEC_SUB = "DECam/NOIRLab instcal exposures (CTIO Blanco, 2012–2026), per-CCD Moffat matched-filter stack, NSC-star flux scale, band {band}"
recs = load(DEC_REC)
fc = {(r["extra"]["set"], r["endpoint_id"], r["role"], r["band"], tuple(r["z_interval_au"])): r for r in recs if r["completeness_kind"] == "final_candidate"}
flow = []
for r in sorted([r for r in recs if r["completeness_kind"] == "threshold"],
                key=lambda r: (r["extra"]["set"] != "confirmatory", r["endpoint_id"], r["role"], "grizY".index(r["band"]), r["z_interval_au"][0])):
    ep, role, band, z, st = r["endpoint_id"], r["role"], r["band"], tuple(r["z_interval_au"]), r["extra"]["set"]
    f = fc[(st, ep, role, band, z)]
    if f["kind"] == "recovery_curve":
        fpm = f["flux_limit"]["per_model_m90"]
        fc_note = "final-candidate curve (after the NSC DR2 flux-consistent static veto): " + ", ".join(f"{k} {v}" for k, v in fpm.items()) + f" ({f['constraint_id']})"
    else:
        fc_note = f"final-candidate curve not_constrainable ({f['extra'].get('not_constrainable_reason')}; {f['constraint_id']})"
    base = dict(target_id=ep, channel="corridor", rung=rung(z), band=band, substrate=DEC_SUB.format(band=band), n_events=1)
    if r["kind"] == "recovery_curve":
        fl = r["flux_limit"]
        pm = fl["per_model_m90"]
        if "persistent" in pm:
            val = pm["persistent"]
            head = f"role={role}; set={st}; threshold-completeness persistent m90 (record {r['constraint_id']})"
        else:
            val = fl["value"]
            head = (f"role={role}; set={st}; NO persistent-model recovery curve — value is the record's worst-model "
                    f"({fl['worst_temporal_model']}) threshold m90 (record {r['constraint_id']})")
        notes = head + "; per-model threshold m90 " + ", ".join(f"{k} {v}" for k, v in pm.items()) + f"; ci_68 {r.get('ci_68')}; n_injections {r['n_injections']}; " + fc_note
        flow.append(row(**base, status="searched", limit_kind="m90", limit_value=val, limit_unit="AB mag",
                        source=f"{DEC_REC} {r['constraint_id']}", notes=notes))
    else:
        reason = r["extra"].get("not_constrainable_reason")
        bl = r["extra"].get("bright_limit_mag")
        if reason == "null_unstable":
            why = "reason null_unstable (void ring null; one of the 2 confirmatory void cells, by band z 1 / Y 1)"
        else:
            why = (f"reason {reason} (report: single-phase, sub-floor epoch counts, or bright limits" +
                   (f"; single-epoch clip limit {bl:.2f} AB" if bl is not None else "") + ")")
        flow.append(row(**base, status="structurally_open", source=f"{DEC_REC} {r['constraint_id']}",
                        notes=f"role={role}; set={st}; record kind not_constrainable, {why}; n_injections {r['n_injections']}; {fc_note}"))

dec_conf = {"g": (23.18, 78, 23.28, 53, "flicker 23.01 / 23.00 (13), visit 23.63 / 23.65 (42), block 23.64 / 23.64 (15)"),
            "r": (22.13, 39, 23.38, 16, "flicker 23.34 / 23.35 (27), visit 23.25 / 23.25 (9), block 22.22 / 21.47 (18)"),
            "i": (22.73, 64, 22.80, 56, "flicker 22.91 / 22.91 (8), visit — (0), block 22.16 / 22.66 (8)"),
            "z": (22.08, 44, 22.08, 44, "flicker 21.89 / 21.89 (16), visit 21.94 / 21.94 (4), block — (0)"),
            "Y": (21.00, 27, 21.00, 27, "flicker 20.99 / 20.99 (6), visit 20.59 / 20.59 (16), block — (0)")}
dec_dev = {"g": (23.54, 57, 23.55, 43), "r": (23.42, 48, 23.44, 41), "i": (22.89, 41, 22.89, 41), "z": (22.56, 5, 22.56, 5)}
block = []
for b, (thr, n, fcv, nfc, others) in dec_conf.items():
    block.append(row(target_id="programme", channel="corridor", rung="550-10000 AU", band=b, substrate=DEC_SUB.format(band=b),
                     status="searched", limit_kind="m90", limit_value=thr, limit_unit="AB mag", n_events=n,
                     source="report/decam_survey.md §Completeness table; surveys/decam/results/report_tables.md §Confirmatory set",
                     notes=f"confirmatory set (14 endpoints / 86 cells over 11 corridors, blind, R̃_FWER 1.564, 11 R>1 vs 11.9 expected): median persistent m90 over cells × z intervals, threshold completeness {thr} AB (n = {n}) / final-candidate {fcv} AB (n = {nfc}); other temporal models threshold / final-candidate (n): {others}"))
for b, (thr, n, fcv, nfc) in dec_dev.items():
    block.append(row(target_id="programme", channel="corridor", rung="550-10000 AU", band=b, substrate=DEC_SUB.format(band=b),
                     status="searched", limit_kind="m90", limit_value=thr, limit_unit="AB mag", n_events=n,
                     source="surveys/decam/results/report_tables.md §Development set",
                     notes=f"development set (5 endpoints / 30 cells: 3 pilot corridors forced + struve2398 drawn; R̃_FWER 1.247, 4 R>1 vs 4.0 expected): median persistent m90 threshold {thr} AB (n = {n}) / final-candidate {fcv} AB (n = {nfc})"))
block.append(row(target_id=["hd-219134", "lalande-21185", "sigma-dra", "struve-2398-a", "struve-2398-b"], channel="corridor",
                 rung="550-10000 AU", band="Y", substrate=DEC_SUB.format(band="Y"), status="structurally_open",
                 source="report/decam_survey.md §Summary; surveys/decam/results/report_tables.md §Development set",
                 notes="'Y is void in the development set entirely' — no Y constraint records exist for the 5 development endpoints (all four temporal models n = 0)"))
block.append(row(target_id="programme", channel="corridor", rung="550-10000 AU", band="any", substrate="DECam grizY confirmatory cells", status="structurally_open", n_events=378,
                 source="report/decam_survey.md §Summary + §Coverage and data notes",
                 notes="378 of 712 confirmatory threshold constraints (53 %) are not_constrainable — single-phase, sub-floor epoch counts, or bright limits — concentrated in Y/z and single-season cells; 2 confirmatory cells void (z 1, Y 1), 0 dev void; per-endpoint rows below carry each cell×interval"))
header = {
    "survey_id": "decam_survey", "report": "report/decam_survey.md", "plan_section": "4.6", "pipeline": "A",
    "archives": ["decam"], "hypothesis_version": "v1.0", "era": "2012–2026 (epoch hold-out annotation at MJD 60400; archive still accumulating)",
    "run_dir": "runs/decam/v2", "records_dir": "runs/decam/v2/records", "candidates": 0,
    "candidate_notes": "0 candidates: confirmatory 86 cells (R̃_FWER 1.564; top cell wolf-1069/tx/g R̃ 1.10, global p 0.695), development 30 cells (R̃_FWER 1.247); mask sensitivity changes nothing; positive control (60000) Miminko recovered in g/r/i/z under no-clip scoring.",
}
write("decam_survey.yaml", header, block, flow, [
    "Harvest of report/decam_survey.md (DECam southern survey, v2 design, blind confirmatory run). Statuses per SCHEMA.md.",
    "Per-endpoint limit_value = THRESHOLD-completeness persistent m90 from runs/decam/v2/records/constraint.jsonl (the report's headline",
    "completeness); the final-candidate (post-veto) curve is quoted in notes. Cells with no persistent-model curve carry the record's",
    "worst-model m90 (flagged in notes). Record kind not_constrainable is mapped to status structurally_open (no usable depth).",
])
