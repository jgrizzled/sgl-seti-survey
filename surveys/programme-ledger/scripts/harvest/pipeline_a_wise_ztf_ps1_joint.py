"""Harvest the four Pipeline A corridor reports into programme-ledger YAML.

Every number is copied from the report, from the cited results table, or from
the constraint.jsonl records the report cites. Nothing is computed except the
min/max of quoted bright limits (labelled as such in the row notes).
"""
import json
import collections
import os

ROOT = "/home/jgreene/code/sgl-seti-survey"
import sys
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "runs/programme-ledger/v2/harvest_generated")

KEYS = ["target_id", "channel", "rung", "band", "substrate", "status", "limit_kind",
        "limit_value", "limit_unit", "limit_range", "power_mw", "power_mw_range",
        "n_events", "n_trials", "epoch_range", "source", "notes"]


def row(**kw):
    r = {k: None for k in KEYS}
    r["channel"] = "corridor"
    r.update(kw)
    unknown = set(kw) - set(KEYS)
    assert not unknown, unknown
    return r


def zstr(iv):
    a, b = iv
    fa = ("%g" % a)
    fb = ("%g" % b)
    return f"{fa}-{fb} AU"


def load(path):
    return [json.loads(l) for l in open(os.path.join(ROOT, path))]


def endpoints(recs):
    return sorted(set(r["endpoint_id"] for r in recs))


def record_rows(recs, records_path, substrate, unit, reason_text, extra_note=None,
                use_kind_priority=("final_candidate", "threshold")):
    """One row per (endpoint, role, band, z interval); the final-candidate record
    is the published depth where it exists, otherwise the threshold record."""
    by = collections.defaultdict(dict)
    for r in recs:
        key = (r["endpoint_id"], r["role"], r["band"], tuple(r["z_interval_au"]))
        by[key][r["completeness_kind"]] = r
    rows = []
    for key in sorted(by, key=lambda k: (k[0], k[1], k[2], k[3][0])):
        v = by[key]
        rec = None
        for ck in use_kind_priority:
            if ck in v:
                rec = v[ck]
                break
        other = v.get("threshold") if rec["completeness_kind"] == "final_candidate" else None
        ep, role, band, iv = key
        ex = rec["extra"]
        notes = [f"role={role}", f"set={ex['set']}", f"kind={rec['kind']}",
                 f"completeness_kind={rec['completeness_kind']}",
                 f"exclusion_claim={str(ex['exclusion_claim']).lower()}"]
        if rec["kind"] == "recovery_curve":
            fl = rec["flux_limit"]
            pm = fl["per_model_m90"]
            notes.append(f"worst_model={fl['worst_temporal_model']}")
            notes.append("m90 persistent/flicker/visit/block="
                         f"{pm.get('persistent','—')}/{pm.get('flicker','—')}/{pm.get('visit','—')}/{pm.get('block','—')}")
            if rec.get("ci_68"):
                notes.append(f"ci68={rec['ci_68'][0]}-{rec['ci_68'][1]}")
            notes.append(f"n_inj={rec['n_injections']}")
            if other is not None:
                if other["kind"] == "recovery_curve":
                    tv = other["flux_limit"]["value"]
                    if tv != fl["value"]:
                        notes.append(f"threshold-kind m90={tv}")
                    else:
                        notes.append("threshold-kind m90 identical")
                else:
                    notes.append("threshold-kind record not_constrainable")
            status, lk, lv, lu = "searched", "m90", fl["value"], unit
        else:
            reason = ex.get("not_constrainable_reason")
            if reason is None and "not_constrainable_reason" not in ex:
                notes.append("reason not recorded in record")
            else:
                notes.append(f"reason={reason} ({reason_text.get(reason, 'see report')})")
            notes.append(f"n_inj={rec['n_injections']}")
            if other is not None and other["kind"] == "recovery_curve":
                notes.append(f"threshold-kind record gives m90={other['flux_limit']['value']} "
                             "(final-candidate fit unsupported)")
            status, lk, lv, lu = "structurally_open", None, None, None
        if "bright_limit_mag" in ex:
            notes.append(f"bright_limit={round(ex['bright_limit_mag'], 3)} AB "
                         "(single-epoch clip; brighter sources belong to the catalogue layer)")
        if "colour_veto" in ex:
            cv = ex["colour_veto"]
            notes.append(f"colour_veto usable/veto={cv['n_usable']}/{cv['n_veto']} (nu={cv['nu']})")
        if extra_note:
            notes.append(extra_note(rec))
        if ex.get("w34_note"):
            notes.append(ex["w34_note"])
        rows.append(row(target_id=ep, rung=zstr(iv), band=band, substrate=substrate,
                        status=status, limit_kind=lk, limit_value=lv, limit_unit=lu,
                        source=f"{records_path} {rec['constraint_id']}",
                        notes="; ".join(notes)))
    return rows


def bright_limit_rows(recs, records_path, substrate, report):
    bl = collections.defaultdict(dict)
    for r in recs:
        bl[r["band"]][(r["endpoint_id"], r["role"])] = r["extra"]["bright_limit_mag"]
    rows = []
    for band in sorted(bl):
        vals = list(bl[band].values())
        lo, hi = min(vals), max(vals)
        rows.append(row(target_id="programme", rung="550-10000 AU", band=band, substrate=substrate,
                        status="structurally_open", limit_kind=None, limit_unit="AB mag",
                        limit_range=f"{lo:.2f}–{hi:.2f}", n_events=len(vals),
                        source=f"{report} (layered-search rule); {records_path} extra.bright_limit_mag",
                        notes="Bright regime not searched by the stack layer: sources brighter than the "
                              "cell's single-epoch 5σ clip limit are clipped from the stack and belong to "
                              "the catalogue-screening layer; limit_range is the min–max of the per-cell "
                              "bright_limit_mag values in the records (one per endpoint × role cell)."))
    return rows


def not_visible_rows(portfolio, subset, report, records_path, scope_text):
    missing = sorted(set(portfolio) - set(subset))
    return [row(target_id=missing, rung="550-10000 AU", band="any", substrate=None,
                status="no_survey", n_events=len(missing),
                source=f"{report} §Scope; set difference of endpoint ids in "
                       f"runs/wise/v2/records/constraint.jsonl (88-endpoint portfolio) and {records_path}",
                notes=scope_text)]


def write_yaml(path, header, rows):
    lines = []
    for k, v in header.items():
        lines.append(f"{k}: {json.dumps(v, ensure_ascii=False)}")
    lines.append("rows:")
    for r in rows:
        lines.append("  - " + json.dumps(r, ensure_ascii=False))
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return len(rows)


# ---------------------------------------------------------------------------
# WISE
# ---------------------------------------------------------------------------
def wise():
    rp = "runs/wise/v2/records/constraint.jsonl"
    recs = load(rp)
    report = "report/wise_survey.md"
    sub = "WISE/NEOWISE L1b frame stack (capped inverse-variance stack on the 64 × 5 × 5 z–µ grid)"
    header = dict(
        survey_id="wise_survey", report=report, plan_section="4.1", pipeline="A",
        archives=["wise"],
        hypothesis_version="v2.1 (wise-hypotheses-v2.1; freeze configs/v2_0_freeze.json sha256:3c83e91e…; physics cell unchanged from v1.0)",
        era="2010 → 2024 (as stated for the positive-control frames, §7; covariance epochs 2010 / 2017 / 2024)",
        run_dir="runs/wise/v2", records_dir="runs/wise/v2/records", candidates=0,
        candidate_notes="Confirmatory (blind) 61 endpoints / 488 cells: 0 candidates at family-wise α = 0.05 "
                        "(R̃_FWER 2.444; highest cell gj-1087/rx/W1 R̃ 1.33, global p 0.960). Development 27 / 216: "
                        "0 (R̃_FWER 2.390; highest lalande-21185/tx/W2 R̃ 1.03). candidate_confirmatory.jsonl empty. "
                        "Supersedes the withdrawn v1–v3 reports.")
    rows = []
    # Aggregate per-band rows from the report §4 table (confirmatory set, Vega, median over cells × intervals).
    agg = {
        "W1": dict(p=13.44, fl=12.71, vi=12.53, bl=12.69, worst=12.5, dev=13.46, n=784),
        "W2": dict(p=12.28, fl=11.44, vi=11.12, bl=11.36, worst=11.1, dev=12.41, n=888),
        "W3": dict(p=9.70, fl=8.76, vi=7.94, bl=8.23, worst=7.9, dev=10.04, n=751),
        "W4": dict(p=7.35, fl=6.30, vi=4.28, bl=5.66, worst=4.3, dev=7.48, n=880),
    }
    for band, a in agg.items():
        notes = (f"Confirmatory-set median 90 %-completeness (cells × intervals), persistent source; "
                 f"flicker {a['fl']} / visit {a['vi']} / long block {a['bl']}; worst of four (duty ≥ 0.5 coverage) {a['worst']}; "
                 f"development-set persistent {a['dev']}. n_events = cell × interval count with a valid fit "
                 f"(surveys/wise/results/report_tables.md). Threshold and final-candidate curves coincide. "
                 f"Report states the calibrated survey is ≈ 1–1.5 mag shallower than the withdrawn v1 claim (v1 depths overstated).")
        if band in ("W3", "W4"):
            notes += (" W3/W4: raw threshold sensitivity for a 300 K blackbody with the empirical PRF "
                      "(completeness_kind = threshold, exclusion_claim = false) — no exclusion of warm structures is claimed; "
                      "confirmation procedure never exercised on a real exceedance.")
        rows.append(row(target_id="programme", rung="550-10000 AU", band=band, substrate=sub, status="searched",
                        limit_kind="m90", limit_value=a["p"], limit_unit="Vega mag", n_events=a["n"],
                        source=f"{report} §4 table; surveys/wise/results/report_tables.md", notes=notes))
    for band in ("W3", "W4"):
        rows.append(row(target_id="programme", rung="550-10000 AU", band=band, substrate=sub,
                        status="structurally_open",
                        source=f"{report} §4, §10 (open)",
                        notes="Long-block temporal model: W3/W4 long-block injections frequently have no epoch on at all "
                              "(cell spans one or two cryo visits); W3/W4 cannot constrain long-block intermittency "
                              "('W3/W4 long-block coverage is nil by cadence')."))
    rows.append(row(target_id="programme", rung="550-10000 AU", band="any", substrate=sub, status="structurally_open",
                    n_events=70, source=f"{report} §1, §2; surveys/wise/results/report_tables.md",
                    notes="Confirmatory set: 70 of 488 cells void (53 heavy-tail; by band W1 24, W2 11, W3 24, W4 11) — "
                          "ring null heavy-tailed (max > 2.5 × q95) or radius-dependent (inner vs outer KS α = 0.01); "
                          "per-cell rows carry reason=null_unstable."))
    rows.append(row(target_id="programme", rung="550-10000 AU", band="any", substrate=sub, status="structurally_open",
                    n_events=27, source=f"{report} §1, §3; surveys/wise/results/report_tables.md",
                    notes="Development set: 27 of 216 cells void (15 heavy-tail; by band W1 8, W2 5, W3 8, W4 6); "
                          "per-cell rows carry reason=null_unstable."))
    reason_text = {
        "null_unstable": "void null: ring heavy-tailed or radius-dependent",
        "insufficient_recovery_fit": "logistic fit extrapolated outside the injected range",
    }
    rows += record_rows(recs, rp, sub, "Vega mag", reason_text)
    return header, rows, recs


# ---------------------------------------------------------------------------
# ZTF
# ---------------------------------------------------------------------------
def ztf(portfolio):
    rp = "runs/ztf/v2/records/constraint.jsonl"
    recs = load(rp)
    report = "report/ztf_survey.md"
    sub = "ZTF difference/science hybrid search-image stack (192 × 5 × 5 grid, T0 = 59800)"
    header = dict(
        survey_id="ztf_survey", report=report, plan_section="4.2", pipeline="A", archives=["ztf"],
        hypothesis_version="v2.0 (ztf-hypotheses-v2.0; freeze configs/v2_freeze.json sha256:9daed5f2…; physics cell as v1.0)",
        era=None, run_dir="runs/ztf/v2", records_dir="runs/ztf/v2/records", candidates=0,
        candidate_notes="Confirmatory (blind, once) 45 endpoints / 228 cells: 0 candidates at family-wise α = 0.05 "
                        "(R̃_FWER 1.658; 19 cells R > 1 vs 26.3 expected). Development 24 / 102: 0 (R̃_FWER 1.872). "
                        "Supersedes the v1 report.")
    rows = []
    agg = {
        "zg": dict(p=23.21, fl=22.44, vi=22.35, bl=22.43, worst=22.4, n=674, fc="23.20 (662)"),
        "zr": dict(p=23.08, fl=22.29, vi=22.14, bl=22.21, worst=22.1, n=690, fc="23.07 (679)"),
        "zi": dict(p=20.96, fl=20.34, vi=20.26, bl=20.32, worst=20.3, n=329, fc="21.11 (285)"),
    }
    for band, a in agg.items():
        rows.append(row(target_id="programme", rung="550-10000 AU", band=band, substrate=sub, status="searched",
                        limit_kind="m90", limit_value=a["p"], limit_unit="AB mag", n_events=a["n"],
                        source=f"{report} §Completeness table; surveys/ztf/results/report_tables.md",
                        notes=(f"Confirmatory-set median 90 %-completeness (cells × intervals), persistent source (threshold curve); "
                               f"flicker {a['fl']} / visit {a['vi']} / block {a['bl']}; worst of four (duty ≥ 0.5 coverage) {a['worst']}; "
                               f"final-candidate persistent median {a['fc']} (results table; within 0.05 mag of threshold). "
                               f"n_events = cell × interval count with a valid fit (results table, threshold). "
                               f"Confirmatory constraints 3,680 of which 3,356 recovery curves; rest not_constrainable "
                               f"(void null, or m90 would sit brighter than the single-epoch clip).")))
    rows += bright_limit_rows(recs, rp, sub, report)
    rows.append(row(target_id="programme", rung="550-10000 AU", band="any", substrate=sub, status="structurally_open",
                    n_events=4, source=f"{report} §Results; surveys/ztf/results/report_tables.md",
                    notes="Confirmatory set: 4 of 228 cells void (0 heavy-tail; zg 2, zr 2); per-cell rows carry reason=null_unstable."))
    rows.append(row(target_id="programme", rung="550-10000 AU", band="any", substrate=sub, status="structurally_open",
                    n_events=1, source=f"{report} §Results; surveys/ztf/results/report_tables.md",
                    notes="Development set: 1 of 102 cells void (1 heavy-tail; zr 1)."))
    rows += not_visible_rows(portfolio, endpoints(recs), report, rp,
                             "Endpoints of the 88-endpoint portfolio outside the '69-endpoint ZTF-visible subset' "
                             "(report §Scope); the report does not name them — this list is the set difference of "
                             "endpoint ids between the WISE records and the ZTF records. The report gives no reason for their absence "
                             "and the list includes northern endpoints (e.g. lalande-21185, 61-cyg-a/b), so it is not simply a declination cut; "
                             "these endpoints were not searched by this survey.")
    reason_text = {
        "null_unstable": "void null: ring heavy-tailed or radius-dependent",
        "insufficient_recovery_fit": "fit whose m90 would sit brighter than the single-epoch clip limit, i.e. in the catalogue layer's regime",
    }
    rows += record_rows(recs, rp, sub, "AB mag", reason_text)
    return header, rows


# ---------------------------------------------------------------------------
# PS1
# ---------------------------------------------------------------------------
def ps1(portfolio):
    rp = "runs/panstarrs/v2/records/constraint.jsonl"
    recs = load(rp)
    report = "report/ps1_survey.md"
    sub = "Pan-STARRS1 DR2 warp cutout stack (~20 warps per band; 360 × 5 × 5 grid at 1″, T0 = 59800)"
    header = dict(
        survey_id="ps1_survey", report=report, plan_section="4.4", pipeline="A", archives=["ps1"],
        hypothesis_version="v2.0 (ps1-hypotheses-v2.0; freeze sha256:be6f3a1c…)",
        era=None, run_dir="runs/panstarrs/v2", records_dir="runs/panstarrs/v2/records", candidates=0,
        candidate_notes="Confirmatory (blind, once) 45 endpoints / 450 cells: 0 candidates at family-wise α = 0.05 "
                        "(R̃_FWER 1.470; 49 cells R > 1 vs 54.9 expected). Development 24 / 240: 0 (R̃_FWER 1.504). "
                        "Supersedes the v1 report.")
    rows = []
    agg = {
        "g": dict(p=21.71, pn=107, fl="21.53 (9)", vi="21.73 (36)", bl="21.49 (24)"),
        "r": dict(p=21.56, pn=64, fl="21.14 (15)", vi="21.50 (25)", bl="21.41 (4)"),
        "i": dict(p=20.99, pn=71, fl="21.01 (8)", vi="21.01 (17)", bl="20.97 (30)"),
        "z": dict(p=20.46, pn=60, fl="20.38 (32)", vi="20.52 (44)", bl="20.03 (10)"),
        "y": dict(p=19.38, pn=104, fl="19.58 (3)", vi="19.29 (42)", bl="19.35 (12)"),
    }
    for band, a in agg.items():
        rows.append(row(target_id="programme", rung="550-10000 AU", band=band, substrate=sub, status="searched",
                        limit_kind="m90", limit_value=a["p"], limit_unit="AB mag", n_events=a["pn"],
                        source=f"{report} §Completeness table; surveys/panstarrs/results/report_tables.md",
                        notes=(f"Confirmatory-set median 90 %-completeness (cells × intervals), persistent source; "
                               f"flicker {a['fl']} / visit {a['vi']} / block {a['bl']} (parentheses: cell × interval count with a valid fit). "
                               f"Applies only where the stack can constrain at all: 1,336 of 7,200 confirmatory records are recovery curves. "
                               f"Report states v1's uniform '~21 AB' stack-depth claim did not account for the single-epoch clip or the "
                               f"contamination threshold (v1 depth overstated / localised to the cells that support it).")))
    rows += bright_limit_rows(recs, rp, sub, report)
    rows.append(row(target_id="programme", rung="550-10000 AU", band="any", substrate=sub, status="structurally_open",
                    source=f"{report} §Completeness (structural result); surveys/panstarrs/results/report_tables.md §Ledger",
                    notes="Structural result: with ~20 warps per band, in most cells the PS1 stack has no 90 %-complete regime "
                          "between the single-epoch clip and the contamination-driven family threshold; only 1,336 of 7,200 "
                          "confirmatory records are recovery curves; elsewhere the constraint belongs to the catalogue-screening "
                          "layer and the joint stage (per-cell rows carry reason=insufficient_recovery_fit)."))
    rows.append(row(target_id="programme", rung="550-10000 AU", band="any", substrate=sub, status="structurally_open",
                    n_events=8, source=f"{report} §Results; surveys/panstarrs/results/report_tables.md",
                    notes="Confirmatory set: 8 of 450 cells void (0 heavy-tail; inner/outer-ring test; z 1, y 2, g 2, i 2, r 1); "
                          "per-cell rows carry reason=null_unstable."))
    rows.append(row(target_id="programme", rung="550-10000 AU", band="any", substrate=sub, status="structurally_open",
                    n_events=4, source=f"{report} §Results; surveys/panstarrs/results/report_tables.md",
                    notes="Development set: 4 of 240 cells void (0 heavy-tail; i 2, z 1, r 1)."))
    rows += not_visible_rows(portfolio, endpoints(recs), report, rp,
                             "Endpoints of the 88-endpoint portfolio outside the '69-endpoint PS1 subset' (report §Scope); "
                             "the report does not name them — this list is the set difference of endpoint ids between the "
                             "WISE records and the PS1 records. The report gives no reason for their absence and the list includes northern "
                             "endpoints (e.g. lalande-21185, 61-cyg-a/b), so it is not simply a declination cut; these endpoints were not "
                             "searched by this survey.")
    reason_text = {
        "null_unstable": "void null: inner/outer-ring test",
        "insufficient_recovery_fit": "fit whose 90 % point would land brighter than the single-epoch clip (extrapolation into the catalogue layer's regime) or logistic fit lacking support",
    }
    rows += record_rows(recs, rp, sub, "AB mag", reason_text)
    return header, rows


# ---------------------------------------------------------------------------
# Joint PS1 + ZTF + WISE colour axis
# ---------------------------------------------------------------------------
def joint(portfolio):
    rp = "runs/joint/v3/records/constraint.jsonl"
    recs = load(rp)
    report = "report/joint_ps1_ztf_wise.md"
    sub = "Joint PS1 + ZTF stack on the PS1 grid (S_joint = (A_ztf + A_ps1)/√(B_ztf + B_ps1)); WISE W1/W2 stacks as colour-axis annotations"
    header = dict(
        survey_id="joint_ps1_ztf_wise", report=report, plan_section="4.5", pipeline="A",
        archives=["ps1", "ztf", "wise"],
        hypothesis_version="v3.0 (joint-ps1-ztf-wise-v3.0; freeze surveys/joint/configs/v3_freeze.json sha256:bb68869c…; "
                           "PS1 v2.0 + ZTF v2.0 + WISE v3.0-joint-conventions freezes)",
        era=None, run_dir="runs/joint/v3", records_dir="runs/joint/v3/records", candidates=0,
        candidate_notes="Confirmatory (blind, once) 45 endpoints / 230 joint cells: 0 candidates at family-wise α = 0.05 "
                        "(R̃_FWER 1.778; 14 cells R > 1 vs 25.2 expected); optical statistic bit-identical to the v2 run. "
                        "Development 24 / 102: 0 (R̃_FWER 1.540). Colour axis: 230/230 cells usable, 0 would-fire vetoes. "
                        "Supersedes the v2 report (2026-08-23).")
    rows = []
    agg = {
        "g": dict(fc=22.83, thr=23.06, fl="22.36 / 22.22", vi="22.16 / 22.09", bl="22.22 / 22.10", n=607, pair="PS1 g + ZTF zg"),
        "r": dict(fc=22.86, thr=22.93, fl="22.19 / 22.15", vi="22.05 / 21.98", bl="22.09 / 22.06", n=639, pair="PS1 r + ZTF zr"),
        "i": dict(fc=21.09, thr=21.07, fl="21.17", vi="21.02", bl="20.97", n=78, pair="PS1 i + ZTF zi"),
    }
    for band, a in agg.items():
        rows.append(row(target_id="programme", rung="550-10000 AU", band=band, substrate=sub, status="searched",
                        limit_kind="m90", limit_value=a["fc"], limit_unit="AB mag", n_events=a["n"],
                        source=f"{report} §Completeness table; surveys/joint/results/report_tables.md",
                        notes=(f"Band pair {a['pair']}. Confirmatory-set median 90 %-completeness (cells × intervals), persistent source, "
                               f"final-candidate curve {a['fc']} (threshold curve {a['thr']}); flicker {a['fl']} / visit {a['vi']} / block {a['bl']} "
                               f"(threshold / final-candidate). n_events = cell × interval count with a valid final-candidate fit (results table). "
                               f"Unchanged from v2; static veto on 2.9 % of threshold-recovered joint injections, colour veto on none in the fitted region; "
                               f"injections brighter than the fainter of the two archives' single-epoch clip limits excluded from the fits.")))
    rows.append(row(target_id="programme", rung="550-10000 AU", band="W1/W2 (colour axis)", substrate="WISE W1/W2 joint-convention stacks (runs/wise/v3; T0 59800, AB ZP 25; Vega→AB W1 +2.699 / W2 +3.339)",
                    status="ledger_only", n_events=230,
                    source=f"{report} §Summary, §Colour veto calibration, §Results",
                    notes="WISE enters only as per-cell W1/W2 stack annotations at every joint peak and a calibrated flat-Fν colour-consistency "
                          "veto (deficit ≥ 5σ in every usable W band; σ_w floored by the 48-ring scatter). Confirmatory: 230/230 cells both W bands "
                          "usable (n ≥ 5); deficit statistic at real peaks median −0.33 to max 3.73 (no cell approaches ν = 5); 0 would-fire. "
                          "Veto fired on 5 / 46,302 threshold-recovered confirmatory injections (1 / 23,117 development), all at mag 16.8–18.8 "
                          "above the optical bright limit (21.0–21.9); false-veto rate in the fitted population 0/46,302 and 0/23,116. "
                          "Not a depth contribution; W1/W2 depth statements stay with wise_survey."))
    for band in ("W3", "W4"):
        rows.append(row(target_id="programme", rung="550-10000 AU", band=band, substrate=None, status="no_survey",
                        source=f"{report} §Scope and provenance notes",
                        notes="W3/W4 take no part in the joint stage (scientific-review rule); threshold statements stay with the standalone WISE survey."))
    rows.append(row(target_id="programme", rung="550-10000 AU", band="any", substrate=sub, status="structurally_open",
                    n_events=5, source=f"{report} §Results; surveys/joint/results/report_tables.md",
                    notes="Confirmatory set: 5 of 230 cells void (0 heavy-tail; r 2, g 2, i 1). The joint records carry no "
                          "not_constrainable_reason field, so void cells are not flagged per row here."))
    rows.append(row(target_id="programme", rung="550-10000 AU", band="any", substrate=sub, status="structurally_open",
                    n_events=1, source=f"{report} §Results; surveys/joint/results/report_tables.md",
                    notes="Development set: 1 of 102 cells void (1 heavy-tail; r 1)."))
    rows += not_visible_rows(portfolio, endpoints(recs), report, rp,
                             "Endpoints of the 88-endpoint portfolio absent from the joint PS1 + ZTF records (the 69-endpoint "
                             "PS1/ZTF subset); the report does not name them — this list is the set difference of endpoint ids "
                             "between the WISE records and the joint records. The report gives no reason for their absence; these endpoints "
                             "were not searched by this stage.")
    rows += record_rows(recs, rp, sub, "AB mag", {})
    return header, rows


def main():
    os.makedirs(OUT, exist_ok=True)
    h, r, wrecs = wise()
    portfolio = endpoints(wrecs)
    assert len(portfolio) == 88
    counts = {}
    counts["wise_survey"] = write_yaml(os.path.join(OUT, "wise_survey.yaml"), h, r)
    h, r = ztf(portfolio)
    counts["ztf_survey"] = write_yaml(os.path.join(OUT, "ztf_survey.yaml"), h, r)
    h, r = ps1(portfolio)
    counts["ps1_survey"] = write_yaml(os.path.join(OUT, "ps1_survey.yaml"), h, r)
    h, r = joint(portfolio)
    counts["joint_ps1_ztf_wise"] = write_yaml(os.path.join(OUT, "joint_ps1_ztf_wise.yaml"), h, r)
    print(counts)


if __name__ == "__main__":
    main()
