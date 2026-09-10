"""Build harvest YAML files for the five sunward Pipeline-B reports.

Every number is taken from the report text or from the results files the
report cites (loaded here by path). Nothing is estimated.
"""
import json
import math
import yaml

ROOT = "/home/jgreene/code/sgl-seti-survey"
import sys
OUT = sys.argv[1] if len(sys.argv) > 1 else f"{ROOT}/surveys/programme-ledger/harvest"


def sig3(x):
    """Round to 3 significant figures (numbers are quoted from the results files)."""
    if x is None:
        return None
    if x == 0:
        return 0.0
    return float(f"{x:.3g}")


def mw(w):
    return sig3(w / 1e6) if w is not None else None


def row(**kw):
    base = dict(
        target_id=None, channel=None, rung=None, band="any", substrate=None,
        status=None, limit_kind=None, limit_value=None, limit_unit=None,
        limit_range=None, power_mw=None, power_mw_range=None, n_events=None,
        n_trials=None, epoch_range=None, source=None, notes=None,
    )
    base.update(kw)
    return base


def dump(doc, stem):
    path = f"{OUT}/{stem}.yaml"
    with open(path, "w") as f:
        yaml.safe_dump(doc, f, sort_keys=False, allow_unicode=True, width=1000)
    print(stem, "rows:", len(doc["rows"]), "->", path)


# ---------------------------------------------------------------------------
# 1. LASCO
# ---------------------------------------------------------------------------
def lasco():
    R = "report/lasco_crossings.md"
    RES = "surveys/heliospheric-crossings/results"
    pl = json.load(open(f"{ROOT}/{RES}/power_limits_v1.json"))
    cov = json.load(open(f"{ROOT}/{RES}/coverage_v1.json"))["aggregates"]
    rows = []
    sub_c2 = "SOHO/LASCO C2 level-0.5 synoptic frames, forced ring-differential photometry"
    sub_c3 = "SOHO/LASCO C3 level-0.5 synoptic frames, forced ring-differential photometry"

    # S1 1.2 Rsun: not_constrainable behind the C2 occulter (per-target coverage counts)
    for t, a in cov["S1_1.2Rsun"].items():
        rows.append(row(
            target_id=t, channel="S1", rung="1.2 Rsun", band="C2 Orange", substrate=sub_c2,
            status="not_constrainable", n_events=a["events"],
            source=f"{R} §1, §6; {RES}/coverage_v1.json aggregates S1_1.2Rsun",
            notes=f"{a['covered']} of {a['events']} windows covered, 0 covered+visible: apparent source stays inside the C2 usable inner radius 2.2 Rsun (occulter); ledger nominal-covered / occulter-unusable",
        ))

    # S1 2.5 Rsun confirmatory (C2): unit table + power_limits
    conf25 = {
        "van-maanen": ("−1.84/0.49", "4.51/1.19 exc", "12.8/20.5", "searched; S_event adjudicated CME-period systematic (see retained_ambiguous row)"),
        "wolf-359": ("0.08/0.31", "0.11/1.20", "6.6/14.2", "searched-null"),
        "teegarden": ("0.62/2.33", "3.70/8.67", "13.3/24.1", "searched-null"),
        "gj-1276": ("0.25/0.60", "0.76/1.70", "18.7/17.4 exc", "searched; S_pulse exceedance vetoed single-frame/cosmic-ray class (see pulse row)"),
    }
    for t, (ss, se, sp, disp) in conf25.items():
        u = pl[f"{t}|S1_2.5Rsun"]
        ev = u["S_event"]
        rows.append(row(
            target_id=t, channel="S1", rung="2.5 Rsun", band="C2 Orange", substrate=sub_c2,
            status="searched", limit_kind="m90", limit_value=u["S_stack"]["m90_V"], limit_unit="V mag",
            power_mw=mw(u["S_stack"]["P_cone_W"]), n_events=cov["S1_2.5Rsun"][t]["covered_and_visible"], n_trials=3,
            source=f"{RES}/confirmatory_v1.md unit table; {RES}/power_limits_v1.json; {R} §3–4",
            notes=f"S_stack {ss}, S_event {se}, S_pulse {sp} (S/T); {disp}; S_stack m90 -> P_cone downlink (W_eff 100 nm C2 Orange); S_event m90 {ev['m90_V']} V -> {mw(ev['P_cone_W'])} MW; wing regime [2.2, 2.5] Rsun; n_events = covered+visible windows (coverage_v1.json); ±~0.3 mag band-conversion systematic declared",
        ))
    # gj-1276 pulse cell (report headline)
    u = pl["gj-1276|S1_2.5Rsun"]["S_pulse"]
    rows.append(row(
        target_id="gj-1276", channel="S1", rung="2.5 Rsun", band="C2 Orange", substrate=sub_c2 + " (full cadence, S_pulse)",
        status="searched", limit_kind="m90", limit_value=u["m90_V"], limit_unit="V mag",
        power_mw=mw(u["P_cone_W"]), n_trials=1,
        source=f"{R} §4 table (pulse cell), §3; {RES}/power_limits_v1.json gj-1276|S1_2.5Rsun S_pulse",
        notes="Pulse cell: >= 1 frame / 12 min, power quoted by the report as 6.7 MW·(12 min); S_pulse 18.7 vs T 17.4 exceedance = one 12-min frame 2014-09-04 07:12 (+2078, neighbours −470/−264), vetoed under the frozen persistence rule (single-frame / cosmic-ray class)",
    ))
    rows.append(row(
        target_id="van-maanen", channel="S1", rung="2.5 Rsun", band="C2 Orange", substrate=sub_c2,
        status="retained_ambiguous", n_events=1, n_trials=1,
        epoch_range="2020-10-05/06",
        source=f"{R} §3; {RES}/confirmatory_v1.md adjudication 2",
        notes="S_event 4.51 vs T 1.19: one window, ~3 h all-position-angle annulus disturbance whose onset follows a CDAW-catalogued C2 CME (20:48) by 48 min; both ±25° rings swing by hundreds in both signs; no recurrence in 27 sibling windows (S_stack −1.84); adjudicated CME-period systematic, retained, non-promotable",
    ))
    # ross-128 S1 2.5 Rsun: dev unit
    rows.append(row(
        target_id="ross-128", channel="S1", rung="2.5 Rsun", band="C2 Orange", substrate=sub_c2,
        status="constraint_only", n_events=29, n_trials=3,
        source=f"{R} §3 (dev stage); {RES}/dev_search_v11.json ross-128|S1_2.5Rsun",
        notes="Dev-stage unit (no injection depth reported): S_stack −1.35/0.39, S_event 1.45/1.20 exc, S_pulse 16.3/49.6; the marginal S_event on the 2015 St. Patrick's Day CME storm window adjudicated (report §3)",
    ))

    # S1 0.1 AU confirmatory (C3)
    conf_s1 = {
        "van-maanen": ("18.3/35.7", "13.1/47.3"),
        "wolf-359": ("12.1/12.9", "22.6/25.7"),
        "teegarden": ("2.0/21.0", "5.9/54.2"),
        "gj-1276": ("3.6/14.7", "5.0/30.9"),
    }
    for t, (ss, se) in conf_s1.items():
        u = pl[f"{t}|S1_0.1AU"]
        rows.append(row(
            target_id=t, channel="S1", rung="0.1 AU", band="C3 Clear", substrate=sub_c3,
            status="searched", limit_kind="m90", limit_value=u["S_stack"]["m90_V"], limit_unit="V mag",
            power_mw=mw(u["S_stack"]["P_cone_W"]), n_events=cov["S1_0.1AU"][t]["covered_and_visible"], n_trials=2,
            source=f"{RES}/confirmatory_v1.md unit table; {RES}/power_limits_v1.json; {R} §3–4",
            notes=f"searched-null: S_stack {ss}, S_event {se} (S/T); S_stack m90 -> P_cone downlink (line-equivalent, W_eff 300 nm C3 Clear; 532 nm in band); S_event m90 {u['S_event']['m90_V']} V -> {mw(u['S_event']['P_cone_W'])} MW; usable annulus [4.4, 29] Rsun; n_events = covered+visible windows",
        ))
    rows.append(row(
        target_id="ross-128", channel="S1", rung="0.1 AU", band="C3 Clear", substrate=sub_c3,
        status="constraint_only", n_events=cov["S1_0.1AU"]["ross-128"]["covered_and_visible"],
        source=f"{R} §3; {RES}/confirmatory_v1.md unit table",
        notes="Permanently blended: the fixed S1 sky position (356.94°, −0.79°) lies 121″ (2.2 px) from a VT 7.2 Tycho star inside the frozen 3 px mask at every epoch; C3's 56″ pixels cannot separate it; ledger nominal-covered / resolution-blended",
    ))

    # S2 0.1 AU confirmatory (C3)
    conf_s2 = {
        "van-maanen": ("8.1/22.1", "10.8/44.7", "stellar template applied"),
        "wolf-359": ("1.9/29.5", "4.9/35.9", None),
        "teegarden": ("4.5/19.3", "5.5/24.1", None),
        "gj-1276": ("11.4/18.3", "14.3/25.8", None),
        "ross-128": ("−0.5/12.6", "10.0/34.8", "stellar template applied"),
    }
    for t, (ss, se, tmpl) in conf_s2.items():
        u = pl[f"{t}|S2_0.1AU"]
        rows.append(row(
            target_id=t, channel="S2", rung="0.1 AU", band="C3 Clear", substrate=sub_c3,
            status="searched", limit_kind="m90", limit_value=u["S_stack"]["m90_V"], limit_unit="V mag",
            power_mw=mw(u["S_stack"]["P_tx_10m_W"]), n_events=cov["S2_0.1AU"][t]["covered_and_visible"], n_trials=2,
            source=f"{RES}/confirmatory_v1.md unit table; {RES}/power_limits_v1.json; {R} §3–4",
            notes=f"searched-null{' (' + tmpl + ')' if tmpl else ''}: S_stack {ss}, S_event {se} (S/T); power_mw = P_tx_10m (10-m-class uplink transmitter at registry distance) from the results file; P_cone {mw(u['S_stack']['P_cone_W'])} MW; S_event m90 {u['S_event']['m90_V']} V -> P_tx_10m {mw(u['S_event']['P_tx_10m_W'])} MW",
        ))

    # dev units on the 0.1 AU rungs (no injection depth)
    dev = json.load(open(f"{ROOT}/{RES}/dev_search_v11.json"))["units"]
    dev_notes = {
        "gj-908|S1_0.1AU": "Dev-stage unit (no injection depth reported)",
        "gj-908|S2_0.1AU": "Dev-stage unit (no injection depth reported); gj-908 (V 8.98) detected window-locked at the predicted S2 position and nulled by the stellar template — the survey's in-situ positive control",
        "ross-154|S1_0.1AU": "Dev-stage unit (no injection depth reported); ross-154's 8-yr Venus-synodic stray-light contamination found in dev and vetoed via the bright-planet in-FOV veto (channel of the contamination not stated in the report)",
        "ross-154|S2_0.1AU": "Dev-stage unit (no injection depth reported); ross-154's 8-yr Venus-synodic stray-light contamination found in dev and vetoed via the bright-planet in-FOV veto (channel of the contamination not stated in the report)",
    }
    for k, note in dev_notes.items():
        t, cell = k.split("|")
        ch = cell.split("_")[0]
        u = dev[k]
        ss = f"{u['S_stack']['S']:.2f}/{u['S_stack']['T']:.2f}"
        se = f"{u['S_event']['S']:.2f}/{u['S_event']['T']:.2f}"
        rows.append(row(
            target_id=t, channel=ch, rung="0.1 AU", band="C3 Clear", substrate=sub_c3,
            status="constraint_only", n_events=u["per_dpa"]["dpa+0"]["n_events"], n_trials=2,
            source=f"{R} §3 (dev stage); {RES}/dev_search_v11.json {k}",
            notes=f"{note}; dev v1.1 S_stack {ss}, S_event {se} (S/T)",
        ))

    rows.append(row(
        target_id="programme", channel="S2", rung="1.0 AU", band="any", substrate="SOHO/LASCO",
        status="no_survey", source=f"{R} §1",
        notes="S2 1.0 AU declared out of scope (non-exclusive twilight skirt); no coverage count given in the report",
    ))

    doc = dict(
        survey_id="lasco_crossings", report=R, plan_section="5.11", pipeline="B", archives=["lasco"],
        hypothesis_version="v1.0 (frozen 2026-08-25) + threshold amendments v1.1/v1.2",
        era="1996-01 → 2026-06", run_dir="runs/heliospheric-crossings", records_dir=None,
        candidates=0,
        candidate_notes="30 trials / 13 units, 2 exceedances vs 3.3 expected: gj-1276 S1 2.5 Rsun S_pulse vetoed (single-frame/cosmic-ray); van-maanen S1 2.5 Rsun S_event adjudicated CME-period systematic, retained non-promotable",
        rows=rows,
    )
    dump(doc, "lasco_crossings")


# ---------------------------------------------------------------------------
# 2. STEREO-A HI-1
# ---------------------------------------------------------------------------
def stereo():
    R = "report/stereo_hi_crossings.md"
    RES = "surveys/stereo-hi-crossings/results"
    pl = json.load(open(f"{ROOT}/{RES}/power_limits_v1.json"))
    sub = "STEREO-A HI-1 level-2 24h1A_br01 frames, ICRS-fixed sky-patch aperture photometry with star-fixed baseline"
    band = "HI-1 630–730 nm"
    rows = []
    for rung in ("1.2 Rsun", "2.5 Rsun"):
        rows.append(row(
            target_id="programme", channel="S1", rung=rung, band=band, substrate=sub,
            status="not_constrainable", source=f"{R} §1, §5",
            notes="Grazing rungs (ε ≤ 0.7°) never enter the HI-1 field (inner edge |HPLN| 3.9°): 84 of 105 ledger events across the grazing rungs are not_constrainable; closed by geometry, not by data",
        ))
    # confirmatory table (report §3)
    conf = {
        ("gj-1276", "S1"): (20, "0.11 / 9.66", "4.06 / 7.75", "5.73 / 14.65", "searched-null"),
        ("gj-1276", "S2"): (18, "−0.21 / 17.64", "8.22 / 7.75 exc", "15.71 / 11.16 exc", "S_pulse 2009-01-14 four consecutive frames = ~10-px diagonal band sweeping 2–3 px/frame over 2.5 h, extended moving front (CME/streamer class), vetoed non-point-source; S_event plateau retained (see retained_ambiguous row)"),
        ("ross-128", "S2"): (19, "−9.16 / 12.41", "3.43 / 37.20", "8.43 / 73.08", "searched-null"),
        ("teegarden", "S1"): (17, "−6.66 / 3.72", "3.58 / 5.78", "19.57 / 7.71 exc", "S_pulse 2010-08-14 last two arc frames (+6.4, +10.0, back to −2.6): whole stamp lifts uniformly ~0.8 DN/s, extended transient front, vetoed non-point-source"),
        ("teegarden", "S2"): (18, "1.54 / 10.94", "6.45 / 10.91", "10.59 / 29.22", "searched-null"),
        ("van-maanen", "S1"): (18, "1.51 / 21.08", "4.81 / 15.73", "9.68 / 106.35", "searched-null"),
        ("van-maanen", "S2"): (18, "−1.44 / 9.83", "3.75 / 12.54", "3.94 / 18.28", "searched-null"),
        ("wolf-359", "S1"): (18, "−3.27 / 6.85", "5.62 / 6.51", "5.53 / 13.31", "searched-null"),
        ("wolf-359", "S2"): (20, "2.88 / 13.94", "6.82 / 14.04", "5.68 / 14.70", "searched-null"),
    }
    for (t, ch), (nev, ss, se, sp, disp) in conf.items():
        u = pl[f"{t}|{ch}_0.1AU"]
        st = u["S_stack"]
        if ch == "S1":
            p = mw(st["P_cone_W"]); pnote = "P_cone downlink through the 0.1 AU cone"
        else:
            p = mw(st["P_tx_10m_W"]); pnote = f"P_tx_10m (10-m-class uplink); P_cone {mw(st['P_cone_W'])} MW"
        extra = []
        if "S_event" in u:
            e = u["S_event"]
            extra.append(f"S_event m90 {e['m90_Veq']} -> {mw(e['P_tx_10m_W'] if ch == 'S2' else e['P_cone_W'])} MW")
        else:
            extra.append("S_event depth degenerate (null exceeds T)")
        if "S_pulse" in u:
            q = u["S_pulse"]
            extra.append(f"S_pulse m90 {q['m90_Veq']} -> {mw(q['P_tx_10m_W'] if ch == 'S2' else q['P_cone_W'])} MW (40-min frames)")
        rows.append(row(
            target_id=t, channel=ch, rung="0.1 AU", band=band, substrate=sub,
            status="searched", limit_kind="m90", limit_value=st["m90_Veq"], limit_unit="V_eq mag (B−V 0.65 source)",
            power_mw=p, n_events=nev, n_trials=3,
            source=f"{R} §3 table; {RES}/power_limits_v1.json {t}|{ch}_0.1AU; {RES}/confirmatory_v1.md",
            notes=f"S_stack {ss}, S_event {se}, S_pulse {sp} (S/T); {disp}; power_mw = {pnote}; {'; '.join(extra)}; line-equivalent through W_eff 100 nm at 680 nm, ±0.3 mag band systematic; 532 nm out of band (broadband leakage cell)",
        ))
    rows.append(row(
        target_id="gj-1276", channel="S2", rung="0.1 AU", band=band, substrate=sub,
        status="retained_ambiguous", n_events=1, n_trials=1, epoch_range="2012-10-22/24",
        source=f"{R} §3; {RES}/confirmatory_v1.md adjudication",
        notes="S_event 8.22 vs T 7.75: a +0.39-unit plateau over 54 epochs (V_eq 12.1) with no point source in the pixels at any epoch (gj-1276 is V 16); patch's own baseline sits at −0.27 units — a regional level-2 residual-background offset between baseline and arc regions the flank controls do not share; null in the unit's other 17 windows; retained, non-promotable",
    ))
    rows.append(row(
        target_id="ross-128", channel="S1", rung="0.1 AU", band=band, substrate=sub,
        status="constraint_only", n_events=0,
        source=f"{R} §3 table, §6; {RES}/confirmatory_v1.md",
        notes="Antipode patch Tycho-masked: VT 7.2 star at 1.7 px — the same permanently blended fixed sky point LASCO found; 0 events searchable",
    ))
    # dev units (report §3 dev stage; results/dev_v1.md)
    devs = {
        ("gj-908", "S1"): (19, "−1.32 / 6.61", "7.80 / 7.46 exc", "4.23 / 6.40", "constraint_only", "Dev-stage unit (no injection depth): S_event exceedance 2022-10 extended-background plateau (same class as gj-1276 S2); gj-908 (V 9.0) repeats to 0.2 % across 18 transits — the in-situ bright-star control"),
        ("gj-908", "S2"): (18, "6.62 / 9.36", "7.77 / 14.66", "6.13 / 21.76", "constraint_only", "Dev-stage unit (no injection depth); gj-908 (V 9.0) repeats to 0.2 % across 18 transits — the in-situ bright-star control"),
        ("ross-154", "S1"): (9, "1.37 / 7.25", "3.97 / 4.14", "4.42 / 8.96", "constraint_only", "Dev-stage unit (no injection depth); ross-154's crowded west-era fields fail the ZP-scatter gate (MAD ≤ 0.12): dev units keep 9–10 of 20 events"),
        ("ross-154", "S2"): (10, "6.67 / 8.50", "4.56 / 10.92", "9.41 / 5.13 exc", "vetoed_known_source", "Dev-stage unit (no injection depth): S_pulse exceedance 2010-10-06 single-frame +6.4-unit flare-class event on V1216 Sgr, vetoed by the persistence rule; crowded west-era fields fail the ZP-scatter gate (9–10 of 20 events kept)"),
    }
    for (t, ch), (nev, ss, se, sp, status, note) in devs.items():
        rows.append(row(
            target_id=t, channel=ch, rung="0.1 AU", band=band, substrate=sub,
            status=status, n_events=nev, n_trials=3,
            source=f"{R} §3 (dev stage); {RES}/dev_v1.md unit table",
            notes=f"{note}; dev S_stack {ss}, S_event {se}, S_pulse {sp} (S/T)",
        ))
    rows.append(row(
        target_id="programme", channel="S2", rung="1.0 AU", band=band, substrate="STEREO-A HI-1 whole-transit in-beam passages",
        status="ledger_only", n_events=289, source=f"{R} §1",
        notes="S2 1 AU (whole HI-1 transit in-beam): 16 targets, 289 transits; out of scope by decision D7 — coverage logged, no statistic",
    ))
    doc = dict(
        survey_id="stereo_hi_crossings", report=R, plan_section="5.16", pipeline="B", archives=["stereo-hi1"],
        hypothesis_version="v1.0 + threshold amendments v1.1/v1.2 (dev, pre-confirmatory)",
        era="2007-01 → 2026-08", run_dir="runs/stereo-hi-crossings", records_dir=None,
        candidates=0,
        candidate_notes="27 trials / 9 units, 3 exceedances vs 3.0 expected, all adjudicated non-point-source systematics (two coronal-transient fronts vetoed, one background-offset plateau retained non-promotable); coverage 272 of 295 event rows (2014-08 → 2015-11 safe-mode/conjunction gap takes 19)",
        rows=rows,
    )
    dump(doc, "stereo_hi_crossings")


# ---------------------------------------------------------------------------
# 3. PSP/WISPR
# ---------------------------------------------------------------------------
def wispr():
    R = "report/wispr_crossings.md"
    RES = "surveys/wispr-crossings/results"
    pl = json.load(open(f"{ROOT}/{RES}/power_limits_v1.json"))
    sub = "PSP/WISPR-I level-3 synoptic frames (E1–E27), ICRS-fixed patch aperture photometry with Tycho/Hipparcos stellar template"
    band = "WISPR-I λ_c 615 nm, W_eff 250 nm (532 nm in band)"
    rows = []
    for rung in ("1.2 Rsun", "2.5 Rsun"):
        rows.append(row(
            target_id="programme", channel="S1", rung=rung, band=band, substrate=sub,
            status="not_constrainable", source=f"{R} §1",
            notes="The grazing cones never enter the WISPR-I field (≤ 6.9° / 13.7° elongation vs the 13.5° inner edge)",
        ))
    # confirmatory unit table (results/confirmatory_v1.md)
    conf = {
        "ez-aqr:S2": (19, "11.07 / 12.54", "21.04 / 26.15", "32.25 / 63.44", "searched", "searched-null"),
        "gj-1002:S1": (11, "16.61 / 18.24", "12.46 / 12.94", "2.85 / 6.22", "searched", "searched-null"),
        "gj-1002:S2": (24, "0.39 / 22.61", "6.55 / 43.77", "20.11 / 89.09", "searched", "searched-null"),
        "gj-1087:S2": (10, "33.34 / 17.87 exc", "72.13 / 8.76 exc", "3.84 / 2.51 exc", "vetoed_known_source", "static content (template incompleteness): z > 0 in all 10 events, Gaia DR3 predicts 0.61 units of static flux in the aperture (two G 10.3–10.8 stars 1.6 px off-centre absent from the Tycho extract) vs template 0.20 — excess 0.42 is the missing static flux; the pulse is an arc-exit ramp (last two frames of a 16-frame arc); no S_stack/S_event depth (degenerate)"),
        "gj-1111:S1": (17, "−7.65 / 13.37", "5.92 / 30.72", "22.11 / 27.24", "searched", "searched-null"),
        "gj-1276:S2": (23, "4.89 / 22.83", "6.34 / 48.27", "18.84 / 80.23", "searched", "searched-null (LASCO/HI-1 family)"),
        "gj-251:S1": (18, "3.77 / 9.16", "7.72 / 26.41", "6.34 / 24.85", "searched", "searched-null"),
        "gj-514:S1": (26, "−11.51 / 30.00", "3.38 / 25.48", "12.62 / 75.27", "searched", "searched-null"),
        "gj-514:S2": (19, "6.06 / 11.74", "3.99 / 13.44", "5.76 / 4.25 exc", "searched", "S_pulse exceedance = arc-exit ramp (4-frame ramp in the last frames of the E26 arc at the 53.5° field exit), curvature-corrected z null, non-recurrent"),
        "gj-518:S1": (24, "19.00 / 27.71", "16.23 / 28.43", "94.60 / 224.04", "searched", "searched-null"),
        "gj-518:S2": (22, "−1.85 / 10.88", "4.33 / 8.80", "5.97 / 12.01", "searched", "searched-null"),
        "gj-54:S2": (26, "8.06 / 12.56", "3.60 / 34.32", "18.70 / 42.00", "searched", "searched-null"),
        "gj-581:S1": (20, "24.99 / 21.17 exc", "20.88 / 51.75", "54.28 / 68.50", "vetoed_known_source", "latitude curvature + static content: z +5.0 median; curvature-corrected S_stack 11 < T; residual = faint Gaia stars (G 12.5) below the Tycho limit; S_stack depth degenerate"),
        "gj-581:S2": (23, "14.56 / 27.16", "8.88 / 50.37", "17.89 / 47.53", "searched", "searched-null"),
        "gj-667-c:S1": (9, "44.20 / 10.25 exc", "21.32 / 5.22 exc", "1.57 / 3.58", "vetoed_known_source", "static content (bright-neighbour template error): z +14.5 in every event; Gaia 4.39 units vs response-scaled template 2.13 (a G 7.0 star 2.7 px off-centre carried at the wrong encircled fraction/magnitude); no depth (S_stack/S_event degenerate; S_pulse not tabulated)"),
        "gj-667-c:S2": (21, "0.60 / 15.10", "19.61 / 16.94 exc", "11.46 / 24.76", "searched", "S_event exceedance = one event (E01, V3-processed first encounter) level offset on a patch dominated by the HD 156384 AB pair (V 6.7, 35″; Gaia 6.1 vs template 14.2, pair over-counted) — bright-neighbour template systematic; other 20 events z median −2.4; S_event depth degenerate"),
        "gj-783:S2": (19, "19.32 / 10.77 exc", "34.00 / 29.73 exc", "17.00 / 8.90 exc", "vetoed_known_source", "bright-star class (declared v1.3 limit): V 5.3 star; template counted the 1.6″/yr star more than once (Tycho positions not PM-propagated); S_event/S_stack from per-encounter calibration scatter of a S/N-100 star; the pulse is a 2-frame 17σ spike of the star itself (E14); no depth tabulated"),
        "gj-876:S2": (21, "15.25 / 21.62", "57.75 / 66.78", "56.71 / 167.99", "searched", "searched-null"),
        "lacaille-8760:S2": (19, "10.77 / 29.33", "13.61 / 29.71", "5.49 / 61.34", "searched", "searched-null"),
        "ross-128:S1": (24, "24.30 / 16.69 exc", "13.35 / 25.75", "22.29 / 102.05", "vetoed_known_source", "static content: z +4.9 median, 96 % of events positive; curvature leaves +4.0; Gaia 6.0 vs template 4.1 (a G 6.8 star 0.8 px from centre); S_stack depth degenerate"),
        "ross-128:S2": (11, "12.03 / 16.25", "6.12 / 15.55", "3.07 / 10.30", "searched", "searched-null (LASCO/HI-1 family)"),
        "teegarden:S1": (23, "9.41 / 19.66", "7.18 / 26.52", "11.67 / 42.66", "searched", "searched-null (LASCO/HI-1 family)"),
        "teegarden:S2": (18, "7.20 / 24.28", "12.78 / 34.67", "39.16 / 44.71", "searched", "searched-null (LASCO/HI-1 family)"),
        "van-maanen:S1": (18, "17.20 / 14.92 exc", "8.04 / 10.58", "7.34 / 5.69 exc", "searched", "S_stack exceedance = latitude curvature (source at orbit latitude −3° on the brightness ridge; curvature-corrected S_stack 5.1 < T), S_stack depth degenerate; S_pulse exceedance (E17) = whole-arc plateau with controls co-elevated 1–4σ, frame-common coronal transient; limit here is the S_event depth"),
        "van-maanen:S2": (22, "−7.21 / 24.30", "12.38 / 39.27", "60.15 / 103.73", "searched", "searched-null (LASCO/HI-1 family)"),
        "wolf-1061:S1": (18, "13.18 / 15.13", "5.50 / 21.65", "6.77 / 5.74 exc", "searched", "S_pulse exceedance (E22, 20 min, 6–7σ) = extended coronal front at the inner edge: 11 × 11 stamp uniformly elevated 0.5–1 unit/px with no point-source peak, within 3° of the inner edge; non-recurrent"),
        "wolf-1061:S2": (23, "7.93 / 9.78", "5.95 / 17.75", "9.27 / 51.21", "searched", "searched-null"),
        "wolf-359:S1": (22, "5.09 / 15.35", "6.20 / 32.73", "37.96 / 108.54", "searched", "searched-null (LASCO/HI-1 family)"),
        "wolf-437:S1": (26, "4.24 / 25.50", "5.56 / 43.11", "63.93 / 123.82", "searched", "searched-null"),
        "wolf-437:S2": (12, "5.04 / 16.38", "4.82 / 8.74", "3.08 / 5.82", "searched", "searched-null"),
    }
    for key, (nev, ss, se, sp, status, disp) in conf.items():
        t, ch = key.split(":")
        u = pl.get(key, {})
        pk = "P_tx_10m_W" if ch == "S2" else "P_cone_W"
        lv = pw = None
        extra = []
        if "S_stack" in u:
            lv = u["S_stack"]["m90_Veq"]; pw = mw(u["S_stack"][pk])
        elif "S_event" in u:
            # no recurrence depth; carry the S_event depth as the row limit and say so
            lv = u["S_event"]["m90_Veq"]; pw = mw(u["S_event"][pk])
            extra.append("limit = S_event m90 (S_stack degenerate)")
        if "S_stack" in u and "S_event" in u:
            extra.append(f"S_event m90 {u['S_event']['m90_Veq']} -> {mw(u['S_event'][pk])} MW")
        if "S_pulse" in u:
            extra.append(f"S_pulse m90 {u['S_pulse']['m90_Veq']} -> {mw(u['S_pulse'][pk])} MW (2 consecutive frames, 5–15 min)")
        if ch == "S2" and "S_stack" in u:
            extra.append(f"P_cone {mw(u['S_stack']['P_cone_W'])} MW")
        pnote = "power_mw = P_tx_10m (10-m-class uplink)" if ch == "S2" else "power_mw = P_cone downlink through the 0.1 AU cone"
        rows.append(row(
            target_id=t, channel=ch, rung="0.1 AU", band=band, substrate=sub,
            status=status, limit_kind="m90" if lv is not None else None, limit_value=lv,
            limit_unit="V_eq mag" if lv is not None else None, power_mw=pw, n_events=nev, n_trials=3,
            source=f"{RES}/confirmatory_v1.md unit + completeness + adjudication tables; {RES}/power_limits_v1.json {key}; {R} §3–4",
            notes=f"S_stack {ss}, S_event {se}, S_pulse {sp} (S/T); {disp}; {pnote}; {'; '.join(extra) if extra else 'no depth tabulated'}; ±0.3 mag band and ±0.3 mag flux-scale systematics declared",
        ))
    # five not-covered S1 units
    for t in ("procyon-a", "procyon-b", "luyten-star", "gj-54", "gj-1087"):
        rows.append(row(
            target_id=t, channel="S1", rung="0.1 AU", band=band, substrate=sub,
            status="not_constrainable", n_events=0,
            source=f"{R} §3; {RES}/confirmatory_v1.md unit table (status constraint_only)",
            notes="Not covered: source at orbit latitude +17.5° to +19.6°, above the real WISPR-I field's +15° north edge (planning model's ±20° was optimistic); results table labels the unit constraint_only",
        ))
    rows.append(row(
        target_id="fomalhaut", channel="S2", rung="0.1 AU", band=band, substrate=sub,
        status="structurally_open", n_events=0,
        source=f"{R} §3; {RES}/confirmatory_v1.md unit table",
        notes="bright_star_unsearchable (V 1.2): saturated target star",
    ))
    # dev units (results/dev_v1.md)
    devs = {
        ("61-vir", "S1"): (23, "−12.75 / 25.28", "4.10 / 32.84", "9.78 / 56.66", "constraint_only", "Dev-stage unit (no injection depth)"),
        ("61-vir", "S2"): (22, "18.70 / 20.05", "31.59 / 7.00 exc", "3.02 / 32.57", "retained_ambiguous", "Dev-stage unit: S_event exceedance = monotonic 0.95 → 0.88 calibration drift of the V 4.7 star; retained, non-promotable"),
        ("gj-908", "S1"): (11, "−7.63 / 42.73", "1.73 / 30.56", "3.22 / 66.87", "constraint_only", "Dev-stage unit (no injection depth)"),
        ("gj-908", "S2"): (24, "4.65 / 18.62", "9.49 / 19.29", "21.01 / 66.58", "constraint_only", "Dev-stage unit (no injection depth): gj-908 (V 9.0) template null, z median +0.3 over 24 events — in-situ control"),
        ("ross-154", "S2"): (17, "7.70 / 12.23", "9.15 / 17.01", "5.94 / 69.45", "constraint_only", "Dev-stage unit (no injection depth)"),
    }
    for (t, ch), (nev, ss, se, sp, status, note) in devs.items():
        rows.append(row(
            target_id=t, channel=ch, rung="0.1 AU", band=band, substrate=sub,
            status=status, n_events=nev, n_trials=3,
            source=f"{R} §3 (dev); {RES}/dev_v1.md unit table",
            notes=f"{note}; dev S_stack {ss}, S_event {se}, S_pulse {sp} (S/T); dev stage 5 units, 15 trials, 1 exceedance vs 1.67 expected",
        ))
    for ch in ("S1", "S2"):
        rows.append(row(
            target_id="programme", channel=ch, rung="0.1 AU", band="WISPR-O (ε 50°–108°)", substrate="PSP/WISPR-O frames (ledger substrate)",
            status="ledger_only", source=f"{R} §6 item 4",
            notes="WISPR-O carried as a ledger substrate; its arcs cover the outer cone only; not searched",
        ))
    rows.append(row(
        target_id="programme", channel="S2", rung="1.0 AU", band=band, substrate=sub,
        status="no_survey", source=f"{R} §1",
        notes="S2 1 AU rung out of scope (D7); no coverage count given in the report",
    ))
    doc = dict(
        survey_id="wispr_crossings", report=R, plan_section="5.18–5.19", pipeline="B", archives=["wispr"],
        hypothesis_version="v1.0 (D0–D10 approved 2026-09-05) + threshold amendments v1.1 (pre-dev), v1.2–v1.5 (dev-driven, pre-confirmatory)",
        era="2018-11 → 2026-03 (encounters E1–E27)", run_dir="runs/wispr-crossings", records_dir=None,
        candidates=0,
        candidate_notes="90 trials / 30 units, 15 exceedances vs 10.0 expected, all adjudicated: 9 static-content / latitude-curvature (Gaia-matched), 3 bright-star class, 2 arc-exit ramps, 1 extended coronal front; dev 61-vir S2 calibration drift retained non-promotable",
        rows=rows,
    )
    dump(doc, "wispr_crossings")


# ---------------------------------------------------------------------------
# 4. Solar Orbiter SoloHI
# ---------------------------------------------------------------------------
def solohi():
    R = "report/solohi_crossings.md"
    RES = "surveys/solohi-crossings/results"
    pl = json.load(open(f"{ROOT}/{RES}/power_limits_v1.json"))
    sub = "Solar Orbiter SoloHI level-2 inner-tile frames (P04–P12), ICRS-fixed patch aperture photometry with star-fixed off-beam baseline"
    band = "SoloHI λ_c 610 nm, W_eff 250 nm (532 nm in band)"
    rows = []
    for rung in ("1.2 Rsun", "2.5 Rsun"):
        rows.append(row(
            target_id="programme", channel="S1", rung=rung, band=band, substrate=sub,
            status="not_constrainable", source=f"{R} §1",
            notes="The grazing cones never approach the 5.2° inner edge of the SoloHI tiles (≤ 1.1° / 2.4°)",
        ))
    conf = {
        "ez-aqr:S1": (5, "2.52 / 4.13", "4.92 / 12.48", "6.99 / 14.80", "searched", "searched-null"),
        "gj-1002:S1": (3, "−0.63 / 6.56", "3.99 / 15.37", "2.78 / 10.11", "searched", "searched-null"),
        "gj-1111:S2": (6, "19.16 / 13.33 exc", "16.93 / 12.81 exc", "38.43 / 30.16 exc", "searched", "all three exceedances = control-interpolation artifact: a bright star in the +1.5° control patch loses its core to the saturation mask on single frames (E′ −47 to −54 units) and the frozen quadratic, extrapolated to the source offset on a one-sided ladder, turns that into +50–60 units on the source; source pixels flat (−3.6 … +0.4 units); diagnostic median interpolation S_stack −3.5 vs 5.6, S_pulse 2.4 vs 5.4; all three statistics degenerate, no depth"),
        "gj-1276:S1": (3, "−7.59 / 14.21", "−3.23 / 19.59", "2.17 / 7.22", "searched", "searched-null (LASCO/HI-1/WISPR family)"),
        "gj-251:S2": (4, "−1.02 / 9.79", "0.84 / 8.56", "5.48 / 5.42 exc", "searched", "S_pulse exceedance (marginal) = control-interpolation artifact, same mechanism on a one-sided ladder at +6.5°; diagnostic median 2.38 / 3.12"),
        "gj-581:S1": (4, "3.44 / 13.91", "3.48 / 12.14", "1.47 / 8.43", "searched", "searched-null"),
        "gj-588:S1": (3, "2.75 / 5.74", "4.30 / 6.41", "2.63 / 10.70", "searched", "searched-null"),
        "gj-667-c:S1": (4, "−1.43 / 10.67", "1.45 / 13.46", "3.56 / 4.39", "searched", "searched-null"),
        "gj-674:S1": (3, "4.71 / 5.31", "3.10 / 5.00", "3.14 / 3.74", "searched", "searched-null; second-deepest unit"),
        "gj-682:S1": (5, "−0.83 / 7.85", "3.56 / 6.49", "2.50 / 7.64", "searched", "searched-null"),
        "gj-783:S1": (7, "2.89 / 3.47", "6.67 / 8.80", "6.79 / 10.13", "searched", "searched-null; deepest unit"),
        "gj-876:S1": (5, "4.60 / 3.71 exc", "7.87 / 4.85 exc", "2.10 / 17.96", "searched", "S_stack/S_event exceedances = extended level offsets of alternating sign on a starless antipode patch: per-orbit z +5.9 (P09), −5.4 (P10), +7.9 (P11); aligned mean stacks uniform ±0.01–0.04 units/px over 15 × 15 px, no point-source peak; marginal against 8 controls spanning −4.2 … +3.7; no consistent-sign recurrence; S_stack/S_event degenerate"),
        "ross-128:S2": (3, "0.82 / 8.36", "4.12 / 15.33", "2.47 / 17.75", "searched", "searched-null (LASCO/HI-1/WISPR family)"),
        "teegarden:S2": (8, "1.78 / 4.37", "7.80 / 6.77 exc", "15.41 / 6.37 exc", "searched", "S_event (P12, 16-epoch arc) and S_pulse exceedances = uncatalogued moving objects through the fixed patch: compact peaks displaced 4 px between consecutive frames, absent before/after (2026-02-23 20:02 → 20:26; 2025-03-24 05:13 → 05:37, diffuse); no H < 7 asteroid within 1°; star's baseline 0.14–0.19 units; S_event/S_pulse degenerate"),
        "van-maanen:S1": (3, "−1.92 / 7.76", "3.86 / 12.07", "4.58 / 5.90", "searched", "searched-null (family)"),
        "van-maanen:S2": (3, "0.75 / 8.89", "1.32 / 9.20", "4.26 / 33.89", "searched", "searched-null (family)"),
        "wolf-1061:S1": (3, "−2.66 / 6.80", "0.81 / 12.73", "1.74 / 8.26", "searched", "searched-null"),
        "wolf-359:S2": (3, "−4.34 / 1.12", "−0.48 / 3.74", "2.44 / 8.22", "searched", "searched-null (family)"),
    }
    for key, (nev, ss, se, sp, status, disp) in conf.items():
        t, ch = key.split(":")
        u = pl.get(key, {})
        pk = "P_tx_10m_W" if ch == "S2" else "P_cone_W"
        lv = pw = None
        extra = []
        if "S_stack" in u:
            lv = u["S_stack"]["m90_Veq"]; pw = mw(u["S_stack"][pk])
        if "S_event" in u:
            extra.append(f"S_event m90 {u['S_event']['m90_Veq']} -> {mw(u['S_event'][pk])} MW")
        if "S_pulse" in u:
            extra.append(f"S_pulse m90 {u['S_pulse']['m90_Veq']} -> {mw(u['S_pulse'][pk])} MW (2 consecutive frames, 12–48 min)")
        if ch == "S2" and "S_stack" in u:
            extra.append(f"P_cone {mw(u['S_stack']['P_cone_W'])} MW")
        pnote = "power_mw = P_tx_10m (10-m-class uplink)" if ch == "S2" else "power_mw = P_cone downlink through the 0.1 AU cone"
        rows.append(row(
            target_id=t, channel=ch, rung="0.1 AU", band=band, substrate=sub,
            status=status, limit_kind="m90" if lv is not None else None, limit_value=lv,
            limit_unit="V_eq mag" if lv is not None else None, power_mw=pw, n_events=nev, n_trials=3,
            source=f"{RES}/confirmatory_v1.md unit + completeness + adjudication tables; {RES}/power_limits_v1.json {key}; {R} §3–4",
            notes=f"S_stack {ss}, S_event {se}, S_pulse {sp} (S/T); {disp}; {pnote}; {'; '.join(extra) if extra else 'no depth tabulated'}; n_events = included events over 3–8 orbits; ±0.3 mag band and ±0.3 mag flux-scale systematics declared",
        ))
    # 9 constraint_only units
    co = {
        ("ez-aqr", "S2"): (1, "included, insufficient_after_veto, insufficient_epochs, no_baseline"),
        ("gj-1002", "S2"): (2, "included, insufficient_after_veto, insufficient_epochs"),
        ("gj-1111", "S1"): (0, "insufficient_epochs"),
        ("gj-1276", "S2"): (2, "included, insufficient_epochs"),
        ("gj-876", "S2"): (2, "included, insufficient_after_veto, insufficient_epochs"),
        ("ross-128", "S1"): (2, "included, no_baseline"),
        ("ross-154", "S2"): (0, "insufficient_epochs"),
        ("teegarden", "S1"): (2, "included, insufficient_epochs"),
        ("wolf-359", "S1"): (1, "included, insufficient_epochs, no_baseline"),
    }
    for (t, ch), (nev, flags) in co.items():
        seam = "; report §6: the seam rows (ross-128 P05–P09, in-plane) are ledger entries" if t == "ross-128" else ""
        rows.append(row(
            target_id=t, channel=ch, rung="0.1 AU", band=band, substrate=sub,
            status="constraint_only", n_events=nev,
            source=f"{R} §3, §6; {RES}/confirmatory_v1.md unit table",
            notes=f"Report state constraint_only: < 3 events survive the edge gates and the baseline rule ({flags}) — arcs lie mostly inside the sunward-edge band (ε ≲ 8°) or lack the same-tile baseline; report calls these ledger entries needing the 16-s perihelion geometry to become searchable{seam}",
        ))
    # dev units
    devs = {
        ("61-vir", "S1"): (4, "−4.63 / 6.59", "2.48 / 8.88", "2.62 / 14.46", "constraint_only", "Dev-stage unit (no injection depth)"),
        ("gj-908", "S1"): (5, "−0.11 / 3.15", "3.21 / 15.58", "3.92 / 22.17", "constraint_only", "Dev-stage unit (no injection depth); gj-908 S1 z median 0.00"),
        ("gj-908", "S2"): (4, "−2.38 / 11.77", "1.77 / 13.70", "3.61 / 28.19", "constraint_only", "Dev-stage unit (no injection depth); gj-908 S2 (V 9.0) baseline null, z median −1.2 over 4 orbits — in-situ control"),
        ("ross-154", "S1"): (6, "10.59 / 16.75", "9.56 / 18.84", "2.61 / 20.65", "constraint_only", "Dev-stage unit (no injection depth)"),
    }
    for (t, ch), (nev, ss, se, sp, status, note) in devs.items():
        rows.append(row(
            target_id=t, channel=ch, rung="0.1 AU", band=band, substrate=sub,
            status=status, n_events=nev, n_trials=3,
            source=f"{R} §3 (dev); {RES}/dev_v1.md unit table",
            notes=f"{note}; dev S_stack {ss}, S_event {se}, S_pulse {sp} (S/T); dev stage 5 units, 12 trials, 0 exceedances vs 1.33 expected",
        ))
    rows.append(row(
        target_id="61-vir", channel="S2", rung="0.1 AU", band=band, substrate=sub,
        status="structurally_open", n_events=1,
        source=f"{R} §3 (dev); {RES}/dev_v1.md unit table",
        notes="bright_star_saturated: the V 4.7 star's aperture fails the mask/validity gates in 5 of 6 events (dev)",
    ))
    for ch in ("S1", "S2"):
        rows.append(row(
            target_id="programme", channel=ch, rung="0.1 AU", band=band, substrate="SoloHI inner-tile sunward-edge band (ε ≲ 8°, inner ~120 px of every arc)",
            status="ledger_only", n_events=49388, source=f"{R} §2, §5, §6 item 2",
            notes="F-corona / sunward-edge strip: saturation plateau in the ≥ 45-s regime and unresolved corona structure (6–11 units/frame) make the frames unusable for point sources; the structure-noise gate dropped 49,388 of 193,902 confirmatory arc records (25 %) — the inner quarter of every arc; a K/F-corona structure model or 16-s perihelion frames would be needed to recover it",
        ))
    rows.append(row(
        target_id="programme", channel="S2", rung="1.0 AU", band=band, substrate="SoloHI outer tiles (ε 25°–45°), used only as baseline",
        status="ledger_only", n_events=214, source=f"{R} §1, §6 item 5",
        notes="S2 1 AU transits in the outer tiles: 36 stars, 214 covered star-orbit pairs; out of scope (D7) — coverage logged, no statistic",
    ))
    doc = dict(
        survey_id="solohi_crossings", report=R, plan_section="5.23", pipeline="B", archives=["solohi"],
        hypothesis_version="v1.0 (D0–D9 approved 2026-09-06) + threshold amendments v1.1 (saturation mask), v1.2 (structure-noise gate, colour calibrators V ≥ 4.5)",
        era="2022-01 → 2026-04 (orbits P04–P12)", run_dir="runs/solohi-crossings", records_dir=None,
        candidates=0,
        candidate_notes="54 trials / 18 units, 8 exceedances vs 6.0 expected, all adjudicated: 4 control-interpolation artifacts (gj-1111 S2 ×3, gj-251 S2), 2 uncatalogued moving objects (teegarden S2), 2 extended level offsets of alternating sign (gj-876 S1); none retained",
        rows=rows,
    )
    dump(doc, "solohi_crossings")


# ---------------------------------------------------------------------------
# 5. S2 1 AU outer skirt (ATLAS)
# ---------------------------------------------------------------------------
def skirt():
    R = "report/skirt_crossings.md"
    RES = "surveys/skirt-crossings/results"
    comp = {(u["target_id"], u["rung_au"]): u for u in json.load(open(f"{ROOT}/{RES}/completeness_v1.json"))["units"]}
    search = {}
    for fn in ("confirmatory_search_v1.json", "dev_search_v1.json"):
        d = json.load(open(f"{ROOT}/{RES}/{fn}"))
        for u in d["units"]:
            search[(u["target_id"], u["rung_au"])] = (u, d["split"])
    adj = {}
    for fn in ("confirmatory_adjudication_v1.json", "dev_adjudication_v1.json"):
        for e in json.load(open(f"{ROOT}/{RES}/{fn}"))["exceedances"]:
            adj[(e["target_id"], e["rung_au"])] = e
    ledger = json.load(open(f"{ROOT}/{RES}/ledger_1au_v1.json"))
    sub = "ATLAS reduced-mode forced photometry (o band), elongation-locked step vs 8-star Gaia DR3 null ensemble"
    rung_str = {"0.9": "0.90 AU", "0.95": "0.95 AU"}
    o_med = {t["target_id"]: t["o_median"] for t in ledger}
    rows = []
    # per-unit rows (28 searched units), in results-file order
    for (t, ru), (u, split) in search.items():
        c = comp[(t, ru)]["statistics"]
        st = u["statistics"]
        stat_bits = []
        for k in ("S_sym", "S_skirt", "S_year"):
            if k in st:
                v = st[k]
                stat_bits.append(f"{k} {v['S']:.2f}/{v['T']:.2f}{' exc' if v['exceedance'] else ''}")
        depth_bits = []
        for k in ("S_sym", "S_skirt", "S_year"):
            if k in c:
                x = c[k]
                if x["constraining"]:
                    depth_bits.append(f"{k} m90 {x['m90']} -> {mw(x['power_W_m90'])} MW")
                else:
                    depth_bits.append(f"{k} no depth (exceedance trial)")
        a = adj.get((t, ru))
        disp = ""
        if a:
            disp = "; dispositions: " + ", ".join(f"{k} {v}" for k, v in a["disposition"].items())
            col = a["diagnostics"].get("colour", {})
            if col.get("unmatched"):
                disp += f" (no control within 1.0 mag of BP−RP; nearest gap {col['min_control_gap_mag']:.2f} mag)"
        n_tr = len(st)
        sym = c.get("S_sym")
        constraining_sym = sym is not None and sym["constraining"]
        retained = a is not None and any(v.startswith("retained_ambiguous") for v in a["disposition"].values())
        target_note = ""
        if t == "gj-13157":
            target_note = "; crowded (§4.2): 14″ from a G 9.4 star, 11″ from a G 12.1 star, tphot chi/N 35,800 vs 13 for controls — depths reported systematics-limited"
        if t == "eps-ind-b":
            target_note = "; one of the three faintest targets (0.6–1.5 MW class)"
        common = dict(
            target_id=t, channel="S2", rung=rung_str[ru], band="o", substrate=sub,
            n_events=u["n_nights_target"], n_trials=n_tr,
            source=f"{RES}/{split}_search_v1.json; {RES}/completeness_v1.json; {RES}/{split}_adjudication_v1.json; {R} §3, §5",
        )
        if retained and not constraining_sym:
            # teegarden: every S_sym trial retained under A4
            rows.append(row(**common, status="retained_ambiguous",
                            notes=f"{split} unit, o median {o_med[t]:.2f}; {', '.join(stat_bits)} (S/T){disp}; {'; '.join(depth_bits)}; teegarden is an M7 flare star 2.7 mag redder than its nearest control: retained_ambiguous_colour_unmatched (A4), not promotable; n_events = nightly units; per-trial depth from completeness_v1.json (power = line power through the beam, 690 nm / 260 nm effective width){target_note}"))
            continue
        if constraining_sym:
            rows.append(row(**common, status="searched", limit_kind="m90", limit_value=sym["m90"],
                            limit_unit="AB mag (ATLAS o)", power_mw=mw(sym["power_W_m90"]),
                            notes=f"{split} unit, o median {o_med[t]:.2f}; {', '.join(stat_bits)} (S/T){disp}; limit = S_sym (symmetric recurrent step) m90 -> line power through the {rung_str[ru]} beam (690 nm / 260 nm effective width, ±0.1 mag scale, ±0.2 % parallax residual); {'; '.join(depth_bits)}; n_events = nightly units; 8 valid controls{target_note}"))
        else:
            rows.append(row(**common, status="searched",
                            notes=f"{split} unit, o median {o_med[t]:.2f}; {', '.join(stat_bits)} (S/T){disp}; S_sym is an exceedance trial — no S_sym depth claimed (hypotheses §12 A5); {'; '.join(depth_bits)}; n_events = nightly units; 8 valid controls{target_note}"))
        if retained:
            # gj-3512 0.95: S_skirt / S_year retained under A4 alongside a constraining S_sym
            rows.append(row(target_id=t, channel="S2", rung=rung_str[ru], band="o", substrate=sub,
                            status="retained_ambiguous", n_trials=sum(1 for v in a["disposition"].values() if v.startswith("retained_ambiguous")),
                            source=f"{RES}/{split}_adjudication_v1.json; {R} §3 table",
                            notes=f"{', '.join(f'{k} {v}' for k, v in a['disposition'].items())}; S_skirt 0.85/−0.12, S_year 2.77/2.14, Δ skirt/opp +0.3 / +0.3 %; colour-unmatched (gap 1.8 mag), not promotable (A4)"))
    # ledger-only units (report §3 Ledger): |β| >= ε_b
    for t, ru, why in (("gj-293", "0.9", "|β| ≥ ε_b (ecliptic latitude −79.1°)"), ("gj-293", "0.95", "|β| ≥ ε_b (ecliptic latitude −79.1°)"),
                       ("gj-13157", "0.9", "|β| ≥ ε_b at the 0.90 AU rung (ecliptic latitude 63.7°)"), ("gj-3112", "0.9", "|β| ≥ ε_b at the 0.90 AU rung (ecliptic latitude −65.9°)"),
                       ("wolf-1069", "0.9", "|β| ≥ ε_b (ecliptic latitude 71.7°); also saturation-excluded (o median 12.38)"), ("wolf-1069", "0.95", "|β| ≥ ε_b (ecliptic latitude 71.7°); also saturation-excluded (o median 12.38)")):
        rows.append(row(target_id=t, channel="S2", rung=rung_str[ru], band="o", substrate=sub,
                        status="ledger_only", source=f"{R} §3 Ledger; {RES}/ledger_1au_v1.json (ecl_lat_deg)",
                        notes=f"Ledger-only unit: {why} — the star never leaves the beam, no in/out step to search"))
    # saturation-excluded target gj-1111
    for ru in ("0.9", "0.95"):
        rows.append(row(target_id="gj-1111", channel="S2", rung=rung_str[ru], band="o", substrate=sub,
                        status="structurally_open", source=f"{R} §2 population; {RES}/ledger_1au_v1.json",
                        notes="Saturation-excluded after the measured cut (o median 12.44); 1 AU coverage still logged"))
    for ru in ("0.9", "0.95"):
        rows.append(row(target_id=["luhman16-a", "luhman16-b"], channel="S2", rung=rung_str[ru], band="o", substrate=sub,
                        status="no_survey", source=f"{R} §2 population",
                        notes="luhman16 excluded at recon (no verified Gaia counterpart)"))
    # 1 AU ledger per target
    for tl in ledger:
        rows.append(row(target_id=tl["target_id"], channel="S2", rung="1.0 AU", band="o", substrate="ATLAS o-band nightly coverage (1 AU coverage ledger)",
                        status="ledger_only", n_events=tl["s2_nights_ge_floor_total"],
                        source=f"{RES}/ledger_1au_v1.json; {R} §5 (1 AU ledger)",
                        notes=f"Coverage record, not a constraint (a 1 AU beam contains Earth all year; no temporal signature): {tl['n_cycles']} conjunction cycles ({tl['n_cycles_s2_covered']} S2-covered), S2-side nights at ε ≥ 50° floor {tl['s2_nights_ge_floor_total']} total / median {tl['s2_nights_median_per_cycle']:g} per cycle, A-side nights {tl['a_side_nights_total']}, minimum observed ε {tl['min_eps_observed']:.1f}°, o median {tl['o_median']:.2f}, saturation screen {tl['saturation']}"))
    doc = dict(
        survey_id="skirt_crossings", report=R, plan_section="5.25", pipeline="B", archives=["atlas"],
        hypothesis_version="v1.0 (D1–D8 approved 2026-09-07) + amendments A1–A4 (pre-confirmatory) and A5 (completeness reporting rule)",
        era="ATLAS full history, 5–12 conjunction cycles per target (hypotheses §3: MJD 57227 → probe end MJD 61252)",
        run_dir="runs/skirt-crossings", records_dir=None,
        candidates=0,
        candidate_notes="81 trials / 28 units, 24 exceedances vs 9.0 expected (dev 7 vs 1.3; confirmatory 17 vs 7.7): all non-promotable (asymmetric / single-cycle / ensemble-shared) or retained_ambiguous_colour_unmatched under A4 (teegarden both rungs, gj-3512 0.95 AU)",
        rows=rows,
    )
    dump(doc, "skirt_crossings")


if __name__ == "__main__":
    lasco(); stereo(); wispr(); solohi(); skirt()
