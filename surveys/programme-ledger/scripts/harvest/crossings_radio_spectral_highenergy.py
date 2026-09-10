#!/usr/bin/env python
"""Emit harvest YAML for the five Pipeline-B reports (radio ×3, spectral, high-energy).

Every number is transcribed from the report or the results file named in `source`.
"""
import json, re, sys
import yaml
from pathlib import Path

ROOT = Path('/home/jgreene/code/sgl-seti-survey')
import sys
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / 'surveys/programme-ledger/harvest'

FIELDS = ['target_id', 'channel', 'rung', 'band', 'substrate', 'status', 'limit_kind',
          'limit_value', 'limit_unit', 'limit_range', 'power_mw', 'power_mw_range',
          'n_events', 'n_trials', 'epoch_range', 'source', 'notes']


def row(**kw):
    unknown = set(kw) - set(FIELDS)
    assert not unknown, unknown
    return {f: kw.get(f) for f in FIELDS}


def dump(name, doc):
    path = OUT / f'{name}.yaml'
    text = yaml.safe_dump(doc, sort_keys=False, allow_unicode=True, width=110, default_flow_style=False)
    path.write_text(text)
    yaml.safe_load(path.read_text())
    print(f'{path}: {len(doc["rows"])} rows, candidates={doc["candidates"]}')


# ---------------------------------------------------------------------------
# 1. report/radio_crossings.md  (plan §5.5)
# ---------------------------------------------------------------------------
R = 'report/radio_crossings.md'
BL = 'Breakthrough Listen Open Data archive pointing metadata (GBT/Parkes/MeerKAT), 5-min bins, in-beam test'
VL = 'VLASS quick-look tile epoch dates (NRAO dynamic-schedule summary), ±3 d tile tolerance'
rows = []
rows.append(row(target_id='programme', channel='A', rung='0.1 AU', band='any', substrate=BL,
                status='structurally_open',
                source=f'{R} §2 table; §3 item 1',
                notes='Zero GBT/Parkes/MeerKAT pointings inside any A 0.1 AU window (93 cones, 0 query failures); '
                      'targeted narrow-rung radio re-analysis moot. Hand-off (c): re-run the intersection if BL '
                      'MeerKAT holdings reach the public archive. Scope: BL public open-data archive as snapshotted 2026-08-24.'))
rows.append(row(target_id='teegarden', channel='A', rung='0.1 AU', band='optical (APF spectra)',
                substrate='BL APF optical spectra of teegarden (as SO0253) inside an A 0.1 AU window',
                status='ledger_only', n_events=1, epoch_range='2016-11-13',
                source=f'{R} §2 table; §3 item 4(b); surveys/radio-crossings/results/bl_inwindow_v1.ecsv',
                notes='3 APF spectra 5.3 d inside an 11.5 d window (non-radio). Hand-off (b): optical laser-line '
                      'inspection → spectral-archive family (plan §3.6 row 10); report/spectral_archives.md §6 lists '
                      'the raw 2D frames as a pending hand-off.'))
for rung in ['1.2 Rsun', '2.5 Rsun', '0.1 AU']:
    rows.append(row(target_id='programme', channel='B', rung=rung, band='any', substrate=BL,
                    status='structurally_open',
                    source=f'{R} §2 table; §3 item 3; §5',
                    notes='Zero BL pointings within 0.5° of any antipode in-window ("nobody has ever looked"); '
                          'future-observation recommendation, not an archival opportunity. Qualifier (§5, 2026-09-04): '
                          'holds for targeted radio (BL) and the grazing rungs; wide-field ASKAP epochs exist at the '
                          '0.1 AU antipode rung (report/radio_crossings_ext.md).'))
rows.append(row(target_id='programme', channel='A', rung='1.0 AU', band='any', substrate=BL,
                status='ledger_only', n_events=22,
                source=f'{R} §2 table; §3 item 2; surveys/radio-crossings/results/coverage_v1_summary.json',
                notes='143 in-beam observations, 22 events, 18 targets (incl. proxima-cen, barnard-era targets, '
                      'teegarden); windows tile ~1/3–2/3 of the calendar; BL published blind nulls cited as the '
                      'constraint for these cells (their own thresholds / drift-range caveats). 221 rows / 39 events / '
                      '28 targets within the 0.5° cone.'))
vl_src = f'{R} §2 table; §3 item 4(a); surveys/radio-crossings/results/vlass_inwindow_v1.ecsv'
rows.append(row(target_id='gj-908', channel='A', rung='0.1 AU', band='any', substrate=VL, status='ledger_only',
                n_events=2, epoch_range='2017-09 (VLASS1.1), 2025-09 (VLASS4.1)', source=vl_src,
                notes='Two strict in-window tile epochs (offsets +2.5 d, +2.4 d). Hand-off (a) executed 2026-09-04 '
                      'as the six-cutout quick-look (report/radio_quicklook.md): nothing above 3σ; broadband only '
                      '(continuum imaging dilutes CW lines by ~1e9).'))
rows.append(row(target_id='ross-154', channel='A', rung='0.1 AU', band='any', substrate=VL, status='ledger_only',
                n_events=1, epoch_range='2019-07 (VLASS1.2)', source=vl_src,
                notes='Strict in-window tile epoch (offset −1.9 d). Hand-off (a) executed: report/radio_quicklook.md.'))
rows.append(row(target_id='teegarden', channel='A', rung='0.1 AU', band='any', substrate=VL, status='ledger_only',
                n_events=1, epoch_range='2021-11 (VLASS2.2)', source=vl_src,
                notes='Strict in-window tile epoch (offset −2.3 d). Hand-off (a) executed: report/radio_quicklook.md.'))
rows.append(row(target_id='gj-1276', channel='A', rung='0.1 AU', band='any', substrate=VL, status='structurally_open',
                n_events=1, epoch_range='2020-09 (VLASS2.1)', source=vl_src + '; §3 item 4(a)',
                notes='Tolerance-edge row only (offset +8.4 d vs ±3 d tile tolerance); exact CADC tile times '
                      '(report/radio_quicklook.md §3) put it outside the strict window — no in-window VLASS epoch.'))
rows.append(row(target_id='wolf-359', channel='B', rung='0.1 AU', band='any', substrate=VL, status='structurally_open',
                n_events=1, epoch_range='2020-09 (VLASS2.1)', source=vl_src + '; §3 item 4(a)',
                notes='Tolerance-edge row only (offset +8.0 d); exact CADC tile times put it outside the strict '
                      'window — no in-window VLASS antipode epoch.'))
for rung in ['1.2 Rsun', '2.5 Rsun']:
    rows.append(row(target_id='programme', channel='B', rung=rung, band='any', substrate=VL,
                    status='structurally_open', source=f'{R} §2 table',
                    notes='Zero VLASS tile epochs in any grazing-rung antipode window.'))
rows.append(row(target_id='programme', channel='A', rung='1.0 AU', band='any', substrate=VL, status='ledger_only',
                n_events=69, source='surveys/radio-crossings/results/coverage_v1_summary.json (vlass A_1.0)',
                notes='69 in-window tile epochs, 48 targets; not discussed in the report text (products file only).'))
dump('radio_crossings', {
    'survey_id': 'radio_crossings', 'report': R, 'plan_section': '5.5', 'pipeline': 'B',
    'archives': ['bl', 'vlass'], 'hypothesis_version': None,
    'era': 'BL: crossing events 2016–2028 (archive snapshot 2026-08-24); VLASS: four epochs 2018–2026',
    'run_dir': 'runs/radio-crossings', 'records_dir': None, 'candidates': 0,
    'candidate_notes': 'Metadata-only coverage intersection executed 2026-08-24; no signal statistic formed. '
                       'Decision (plan §11.2 step 7): geometry products stay geometry-only; no radio signal search in this repository.',
    'rows': rows})

# ---------------------------------------------------------------------------
# 2. report/radio_crossings_ext.md  (plan §5.7)
# ---------------------------------------------------------------------------
R = 'report/radio_crossings_ext.md'
AS = 'CASDA ObsCore Stokes-I continuum images (RACS/VAST/EMU/FLASH/…): SBID start/stop vs flat-chord window; in_footprint = sep ≤ 2.25° + HPBW/2'
LO = 'LoTSS DR3 HBA pointing observation mid-MJDs (ASTRON VO); in_beam ≤ 1.98°'
ECSV = 'surveys/radio-crossings/results/askap_inwindow_v2.ecsv'
rows = []
rows.append(row(target_id='programme', channel='B', rung='1.2 Rsun', band='any', substrate=AS, status='structurally_open',
                n_events=29, source=f'{R} §2 table; §3 item 3',
                notes='29 in-span events (4 targets), 0 covered. ±0.3–0.7 d windows vs 12–15-min visits at cadences of '
                      'weeks (VAST) to years (RACS); grazing rungs stay open in radio — scheduled observation only.'))
rows.append(row(target_id='programme', channel='B', rung='2.5 Rsun', band='any', substrate=AS, status='structurally_open',
                n_events=36, source=f'{R} §2 table; §3 item 3',
                notes='36 in-span events (5 targets), 0 covered (65 in-span ASKAP grazing events in all, none touched).'))
rows.append(row(target_id='programme', channel='B', rung='0.1 AU', band='any', substrate=AS, status='ledger_only',
                n_events=50, source=f'{R} §2 table; §2.1; §3 item 2',
                notes='50 in-span events (7 targets): 7 covered in footprint, 6 validated (5 targets) — the antipode '
                      'channel is no longer archivally virgin at 0.1 AU. Hand-off (a) antipode cutout quick-look '
                      'executed 2026-09-04 (report/radio_quicklook.md).'))
rows.append(row(target_id='wolf-359', channel='B', rung='0.1 AU', band='888 MHz', substrate=AS, status='ledger_only',
                n_events=2, epoch_range='2023-09-03, 2024-09-10', source=f'{R} §2.1 table; {ECSV}',
                notes='b 0.71 R☉. VAST SB52549 VAST_2257-06 (UNCERTAIN; offset −2.1 of ±5.9 d; sep 0.80° / 3.13°) and '
                      'VAST SB65727 VAST_2257-06 (GOOD; +5.6 / 5.9 d; 0.80° / 3.13°): the two clean primary-beam-core '
                      'cases, the September wolf-359 B family also covered by ZTF/PGIR. 44 dated ASKAP antipode visits over all dates.'))
rows.append(row(target_id='ross-154', channel='B', rung='0.1 AU', band='1656 MHz', substrate=AS, status='ledger_only',
                n_events=1, epoch_range='2021-12-29', source=f'{R} §2.1 table; {ECSV}',
                notes='b 3.30 R☉. RACS-high SB34957 RACS_0651+23 (GOOD; −3.4 / 5.6 d; sep 0.63° / 2.72°). 13 dated antipode visits over all dates.'))
rows.append(row(target_id='ross-128', channel='B', rung='0.1 AU', band='1368 MHz', substrate=AS, status='ledger_only',
                n_events=1, epoch_range='2021-09-21', source=f'{R} §2.1 table; {ECSV}',
                notes='b 1.85 R☉. VAST pilot SB32330 VAST_2338+00 (GOOD; +1.7 / 5.8 d; sep 2.45° / 2.82°); plus one '
                      'in-image-only row (second pilot field at 3.16°). 57 dated antipode visits over all dates.'))
rows.append(row(target_id='van-maanen', channel='B', rung='0.1 AU', band='888 MHz', substrate=AS, status='ledger_only',
                n_events=1, epoch_range='2022-03-29', source=f'{R} §2.1 table; {ECSV}',
                notes='b 0.25 R☉. RACS-low SB38682 RACS_1237-06 (GOOD; −4.8 / 5.8 d; sep 3.12° / 3.13°, footprint edge); '
                      'plus one in-image-only row (second RACS field at 3.32°). 79 dated antipode visits over all dates.'))
rows.append(row(target_id='van-maanen', channel='B', rung='0.1 AU', band='856 MHz', substrate=AS, status='structurally_open',
                n_events=1, epoch_range='2026-03-29', source=f'{R} §2.1 table; {ECSV}',
                notes='b 0.24 R☉. FLASH SB83234 FLASH_497 (−5.0 / 5.8 d; sep 3.15° / 3.16°, edge) is REJECTED — failed '
                      'the survey validation, not released; excluded from the validated count.'))
rows.append(row(target_id='teegarden', channel='B', rung='0.1 AU', band='856 MHz', substrate=AS, status='ledger_only',
                n_events=1, epoch_range='2026-05-03', source=f'{R} §2.1 table; {ECSV}',
                notes='b 0.98 R☉. FLASH SB84179 FLASH_389 (UNCERTAIN; −3.2 / 5.9 d; sep 2.18° / 3.16°). 21 dated antipode visits over all dates.'))
rows.append(row(target_id=['gj-908', 'gj-1276'], channel='B', rung='0.1 AU', band='any', substrate=AS,
                status='structurally_open', source=f'{R} §2.1; {ECSV}',
                notes='No in-footprint in-window ASKAP observation for these two of the seven in-span targets; gj-908 '
                      '2022 RACS-low in-image only at 3.92°. Antipode visits over all dates: gj-908 89, gj-1276 45.'))
rows.append(row(target_id='programme', channel='A', rung='0.1 AU', band='any', substrate=AS, status='ledger_only',
                n_events=50, source=f'{R} §2 table; §2.2',
                notes='50 in-span events (7 targets): 3 covered / 3 validated. Wide-field continuum epochs only — none '
                      'targeted or high-time-resolution. Hand-off (b) executed 2026-09-04 (report/radio_quicklook.md).'))
rows.append(row(target_id='gj-1276', channel='A', rung='0.1 AU', band='888 MHz', substrate=AS, status='ledger_only',
                n_events=1, epoch_range='2023-09-03', source=f'{R} §2.2 table; {ECSV}',
                notes='b 0.83 R☉. VAST SB52549 VAST_2257-06 (UNCERTAIN; −1.6 / 5.9 d; sep 1.10°).'))
rows.append(row(target_id='teegarden', channel='A', rung='0.1 AU', band='1368 MHz', substrate=AS, status='ledger_only',
                n_events=1, epoch_range='2024-11-11', source=f'{R} §2.2 table; {ECSV}',
                notes='b 1.08 R☉. RACS-mid SB67840 RACS_0248+18 (+3.6 / 5.8 d; sep 2.09°); two further RACS-mid fields in-image only at 3.7–4.2°.'))
rows.append(row(target_id='van-maanen', channel='A', rung='0.1 AU', band='888 MHz', substrate=AS, status='ledger_only',
                n_events=1, epoch_range='2024-10-07', source=f'{R} §2.2 table; {ECSV}',
                notes='b 0.33 R☉. VAST SB66449 VAST_0037+06 (+1.2 / 5.8 d; sep 3.11°, footprint edge).'))
rows.append(row(target_id=['wolf-359', 'ross-128', 'ross-154', 'gj-908'], channel='A', rung='0.1 AU', band='any',
                substrate=AS, status='structurally_open', source=f'{R} §2.2; {ECSV}',
                notes='No in-footprint in-window ASKAP observation; gj-908 2021-09 in two VAST-pilot fields in-image only at 3.6°.'))
rows.append(row(target_id='programme', channel='A', rung='1.0 AU', band='any', substrate=AS, status='ledger_only',
                n_events=606, source=f'{R} §2 table; §2.3',
                notes='606 in-span events (86 targets): 211 covered in footprint / 206 validated (64 targets); over all '
                      'in-era events VAST covers 124, RACS 96. Constraint = the surveys\' own published transient/variable searches; nothing added.'))
rows.append(row(target_id='programme', channel='B', rung='1.2 Rsun', band='144 MHz', substrate=LO, status='structurally_open',
                n_events=10, source=f'{R} §2 table', notes='10 in-span events (1 target), 0 covered.'))
rows.append(row(target_id='programme', channel='B', rung='2.5 Rsun', band='144 MHz', substrate=LO, status='structurally_open',
                n_events=20, source=f'{R} §2 table', notes='20 in-span events (2 targets), 0 covered (30 LoTSS grazing events in all, zero observations).'))
rows.append(row(target_id='programme', channel='B', rung='0.1 AU', band='144 MHz', substrate=LO, status='structurally_open',
                n_events=30, source=f'{R} §2 table; §2.1',
                notes='30 in-span events (3 targets), 0 covered. LoTSS reaches only gj-1276 B (18), ross-128 B (4), ross-154 B (1) antipode visits over all dates, none in-window.'))
rows.append(row(target_id='teegarden', channel='A', rung='0.1 AU', band='144 MHz', substrate=LO, status='ledger_only',
                n_events=1, epoch_range='2023-11-04 22:09', source=f'{R} §2.2 table; surveys/radio-crossings/results/lotss_inwindow_v2.ecsv',
                notes='b 1.09 R☉. LoTSS pointing P043+19 (−3.9 / 5.8 d; sep 1.89°, in beam); 8-h run. Hand-off (b) executed (report/radio_quicklook.md §4).'))
rows.append(row(target_id='programme', channel='A', rung='0.1 AU', band='144 MHz', substrate=LO, status='ledger_only',
                n_events=50, source=f'{R} §2 table',
                notes='50 in-span events (5 targets at δ > −10°), 1 covered (teegarden); the other four targets have no in-beam in-window LoTSS observation.'))
rows.append(row(target_id='programme', channel='A', rung='1.0 AU', band='144 MHz', substrate=LO, status='ledger_only',
                n_events=300, source=f'{R} §2 table; §2.3',
                notes='300 in-span events (30 targets): 13 covered (10 targets). Constraint = LoTSS published transient searches; nothing added.'))
dump('radio_crossings_ext', {
    'survey_id': 'radio_crossings_ext', 'report': R, 'plan_section': '5.7', 'pipeline': 'B',
    'archives': ['casda', 'racs', 'vast', 'emu', 'flash', 'lotss'], 'hypothesis_version': None,
    'era': 'crossing events 2014–2028; CASDA archive span 2019-04-21 → 2026-08-28; LoTSS DR3 observations 2014-05-23 → 2024-08-18',
    'run_dir': 'runs/radio-crossings/v2', 'records_dir': None, 'candidates': 0,
    'candidate_notes': 'Metadata-only coverage intersection executed 2026-09-04; no pixel touched, no statistic formed. '
                       'Radio decision (§5.5, geometry-only) unchanged. Hand-offs (c) VLASS refresh at the yearly crossing-list '
                       'refresh and (d) BL-MeerKAT re-check remain open.',
    'rows': rows})

# ---------------------------------------------------------------------------
# 3. report/radio_quicklook.md  (plan §5.8)
# ---------------------------------------------------------------------------
R = 'report/radio_quicklook.md'
QL = 'surveys/radio-crossings/results/quicklook_v1.json'
CUT = 'ASKAP continuum image cutout (3′) + per-epoch/release component catalogue cone (20″ coincidence)'
VCUT = 'VLASS quick-look tile cutout (3 GHz, 1″ pixels; CADC CAOM2 exact tile times)'
rows = []


def ql(target, ch, band, sub, status, rms, n_events, epoch, notes, src=f'{R} §2 table; ' + QL):
    return row(target_id=target, channel=ch, rung='0.1 AU', band=band, substrate=sub, status=status,
               limit_kind='flux' if rms is not None else None, limit_value=rms,
               limit_unit='mJy/beam rms at the position (cutout annulus MAD; noise floor, not a frozen limit)' if rms is not None else None,
               n_events=n_events, epoch_range=epoch, source=src, notes=notes)


rows.append(ql('wolf-359', 'B', '888 MHz', CUT, 'ledger_only', 0.22, 1, '2023-09-03',
               'b 0.71 R☉; VAST SB52549 VAST_2257-06, Δt −2.1 d; per-SBID epoch-resolved catalogue, 77 in 0.5°, nearest 288″ (1.4 mJy); '
               'cutout value +0.03 mJy (0.1σ), peak 0.29, beam 12.8″ — empty.'))
rows.append(ql('wolf-359', 'B', '888 MHz', CUT, 'ledger_only', 0.17, 1, '2024-09-10',
               'b 0.71 R☉; VAST SB65727 VAST_2257-06, Δt +5.6 d; per-SBID catalogue, 107, nearest 292″ (1.8 mJy); '
               'value −0.12 mJy (−0.7σ), peak 0.13, beam 13.1″ — empty.'))
rows.append(ql('ross-154', 'B', '1656 MHz', CUT, 'ledger_only', 0.16, 1, '2021-12-29',
               'b 3.30 R☉; RACS-high SB34957 RACS_0651+23, Δt −3.4 d; RACS-high release catalogue (single epoch), 77, nearest 183″ (8.1 mJy); '
               'value −0.06 (−0.4σ), peak 0.20, beam 12.3″ — empty.'))
rows.append(ql('ross-128', 'B', '1368 MHz', CUT, 'ledger_only', 0.26, 1, '2021-09-21',
               'b 1.85 R☉; VAST pilot SB32330 VAST_2338+00, Δt +1.7 d; Selavy file (v2 reprocessing; v1 no permission), 49, nearest 219″ (9.1 mJy); '
               'value −0.17 (−0.7σ), peak 0.22, beam 12.7″ — empty.'))
rows.append(ql('van-maanen', 'B', '888 MHz', CUT, 'ledger_only', 0.23, 1, '2022-03-29',
               'b 0.25 R☉; RACS-low SB38682 RACS_1237-06, Δt −4.8 d; RACS-low2 release catalogue, 83 (Stokes V table empty), nearest 29″ (3.1 mJy) — '
               'outside the 20″ radius, a fixed-sky background source (deep-graze antipode is a fixed sky point; §4); value +0.18 (0.8σ), peak 0.34, beam 13.5″ — empty at the position.'))
rows.append(ql('teegarden', 'B', '856 MHz', CUT, 'ledger_only', 0.07, 1, '2026-05-03',
               'b 0.98 R☉; FLASH SB84179 FLASH_389, Δt −3.2 d; Selavy file, 234 (2-h FLASH integration), nearest 163″ (0.7 mJy); '
               'value +0.03 (0.4σ), peak 0.10, beam 24.3″ — empty; the deepest antipode epoch (5σ ≈ 0.36 mJy).'))
rows.append(ql('van-maanen', 'B', '856 MHz', CUT, 'structurally_open', None, 1, '2026-03-29',
               'b 0.24 R☉; FLASH SB83234 FLASH_497, Δt −5.0 d, REJECTED: image and catalogue return "no permission" (failed the survey validation) — not checked.'))
rows.append(ql('gj-1276', 'A', '888 MHz', CUT, 'ledger_only', 0.24, 1, '2023-09-03',
               'b 0.83 R☉; VAST SB52549 VAST_2257-06, Δt −1.6 d; per-SBID catalogue, 70, nearest 219″ (10.7 mJy); '
               'value +0.64 mJy (2.7σ), peak 0.75, beam 12.8″ — the largest value at any position; §4: maximum of the 16 VAST_2257-06 epochs '
               '(2023-07 → 2026-06; mean +0.09 mJy, scatter 0.25, 16-epoch mean 0.09 ± 0.06 mJy = 1.5σ; no catalogue component within 60″ in any of 10 epochs; '
               'max of 16 Gaussian samples reaches 2.7σ ~5 % of the time) — not a detection, not a flare candidate (12-min epoch needs ≳ 1 mJy for 5σ). '
               'Source: surveys/radio-crossings/results/quicklook_gj1276_epochs.json.'))
rows.append(ql('van-maanen', 'A', '888 MHz', CUT, 'ledger_only', 0.38, 1, '2024-10-07',
               'b 0.33 R☉; VAST SB66449 VAST_0037+06, Δt +1.2 d; per-SBID catalogue, 59, nearest 266″ (27.8 mJy); value +0.33 (0.9σ), peak 0.63, beam 14.6″ — empty.'))
rows.append(ql('teegarden', 'A', '1368 MHz', CUT, 'ledger_only', 0.18, 1, '2024-11-11',
               'b 1.08 R☉; RACS-mid SB67840 RACS_0248+18, Δt +3.6 d; own Selavy catalogue (SB not in the release), 87, nearest 156″ (2.5 mJy); '
               'value +0.12 (0.7σ), peak 0.37, beam 11.9″ — empty.'))
rows.append(ql('teegarden', 'A', '144 MHz', 'LoTSS DR3 mosaic source catalogue (lotss_dr3.main_sources; all runs on the pointing, not epoch-resolved)',
               'ledger_only', 0.08, 1, '2023-11-04',
               'b 1.09 R☉; P043+19, Δt −3.9 d; 503 sources in 0.5°, nearest 71″ at 0.28 mJy — empty, but a mosaic of all runs, never epoch-resolved.',
               src=f'{R} §4 (LoTSS bullet); ' + QL))
rows.append(ql('gj-908', 'A', '3 GHz', VCUT, 'ledger_only', 0.13, 1, '2017-09-24',
               'b 12.3 R☉; VLASS1.1 T11t36.J235001+023000, Δt +2.7 d (in window); value +0.29 mJy (2.3σ), peak 0.30 — sub-threshold, '
               'the largest pixel value at any VLASS position; 1.2σ at the same star in epoch 4.1.', src=f'{R} §3 table; ' + QL))
rows.append(ql('gj-908', 'A', '3 GHz', VCUT, 'ledger_only', 0.13, 1, '2025-09-24',
               'b 12.3 R☉; VLASS4.1 T11t36.J235001+023000, Δt +2.8 d (in window); value +0.15 (1.2σ), peak 0.30 — noise.', src=f'{R} §3 table; ' + QL))
rows.append(ql('ross-154', 'A', '3 GHz', VCUT, 'ledger_only', 0.16, 1, '2019-07-02',
               'b 3.4 R☉; VLASS1.2 T05t29.J185051-233000, Δt −1.5 d (in window); value +0.10 (0.6σ), peak 0.22 — noise.', src=f'{R} §3 table; ' + QL))
rows.append(ql('teegarden', 'A', '3 GHz', VCUT, 'ledger_only', 0.11, 1, '2021-11-06',
               'b 1.1 R☉; VLASS2.2 T15t04.J025431+163000, Δt −2.0 d (in window); value −0.06, peak −0.02 — noise.', src=f'{R} §3 table; ' + QL))
rows.append(ql('gj-1276', 'A', '3 GHz', VCUT, 'structurally_open', 0.14, 1, '2020-09-13',
               'b 0.83 R☉; VLASS2.1 T09t35.J225404-063000, Δt +8.8 d (exact) — OUTSIDE the strict window (half-window 5.9 d; v1 tolerance edge); '
               'value +0.07, peak 0.30 — noise, out of window. Largest one-beam peak at any VLASS position is 2.9σ (gj-1276, VLASS3.1, 2023, out of window).',
               src=f'{R} §3 table; ' + QL))
rows.append(ql('wolf-359', 'B', '3 GHz', VCUT, 'structurally_open', 0.14, 1, '2020-09-13',
               'b 0.72 R☉; VLASS2.1 T09t35 (two overlapping tiles), Δt +8.3 d (exact) — OUTSIDE the strict window (tolerance edge); '
               'values −0.15 / −0.16, peaks 0.19 / 0.22 — noise, out of window. Strict VLASS in-window set is four channel-A epochs, none at the antipode.',
               src=f'{R} §3 table; ' + QL))
dump('radio_quicklook', {
    'survey_id': 'radio_quicklook', 'report': R, 'plan_section': '5.8', 'pipeline': 'B',
    'archives': ['casda', 'vast', 'racs', 'flash', 'lotss', 'vlass'], 'hypothesis_version': None,
    'era': 'in-window epochs 2017-09-24 → 2026-05-03 (§2–3 tables)',
    'run_dir': 'runs/radio-crossings/quicklook', 'records_dir': None, 'candidates': 0,
    'candidate_notes': '"A look, not a search": no threshold frozen, no statistic claimed, no power constraint derived; '
                       '17 epochs (10 ASKAP, 1 LoTSS, 6 VLASS): 0 coincident catalogue components within 20″, 0 cutout pixels above 3σ '
                       '(12 ASKAP + 26 VLASS cutouts); largest value 2.7σ (gj-1276 star, VAST 2023-09) = noise-series maximum. '
                       'A 3σ point source is 0.2–0.7 mJy at 856–1656 MHz (ASKAP) and 0.35–0.5 mJy at 3 GHz (VLASS).',
    'rows': rows})

# ---------------------------------------------------------------------------
# 4. report/spectral_archives.md  (plan §5.17)
# ---------------------------------------------------------------------------
R = 'report/spectral_archives.md'
CV = 'surveys/spectral-archives/results/confirmatory_v1.md'
DV = 'surveys/spectral-archives/results/dev_v1.md'
rows = []
SUB = {
    'NIRPS': 'ESO phase-3 NIRPS 1D spectra (telluric-corrected, sky fibre)',
    'SPIRou': 'CFHT SPIRou APERO 1D spectra via CADC',
    'HARPS': 'ESO phase-3 HARPS s1d spectra',
    'ESPRESSO': 'ESO phase-3 ESPRESSO s1d spectra (vacuum WAVE, A9)',
    'CARMENES-VIS': 'CARMENES DR1 VIS spectra (no sky fibre)',
}
BAND = {'532': '532 nm (515–545 nm)', '1064': '1064 nm (1030–1090 nm)', '1550': '1550 nm (1530–1570 nm)',
        'generic': 'generic (whole band)'}
RUNG = {'1.2Rsun': '1.2 Rsun', '2.5Rsun': '2.5 Rsun', '0.1AU': '0.1 AU'}
PW = 'W (P = 1.0645 A90 F_λ FWHM_λ π b²; continuous line through the rung cone)'


def sp(target, inst, date, b, cell, rung, watts, a90, n_spec, n_trials, exc, note, dev=False):
    rung_s = RUNG[rung]
    if watts is None:
        status, lk, lu, pmw = 'constraint_only', None, None, None
        note = ('Unconstrained at v1: A90 not recovered ≤ 2.0 (200 % injection cap; systematics-limited). ' + note).strip()
    else:
        status, lk, lu, pmw = 'searched', 'power', PW, float(f'{watts / 1e6:.4g}')
        note = (note + ' power_mw converted from W.').strip()
    return row(target_id=target, channel='A', rung=rung_s, band=BAND[cell], substrate=SUB[inst], status=status,
               limit_kind=lk, limit_value=watts, limit_unit=lu, power_mw=pmw, n_events=n_spec, n_trials=n_trials,
               epoch_range=date,
               source=f'{R} §3–4 tables; {DV if dev else CV} (unit {target}|{inst}|{date})',
               notes=(f'b {b} R☉; A90 {a90}; ' if a90 else f'b {b} R☉; ') +
                     (f'exceedances {exc}; ' if exc else 'no exceedance; ') +
                     ('dev unit (not blind); ' if dev else '') + note)


# wolf-359 NIRPS 2025-03-03 (b 0.76): 1550, 1064, generic × three rungs
for cell, a90, p in [('1550', 0.066, {'1.2Rsun': 12.3, '2.5Rsun': 53.2, '0.1AU': 3.94e3}),
                     ('1064', 0.100, {'1.2Rsun': 21.8, '2.5Rsun': 94.6, '0.1AU': 7.0e3}),
                     ('generic', 1.622, {'1.2Rsun': 203, '2.5Rsun': 882, '0.1AU': 6.52e4})]:
    for rung, w in p.items():
        exc = {'1550': '2 (S_line 8.4/7.5, S_coadd 17.4/9.0 at 1569.499 nm → telluric: CO₂ line, NIRPS ATM_TRANSM 0.79; absent in SPIRou 2019/2021 in-window spectra)',
               '1064': None, 'generic': '1 (S_coadd 62/58 at 1807 nm → telluric, 1.8 µm band)'}[cell]
        note = {'1550': 'the first 1550 nm line-SED constraint of the programme (deepest: 12 W through 1.2 R☉).',
                '1064': 'the first 1064 nm constraint (22 W through 1.2 R☉); cell null (5.5/7.3, 10.3/10.8).',
                'generic': 'generic cell 966–1923 nm; F_λ from the measured NIRPS absolute SED.'}[cell]
        rows.append(sp('wolf-359', 'NIRPS', '2025-03-03', 0.76, cell, rung, w, a90, 8, 2, exc, note))
# wolf-359 SPIRou 2021-03-03 (0.77)
for cell, a90, p in [('1064', 0.93, {'1.2Rsun': 193, '2.5Rsun': 839, '0.1AU': 6.21e4}),
                     ('1550', 1.24, {'1.2Rsun': 221, '2.5Rsun': 960, '0.1AU': 7.1e4}),
                     ('generic', None, {'1.2Rsun': None, '2.5Rsun': None, '0.1AU': None})]:
    for rung, w in p.items():
        exc = {'1064': '1 (S_coadd 37.7/28.6 at 1080.088 nm → telluric: OH-subtraction plateau, 2021-03-02 night)',
               '1550': None, 'generic': '1 (S_coadd 281/71 at 1355 nm → telluric, 1.35 µm band)'}[cell]
        rows.append(sp('wolf-359', 'SPIRou', '2021-03-03', 0.77, cell, rung, w, a90, 20, 2, exc, ''))
# wolf-359 SPIRou 2019-03-03 (0.78), 0.1 AU only
rows.append(sp('wolf-359', 'SPIRou', '2019-03-03', 0.78, '1064', '0.1AU', 8.34e4, 1.25, 4, 2, None, 'null (all cells).'))
rows.append(sp('wolf-359', 'SPIRou', '2019-03-03', 0.78, '1550', '0.1AU', 7.23e4, 1.26, 4, 2, None, 'null (all cells).'))
rows.append(sp('wolf-359', 'SPIRou', '2019-03-03', 0.78, 'generic', '0.1AU', None, None, 4, 2, None, ''))
# wolf-359 CARMENES-VIS 2016/17/18 (0.79), generic, 0.1 AU
for date, n, nt, exc in [('2016-03-03', 1, 1, None),
                         ('2017-03-03', 7, 2, '1 (S_coadd 97.9/90.3 at 760.024 nm → telluric: O₂ A band)'),
                         ('2018-03-03', 1, 1, None)]:
    rows.append(sp('wolf-359', 'CARMENES-VIS', date, 0.79, 'generic', '0.1AU', None, None, n, nt, exc,
                   'CARMENES generic cell declared unconstrained (airglow-dominated statistics; no sky fibre).'))
# ross-128 HARPS 2021-03-17 (1.81): 532 at 2.5 Rsun + 0.1 AU; generic unconstrained
rows.append(sp('ross-128', 'HARPS', '2021-03-17', 1.81, '532', '2.5Rsun', 37.2, 0.556, 18, 2, None,
               'first 532 nm constraint through a 2.5 R☉ cone (37 W); cell null (S_line 14.7/16.0; coadd 4.3/9.1). F_λ from BP–G interpolation, ±30 % systematic.'))
rows.append(sp('ross-128', 'HARPS', '2021-03-17', 1.81, '532', '0.1AU', 2.75e3, 0.556, 18, 2, None, 'cell null.'))
for rung in ['2.5Rsun', '0.1AU']:
    rows.append(sp('ross-128', 'HARPS', '2021-03-17', 1.81, 'generic', rung, None, None, 18, 2,
                   '1 (S_coadd 69.5/35.9 at 557.891 nm → telluric: [O I] 5577 airglow, in all 18 spectra)', ''))
# ross-128 SPIRou 2022-03-17 (1.81)
for cell, a90, p in [('1550', 0.235, {'2.5Rsun': 404, '0.1AU': 2.99e4}),
                     ('1064', 1.024, {'2.5Rsun': 2.9e3, '0.1AU': 2.14e5}),
                     ('generic', 1.822, {'2.5Rsun': 4.18e3, '0.1AU': 3.09e5})]:
    for rung, w in p.items():
        exc = '1 (S_coadd 139/127 at 1456.410 nm → telluric: OH sky line, model 208)' if cell == 'generic' else None
        rows.append(sp('ross-128', 'SPIRou', '2022-03-17', 1.81, cell, rung, w, a90, 40, 2, exc,
                       'NIR F_λ = photometric interpolation × wolf-359 ratio (conservative, limit-weakening).'))
# ross-128 HARPS 2006/13/14/15/19 (1.80–1.81), 0.1 AU
for date, b, n, nt, w, a90, exc in [('2006-03-17', 1.80, 2, 2, 2.74e3, 0.554, None),
                                    ('2013-03-17', 1.81, 2, 2, 4.05e3, 0.818, None),
                                    ('2014-03-17', 1.81, 2, 2, 3.73e3, 0.753, None),
                                    ('2015-03-18', 1.81, 3, 2, 4.59e3, 0.928, None),
                                    ('2019-03-18', 1.81, 1, 1, 2.79e3, 0.564, None)]:
    rows.append(sp('ross-128', 'HARPS', date, b, '532', '0.1AU', w, a90, n, nt, None, 'cell null.'))
    gexc = '1 (S_coadd 37.6/28.8 at 557.888 nm → telluric: [O I] 5577 airglow)' if date == '2015-03-18' else None
    rows.append(sp('ross-128', 'HARPS', date, b, 'generic', '0.1AU', None, None, n, nt, gexc, ''))
# ross-128 ESPRESSO 2019/20/21/22 (1.81), 0.1 AU
for date, n, nt, w, a90 in [('2019-03-18', 2, 2, 872, 0.214), ('2020-03-17', 2, 2, 922, 0.227),
                            ('2021-03-17', 3, 2, 944, 0.232), ('2022-03-17', 1, 1, 1.4e3, 0.345)]:
    exc532 = exc_gen = None
    note = 'cell null; ESPRESSO grid re-run in the corrected vacuum frame (A9).'
    if date == '2020-03-17':
        exc532 = '2 (S_line 20.3/14.9, S_coadd 17.6/10.3 → defect/contamination: Hg I lamp light in the science fibre of the 2020-03-16 spectrum, 31 peaks above threshold; the 2020-03-22 spectrum is clean, 532 S 6.0/14.9)'
        exc_gen = '2 (S_line 465/127, S_coadd 458/110 at 546.228 nm → defect/contamination, same Hg-lamp spectrum; clean spectrum generic 44/127)'
        note = 'archive fact: ESPRESSO ross-128 2020-03-16 (ADP.2021-04-27T19:22:55.144) carries Hg I 404.8/436.0/546.2/577.1/579.2 nm lamp contamination.'
    rows.append(sp('ross-128', 'ESPRESSO', date, 1.81, '532', '0.1AU', w, a90, n, nt, exc532, note))
    rows.append(sp('ross-128', 'ESPRESSO', date, 1.81, 'generic', '0.1AU', None, None, n, nt, exc_gen, ''))
# ross-128 CARMENES-VIS 2016/17
for date, n, nt in [('2016-03-17', 1, 1), ('2017-03-17', 2, 2)]:
    rows.append(sp('ross-128', 'CARMENES-VIS', date, 1.81, 'generic', '0.1AU', None, None, n, nt, None,
                   'CARMENES generic cell declared unconstrained.'))
# teegarden CARMENES-VIS 2017-11-08 (1.12): 0.1 AU + 2.5 Rsun
for rung in ['2.5Rsun', '0.1AU']:
    rows.append(sp('teegarden', 'CARMENES-VIS', '2017-11-08', 1.12, 'generic', rung, None, None, 7, 2,
                   '2 (S_line 83/67 at 882.950 nm and S_coadd 93/73 at 775.276 nm → telluric: OH 7-3 8827.10 Å and OH 9-4 7750.64 Å airglow; OH forest at z 9–50; amplitude anticorrelates with SNR)',
                   'CARMENES generic cell declared unconstrained; the only grazing-rung teegarden unit at v1.'))
for date, b, n, nt, exc in [('2016-11-08', 1.12, 1, 1, None), ('2018-11-08', 1.11, 3, 2, None),
                            ('2019-11-08', 1.11, 1, 1, '1 (S_line 82/67 at 882.956 nm → telluric: the same OH 8827 line, same observer wavelength)')]:
    rows.append(sp('teegarden', 'CARMENES-VIS', date, b, 'generic', '0.1AU', None, None, n, nt, exc,
                   'CARMENES generic cell declared unconstrained.'))
# teegarden ESPRESSO 2024-11-08 (1.08)
rows.append(sp('teegarden', 'ESPRESSO', '2024-11-08', 1.08, '532', '0.1AU', 129, 0.976, 4, 2,
               '2 (S_line 37.9/8.5, S_coadd 31.9/7.7 at 532.100 nm → defect/cosmic: 2-px spike z 23/29/8 in one spectrum, 2024-11-06T03:06, absent in the other three)',
               'deepest 0.1 AU line limit of the survey (faint star: F_λ 4.67e-15; report §4 quotes 107 W in the reading paragraph vs 129 W in the table — table value transcribed). 27 usable nulls.'))
rows.append(sp('teegarden', 'ESPRESSO', '2024-11-08', 1.08, 'generic', '0.1AU', None, None, 4, 2, None, 'generic null (134/171, 50/221).'))
# dev units
rows.append(sp('gj-908', 'HARPS', '2011-09-21', 12.32, '532', '0.1AU', 8.5e3, 0.268, 3, 2, None,
               'dev pass v1.1: 0 exceedances (532 4.9/19.3; generic 34.3/48.6). Bright star, photometric F_λ.', dev=True))
rows.append(row(target_id='gj-908', channel='A', rung='0.1 AU', band=BAND['generic'], substrate=SUB['HARPS'],
                status='constraint_only', n_events=3, n_trials=2, epoch_range='2011-09-21',
                source=f'{R} §2, §4 (dev completeness); {DV}',
                notes='b 12.32 R☉; dev unit; A90 2.000 / 2.000 sits at the 200 % injection cap (dev_v1.md lists 5.59e+04 W); '
                      'report §4 declares generic cells beyond 200 % of the local continuum unconstrained — limit not carried.'))
rows.append(sp('ross-154', 'HARPS', '2017-07-03', 3.40, '532', '0.1AU', 4.41e3, 0.513, 6, 2, None,
               'dev pass v1.1: 0 exceedances (532 14.5/15.6). Bright star, photometric F_λ.', dev=True))
rows.append(row(target_id='ross-154', channel='A', rung='0.1 AU', band=BAND['generic'], substrate=SUB['HARPS'],
                status='searched', limit_kind='power', limit_value=1.29e4, limit_unit=PW, power_mw=0.0129,
                n_events=6, n_trials=2, epoch_range='2017-07-03', source=f'{DV} (unit ross-154|HARPS|2017-07-03, A90 1.865 / 2.000)',
                notes='b 3.40 R☉; dev unit (not blind); power_mw converted from W; generic 0.1 AU value from dev_v1.md only — the report §4 quotes only the 532 cell for dev units.'))
rows.append(sp('ross-154', 'CARMENES-VIS', '2018-07-03', 3.40, 'generic', '0.1AU', None, None, 3, 2,
               '1 (S_coadd 98.7/92.7 at 760.026 nm → telluric: O₂ A band)', 'dev outcome: 10 trials, 1 exceedance vs 0.6 expected over the three v1.1 dev units.', dev=True))
# not searchable / excluded / uncovered
rows.append(row(target_id='gj-908', channel='A', rung='0.1 AU', band='any', substrate='ESO phase-3 X-shooter NIR spectra',
                status='structurally_open', n_events=36, epoch_range='2010-09-21',
                source=f'{R} §2 (amendment A8); {DV}',
                notes='b 12.32 R☉; dev unit excluded by amendment A8 (4 product variants per exposure, single-night nulls, no BERV); '
                      'its 5 dev-pass exceedances are not carried (dev_v1.md).'))
rows.append(row(target_id='wolf-359', channel='A', rung='0.1 AU', band='any', substrate=SUB['HARPS'],
                status='structurally_open', epoch_range='2019-03-03', source=f'{R} §3 (combos); {CV}',
                notes='wolf-359|HARPS combo has only 18 usable out-of-window nulls (V 13.5, short exposures; ≥ 20 required) → not searchable.'))
rows.append(row(target_id='wolf-359', channel='A', rung='2.5 Rsun', band='any', substrate='OHP SOPHIE public s1d spectra',
                status='structurally_open', epoch_range='2022-03-04', source=f'{R} §2; §6 hand-offs',
                notes='SOPHIE wolf-359 2022 unit not admitted: public s1d header is date-stripped (BJD rounded to the day), '
                      'so the frozen "header UT inside the window" condition cannot be met; hand-off (hour requires the observatory).'))
for rung in ['1.2 Rsun', '2.5 Rsun', '0.1 AU']:
    rows.append(row(target_id='van-maanen', channel='A', rung=rung, band='any',
                    substrate='11 spectral archives (recon 2026-09-05)', status='structurally_open',
                    source=f'{R} §1 (units); §6',
                    notes='No unit: van-maanen has no échelle spectrum in any grazing window in any archive (deepest graze, b 0.32–0.46 R☉); '
                          'future-observation recommendation: one 20-min HARPS/ESPRESSO/NIRPS exposure at t_ca ± 7 h on Oct 6–7 of any year.'
                          if rung != '0.1 AU' else 'No unit fixed for van-maanen at any rung (§1: "van-maanen and gj-1276 have none").'))
    rows.append(row(target_id='gj-1276', channel='A', rung=rung, band='any',
                    substrate='11 spectral archives (recon 2026-09-05)', status='structurally_open',
                    source=f'{R} §1 (units); §6', notes='No unit: gj-1276 is uncovered in every archive.'))
# hand-offs (raw / non-public frames in windows)
HO = f'{R} §5–6 hand-offs'
for tgt, rung, band, sub, ep, note in [
        ('wolf-359', '1.2 Rsun', 'any', 'KOA HIRES raw frames (I₂ cell in)', '2010-03-03',
         'Hand-off: KOA HIRES wolf-359 2010-03-03 (1.2 R☉); Tellis & Marcy 2017 may already carry a per-spectrum null — cross-match their table.'),
        ('teegarden', '2.5 Rsun', 'any', 'Gemini MAROON-X raw frames', '2021-11-08',
         'Hand-off: MAROON-X teegarden 2021-11-08 at +3.3 h (2.5 R☉); non-public/raw.'),
        ('ross-128', '2.5 Rsun', 'any', 'Gemini MAROON-X raw frames', '2024-03-17',
         'Hand-off: MAROON-X ross-128 2024-03-17 at +1.4 h (2.5 R☉); non-public/raw.'),
        ('teegarden', '0.1 AU', 'generic (K band)', 'ESO CRIRES raw frames', '2009-11-08',
         'Hand-off: CRIRES teegarden 2009-11-08 (K band, generic cell).'),
        ('teegarden', '0.1 AU', 'any', 'BL APF raw 2D frames', '2016-11-13',
         'Hand-off: BL APF teegarden 2016 (raw 2D, 0.1 AU) — the 3 spectra flagged in report/radio_crossings.md §2.'),
        ('teegarden', '1.2 Rsun', 'any', 'CFHT SPIRou (proprietary to ~2026-11)', '2025-11-08',
         'Hand-off: SPIRou teegarden 2025-11-08 ×4 at +0.2 h (1.2 R☉); pre-registered re-run when public.'),
        ('teegarden', '0.1 AU', 'generic (810–1280 nm)', 'HPF spectra (data ask)', '2019–2020',
         'Hand-off: HPF teegarden 2019/2020, 8 epochs at 0.1 AU; data ask.'),
        ('teegarden', '2.5 Rsun', 'any', 'CARMENES NIR channel (Calar Alto form)', '2017-11-08',
         'Hand-off: CARMENES NIR twins of the 2017 grazing pair (OH-corrected reduction).')]:
    rows.append(row(target_id=tgt, channel='A', rung=rung, band=band, substrate=sub, status='structurally_open',
                    epoch_range=ep, source=HO, notes=note))
dump('spectral_archives', {
    'survey_id': 'spectral_archives', 'report': R, 'plan_section': '5.17', 'pipeline': 'B',
    'archives': ['eso', 'harps', 'espresso', 'nirps', 'cadc', 'spirou', 'carmenes'],
    'hypothesis_version': 'v1.0 + A1–A9 (amendment v1.1; A9 ESPRESSO vacuum frame during adjudication)',
    'era': '2006→2025', 'run_dir': 'runs/spectral-archives/v1', 'records_dir': None, 'candidates': 0,
    'candidate_notes': 'Blind confirmatory: 24 units, 77 trials, 18 exceedances vs 4.7 expected — 12 sky/telluric, 4 Hg-lamp contamination '
                       '(one ESPRESSO spectrum), 2 cosmic (one 2-px feature); 0 retained-ambiguous (results/adjudication_confirmatory_v1.json). '
                       'n_trials is per (unit, cell): rows of one unit at different rungs repeat the same trials. Dev: 10 trials, 1 exceedance (telluric). Recon: 7,644 on-star spectra; 633/30/14 in 0.1 AU / 2.5 R☉ / 1.2 R☉ windows; 29 units fixed, 24 confirmatory searched.',
    'rows': rows})

# ---------------------------------------------------------------------------
# 5. report/highenergy_crossings.md  (plan §5.26)
# ---------------------------------------------------------------------------
R = 'report/highenergy_crossings.md'
HE = ROOT / 'surveys/highenergy-crossings/results'
comp = {u['unit']: u for u in json.load(open(HE / 'completeness_v1.json'))['units']}
conf = {u['unit']: u for u in json.load(open(HE / 'confirmatory_v1.json'))['units']}
cov = {}
for line in (HE / 'coverage_v1.md').read_text().splitlines():
    m = re.match(r'\| ([AB] \S+ \S+) \| (\d+) \| (\d+) \| ([\d.]+) \|', line)
    if m:
        cov[m.group(1)] = (int(m.group(2)), int(m.group(3)), float(m.group(4)))
assert sum(v[0] for v in cov.values()) == 580, cov
unc = {}
for line in (HE / 'coverage_v1.md').read_text().splitlines():
    m = re.match(r'- ([AB]) (\S+) (\S+) (\d{4}-\d\d-\d\d):', line)
    if m:
        unc[(m.group(1), m.group(3))] = unc.get((m.group(1), m.group(3)), 0) + 1
assert sum(unc.values()) == 62, unc
HRUNG = {'1.2Rsun': '1.2 Rsun', '2.5Rsun': '2.5 Rsun', '0.1AU': '0.1 AU'}
LAT_L = 'Fermi-LAT P8R3 SOURCE photons, lane L (100 MeV–1 GeV, r = 2.0°), gated (θ ≤ 60°, DATA_QUAL > 0, LAT_CONFIG = 1, Moon/Sun > 8°); Poisson vs exposure-matched pseudo-window rate'
LAT_H = 'Fermi-LAT P8R3 SOURCE photons, lane H (1–300 GeV, r = 0.7°), same gating; Poisson vs exposure-matched pseudo-window rate'
LAT_U = 'Fermi-LAT P8R3 SOURCE photons, union band (0.1–300 GeV, r = 1.0°), max over 1-ks boxcars'
COMP = 'surveys/highenergy-crossings/results/completeness_v1.md'
COMPJ = 'surveys/highenergy-crossings/results/completeness_v1.json'
CONF = 'surveys/highenergy-crossings/results/confirmatory_v1.md'
COVM = 'surveys/highenergy-crossings/results/coverage_v1.md'
rows = []


def f(x, d=2):
    return f'{x:.{d}e}'


order = list(comp.keys())
for unit in order:
    ch, tgt, rg = unit.split()
    rung = HRUNG[rg]
    c, k, cv = comp[unit], conf[unit], cov[unit]
    st, dg = k['stats'], k['diag']
    n_win = k['n_searchable']
    common = f"{cv[0]} windows / {cv[1]} searchable, Σ gated livetime {cv[2]} ks ({COVM})"
    vetoed = (ch == 'B' and tgt == 'van-maanen')
    # lane L (+ burst statistic for this unit)
    if vetoed:
        rows.append(row(target_id=tgt, channel=ch, rung=rung, band='100 MeV–1 GeV (lane L) + 0.1–300 GeV burst', substrate=LAT_L,
                        status='vetoed_known_source', n_events=n_win, n_trials=3, epoch_range='2014-04-03 window (3C 279 April-2014 flare)',
                        source=f'{R} §2 (blind confirmatory); {CONF}; {COMP}',
                        notes=(f"Lane L and burst trials constraint_only at the threshold freeze (amendment v1.2: 3C 279 = 4FGL J1256.1−0547, variability index 33,299, "
                               f"1.8° from the antipode inside the 2° aperture; 4FGL J1249.3−0545 at 0.4°, 4FGL J1245.4−0701 at 1.9°); S_stack_L {st['S_stack_L']:.2f}, "
                               f"S_event_L {st['S_event_L']:.2f}, S_burst {st['S_burst']:.2f} (T 3.451); n_L/λ_L {dg['n_tot_L']} / {dg['lam_tot_L']:.1f}; "
                               f"excess-photon centroids 0.78–0.84° from 3C 279; strict re-run (θ ≤ 50°, zenith ≤ 90°, Moon/Sun > 16°) leaves them in place; "
                               f"disposition vetoed_known_source (ladder step 1). Completeness: not_constrainable in lane L and burst (real statistic exceeds T); "
                               f"pseudo-window over-dispersion 7–20 %. {common}")))
    else:
        rows.append(row(target_id=tgt, channel=ch, rung=rung, band='100 MeV–1 GeV (lane L)', substrate=LAT_L, status='searched',
                        limit_kind='flux', limit_value=float(f(c['F90_stack_L'])),
                        limit_unit='ph cm^-2 s^-1 (F90, recurrence-stacked over all searchable windows; ⟨E⟩ = 256 MeV, Γ = 2)',
                        power_mw=round(c['P90_stack_L'] / 1e6, 4), n_events=n_win, n_trials=2,
                        source=f'{R} §2; {COMP}; {CONF}; {COVM}',
                        notes=(f"S_stack_L {st['S_stack_L']:.2f}, S_event_L {st['S_event_L']:.2f} (T 3.451); n_L/λ_L {dg['n_tot_L']} / {dg['lam_tot_L']:.1f}; "
                               f"F90 event L best {f(c['F90_event_L_best'])} ({c['event_L_best_window']}) / median {f(c['F90_event_L_median'])}; "
                               f"P90 stack L {f(c['P90_stack_L'])} W (power_mw converted from W), P90 event L best {f(c['P90_event_L_best'])} W; "
                               f"burst Φ90 / Φ50 {f(c['Phi90_burst'])} / {f(c['Phi50_burst'])} ph cm^-2 (1 ks; S_burst {st['S_burst']:.2f}). {common}")))
    # lane H
    rows.append(row(target_id=tgt, channel=ch, rung=rung, band='1–300 GeV (lane H)', substrate=LAT_H, status='searched',
                    limit_kind='flux', limit_value=float(f(c['F90_stack_H'])),
                    limit_unit='ph cm^-2 s^-1 (F90, recurrence-stacked; ⟨E⟩ = 5.7 GeV, Γ = 2)',
                    power_mw=round(c['P90_stack_H'] / 1e6, 4), n_events=n_win, n_trials=2,
                    source=f'{R} §2; {COMPJ} (P90_stack_H); {COMP}; {CONF}',
                    notes=(f"S_stack_H {st['S_stack_H']:.2f}, S_event_H {st['S_event_H']:.2f} (T 3.451); n_H/λ_H {dg['n_tot_H']} / {dg['lam_tot_H']:.2f}; "
                           f"F90 event H best {f(c['F90_event_H_best'])}; P90 stack H {f(c['P90_stack_H'])} W (power_mw converted from W), P90 event H best {f(c['P90_event_H_best'])} W. "
                           + ("Lane H calibrated for this unit (the 3C 279 veto applies to lane L and burst only)." if vetoed else ''))))
# burst per unit (26 calibrated S_burst units; report per-rung aggregates in notes)
BURST_AGG = {'1.2Rsun': ('2.7e-3 / 1.2e-2', 6), '2.5Rsun': ('3.2e-3 / 1.0e-2', 9), '0.1AU': ('3.7e-3 / 1.5e-2', 11)}
for unit in order:
    ch, tgt, rg = unit.split()
    if ch == 'B' and tgt == 'van-maanen':
        continue
    c, k = comp[unit], conf[unit]
    agg, nu = BURST_AGG[rg]
    rows.append(row(target_id=tgt, channel=ch, rung=HRUNG[rg], band='0.1–300 GeV (burst, 1 ks)', substrate=LAT_U,
                    status='searched', limit_kind='fluence', limit_value=float(f(c['Phi90_burst'])),
                    limit_unit='ph cm^-2 (Φ90, single 1-ks pulse, 90 % recovery in the best-exposure window)',
                    n_events=k['n_searchable'], n_trials=1,
                    source=f'{R} §2 90 %-recovery table; {COMP}; {CONF}',
                    notes=(f"S_burst {k['stats']['S_burst']:.2f} (T 3.451); Φ50 {f(c['Phi50_burst'])} ph cm^-2; burst injection window {c['burst_window']}. "
                           f"Report per-rung Φ50 / Φ90 over the {nu} units of this rung: {agg} ph cm^-2 — Φ90 reached only at 1–2e-2 because a 1-ks pulse near the end of a "
                           f"good-time run is mostly lost. 0 exceedances in 26 burst trials (0.009 expected). Pulses shorter than 1 ks and coherent trains declared unconstrained.")))
# dev units (constraint-only; no completeness)
for unit, note in [('A gj-908 0.1AU', 'dev unit'), ('B gj-908 0.1AU', 'dev unit'),
                   ('A gj-1276 1.2Rsun', 'demoted to dev by amendment v1.1 (statistics printed in a machinery test)')]:
    ch, tgt, rg = unit.split()
    cv = cov[unit]
    rows.append(row(target_id=tgt, channel=ch, rung=HRUNG[rg], band='0.1–300 GeV (lanes L, H, burst)', substrate=LAT_L,
                    status='constraint_only', n_events=cv[1], n_trials=5,
                    source=f'{R} §2 (dev); surveys/highenergy-crossings/results/dev_v1.md; {COVM}',
                    notes=f'{note}; all 15 dev statistics below T (max S_event_L 3.16 over the three dev units); not blind, no completeness/limit derived. '
                          f'{cv[0]} windows / {cv[1]} searchable, Σ gated livetime {cv[2]} ks.'))
# uncovered windows per (channel, rung) — counted from the coverage-without-statistic bullet list
for (ch, rg), n in sorted(unc.items()):
    if n == 0:
        continue
    rows.append(row(target_id='programme', channel=ch, rung=HRUNG[rg], band='0.1–300 GeV', substrate=LAT_L, status='structurally_open',
                    n_events=n, source=f'{COVM} (coverage-without-statistic list, counted per channel and rung); {R} §2 coverage',
                    notes=f'{n} windows with < 1 ks gated lane-L livetime (coverage-without-statistic), all post-2018-03 (solar-array-anomaly survey profile); '
                          f'rung totals in the report: 1.2 R☉ 27, 2.5 R☉ 27, 0.1 AU 8 uncovered of 580 windows. Moon/Sun 8° exclusion removed 8.0 of 41.6 Ms (137 of 580 windows touched).'))
# XMM / eROSITA not constrainable; Chandra none
for ch in 'AB':
    for rg in ['1.2Rsun', '2.5Rsun', '0.1AU']:
        rows.append(row(target_id='programme', channel=ch, rung=HRUNG[rg], band='any', substrate='XMM-Newton pointed archive',
                        status='not_constrainable', source=f'{R} §2 (geometry); §4',
                        notes='Both channels sit at solar elongation 180° in every window; XMM-Newton solar aspect 70–110° never observes a crossing (Gate I of §5.3, now four members). 7 deep-family targets.'))
        rows.append(row(target_id='programme', channel=ch, rung=HRUNG[rg], band='any', substrate='eROSITA all-sky scan (eRASS)',
                        status='not_constrainable', source=f'{R} §2 (geometry); §4',
                        notes='eROSITA scans at 90° elongation; never observes a crossing window (Gate I). 7 deep-family targets.'))
        rows.append(row(target_id='programme', channel=ch, rung=HRUNG[rg], band='any', substrate='Chandra pointed archive',
                        status='structurally_open', source=f'{R} §2 (geometry)',
                        notes='Chandra has no in-window observation in 27 years. 7 deep-family targets.'))
# Swift/BAT arm
BATM = 'surveys/highenergy-crossings/results/bat_v1.md'
for line in (HE / 'bat_v1.md').read_text().splitlines():
    m = re.match(r'\| ([AB]) \| (\S+) \| (\S+) \| (\d+) \| (\d+) \| (\d+) \| (\d+) \| ([\d.e+]+) \| (\S+) \|', line)
    if not m:
        continue
    ch, tgt, rg, nwin, ncov, med, mn, pw, snr = m.groups()
    note = (f'Constraint-only arm (no threshold; BAT survey few-σ systematics floor). {nwin} grazing windows, {ncov} with ≥ 1 usable pointing '
            f'(offset ≤ 20°, partial coding ≥ 0.1); min 3σ UL {mn} mCrab; max pointing SNR {snr}. Crab measured in 24 pointings (0.0408 units) — empirical mCrab scale; '
            f'Crab 14–195 keV = 2.31e-08 erg cm^-2 s^-1; power_mw = median 3σ power {pw} W converted from W. batsurvey (HEASoft 6.24, chbrandt/heasoft container, local CALDB), 238 pointings.')
    if ch == 'A' and tgt == 'gj-1276' and rg == '2.5Rsun':
        note += ' Two windows at combined SNR 3.1–3.2 (2021-09-04, 2 pointings, rate 9.28e-04 ± 2.93e-04; 2025-09-04, 4 pointings, 9.78e-04 ± 3.17e-04) recorded, not adjudicated.'
    rows.append(row(target_id=tgt, channel=ch, rung=HRUNG[rg], band='14–195 keV', substrate='Swift/BAT survey pointings, batsurvey per-pointing rates at the channel position',
                    status='constraint_only', limit_kind='flux', limit_value=int(med),
                    limit_unit='mCrab (median 3σ upper limit over covered windows, 14–195 keV)',
                    power_mw=round(float(pw) / 1e6, 3), n_events=int(ncov),
                    source=f'{R} §2 Swift/BAT arm; {BATM}', notes=note))
# Swift/XRT look
rows.append(row(target_id='wolf-359', channel='A', rung='0.1 AU', band='0.3–10 keV', substrate='Swift/XRT PC-mode cleaned events, obsid 00010119025 (47″ aperture, 100-s bins)',
                status='ledger_only', n_events=1, epoch_range='2018-03-06 (+2.9 d)',
                source=f'{R} §2 Swift/XRT look; surveys/highenergy-crossings/results/xrt_look_v1.md',
                notes='A look, not a trial (hypotheses §9): 471 net counts in 4,160 s, mean 0.112 ct/s, 41 bins, median 0.108, max 0.219, no excursion above 3× median; '
                      'elevated ~3.6× quiescent (2017 campaign 23 visits median 0.031; 2021-12 campaign 80 visits median 0.029; flare visits up to 0.67) — smooth, nothing pulse-like; no constraint. '
                      'The only in-window Swift XRT observation of the survey. UVOT UV-grism spectrum of the same visit handed to §5.17.'))
# corridor catalogue screen (Pipeline A step 1)
cs = json.load(open(HE / 'corridor_screen_v1.json'))
CSM = 'surveys/highenergy-crossings/results/corridor_screen_v1.md'
CSJ = 'surveys/highenergy-crossings/results/corridor_screen_v1.json'
for c in cs['corridors']:
    hits = c['hits']
    nfix = sum(h['disposition'] == 'fixed_sky' for h in hits)
    nunr = sum(h['disposition'] == 'single_epoch_unresolved' for h in hits)
    closest = min((h['sep_arcmin'] for h in hits), default=None)
    ul1 = next((u for u in c['upper_limits'] if u['dr_survey'] == 'DR1_eRASS1' and u['band'] == '024'), None)
    ul3 = next((u for u in c['upper_limits'] if u['dr_survey'] == 'DR2_eRASSc3' and u['band'] == '024'), None)
    de = c['de_sky']
    note = (f"Pipeline A step 1 catalogue screen (report label catalogue_screen_only); corridor status {c['status']}: {c['n_hits']} catalogued X-ray/γ-ray sources within 7′, "
            f"{nfix} fixed_sky (SGL-track test: positions agree to < 10″ across epochs with predicted relay displacement ≥ 15″ at z = 10,000 AU, or compact multi-scan stack), "
            f"{nunr} single_epoch_unresolved" + (f"; closest {closest:.1f}′" if closest is not None else '') + '. ')
    if de:
        note += (f"eROSITA-DE upper limits (0.2–2.3 keV): eRASS1 {ul1['UL_B']:.1e} erg cm^-2 s^-1 (exposure {ul1['Exposure']:.0f} s), "
                 f"eRASS:3 {ul3['UL_B']:.1e} (exposure {ul3['Exposure']:.0f} s).")
        status, lk, lv, lu = 'constraint_only', 'flux', float(f'{ul1["UL_B"]:.2e}'), 'erg cm^-2 s^-1 (eROSITA-DE DR1 eRASS1 0.2–2.3 keV upper limit at the corridor position)'
    else:
        note += 'Eastern (non-DE) sky: eROSITA upper limits no_public_data.'
        status, lk, lv, lu = 'ledger_only', None, None, None
    if c['target'] == 'teegarden':
        note += ' The teegarden-antipode source at 1.0′ (LSXPS J145306.5−165245 = 3eRASS J145306.7−165249) is fixed_sky (LSXPS span 2005–2025 and the compact eRASS:3 stack).'
    if c['target'] in ('gj-687', 'gj-1221'):
        note += ' LMC direction.'
    if c['target'].startswith('luhman16'):
        note += ' Unresolved sources are CSC-only detections.'
    rows.append(row(target_id=c['target'], channel='corridor', rung='550-10000 AU',
                    band='0.2–2.3 keV (eROSITA UL) + X-ray/γ-ray catalogue cones',
                    substrate='5XMM-DR15 / CSC 2.1 / 2SXPS+LSXPS / eRASS1+eRASS:3 / BAT-157m catalogue cones (7′) at the anti-star corridor position; eROSITA-DE UL service',
                    status=status, limit_kind=lk, limit_value=lv, limit_unit=lu, n_events=c['n_hits'],
                    source=f'{R} §2 corridor catalogue screen; {CSM}; {CSJ}', notes=note))
dump('highenergy_crossings', {
    'survey_id': 'highenergy_crossings', 'report': R, 'plan_section': '5.26', 'pipeline': 'B',
    'archives': ['fermi-lat', 'swift-bat', 'swift-xrt', 'xmm', 'erosita', 'chandra'],
    'hypothesis_version': 'v1.0 + amendments v1.1, v1.2',
    'era': '2008–2026 (LAT era 2008-08-04 → 2026-09-07 UTC, weekly files w009–w953; hypotheses.md D2)',
    'run_dir': 'runs/highenergy-crossings/v1', 'records_dir': None, 'candidates': 0,
    'candidate_notes': 'Blind confirmatory 29 units / 145 frozen trials (136 calibrated, T = 3.451): 9 exceedances, all on the three B van-maanen units '
                       '(lane L + burst) whose aperture contains 3C 279 — vetoed_known_source, already constraint_only at the freeze (amendment v1.2); '
                       '136 calibrated trials 0 exceedances vs 0.05 expected (lane L 10,778 photons vs 10,516.5 expected; lane H 212 vs 198.2). '
                       'No retained-ambiguous exceedance. BAT arm: two A gj-1276 2.5 R☉ windows at combined SNR 3.1–3.2 recorded, not adjudicated (no threshold in that arm).',
    'rows': rows})
