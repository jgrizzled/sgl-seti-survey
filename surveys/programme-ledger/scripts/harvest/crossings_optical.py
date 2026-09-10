"""Generate harvest YAMLs for the five Pipeline-B reports.

Every number is transcribed from the report text or from the results
files the report cites; nothing is estimated."""
import json, collections, math
import numpy as np, yaml
from astropy.table import Table

ROOT = '/home/jgreene/code/sgl-seti-survey/'
import sys
OUT = (sys.argv[1].rstrip('/') + '/') if len(sys.argv) > 1 else ROOT + 'surveys/programme-ledger/harvest/'

RUNG = {0.0056: "1.2 Rsun", 0.0116: "2.5 Rsun", 0.1: "0.1 AU", 1.0: "1.0 AU"}


def rung_of(x):
    for k, v in RUNG.items():
        if abs(float(x) - k) < 1e-6:
            return v
    raise ValueError(x)


def r2(x):
    if x is None:
        return None
    x = float(x)
    if math.isnan(x):
        return None
    return round(x, 2)


def row(target, channel, rung, band, substrate, status, limit_kind=None, limit_value=None,
        limit_unit=None, limit_range=None, power_mw=None, power_mw_range=None, n_events=None,
        n_trials=None, epoch_range=None, source=None, notes=None):
    return dict(target_id=target, channel=channel, rung=rung, band=band, substrate=substrate,
                status=status, limit_kind=limit_kind, limit_value=limit_value,
                limit_unit=limit_unit, limit_range=limit_range, power_mw=power_mw,
                power_mw_range=power_mw_range, n_events=n_events, n_trials=n_trials,
                epoch_range=epoch_range, source=source, notes=notes)


def rng(vals, nd=2):
    vals = [float(v) for v in vals if v is not None and not (isinstance(v, float) and math.isnan(v))]
    if not vals:
        return None
    lo, hi = min(vals), max(vals)
    if abs(lo - hi) < 1e-9:
        return None
    return f"{round(lo, nd)}–{round(hi, nd)}"


def m90stats(m90s, censored_flags=None, cens_value=22.0):
    vals = [float(v) for v in m90s if v is not None and not (isinstance(v, float) and math.isnan(v))]
    if not vals:
        return None, None, None, 0
    med = float(np.median(vals)); best = max(vals)
    if censored_flags is None:
        ncens = sum(1 for v in vals if abs(v - cens_value) < 1e-6)
    else:
        ncens = int(sum(bool(c) for c in censored_flags))
    return r2(med), r2(best), rng(vals), ncens


def dump(name, doc):
    path = OUT + name
    with open(path, 'w') as f:
        yaml.safe_dump(doc, f, sort_keys=False, allow_unicode=True, width=100)
    n = len(doc['rows'])
    print(name, 'rows', n)


# ---------------------------------------------------------------- ZTF
def gen_ztf():
    R = 'report/ztf_crossings.md'
    conf = json.load(open(ROOT + 'surveys/ztf-crossings/results/confirmatory_v1.json'))
    comp = json.load(open(ROOT + 'surveys/ztf-crossings/results/completeness_v1.json'))
    cov = Table.read(ROOT + 'surveys/ztf-crossings/results/coverage_v1_events.ecsv')
    tca = {str(r['event_id']): float(r['t_ca_mjd']) for r in cov}
    compB = {(r['target_id'], r['event_id'], r['band'], r['radius_au']): r for r in comp['B']}
    compA = {(r['target_id'], r['band'], r['radius_au']): r for r in comp['A']}
    SRC_B = ("surveys/ztf-crossings/results/confirmatory_v1.json (B) + "
             "surveys/ztf-crossings/results/completeness_v1.json (B m90)")
    SRC_A = ("surveys/ztf-crossings/results/confirmatory_v1.json (A) + "
             "surveys/ztf-crossings/results/completeness_v1.json (A m90)")
    SUB_B = "ZTF difference-image forced photometry along the antipode relay track (z family 550–10000 AU)"
    SUB_A = "ZTF difference-image forced photometry at the propagated star position (blended)"
    rows = []

    # --- headline aggregate rows (report §2 table + §3) ---
    rows.append(row('programme', 'B', '0.1 AU', 'g/r', SUB_B, 'searched', 'm90', 21.8, 'AB mag', None,
                    0.00013, None, 31, 33, None, f'{R} §2 depth table; §3',
                    'report headline aggregate (per-target rows below enumerate the same cells): 31 searchable '
                    'units with m90 (14 grid-censored >=22); 33 searchable B units in total on 5 targets; '
                    '~130 W = 0.00013 MW transmitter power at 90% for a 550 AU relay through the 2.5 Rsun '
                    'cone (gain ~9e9), 532 nm line in g, EIRP 1.1e12 W at 549 AU; per-event statement only'))
    rows.append(row('programme', 'B', '0.1 AU', 'g/r', SUB_B, 'structurally_open', 'm90', 21.2, 'AB mag', None,
                    None, None, 8, None, None, f'{R} §2 depth table',
                    'single_epoch class (report headline aggregate; 8 units with m90); depth tabulated, no trial'))
    rows.append(row(['gj-1276', 'teegarden'], 'B', '2.5 Rsun', 'g/r', SUB_B, 'structurally_open', 'm90', 21.6,
                    'AB mag', None, 0.00015, None, 9, None, None, f'{R} §2 depth table; §3',
                    'single_epoch class (report headline aggregate, 9 units); grazing single-epoch median 21.6 '
                    '-> ~150 W = 0.00015 MW through the 2.5 Rsun cone; 5 grazing windows on 2 targets covered'))
    rows.append(row('gj-1276', 'B', '1.2 Rsun', 'g/r', SUB_B, 'structurally_open', 'm90', 21.9, 'AB mag', None,
                    None, None, 2, None, None, f'{R} §2 depth table',
                    'single_epoch class (report headline aggregate, 2 units)'))
    rows.append(row(['gj-1276', 'wolf-359', 'wolf-1069'], 'A', 'any', 'g', SUB_A, 'searched', 'm90',
                    16.6, 'AB mag', None, None, None, 3, 3, None, f'{R} §2 depth table',
                    'report headline aggregate (rungs 0.1 AU: gj-1276, wolf-359; 1.0 AU: wolf-1069): 3 searchable A units, '
                    'contrast-limited, median m90 16.6 (g); per-target rows below'))
    rows.append(row('programme', 'A', 'any', 'g/r', SUB_A, 'constraint_only', 'm90', None, 'AB mag',
                    '14.4–18.4', None, None, 18, None, None, f'{R} §2 depth table',
                    'report headline aggregate (rungs 0.1 AU and 1.0 AU): 18 constraint-only A units with parallax-factor template; '
                    'per-target rows below'))
    rows.append(row('wise-0855', 'A', '1.0 AU', 'g/r', SUB_A, 'constraint_only', 'm90', 22.0, 'AB mag', None,
                    0.0002, None, 2, None, None, f'{R} §2 depth table; §3',
                    'empty-field limit (m90 >= 22, grid-censored); 10-m-class uplink transmitter >~200 W = '
                    '0.0002 MW at 2.28 pc; deepest uplink limit in the survey'))

    # --- per-target B rows ---
    groups = collections.OrderedDict()
    for r in conf['B']:
        k = (r['target_id'], r['event_id'], r['band'], r['radius_au'])
        c = compB.get(k)
        st = r['status']
        if st == 'searchable' and r.get('exceedance'):
            cls = 'exceedance_adjudicated'
        elif st == 'searchable':
            cls = 'searched'
        elif st == 'single_epoch':
            cls = 'single_epoch' if r['n_epochs'] > 0 else 'single_epoch_0ep'
        elif st == 'unit' and r['n_epochs'] == 0:
            cls = 'no_usable_data'
        else:
            cls = 'unit_invalid'
        groups.setdefault((r['target_id'], r['radius_au'], r['band'], cls), []).append((r, c))
    for (tgt, rad, band, cls), items in groups.items():
        m90s = [c['m90'] for _, c in items if c is not None]
        med, best, mrange, ncens = m90stats(m90s)
        Ss = [x['S'] for x, _ in items]; Ts = [x['T'] for x, _ in items]
        evs = [x['event_id'] for x, _ in items]
        ts = [tca[e] for e in evs if e in tca]
        er = f"MJD {min(ts):.1f}–{max(ts):.1f}" if len(ts) > 1 else (f"MJD {ts[0]:.1f}" if ts else None)
        neps = sorted(set(x['n_epochs'] for x, _ in items))
        base = dict(target=tgt, channel='B', rung=rung_of(rad), band=band, substrate=SUB_B, source=SRC_B,
                    n_events=len(items), epoch_range=er)
        stat = f"S {rng(Ss, 3) or (r2(Ss[0]) if Ss[0] is not None else 'null')}; T {rng(Ts, 3) or (r2(Ts[0]) if Ts[0] is not None else 'null')}; n_epochs {neps}"
        if cls == 'searched':
            rows.append(row(**base, status='searched', limit_kind='m90' if med is not None else None,
                            limit_value=med, limit_unit='AB mag' if med is not None else None, limit_range=mrange,
                            n_trials=len(items),
                            notes=f"searchable units, null; {stat}; best m90 {best}; {ncens} grid-censored at 22.0"
                                  + ("; m90 not tabulated in completeness_v1.json" if med is None else "")))
        elif cls == 'exceedance_adjudicated':
            adj = {'evt-1ff54c0eec27': 'amplitude 1.3 sigma; no recurrence (4 other events, max R = 0.48) -> control crossing',
                   'evt-57e0c494a587': 'sub-1 sigma amplitude; exceeds only via freakishly low T -> control underdispersion',
                   'evt-9c6a08c6df77': '2.2 sigma; no amplitude-consistent recurrence (other 4 events R <= 0.77) -> control crossing'}
            notes = '; '.join(f"{x['event_id']}: S {x['S']} vs T {x['T']}, {adj.get(x['event_id'], '')}" for x, _ in items)
            rows.append(row(**base, status='searched', limit_kind='m90', limit_value=med, limit_unit='AB mag',
                            limit_range=mrange, n_trials=len(items),
                            notes=f"exceedance_adjudicated (control-crossing statistics, vetoed, not candidates; "
                                  f"surveys/ztf-crossings/results/confirmatory_v1.md): {notes}"))
        elif cls == 'single_epoch':
            rows.append(row(**base, status='structurally_open', limit_kind='m90' if med is not None else None,
                            limit_value=med, limit_unit='AB mag' if med is not None else None, limit_range=mrange,
                            notes=f"single_epoch class (constraint-only statement, no trial); {stat}; "
                                  f"{ncens} grid-censored at 22.0"
                                  + ("; m90 not tabulated" if med is None else "")))
        elif cls == 'single_epoch_0ep':
            rows.append(row(**base, status='structurally_open',
                            notes="single_epoch class with 0 usable epochs (no statistic)"))
        elif cls == 'no_usable_data':
            rows.append(row(**base, status='structurally_open',
                            notes="no_usable_data: unit with 0 usable epochs (status 'unit', n_epochs 0)"))
        else:
            x = items[0][0]
            rows.append(row(**base, status='structurally_open', limit_kind='m90' if med is not None else None,
                            limit_value=med, limit_unit='AB mag' if med is not None else None,
                            notes=f"status 'unit' but not searchable under the freeze: {x['event_id']} n_epochs {x['n_epochs']}, "
                                  f"{x['n_controls']}/8 ring controls, S {x['S']} vs T {x['T']}; "
                                  + ("reported as controls_incomplete (5/8 ring controls), anomaly retained, not searched" if x['n_controls'] == 5 else "not searched")))

    # --- per-target A rows ---
    power_A = {('teegarden', 0.1, 'g'): (0.056, '~56 kW = 0.056 MW, 10-m diffraction-limited uplink transmitter (constraint-only)'),
               ('wolf-359', 0.1, 'g'): (0.066, '~66 kW = 0.066 MW, 10-m-class uplink transmitter'),
               ('gj-1276', 0.1, 'g'): (0.4, '~400 kW = 0.4 MW, 10-m-class uplink transmitter')}
    for r in conf['A']:
        k = (r['target_id'], r['band'], r['radius_au'])
        c = compA.get(k)
        m90 = r2(c['m90']) if c and c.get('m90') is not None else None
        pw = power_A.get((r['target_id'], r['radius_au'], r['band']))
        status = 'searched' if r['status'] == 'searchable' else 'constraint_only'
        notes = (f"{r['kind']}; S {r['S']} vs T {r['T']}; valid offsets {r.get('n_valid_offsets')}/8; "
                 f"windows with data {r.get('n_windows_with_data')}; k {r['gate'].get('k_scale')}")
        if r['kind'] == 'single_window':
            notes += '; single-window unit (constraint-only by construction; raw unstandardized S where no template)'
        if pw:
            notes += '; ' + pw[1] + ' (1-m-class x100)'
        if k[0] == 'wise-0855' and r['radius_au'] == 1.0 and r['band'] in ('g', 'r'):
            notes += '; empty-field limit, m90 >= 22 grid-censored'
        rows.append(row(r['target_id'], 'A', rung_of(r['radius_au']), r['band'], SUB_A, status,
                        'm90' if m90 is not None else None, m90, 'AB mag' if m90 is not None else None, None,
                        pw[0] if pw else None, None, r.get('n_windows_with_data'),
                        1 if status == 'searched' else None, None, SRC_A, notes))

    # --- coverage / open-cell census (report §1 item 3 + coverage_v1_summary.json) ---
    CS = 'surveys/ztf-crossings/results/coverage_v1_summary.json; ' + R + ' §1 item 3'
    rows.append(row('programme', 'A', '1.0 AU', 'any', SUB_A, 'not_constrainable', n_events=318, source=CS,
                    notes='402/720 events covered; all 29 uncovered targets have dec < -30 deg (outside the ZTF footprint)'))
    rows.append(row('programme', 'A', '0.1 AU', 'any', SUB_A, 'structurally_open', n_events=15, source=CS,
                    notes='44/59 events covered (coverage_v1_summary.json); uncovered remainder'))
    rows.append(row('programme', 'B', '1.2 Rsun', 'any', SUB_B, 'structurally_open', n_events=34, source=CS,
                    notes='1/35 events covered (grazing rung epoch-starved as the cadence predicts)'))
    rows.append(row('programme', 'B', '2.5 Rsun', 'any', SUB_B, 'structurally_open', n_events=38, source=CS,
                    notes='5/43 events covered (gj-1276 + teegarden only); grazing rung epoch-starved'))
    rows.append(row('programme', 'B', '0.1 AU', 'any', SUB_B, 'structurally_open', n_events=27, source=CS,
                    notes='33/60 events covered; uncovered remainder'))
    rows.append(row('programme', 'A', 'any', 'g/r/i', SUB_A, 'structurally_open', source=R + ' §1 item 4',
                    notes='channel-A saturation cut excludes targets < 13.0 mag; searchable populations g 27 / r 20 / i 14 '
                          'targets (M dwarfs, white dwarfs, late-T/Y); brighter targets saturated (surveys/ztf-crossings/results/saturation_cut_v1.*)'))
    rows.append(row('programme', 'A', '1.0 AU', 'g/r', SUB_A, 'constraint_only', source=R + ' §1 item 6; §4 lesson 1',
                    notes='1.0 AU rung constraint-only by geometry: windows ~ observing season, temporal pseudo-window controls fail by construction'))
    rows.append(row('programme', ['A', 'B'], 'any', '1064/1550 nm', 'ZTF g/r/i', 'not_constrainable', source=R + ' §3',
                    notes='1064/1550 nm links outside the ZTF response; 400–900 nm only'))
    rows.append(row('programme', 'S1', 'any', 'any', None, 'no_survey', source=R + ' §1 item 2; §3',
                    notes='post-lens downlink arrives sunward; sunward-arriving combinations declared out of scope at freeze v1.0'))
    rows.append(row('programme', 'S2', 'any', 'any', None, 'no_survey', source=R + ' §1 item 2; §3',
                    notes='uplink past the Sun (sunward); declared out of scope at freeze v1.0'))
    rows.append(row('programme', 'B', 'any', 'any', SUB_B, 'structurally_open', source=R + ' §3',
                    notes='not constrained: transmitters scheduled to avoid Earth crossings; pulse periods between 30 s and the '
                          'window length; the anti-target downlink outside covered windows'))

    doc = dict(survey_id='ztf_crossings', report=R, plan_section="5.1", pipeline='B', archives=['ztf'],
               hypothesis_version="v1.0 (threshold freeze v1.0 + amendments v1.1, v1.2)",
               era="MJD 58178–61275 (2018–2026)", run_dir='runs/ztf-crossings', records_dir=None,
               candidates=0,
               candidate_notes='3 exceedances vs 4.0 expected control crossings (gj-908 g x1, gj-1276 g x2), all adjudicated as '
                               'control-crossing statistics, none promoted; 772 archive-404 epochs counted',
               rows=rows)
    dump('ztf_crossings.yaml', doc)


# ---------------------------------------------------------------- PS1
def gen_ps1():
    R = 'report/ps1_crossings.md'
    conf = json.load(open(ROOT + 'surveys/ps1-crossings/results/confirmatory_v1.json'))
    dev = json.load(open(ROOT + 'surveys/ps1-crossings/results/dev_search_v1.json'))
    for r in dev['B']:
        r['split'] = 'dev'
        r.setdefault('radius_au', 0.1)  # dev split is the B wide rung (report §1 item 5: 'B-wide teegarden + wolf-359')
    conf['B'] = conf['B'] + dev['B']
    comp = json.load(open(ROOT + 'surveys/ps1-crossings/results/completeness_v1.json'))
    compB = {(r['target_id'], r['event_id'], r['band'], r['radius_au']): r for r in comp['B']}
    compA = {(r['target_id'], r['band']): r for r in comp['A']}
    SRC_B = ("surveys/ps1-crossings/results/confirmatory_v1.json + dev_search_v1.json (B) + "
             "surveys/ps1-crossings/results/completeness_v1.json (B m90)")
    SRC_A = ("surveys/ps1-crossings/results/confirmatory_v1.json (A) + "
             "surveys/ps1-crossings/results/completeness_v1.json (A m90)")
    SUB_B = "PS1 DR2 single-epoch warps, star-calibrated (warp-direct, no differencing), antipode relay track"
    SUB_A = "PS1 DR2 single-epoch warps, star-calibrated, at the propagated star position (blended)"
    rows = []
    # headline aggregates
    rows.append(row('programme', 'B', '0.1 AU', 'g/r/i', SUB_B, 'searched', 'm90', 21.6, 'AB mag', None,
                    0.00017, None, 11, 11, '2010–2014', f'{R} §2 depth table; §3',
                    'report headline aggregate (per-target rows below): 11 searched B units, all bands, median m90 21.6 '
                    '(4 grid-censored >= 22); deepest units ~170 W = 0.00017 MW (all-band median) through the 2.5 Rsun cone '
                    'at 550 AU; 9 searched events on 6 targets; per-event statement'))
    rows.append(row('programme', 'B', '0.1 AU', 'g', SUB_B, 'searched', 'm90', 21.1, 'AB mag', '19.2–22.0',
                    0.00028, None, 4, 4, None, f'{R} §2 depth table; §3',
                    'g units: 19.2 / 20.6 / 21.6 / >=22 (median 21.1); 532 nm line hypothesis; ~280 W = 0.00028 MW '
                    '(EIRP 2.5e12 W at 549 AU)'))
    rows.append(row('programme', 'B', '0.1 AU', 'r', SUB_B, 'searched', 'm90', 21.1, 'AB mag', '20.5–21.4',
                    None, None, 3, 3, None, f'{R} §2 depth table', 'r units: 20.5 / 21.1 / 21.4'))
    rows.append(row('programme', 'B', '0.1 AU', 'i', SUB_B, 'searched', 'm90', 22.0, 'AB mag', '21.7–22.0',
                    0.000037, None, 4, 4, None, f'{R} §2 depth table; §3',
                    'i units: 21.7 / >=22 / >=22 / >=22 (3 grid-censored); grid-censored i units <~ 37 W = 0.000037 MW '
                    '(752 nm leakage)'))
    rows.append(row(['gj-1276', 'teegarden'], 'A', '0.1 AU', 'g/i', SUB_A, 'constraint_only', 'm90', None, 'AB mag',
                    '17.9–21.3', None, '0.001–0.032', 3, None, None, f'{R} §2 depth table; §3',
                    'reference depths, threshold = max(T,0), no calibrated false-alarm control: gj-1276 g 17.9 / i 18.0; '
                    'teegarden g 21.3; 10-m-class 532 nm uplink transmitter >~ ~1 kW = 0.001 MW (teegarden), '
                    '~23–32 kW = 0.023–0.032 MW (gj-1276); qualified sensitivity statements, not calibrated exclusions'))

    # per-target B rows
    groups = collections.OrderedDict()
    for r in conf['B']:
        k = (r['target_id'], r['event_id'], r['band'], r['radius_au'])
        c = compB.get(k)
        st = r['status']
        if st == 'controls_partial_7of8':
            st = 'ok'
        if st == 'ok':
            if r['event_id'] == 'evt-3f7ad50e70e1':
                cls = 'exceedance_vetoed'
            elif r['event_id'] == 'evt-344c12d32a8d':
                cls = 'retained_ambiguous'
            else:
                cls = 'searched'
        else:
            cls = st
        groups.setdefault((r['target_id'], r['radius_au'], r['band'], cls), []).append((r, c))
    for (tgt, rad, band, cls), items in groups.items():
        m90s = [c['m90'] for _, c in items if c is not None]
        cens = [c.get('grid_censored') for _, c in items if c is not None]
        med, best, mrange, ncens = m90stats(m90s, cens)
        ts = [x['t_ca_mjd'] for x, _ in items]
        er = f"MJD {min(ts):.1f}–{max(ts):.1f}" if len(ts) > 1 else f"MJD {ts[0]:.1f}"
        evs = ', '.join(x['event_id'] for x, _ in items)
        zdead = sorted({z for _, c in items if c for z in c.get('z_dead', [])})
        base = dict(target=tgt, channel='B', rung=rung_of(rad), band=band, substrate=SUB_B, source=SRC_B,
                    n_events=len(items), epoch_range=er)
        if cls == 'searched':
            Ss = [x['S'] for x, _ in items]; Ts = [x['T'] for x, _ in items]
            n = f"searched units, null: {evs}; S {rng(Ss,3) or Ss[0]}; T {rng(Ts,3) or Ts[0]}; {ncens} grid-censored at 22.0"
            if zdead:
                n += f"; z = {zdead} node dead (masked in every usable epoch), excluded from the depth statement"
            if any(x.get('split') == 'dev' for x, _ in items):
                n += '; dev split (dev_search_v1.json)'
            if any(x['status'] == 'controls_partial_7of8' for x, _ in items):
                n += '; one unit controls_partial_7of8 (recorded as anomaly P2)'
            rows.append(row(**base, status='searched', limit_kind='m90', limit_value=med, limit_unit='AB mag',
                            limit_range=mrange, n_trials=len(items), notes=n))
        elif cls == 'exceedance_vetoed':
            x = items[0][0]
            rows.append(row(**base, status='vetoed_known_source', limit_kind='m90', limit_value=med, limit_unit='AB mag',
                            n_trials=1,
                            notes=f"{x['event_id']}: S {x['S']} vs T {x['T']}; vetoed — flux-consistent catalogued static "
                                  "(DR2 stack sources 1.46 arcsec r 22.5–22.8 and 1.99 arcsec r 23.5 from the track nodes; "
                                  "node fluxes r ~21.6–22.1); SkyBoT clean; TTI-pair excludes movers"))
        elif cls == 'retained_ambiguous':
            x = items[0][0]
            rows.append(row(**base, status='retained_ambiguous', limit_kind='m90', limit_value=med, limit_unit='AB mag',
                            n_trials=1,
                            notes=f"{x['event_id']}: S {x['S']} vs T {x['T']}; night-consistent i ~22.7–23.4 signal "
                                  "(1.1–2.3 sigma per epoch, 3 usable exposures over ~35 min) at the z = 550 node; census clean, "
                                  "TTI-consistent, no static counterpart to ~23+; non-promotable; recurrence test designated "
                                  "to the joint stage (see joint_crossings.yaml); m90 grid-censored >= 22"))
        elif cls == 'track_masked':
            x = items[0][0]
            extra = ''
            if tgt == 'van-maanen':
                extra = '; van-maanen b = 0.28 Rsun graze (2010-04-03, TTI i-band pair) — deepest graze in the programme, mask-lost'
            if any(x.get('split') == 'dev' for x, _ in items):
                extra += '; dev split (dev_search_v1.json)'
            rows.append(row(**base, status='structurally_open',
                            notes=f"track_masked (nominal-covered / mask-unusable): {evs}; every in-window sample lost to CONV.BAD "
                                  f"chip-gap bands; n_in_epochs {[x['n_in_epochs_raw'] for x,_ in items]}{extra}"))
        else:
            rows.append(row(**base, status='structurally_open', notes=f"status {cls}: {evs}"))

    # per-target A rows
    for r in conf['A']:
        c = compA.get((r['target_id'], r['band']))
        m90 = r2(c['m90']) if c and c.get('m90') is not None else None
        pw = None; pn = ''
        if r['target_id'] == 'teegarden' and r['band'] == 'g':
            pw = 0.001; pn = '; ~1 kW = 0.001 MW 10-m-class 532 nm uplink transmitter (reference sensitivity, threshold 0)'
        if r['target_id'] == 'gj-1276' and r['band'] in ('g', 'i'):
            pn = '; gj-1276 reference ~23–32 kW = 0.023–0.032 MW (report §3, g/i together)'
        g = r['gates']
        rows.append(row(r['target_id'], 'A', '0.1 AU', r['band'], SUB_A, 'constraint_only',
                        'm90' if m90 is not None else None, m90, 'AB mag' if m90 is not None else None, None,
                        pw, '0.023–0.032' if (r['target_id'] == 'gj-1276' and r['band'] in ('g', 'i')) else None,
                        g.get('n_windows_with_data'), None, None, SRC_A,
                        f"{r['class']}; reference_constraint_only (threshold max(T,0), no calibrated control): S {r['S']}, "
                        f"T {r['T']}; windows with data {g.get('n_windows_with_data')}, valid offsets {g.get('n_valid_offsets')}/8"
                        + ('; completeness insufficient_data' if c and c.get('status') == 'insufficient_data' else '') + pn))

    # coverage / open-cell census
    CS = 'surveys/ps1-crossings/results/coverage_v1_summary.json; ' + R + ' §1 item 3'
    rows.append(row('programme', 'B', '0.1 AU', 'any', SUB_B, 'structurally_open', n_events=33, source=CS,
                    notes='14/47 wide-rung events covered on all seven narrow-rung targets; uncovered remainder (report §3 '
                          'lists "the 5 remaining uncovered wide-rung events" under not constrained)'))
    rows.append(row('programme', 'B', '1.2 Rsun', 'any', SUB_B, 'structurally_open', n_events=26, source=CS,
                    notes='1/27 events covered (epoch-starved); the one covered event (van-maanen b = 0.28 Rsun) is mask-unusable'))
    rows.append(row('programme', 'B', '2.5 Rsun', 'any', SUB_B, 'structurally_open', n_events=33, source=CS,
                    notes='1/34 events covered (epoch-starved); the one covered event (van-maanen) is mask-unusable'))
    rows.append(row('programme', 'A', '0.1 AU', 'any', SUB_A, 'structurally_open', n_events=28, source=CS,
                    notes='18/46 events covered; uncovered remainder'))
    rows.append(row('programme', 'A', '1.0 AU', 'any', SUB_A, 'ledger_only', n_events=294, source=CS,
                    notes='294/563 events covered; rung declared constraint-only at freeze (window ~ observing season, ZTF lesson 1); '
                          'coverage recorded, per-window depths deferred, zero trials'))
    rows.append(row('programme', 'A', '1.0 AU', 'any', SUB_A, 'structurally_open', n_events=269, source=CS,
                    notes='uncovered remainder of the 1.0 AU rung (563 - 294)'))
    rows.append(row('programme', 'A', '0.1 AU', 'grizy', SUB_A, 'structurally_open', source=R + ' §1 item 4',
                    notes='channel-A saturation cut: survivors gj-1276 (all bands), teegarden (g + marginal r), van-maanen (y); '
                          'all other targets saturated in PS1 warps (surveys/ps1-crossings/results/saturation_cut_v1.*)'))
    rows.append(row('programme', 'B', 'any', 'any', SUB_B, 'structurally_open', n_events=7, source=R + ' §2',
                    notes='mask attrition: 7 of 18 covered-window unit-rows lost every in-window sample to CONV.BAD chip-gap bands '
                          '(~40% unit-level attrition); enumerated per target above'))
    rows.append(row('programme', ['A', 'B'], 'any', '1064/1550 nm', 'PS1 grizy', 'not_constrainable', source=R + ' §1 item 2; §3',
                    notes='1064/1550 nm outside grizy; unconstrained'))
    rows.append(row('programme', 'S1', 'any', 'any', None, 'no_survey', source='surveys/ps1-crossings/hypotheses.md §1',
                    notes='sunward-arriving (post-lens downlink) declared out of scope at freeze'))
    rows.append(row('programme', 'S2', 'any', 'any', None, 'no_survey', source='surveys/ps1-crossings/hypotheses.md §1',
                    notes='sunward-arriving (uplink past the Sun) declared out of scope at freeze'))
    rows.append(row('programme', 'B', 'any', 'any', SUB_B, 'structurally_open', source=R + ' §3',
                    notes='not constrained: transmitters scheduled to avoid Earth crossings; pulse periods between ~30 s and the window length'))

    doc = dict(survey_id='ps1_crossings', report=R, plan_section="5.2", pipeline='B', archives=['ps1'],
               hypothesis_version="v1.0 (ZTF v1.1/v1.2 amendments adopted at freeze; no dev amendment)",
               era="MJD 54900–57300 (warps span 54985–57067; 2009–2015)", run_dir='runs/ps1-crossings', records_dir=None,
               candidates=0,
               candidate_notes='2 exceedances vs 2.1 expected over the frozen 19-unit family: ross-128 r evt-3f7ad50e70 vetoed '
                               '(stack-catalogued static); gj-1276 i evt-344c12d32a retained-ambiguous, non-promotable, '
                               'handed to the joint recurrence test',
               rows=rows)
    dump('ps1_crossings.yaml', doc)


# ---------------------------------------------------------------- WISE
def gen_wise():
    R = 'report/wise_crossings.md'
    frz = json.load(open(ROOT + 'surveys/wise-crossings/configs/threshold_freeze_v1.json'))
    gate = json.load(open(ROOT + 'surveys/wise-crossings/configs/elongation_gate_v1.json'))
    cov = Table.read(ROOT + 'surveys/wise-crossings/results/coverage_v1_events.ecsv')
    sat = Table.read(ROOT + 'surveys/wise-crossings/results/saturation_cut_v1.ecsv')
    satst = {str(r['target_id']): {b: str(r[f'status_{b}']) for b in ('W1', 'W2', 'W3', 'W4')} for r in sat}
    satest = {str(r['target_id']): {b: float(r[f'{b}_est']) for b in ('W1', 'W2', 'W3', 'W4')} for r in sat}
    SUB = "WISE/NEOWISE L1b single exposures (no image searched; TAP frame inventory only)"
    rows = []
    GS = 'surveys/wise-crossings/configs/elongation_gate_v1.json; ' + R + ' §2 table'
    for key, rung, ch in (('B_0.0056', '1.2 Rsun', 'B'), ('B_0.0116', '2.5 Rsun', 'B'), ('B_0.1', '0.1 AU', 'B'),
                          ('A_0.1', '0.1 AU', 'A')):
        g = gate['rungs'][key]
        rows.append(row('programme', ch, rung, 'W1/W2/W3/W4', SUB, 'not_constrainable', n_events=g['n_events'], source=GS,
                        notes=f"Gate I (elongation theorem): {g['n_events']} in-era events, {g['n_viable']} observable in-window "
                              f"(WISE observes only at solar elongation ~90 deg; beam source sits near the observer–Sun axis); "
                              f"geometric null, carries no ledger rows"))
    g = gate['rungs']['A_1.0']
    rows.append(row('programme', 'A', '1.0 AU', 'W1/W2', SUB, 'not_constrainable', n_events=g['n_events'], n_trials=0, source=GS + '; §3',
                    notes=f"Gate I passes {g['n_viable']} of {g['n_events']} in-era events (45 targets) — but Gate II (control-geometry "
                          "theorem) leaves 0 searchable units: ~116 d windows recurring annually cannot be cleared by the 8 designated "
                          "temporal pseudo-window controls (v1.0 pool ±23–97 d: 0 units; v1.1 pool ±130–245 d: 2–6 valid offsets, still 0); "
                          "0 trials, 0 exceedance budget, no pixel searched"))
    # 30 unit-rows from the freeze
    units = {}
    for cls in ('single_window', 'gate_blocked'):
        for u in frz['search_units'][cls]:
            units[(u['target_id'], u['band'])] = (cls, u)
    # per-target coverage
    per = collections.OrderedDict()
    for r in cov:
        if int(r['n_W1']) == 0 and int(r['n_W2']) == 0 and int(r['n_W3']) == 0 and int(r['n_W4']) == 0:
            continue
        for b in ('W1', 'W2', 'W3', 'W4'):
            if int(r[f'n_{b}']) > 0:
                per.setdefault((str(r['target_id']), b), []).append(r)
    ntargets = len({k[0] for k in per})
    nev = len({str(r['event_id']) for rs in per.values() for r in rs})
    print('wise covered events', nev, 'targets', ntargets)
    SRC_U = ('surveys/wise-crossings/configs/threshold_freeze_v1.json search_units + '
             'surveys/wise-crossings/results/coverage_v1_events.ecsv')
    SRC_C = ('surveys/wise-crossings/results/coverage_v1_events.ecsv + '
             'surveys/wise-crossings/results/saturation_cut_v1.ecsv')
    for (tgt, b), rs in sorted(per.items()):
        ts = [float(r['t_ca_mjd']) for r in rs]
        er = f"MJD {min(ts):.1f}–{max(ts):.1f}" if len(ts) > 1 else f"MJD {ts[0]:.1f}"
        neps = [int(r[f'n_{b}']) for r in rs]
        u = units.get((tgt, b))
        s = satst.get(tgt, {}).get(b)
        if u:
            cls, uu = u
            rows.append(row(tgt, 'A', '1.0 AU', b, SUB, 'not_constrainable', n_events=len(rs), n_trials=0, epoch_range=er,
                            source=SRC_U,
                            notes=f"Gate II unit class {cls} (report: constraint-only, never a candidate): windows {uu['n_windows']}, "
                                  f"off-window epochs {uu['n_off_epochs']}, valid offsets {uu['n_valid_offsets']}/8; in-window {b} "
                                  f"epochs per event {neps}; coverage-without-statistic entry in the joint ledger"))
        else:
            rows.append(row(tgt, 'A', '1.0 AU', b, SUB, 'ledger_only', n_events=len(rs), n_trials=0, epoch_range=er,
                            source=SRC_C,
                            notes=f"covered window(s), in-window {b} epochs per event {neps}; not a freeze unit — saturation cut "
                                  f"status_{b} = {s} ({b}_est {satest.get(tgt, {}).get(b)} mag); coverage-without-statistic "
                                  f"entry in the joint ledger"))
    rows.append(row('programme', 'A', '1.0 AU', 'W1/W2', SUB, 'ledger_only', n_events=91, source=R + ' §4',
                    notes='covered-window inventory: 91 events on 43 targets with in-window primary W1+W2 frames (W3/W4: 5/4 cryo-era '
                          'events), 211–343-epoch off-window baselines per band; enumerated per target above'))
    rows.append(row('programme', 'A', '1.0 AU', 'W1/W2', SUB, 'structurally_open', source=R + ' §1 item 5; §5 lesson 2',
                    notes='saturation cut: W1/W2 searchable population = 11 faint targets (white dwarfs, T/Y dwarfs, latest-M); every '
                          'ordinary M dwarf within 10 pc is brighter than the W1/W2 single-frame saturation limits'))
    rows.append(row('programme', ['A', 'B'], 'any', '3–5 um', SUB, 'not_constrainable', source=R + ' §4',
                    notes='"Not constrained: everything" — no calibrated sensitivity, exclusion or reference depth for any cell; the '
                          'mid-IR crossing regime remains untested'))
    doc = dict(survey_id='wise_crossings', report=R, plan_section="5.3", pipeline='B', archives=['wise'],
               hypothesis_version="v1.0 (+ pre-search amendments: qa_status quality value; threshold freeze v1.1 control-pool rescale)",
               era="2010–2024 (crossings/wise_v1; coverage era MJD 55203–60524)", run_dir='runs/wise-crossings', records_dir=None,
               candidates=0,
               candidate_notes='structural null: 0 searchable units, 0 trials, no signal statistic formed; 30 constraint-only unit-rows '
                               '(18 single-window + 12 gate-blocked)',
               rows=rows)
    dump('wise_crossings.yaml', doc)


# ---------------------------------------------------------------- JOINT
STATUS_MAP = {'searched_null': 'searched', 'coverage_only': 'ledger_only', 'constraint_only': 'constraint_only',
              'single_epoch': 'structurally_open', 'track_masked': 'structurally_open',
              'no_usable_data': 'structurally_open', 'exceedance_vetoed': 'vetoed_known_source',
              'retained_ambiguous': 'retained_ambiguous', 'exceedance_adjudicated': 'searched'}
SUBSTRATE = {'ztf': 'ZTF difference-image forced photometry (covered-window ledger v1)',
             'ps1': 'PS1 DR2 warp-direct star-calibrated photometry (covered-window ledger v1)',
             'wise': 'WISE/NEOWISE L1b frame inventory (covered-window ledger v1; no pixel searched)'}


def gen_joint():
    R = 'report/joint_crossings.md'
    LED = 'surveys/joint-crossings/results/covered_window_ledger_v1.ecsv'
    t = Table.read(ROOT + LED).filled()
    groups = collections.OrderedDict()
    for r in t:
        k = (str(r['archive']), str(r['channel']), float(r['rung_au']), str(r['target_id']), str(r['band']), str(r['status']))
        groups.setdefault(k, []).append(r)
    rows = []
    rec = json.load(open(ROOT + 'surveys/joint-crossings/results/gj1276_recurrence_v1.json'))
    rows.append(row('gj-1276', 'B', '0.1 AU', 'g+r', 'ZTF channel-B stack at z fixed to 550 AU over 11 event-bands (7 events x g/r), '
                    'ring controls paired across event-bands; single joint trial (Arm R1)', 'searched', 'm90', rec['joint_m90'],
                    'AB mag', None, None, None, rec['n_event_bands'], 1, '2019–2026',
                    'surveys/joint-crossings/results/gj1276_recurrence_v1.json; ' + R + ' §2',
                    f"gj-1276 recurrence test of the PS1 retained-ambiguous exceedance (evt-344c12d32a8d, i ~22.8 at z = 550, 2014-03-06): "
                    f"J = {rec['J']}, T_J = {rec['T_J']}, joint exceedance no; joint recovery at 22.8 AB (flat F_nu) = "
                    f"{rec['joint_recovery_at_22p8']}; verdict {rec['verdict']} (flat-SED persistent relay refuted at 90%); amendment v1.1 "
                    f"single-epoch clip |f|/sigma > 20 applied (pre-amendment J = {rec['pre_amendment_v1_0']['J']}, recovery "
                    f"{rec['pre_amendment_v1_0']['joint_recovery_at_22p8']}); anomaly retained, narrowed, not a candidate"))
    rows.append(row('gj-1276', 'B', '0.1 AU', 'i', 'ZTF i-band in-window epochs across all gj-1276 wide-rung windows', 'structurally_open',
                    'm90', 21.1, 'AB mag', None, None, None, None, None, None,
                    'surveys/joint-crossings/results/gj1276_recurrence_v1.json arm_R2_i_band; ' + R + ' §2',
                    'Arm R2 (752 nm line / red SED) untestable: exactly 2 in-window i epochs archive-wide, 2-epoch depth ~21.1 vs '
                    'required 22.8; static arm settled by the PS1 DR2 stack catalog (no counterpart within 6 arcsec to ~23+)'))
    for (arch, ch, rad, tgt, band, st), rs in groups.items():
        m90 = [float(r['m90']) for r in rs]
        cens = [bool(r['m90_censored']) for r in rs]
        med, best, mrange, ncens = m90stats(m90, cens)
        ts = [float(r['t_ca_mjd']) for r in rs if not math.isnan(float(r['t_ca_mjd']))]
        er = None
        if ts:
            er = f"MJD {min(ts):.1f}–{max(ts):.1f}" if len(ts) > 1 and abs(max(ts) - min(ts)) > 0.05 else f"MJD {ts[0]:.1f}"
        Ss = [float(r['S']) for r in rs if not math.isnan(float(r['S']))]
        Ts = [float(r['T']) for r in rs if not math.isnan(float(r['T']))]
        notes_led = sorted({str(r['notes']) for r in rs if str(r['notes']) not in ('', 'nan', 'N/A', '--')})
        neps = sorted({int(r['n_epochs']) for r in rs})
        evs = sorted({str(r['event_id']) for r in rs})
        status = STATUS_MAP[st]
        n = f"ledger status {st}; {len(rs)} ledger row(s), {len(evs)} event(s)"
        if med is not None:
            n += f"; m90 median {med}, best {best}, n with m90 {len([v for v in m90 if not math.isnan(v)])}, {ncens} grid-censored"
        else:
            n += "; no m90 in ledger"
        if Ss:
            n += f"; S {rng(Ss,3) or r2(Ss[0])}"
        if Ts:
            n += f"; T {rng(Ts,3) or r2(Ts[0])}"
        n += f"; n_epochs {neps}"
        if notes_led:
            n += "; ledger notes: " + ' | '.join(notes_led)
        if st == 'exceedance_adjudicated':
            n += "; exceedance_adjudicated (ZTF control-crossing statistics, not candidates)"
        if st == 'exceedance_vetoed':
            n += "; PS1 ross-128 stack-catalogued static"
        if st == 'retained_ambiguous':
            n += "; PS1 gj-1276 evt-344c12d32a8d — resolved by the recurrence test row above (flat-SED refuted at 90%)"
        if st == 'track_masked' and tgt == 'van-maanen':
            n += "; van-maanen b = 0.28 Rsun graze, nominal-covered / mask-unusable"
        if len(evs) <= 3:
            n += "; events " + ', '.join(evs)
        ntr = len(rs) if st in ('searched_null', 'exceedance_adjudicated', 'exceedance_vetoed', 'retained_ambiguous') else (0 if arch == 'wise' else None)
        rows.append(row(tgt, ch, rung_of(rad), band, SUBSTRATE[arch], status,
                        'm90' if med is not None else None, med, 'AB mag' if med is not None else None, mrange,
                        None, None, len(evs), ntr, er, LED + f' (aggregated per archive={arch}, channel, rung, target, band, status)', n))
    # programme tally rows (report §3)
    rows.append(row('programme', 'B', '0.1 AU', 'g/r', 'ZTF (ledger)', 'searched', 'm90', 21.8, 'AB mag', None, 0.00013, None, None, 36, '2018–2026',
                    R + ' §3', 'programme tally: ZTF 36 searched units, B median m90 ~21.8 (relay >~ ~130 W = 0.00013 MW through the 2.5 Rsun cone)'))
    rows.append(row('programme', 'B', '0.1 AU', 'g/r/i', 'PS1 (ledger)', 'searched', 'm90', 21.6, 'AB mag', None, None, '0.00017–0.00028', None, 11, '2010–2014',
                    R + ' §3', 'programme tally: PS1 11 searched units, B median 21.6 (>~ ~170–280 W = 0.00017–0.00028 MW)'))
    rows.append(row('programme', ['A', 'B'], 'any', 'any', 'crossings programme steps 1–6', 'searched', None, None, None, None, None, None, None, 47, None,
                    R + ' §3', '47 searched trials programme-wide (ZTF 36, PS1 11, WISE 0); 5 exceedances vs ~6.1 expected control crossings; 0 candidates'))
    rows.append(row('programme', 'B', 'any', 'any', 'optical archives', 'structurally_open', source=R + ' §3',
                    notes='structurally open: photosphere/coronal grazing rungs (ZTF: 5 windows on 2 targets searched [single-epoch]; PS1: the b = 0.28 Rsun event mask-lost)'))
    rows.append(row('programme', ['A', 'B'], 'any', '1064/1550 nm and line SEDs redward of ~900 nm', 'optical archives', 'not_constrainable', source=R + ' §3',
                    notes='structurally open cell stated for the record'))
    rows.append(row('programme', ['A', 'B'], 'any', 'mid-IR', 'WISE', 'not_constrainable', source=R + ' §3',
                    notes='mid-IR crossings regime closed by the WISE elongation + control-geometry theorems'))
    rows.append(row('programme', 'A', '1.0 AU', 'any', 'all archives', 'constraint_only', source=R + ' §3',
                    notes='the d = 1 wide-beam (1 AU) rung as a controlled search anywhere: three independent proofs it defeats the frozen '
                          'temporal-control family (ZTF, PS1, WISE); designated future route a v2-style null ensemble'))
    rows.append(row('programme', ['A', 'B'], 'any', 'any', 'all archives', 'structurally_open', source=R + ' §3',
                    notes='not constrained: pulse periods between exposure length and window length; transmitters scheduled to avoid Earth crossings'))
    doc = dict(survey_id='joint_crossings', report=R, plan_section="5.4", pipeline='B', archives=['ztf', 'ps1', 'wise'],
               hypothesis_version="joint plan v1.0 + amendment v1.1 (single-epoch clip)",
               era="ZTF 2018–2026, PS1 2009–2015, WISE 2010–2024 (eras disjoint per archive)", run_dir='runs/joint-crossings', records_dir=None,
               candidates=0,
               candidate_notes='PS1 gj-1276 retained-ambiguous exceedance: flat-SED persistent-relay recurrence refuted at 90% (joint recovery 0.94, '
                               'joint m90 22.93); retained in the record, narrowed, not promotable; ledger 593 rows (ZTF 86, PS1 316, WISE 191)',
               rows=rows)
    dump('joint_crossings.yaml', doc)


# ---------------------------------------------------------------- PTF
def gen_ptf():
    R = 'report/ptf_crossings.md'
    SUB_B = "PTF level-1 epochal images, scie-direct (no differencing), per-frame field-star calibration, antipode relay track"
    SUB_A = "PTF level-1 epochal images, scie-direct, at the star position (blended)"
    SRC = R + ' §2 tables (S/T, m90) + §3; surveys/ptf-crossings/results/confirmatory_v1.json, completeness_v1.json'
    rows = []
    covt = Table.read(ROOT + 'surveys/ptf-crossings/results/coverage_v1_events.ecsv')
    tca = {str(r['event_id']): float(r['t_ca_mjd']) for r in covt}
    conf = json.load(open(ROOT + 'surveys/ptf-crossings/results/confirmatory_v1.json'))
    uev = {(u['target_id'], rung_of(u['radius_au']), u['band']): sorted(u['events_usable_epochs']) for u in conf['units']}
    units = [
        ('van-maanen', '2.5 Rsun', 'g', 1, '2 (same-night pair)', '1.45 / 2.16', None, 21.37, False, 0.00019,
         'b = 0.28 Rsun, 2011 (evt-5c45890fa28c) — photosphere-grazing cell, pre-2015 era, first constraint: relay through the 2.5 Rsun cone '
         '>~ ~190 W = 0.00019 MW (g) at 90%; 1 trial (S_event)', None),
        ('van-maanen', '0.1 AU', 'g', 2, '4 + 12', '1.19 / 3.16', '0.95 / 2.98', 22.0, True, 0.008,
         'April recurrence family; ~8 kW = 0.008 MW through the 0.1 AU cone (censored depth); 2 trials (S_event + S_stack); '
         'recurrence-stack cell closes clean', None),
        ('van-maanen', '0.1 AU', 'R', 3, '2 + 4 + 3', '2.02 / 17.81', '1.51 / 8.14', 20.58, False, 0.03,
         'includes the PS1-lost 2010 event evt-4c4ea2c39734 (2 independent R epochs); threshold 17.8 confusion-dominated (bright static '
         'source on one ring trajectory; track passes 3.0 arcsec from nearest catalogued static); ~30 kW = 0.03 MW (confusion-limited); '
         '2 trials; 658 nm generic leakage in R', None),
        ('ross-128', '2.5 Rsun', 'R', 2, '7 + 2', '1.58 / 2.61', '1.25 / 2.24', 22.0, True, 0.00011,
         'ross-128 2.5 Rsun windows (2010–2012) constrain <~ 110 W = 0.00011 MW (R, censored depth); 2 trials; recurrence stack closes clean', None),
        ('ross-128', '0.1 AU', 'R', 2, '18 + 2', '2.76 / 5.48', '2.56 / 5.26', 21.61, False, 0.011,
         '~11 kW = 0.011 MW through the 0.1 AU cone; 2 trials; recurrence stack closes clean', None),
        ('ross-128', '0.1 AU', 'g', 1, '44', '3.86 / 5.12', None, 22.0, True, 0.008,
         '~8 kW = 0.008 MW (censored depth); 1 trial', None),
    ]
    for tgt, rung, band, nev, eps, se, ss, m90, cens, pw, note, er in units:
        ntr = 1 if ss is None else 2
        evs = uev[(tgt, rung, band)]
        ts = [tca[e] for e in evs]
        er = f"MJD {min(ts):.1f}–{max(ts):.1f}" if len(ts) > 1 else f"MJD {ts[0]:.1f}"
        note = note + '; events ' + ', '.join(evs)
        n = f"searched, null: S_event / T = {se}" + (f"; S_stack / T = {ss}" if ss else '') + f"; events (epochs) {nev} ({eps}); "
        n += ('m90 >= 22.0 grid-censored; ' if cens else '') + note
        rows.append(row(tgt, 'B', rung, band, SUB_B, 'searched', 'm90', m90, 'AB mag', None, pw, None, nev, ntr,
                        er, SRC, n))
    rows.append(row('programme', 'B', 'any', 'g/R', SUB_B, 'searched', 'm90', None, 'AB mag', '20.58–22.0', None, '0.00011–0.03',
                    None, 10, '2009–2013', R + ' §2',
                    'blind confirmatory run: 10 trials over 6 units (7 searched units frozen incl. gj-1276 A), 0 exceedances vs 1.11 expected '
                    'control crossings (P(0) ~ 0.31); 600/600 cutouts, zero 404s; k rescales 1.00–1.10; powers rescaled from the ZTF anchor '
                    '(m90 21.8 -> 130 W through the 2.5 Rsun cone at 550 AU, gain ~9e9)'))
    # constraint-only and ledger classes
    rows.append(row('gj-1276', 'A', '0.1 AU', 'R', SUB_A, 'constraint_only', None, None, None, None, None, None, None, 0, None,
                    R + ' §1 item 6; surveys/ptf-crossings/results/dev_search_v1.md',
                    'dev: resolved constraint-only at the frozen offset-validity gate (1/12 valid pseudo-window offsets — campaign cadence); '
                    'S = -3.84 reported (dev_search_v1.md); the only A unit surviving the saturation cut; 0 searched trials'))
    rows.append(row('wolf-359', 'B', '0.1 AU', 'R', SUB_B, 'structurally_open', None, None, None, None, None, None, 1, 0, None,
                    R + ' §1 item 6; surveys/ptf-crossings/results/dev_search_v1.md',
                    'single_epoch class (dev, 0 trials): 1 in-window epoch, S 0.10 vs T 2.10 (8/8 ring controls), clean; nearest catalogued static 8.2 arcsec'))
    for tgt, band, mag in (('teegarden', 'R', 13.7), ('van-maanen', 'g', 12.4), ('ross-128', 'R', 9.9)):
        rows.append(row(tgt, 'A', '0.1 AU', band, SUB_A, 'structurally_open', None, None, None, None, None, None, 1, None, None,
                        R + ' §1 item 4; surveys/ptf-crossings/results/saturation_cut_v1.*',
                        f'saturation-excluded A unit: star {mag} mag vs frozen level E_R 14.5 / E_g 15.0 (Mould-R via Jordi 2006); '
                        + ('ross-128 star core fires dmask bits 8+6 at dev (exclusion confirmed)' if tgt == 'ross-128' else 'excluded-class')))
    CS = 'surveys/ptf-crossings/results/coverage_v1_summary.json; ' + R + ' §1 item 3'
    rows.append(row('programme', 'B', '1.2 Rsun', 'g/R', SUB_B, 'structurally_open', n_events=24, source=CS,
                    notes='0/24 events covered — structurally uncovered: the 0.3–0.6 d windows fall between epochs (closest 0.36 d vs a ±0.3 d edge); '
                          'coverage-without-statistic'))
    rows.append(row('programme', 'B', '2.5 Rsun', 'g/R', SUB_B, 'structurally_open', n_events=27, source=CS,
                    notes='3/30 events covered (van-maanen 2011 b = 0.28 Rsun; ross-128 x2); uncovered remainder'))
    rows.append(row('programme', 'B', '0.1 AU', 'g/R', SUB_B, 'structurally_open', n_events=35, source=CS,
                    notes='7/42 events covered (van-maanen April family 2009/2010/2011/2013, ross-128 x2, wolf-359 x1); uncovered remainder'))
    rows.append(row(['teegarden', 'gj-908', 'ross-154'], 'B', 'any', 'g/R', SUB_B, 'structurally_open', source=R + ' §1 item 3; §3',
                    notes='campaign cadence: these antipodes have zero PTF epochs (zero-coverage positions; declared unconstrained)'))
    rows.append(row('programme', 'A', '0.1 AU', 'g/R', SUB_A, 'structurally_open', n_events=38, source=CS,
                    notes='4/42 events covered; uncovered remainder'))
    rows.append(row('programme', 'A', '1.0 AU', 'g/R', SUB_A, 'ledger_only', n_events=76, source=CS,
                    notes='76/513 events covered; rung constraint-only by the standing window~season theorem (freeze v1.0); zero trials, coverage recorded'))
    rows.append(row('programme', 'A', '1.0 AU', 'g/R', SUB_A, 'structurally_open', n_events=437, source=CS,
                    notes='uncovered remainder of the 1.0 AU rung (513 - 76)'))
    rows.append(row('programme', ['A', 'B'], 'any', '1064/1550 nm', 'PTF g + Mould R', 'not_constrainable', source=R + ' §1 item 2; §3',
                    notes='1064/1550 nm outside both bands; unconstrained'))
    rows.append(row('programme', 'S1', 'any', 'any', None, 'no_survey', source='surveys/ptf-crossings/hypotheses.md §1',
                    notes='sunward combinations out of scope (identical to the PS1 freeze; plan §5.8 item 9 owns them)'))
    rows.append(row('programme', 'S2', 'any', 'any', None, 'no_survey', source='surveys/ptf-crossings/hypotheses.md §1',
                    notes='sunward combinations out of scope (identical to the PS1 freeze; plan §5.8 item 9 owns them)'))
    rows.append(row('programme', 'B', 'any', 'g/R', SUB_B, 'structurally_open', source=R + ' §3',
                    notes='declared unconstrained: sub-exposure pulse schedules beyond x d scaling; transmitters avoiding Earth-crossing windows'))
    doc = dict(survey_id='ptf_crossings', report=R, plan_section="5.10", pipeline='B', archives=['ptf'],
               hypothesis_version="v1.0 (D1–D8) + amendment v1.1 (confirmatory cutouts 256 -> 384 px; no statistic/threshold/gate changed)",
               era="MJD 54891–57051 (2009–2015)", run_dir='runs/ptf-crossings', records_dir=None, candidates=0,
               candidate_notes='10 blind confirmatory trials, 0 exceedances vs 1.11 expected control crossings; no retained or ambiguous items',
               rows=rows)
    dump('ptf_crossings.yaml', doc)


gen_ztf(); gen_ps1(); gen_wise(); gen_joint(); gen_ptf()
