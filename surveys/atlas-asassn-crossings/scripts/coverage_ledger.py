"""ATLAS coverage ledger v1 (hypotheses §6/§7 + amendments §12).

Coverage counting only - no window-locked flux statistic is formed.
For every pulled windowed series (runs/atlas-asassn-crossings/lc/,
sidecar = done-list) apply the frozen primary quality mask (FAQ recipe,
§7) and count, per event x rung x band:
  in-window exposures and distinct nights (window = t_ca +/- half,
  flat chord, §6), off-window baseline exposures inside the +/-110 d
  pull (the D5 detrend / variance-rescale material), and how many of
  the 8 pseudo-windows (+/-23/47/71/97 d, same width) hold >= 1
  exposure (the control family's raw availability).
Channel B 0.1 AU rung: counted per mini-track position (tca / ingress /
egress) plus the union of distinct exposures (Obs ids) over the three.
H-filter rows are counted separately and excluded from o/c.

Runs on whatever series exist (--partial); the summary records the
task completeness so a partial ledger can never pass for the freeze
input. Output: results/coverage_v1_events.ecsv + coverage_v1_summary.json.
"""
import argparse
import json
import os
import sys

import numpy as np
from astropy.table import Table

sys.path.insert(0, os.path.dirname(__file__))
import atlas_api  # noqa: E402

REPO = atlas_api.REPO
RES = os.path.join(os.path.dirname(__file__), '..', 'results')
LC = os.path.join(REPO, 'runs', 'atlas-asassn-crossings', 'lc')
TASKS = os.path.join(RES, 'task_list_v1.ecsv')
BANDS = ('o', 'c')
RUNGS = {'B': (('1.2Rsun', 'half_1p2Rsun_d'), ('2.5Rsun', 'half_2p5Rsun_d'),
               ('0.1AU', 'half_0p1AU_d')),
         'A': (('0.1AU', 'half_0p1AU_d'),)}
PSEUDO_OFFSETS = (-97, -71, -47, -23, 23, 47, 71, 97)


def night(mjd):
    # ATLAS units span Hawaii/Chile/S.Africa; MJD+0.5 floor groups a
    # local night for Hawaii (UTC-10) well enough for counting.
    return np.floor(np.asarray(mjd) + 0.5).astype(int)


def load_series(key):
    path = os.path.join(LC, f'{key}.txt')
    sc = path + '.json'
    if not os.path.exists(sc):
        return None, 'missing'
    with open(sc) as f:
        meta = json.load(f)
    if 'error' in meta:
        return None, 'failed'
    try:
        df = atlas_api.read_result(path)
    except Exception as e:  # noqa: BLE001
        return None, f'parse:{e}'
    if len(df) == 0:
        return df, 'empty'
    return df, 'ok'


def count_block(df, tca, half):
    """Per-band counts for one (series, window) pair."""
    out = {}
    dt = np.asarray(df.MJD) - tca
    inw = np.abs(dt) <= half
    for band in BANDS:
        mb = np.asarray(df.F == band)
        sel = mb & inw
        out[f'n_{band}'] = int(sel.sum())
        out[f'nights_{band}'] = int(len(np.unique(night(df.MJD[sel])))) if sel.any() else 0
        out[f'base_{band}'] = int((mb & ~inw).sum())
        npw = 0
        for i, off in enumerate(PSEUDO_OFFSETS):
            has = bool((mb & (np.abs(dt - off) <= half)).any())
            out[f'pw{i}_{band}'] = int(has)
            npw += has
        out[f'pseudo_{band}'] = npw
    out['n_total'] = out['n_o'] + out['n_c']
    out['nights_total'] = int(len(np.unique(night(df.MJD[inw & np.isin(df.F, BANDS)])))) if inw.any() else 0
    out['obs_ids'] = set(df.Obs[inw & np.isin(df.F, BANDS)])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--partial', action='store_true')
    args = ap.parse_args()
    tasks = Table.read(TASKS)
    status = {}
    series = {}
    for r in tasks:
        k = str(r['task_key'])
        df, st = load_series(k)
        status[k] = st
        if df is not None:
            series[k] = df
    n_ok = sum(1 for s in status.values() if s in ('ok', 'empty'))
    if n_ok < len(tasks) and not args.partial:
        sys.exit(f'{len(tasks) - n_ok} of {len(tasks)} series not available; use --partial')

    rows = []
    masked = {}
    for r in tasks:
        k = str(r['task_key'])
        if k not in series:
            continue
        df = series[k]
        m = atlas_api.faq_quality_mask(df) if len(df) else np.zeros(0, bool)
        masked[k] = {'rows': int(len(df)), 'kept': int(m.sum()),
                     'n_H': int((df.F == 'H').sum()) if len(df) else 0}
        dfm = df[m]
        ch = str(r['channel'])
        tca = float(r['t_ca_mjd'])
        for rung, col in RUNGS[ch]:
            half = float(r[col])
            if not np.isfinite(half):
                continue
            if str(r['pos_role']) != 'tca' and rung != '0.1AU':
                continue
            c = count_block(dfm, tca, half)
            obs = c.pop('obs_ids')
            rows.append({'channel': ch, 'target_id': str(r['target_id']),
                         'event_id': str(r['event_id']), 'task_key': k,
                         'pos_role': str(r['pos_role']), 'rung': rung,
                         't_ca_mjd': tca, 'b_min_rsun': float(r['b_min_rsun']),
                         'half_window_d': half, 'window_d': 2 * half,
                         'rows_raw': masked[k]['rows'], 'rows_kept': masked[k]['kept'],
                         **c, 'n_obs_ids': len(obs)})
            if ch == 'B' and rung == '0.1AU':
                r_union = rows[-1]
                r_union['_obs'] = obs
    # mini-track union for B 0.1AU
    union = {}
    for row in rows:
        if row['channel'] == 'B' and row['rung'] == '0.1AU':
            union.setdefault(row['event_id'], set()).update(row.pop('_obs'))
    for row in rows:
        row['n_union_track'] = len(union[row['event_id']]) if (
            row['channel'] == 'B' and row['rung'] == '0.1AU') else row['n_obs_ids']
    out = Table(rows=rows)
    os.makedirs(RES, exist_ok=True)
    out.write(os.path.join(RES, 'coverage_v1_events.ecsv'), format='ascii.ecsv', overwrite=True)

    summary = {'tasks_total': len(tasks), 'series_available': n_ok,
               'series_failed': sum(1 for s in status.values() if s == 'failed'),
               'partial': n_ok < len(tasks), 'mask': 'FAQ primary (§7)',
               'faq_keep_fraction_median': float(np.median(
                   [v['kept'] / v['rows'] for v in masked.values() if v['rows']])) if masked else None,
               'pseudo_offsets_d': list(PSEUDO_OFFSETS), 'channels': {}}
    for ch in ('B', 'A'):
        summary['channels'][ch] = {}
        for rung, _ in RUNGS[ch]:
            s = out[(out['channel'] == ch) & (out['rung'] == rung) & (out['pos_role'] == 'tca')]
            if len(s) == 0:
                continue
            cov = np.asarray(s['n_total']) > 0
            per_t = {}
            for tid in np.unique(s['target_id']):
                st = s[s['target_id'] == tid]
                per_t[str(tid)] = {'events': len(st),
                                   'covered': int((np.asarray(st['n_total']) > 0).sum()),
                                   'nights_total': int(np.sum(st['nights_total'])),
                                   'exposures_total': int(np.sum(st['n_total']))}
            summary['channels'][ch][rung] = {
                'events_available': len(s), 'events_covered': int(cov.sum()),
                'events_covered_o': int((np.asarray(s['n_o']) > 0).sum()),
                'events_covered_c': int((np.asarray(s['n_c']) > 0).sum()),
                'median_in_window_exposures_covered': float(np.median(s['n_total'][cov])) if cov.any() else 0.0,
                'median_window_d': float(np.median(s['window_d'])),
                'pseudo_windows_with_data_mean_o': float(np.mean(s['pseudo_o'])),
                'per_target': per_t}
            if ch == 'B' and rung == '0.1AU':
                summary['channels'][ch][rung]['events_covered_union_track'] = int(
                    (np.asarray(s['n_union_track']) > 0).sum())
    with open(os.path.join(RES, 'coverage_v1_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    print(json.dumps({k: v for k, v in summary.items() if k != 'channels'}, indent=1))
    for ch, d in summary['channels'].items():
        for rung, v in d.items():
            print(ch, rung, {k: v[k] for k in v if k != 'per_target'})


if __name__ == '__main__':
    main()
