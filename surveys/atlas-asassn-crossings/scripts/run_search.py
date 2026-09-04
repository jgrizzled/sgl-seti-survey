"""Run the frozen statistics on one split (threshold_freeze_v1.json).

  python run_search.py --split dev            -> results/dev_search_v1.json
  python run_search.py --split confirmatory   -> results/confirmatory_search_v1.json

Dev: machinery validation on the D8 dev targets (wolf-359, gj-908) -
every unit incl. constraint-only ones is computed and reported with
diagnostics (k values, baseline sizes, control distributions, KS of
standardized baseline residuals). Confirmatory: the blind run; only the
searched statistics of searched units are trials, everything else is
reported as constraint-only. Exceedances are listed for adjudication -
never dropped.
"""
import argparse
import json
import os
import sys
from datetime import datetime, timezone

import numpy as np
from scipy import stats as sps

sys.path.insert(0, os.path.dirname(__file__))
import search_core as sc  # noqa: E402

RES = os.path.join(sc.D, 'results')


def ks_diag(ch, tid, eid, band):
    """KS of standardized baseline residuals (real-window geometry)."""
    es = sc.EventSeries(ch, tid, eid, band)
    zs = []
    for pos in es.positions:
        det, k, nb = es.detrended(pos, es.t_ca, es.half_widest)
        if len(det):
            base = es.baseline_mask(det.MJD.values, es.t_ca, es.half_widest)
            zs.append(det.r.values[base] / np.sqrt(det.v.values[base]))
    if not zs:
        return None
    z = np.concatenate(zs)
    d, p = sps.kstest(z, 'norm')
    return {'n': int(len(z)), 'ks_D': round(float(d), 4), 'ks_p': float(p),
            'std': round(float(np.std(z)), 3), 'frac_gt3': round(float(np.mean(z > 3)), 4)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--split', choices=('dev', 'confirmatory'), required=True)
    ap.add_argument('--only', nargs='*')
    args = ap.parse_args()
    fz = sc.FREEZE
    min_stack = int(fz['statistics']['S_stack'].split('>= ')[1].split()[0])
    units = [u for u in fz['search_units'] if u['split'] == args.split]
    if args.only:
        units = [u for u in units if u['target_id'] in args.only]
    results, exceedances = [], []
    for u in units:
        key = f"{u['channel']} {u['rung']} {u['target_id']} {u['band']}"
        if u['class'] in ('uncovered', 'saturation_excluded'):
            results.append({**{k: u[k] for k in ('target_id', 'channel', 'rung', 'band', 'class', 'split')},
                            'computed': False})
            print(key, u['class'], flush=True)
            continue
        st = sc.unit_stats(u['channel'], u['target_id'], u['rung'], u['band'], u['covered_event_ids'])
        th = sc.thresholds(st, min_stack)
        rec = {**{k: u[k] for k in ('target_id', 'channel', 'rung', 'band', 'class', 'split')},
               'computed': True, 'searched_statistics': list(u['statistics']),
               'thresholds': th,
               'real': st['real'],
               'ks': {eid: ks_diag(u['channel'], u['target_id'], eid, u['band'])
                      for eid in u['covered_event_ids']} if args.split == 'dev' else None}
        for stat, r in th.items():
            trial = stat in u['statistics']
            r['trial'] = trial
            # sanity: the freeze's validity count must match the engine's
            if trial and r['n_valid_controls'] != u['statistics'][stat]['n_valid_controls']:
                r['validity_mismatch'] = u['statistics'][stat]['n_valid_controls']
            if r['exceedance'] and trial:
                exceedances.append({'unit': key, 'statistic': stat, 'S': round(r['S'], 3),
                                    'T': round(r['T'], 3), 'z': st['real'][stat]['z']})
        results.append(rec)
        print(key, u['class'], {s: (None if r['S'] is None else round(r['S'], 2),
                                    None if r['T'] is None else round(r['T'], 2),
                                    r['n_valid_controls'], 'EXC' if r['exceedance'] else '')
                                for s, r in th.items()}, flush=True)
    n_trials = sum(len(u['statistics']) for u in units if u['class'] == 'searched')
    exp = round(sum(u['expected_control_crossings'] for u in units if u['class'] == 'searched'), 3)
    out = {'split': args.split, 'run_utc': datetime.now(timezone.utc).isoformat(),
           'freeze_hash': fz['freeze_content_hash'], 'n_units': len(units), 'n_trials': n_trials,
           'expected_control_crossings': exp, 'n_exceedances': len(exceedances),
           'exceedances': exceedances, 'units': results}
    path = os.path.join(RES, f'{args.split}_search_v1.json')
    with open(path, 'w') as f:
        json.dump(out, f, indent=1, default=float)
    print(f'trials {n_trials}, expected control crossings {exp}, exceedances {len(exceedances)} -> {os.path.normpath(path)}')


if __name__ == '__main__':
    main()
