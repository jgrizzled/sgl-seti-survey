"""Run the frozen search on one split (dev or confirmatory) for every
searched unit of configs/threshold_freeze_v1.json. Output:
results/<split>_search_v1.json with S, T, exceedances, S_opp
annotation, the c-band annotation for exceedances, and the airmass
fits. The confirmatory run is blind: executed once, unchanged.
"""
import argparse
import json
import os
import sys
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(__file__))
import search_core as sc  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--split', required=True, choices=['dev', 'confirmatory'])
    ap.add_argument('--freeze', default='threshold_freeze_v1.json')
    args = ap.parse_args()
    fz = json.load(open(os.path.join(sc.D, 'configs', args.freeze)))
    out = {'split': args.split, 'freeze_content_hash': fz['freeze_content_hash'], 'run_utc': datetime.now(timezone.utc).isoformat(),
           'units': [], 'n_trials': 0, 'n_exceedances': 0, 'expected_control_crossings': 0.0}
    for u in fz['units']:
        if u['split'] != args.split or u['class'] != 'searched':
            continue
        tid, rb = u['target_id'], u['rung_au']
        us = sc.unit_stats(tid, 'o', rb)
        th = sc.thresholds(us)
        rec = {'target_id': tid, 'rung_au': rb, 'band': 'o', 'airmass_fit_target': us['real']['airmass_fit'],
               'airmass_fit_controls': {k: c['airmass_fit'] for k, c in us['controls'].items()},
               'n_nights_target': us['real']['n_nights'], 'cycles_target': us['real']['cycles'], 'statistics': {}}
        for stat in u['statistics']:
            t = th[stat]
            rec['statistics'][stat] = t
            out['n_trials'] += 1
            out['expected_control_crossings'] += 1 / (t['n_valid_controls'] + 1) if t['searched'] else 0
            if t['exceedance']:
                out['n_exceedances'] += 1
                # c-band annotation
                try:
                    usc = sc.unit_stats(tid, 'c', rb)
                    thc = sc.thresholds(usc)
                    t['c_band'] = {k: {kk: thc[k][kk] for kk in ('S', 'T', 'n_valid_controls', 'S_opp')} for k in sc.STATS}
                except Exception as e:  # noqa: BLE001
                    t['c_band'] = {'error': str(e)[:200]}
            flag = 'EXCEEDANCE' if t['exceedance'] else ('searched' if t['searched'] else 'constraint-only(at-search)')
            print(f"{tid:12s} {rb:5s} {stat:8s} S={t['S'] if t['S'] is None else round(t['S'], 3)} T={t['T'] if t['T'] is None else round(t['T'], 3)} "
                  f"n_ctrl={t['n_valid_controls']} S_opp={None if t['S_opp'] is None else round(t['S_opp'], 3)} {flag}", flush=True)
        out['units'].append(rec)
    out['expected_control_crossings'] = round(out['expected_control_crossings'], 3)
    json.dump(out, open(os.path.join(sc.RES, f'{args.split}_search_v1.json'), 'w'), indent=1)
    print(json.dumps({k: out[k] for k in ('split', 'n_trials', 'n_exceedances', 'expected_control_crossings')}))


if __name__ == '__main__':
    main()
