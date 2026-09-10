"""Adjudication diagnostics for every exceedance of a search split
(hypotheses §6 veto ladder, rules fixed before the confirmatory run):

  (1) symmetry   : a beam gives equal Delta on both sides; report
                   S_skirt, S_opp and the per-side weighted-mean
                   Delta; 'asymmetric' when the two sides differ in
                   sign or the weaker side is < 25 % of the stronger
  (2) chromatic  : c-band S and Delta (annotation)
  (3) per-cycle  : leave-one-cycle-out S; 'single_cycle' when
                   dropping one cycle takes S below max(T, 0)
  (4) ensemble   : the target's per-cycle Delta against the mean and
                   scatter of the 8 controls' Delta in the same cycles
                   ('ensemble_shared' when |Delta_t - mean_c| < 2 sd_c
                   in every valid cycle -- the excess is in the field)
  (5) neighbours : Gaia DR3 neighbours within 10" (recon control_pool
                   target_neighbours_10as) reported
Disposition: non_promotable_<rule> for the first failing rule in the
order (1) (3) (4); else retained_ambiguous_colour_unmatched when no
control lies within 1.0 mag of the target's BP-RP (amendment A4);
otherwise retained_ambiguous (hand-off to the
enumerated Sky Patrol v1 follow-up). Output:
results/<split>_adjudication_v1.json.
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import search_core as sc  # noqa: E402


def loo(cs, name, min_in):
    """Leave-one-cycle-out S_sym-type statistic."""
    ds, vs, ids = [], [], []
    for i, rec in cs.items():
        r = rec[name]
        if r['n_in'] >= min_in and r['n_out'] >= sc.MIN_OUT and r['v']:
            ds.append(r['delta']); vs.append(r['v']); ids.append(i)
    ds, vs = np.array(ds), np.array(vs)
    out = {}
    for j in range(len(ds)):
        m = np.ones(len(ds), bool); m[j] = False
        out[str(ids[j])] = float(np.sum(ds[m] / vs[m]) / np.sqrt(np.sum(1 / vs[m]))) if m.sum() else None
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--split', required=True)
    args = ap.parse_args()
    srch = json.load(open(os.path.join(sc.RES, f'{args.split}_search_v1.json')))
    pools = json.load(open(os.path.join(sc.RES, 'control_pool_v1.json')))['targets']
    out = {'split': args.split, 'exceedances': []}
    for u in srch['units']:
        tid, rb = u['target_id'], u['rung_au']
        exc = [s for s, t in u['statistics'].items() if t['exceedance']]
        if not exc:
            continue
        us = sc.unit_stats(tid, 'o', rb)
        real = us['real']
        rec = {'target_id': tid, 'rung_au': rb, 'statistics': exc, 'diagnostics': {}, 'disposition': {}}
        # (1) symmetry
        sk, op = real['S_skirt'], real['S_opp']
        dsk, dop = sk.get('delta_mean'), op.get('delta_mean')
        asym = None
        if dsk is not None and dop is not None:
            asym = (np.sign(dsk) != np.sign(dop)) or (min(abs(dsk), abs(dop)) < 0.25 * max(abs(dsk), abs(dop)))
        rec['diagnostics']['symmetry'] = {'S_skirt': sk['S'], 'S_opp': op['S'], 'delta_skirt': dsk, 'delta_opp': dop, 'asymmetric': None if asym is None else bool(asym)}
        # (2) chromatic
        try:
            usc = sc.unit_stats(tid, 'c', rb, controls=[])
            rec['diagnostics']['c_band'] = {s: {'S': usc['real'][s]['S'], 'delta_mean': usc['real'][s].get('delta_mean')} for s in ('S_sym', 'S_skirt', 'S_opp')}
        except Exception as e:  # noqa: BLE001
            rec['diagnostics']['c_band'] = {'error': str(e)[:100]}
        # (3) per-cycle
        rec['diagnostics']['leave_one_cycle_out'] = {}
        for stat in exc:
            T = max(u['statistics'][stat]['T'], 0.0)
            name, min_in = ('skirt', sc.MIN_SKIRT) if stat == 'S_skirt' else ('sym', sc.MIN_IN)
            l = loo(real['cycles'], name, min_in) if stat != 'S_year' else {}
            rec['diagnostics']['leave_one_cycle_out'][stat] = {'T': T, 'S_without_cycle': l, 'single_cycle': bool(any(v is not None and v <= T for v in l.values())) if l else None}
        # (4) ensemble
        ens = {}
        for i, cyc in real['cycles'].items():
            dt = cyc['sym']['delta']
            dc = [c['cycles'].get(i, {}).get('sym', {}).get('delta') for c in us['controls'].values()]
            dc = [x for x in dc if x is not None]
            if dt is not None and len(dc) >= 3:
                ens[str(i)] = {'delta_target': dt, 'mean_controls': float(np.mean(dc)), 'sd_controls': float(np.std(dc, ddof=1)), 'n_controls': len(dc),
                               'within_2sd': bool(abs(dt - np.mean(dc)) < 2 * np.std(dc, ddof=1))}
        rec['diagnostics']['ensemble'] = ens
        shared = bool(ens) and all(v['within_2sd'] for v in ens.values())
        rec['diagnostics']['neighbours_10as'] = pools.get(tid, {}).get('target_neighbours_10as', [])
        # A4: colour match of the ensemble + control colour trend (annotation)
        t = sc.tasks()
        tcol = float(sc.task_row(f'T-{tid}')['bp_rp'])
        ccols = {str(k): float(sc.task_row(str(k))['bp_rp']) for k in t['task_key'][(t['target_id'] == tid) & (t['role'] == 'control')]}
        gap = min(abs(c - tcol) for c in ccols.values())
        trend = None
        try:
            xs = np.array([ccols[k] for k in us['controls']]); ys = np.array([us['controls'][k]['S_sym'].get('delta_mean') or np.nan for k in us['controls']])
            m = np.isfinite(ys)
            if m.sum() >= 4 and np.ptp(xs[m]) > 0.05:
                pf = np.polyfit(xs[m], ys[m], 1)
                trend = {'slope_per_mag': float(pf[0]), 'extrapolated_delta_at_target_colour': float(np.polyval(pf, tcol)),
                         'target_delta': real['S_sym'].get('delta_mean'), 'control_colour_span': [float(xs[m].min()), float(xs[m].max())]}
        except Exception as e:  # noqa: BLE001
            trend = {'error': str(e)[:100]}
        rec['diagnostics']['colour'] = {'target_bp_rp': tcol, 'min_control_gap_mag': gap, 'unmatched': bool(gap > 1.0), 'control_trend': trend}
        for stat in exc:
            if asym:
                d = 'non_promotable_asymmetric'
            elif rec['diagnostics']['leave_one_cycle_out'].get(stat, {}).get('single_cycle'):
                d = 'non_promotable_single_cycle'
            elif shared:
                d = 'non_promotable_ensemble_shared'
            elif gap > 1.0:
                d = 'retained_ambiguous_colour_unmatched'
            else:
                d = 'retained_ambiguous'
            rec['disposition'][stat] = d
            print(tid, rb, stat, d, rec['diagnostics']['symmetry'], flush=True)
        out['exceedances'].append(rec)
    json.dump(out, open(os.path.join(sc.RES, f'{args.split}_adjudication_v1.json'), 'w'), indent=1)
    print('exceedances adjudicated:', len(out['exceedances']))


if __name__ == '__main__':
    main()
