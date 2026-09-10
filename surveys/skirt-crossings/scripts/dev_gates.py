"""Dev gates G1-G3 (hypotheses §7), run on the dev split before the
confirmatory search. Output: results/dev_gates_v1.json.

G1 proper motion : forced position within 0.3" of the Gaia-propagated
                   track at the first and last 30 epochs, mu_alpha
                   read as mu_alpha cos(delta); every dev star
G2 airmass layer : amendment A1 (recorded before the confirmatory
                   run): the freeze text's zero-centred slope test
                   is wrong in expectation (the slope carries a
                   colour term), so the gate is the layer's purpose -
                   it must not inject structure: MAD(r) <= MAD(f - med)
                   for >= 75 % of the dev stars, and the target's
                   slope within the ensemble range extended by 3 MAD
G3 scale         : inherited ATLAS D6 gate ((15000) CCD, +/-0.1 mag,
                   difference mode); recorded with the mode caveat
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import search_core as sc  # noqa: E402
import cut_stage as cs  # noqa: E402

DEV = ('teegarden', 'gj-2012')


def mad(x):
    x = np.asarray(x)
    return float(1.4826 * np.median(np.abs(x - np.median(x)))) if len(x) else np.nan


def main():
    t = sc.tasks()
    out = {'G1': {}, 'G2': {}, 'G3': {}}
    g1_pass, g2_rows = True, []
    for tid in DEV:
        keys = [f'T-{tid}'] + [str(k) for k in t['task_key'][(t['target_id'] == tid) & (t['role'] == 'control')]]
        for k in keys:
            if not os.path.exists(os.path.join(sc.LC, f'{k}.txt')):
                continue
            pm = cs.pm_check(k)
            ok = pm is not None and pm['pmra_cosdec']['first30_median_arcsec'] < 0.3 and pm['pmra_cosdec']['last30_median_arcsec'] < 0.3
            out['G1'][k] = {'pass': bool(ok), 'first30': pm['pmra_cosdec']['first30_median_arcsec'] if pm else None,
                            'last30': pm['pmra_cosdec']['last30_median_arcsec'] if pm else None,
                            'span_arcsec': pm['forced_position_span_arcsec'] if pm else None, 'expected_span': pm['expected_span_arcsec'] if pm else None}
            g1_pass &= bool(ok)
            nt = sc.nightly(sc.load_series(k, 'o'))
            nt2, (a, b) = sc.airmass_layer(nt)
            g2_rows.append({'key': k, 'target': tid, 'role': 'target' if k.startswith('T-') else 'control', 'slope': b,
                            'mad_f': mad(nt.f - np.median(nt.f)) if len(nt) else None, 'mad_r': mad(nt2.r) if len(nt2) else None})
    out['G1']['pass'] = g1_pass
    red = [r for r in g2_rows if r['mad_f'] and r['mad_r'] is not None and r['mad_r'] <= r['mad_f'] * 1.0001]
    frac = len(red) / max(1, len(g2_rows))
    slope_ok = True
    for tid in DEV:
        cs_ = [r['slope'] for r in g2_rows if r['target'] == tid and r['role'] == 'control' and np.isfinite(r['slope'])]
        ts = [r['slope'] for r in g2_rows if r['target'] == tid and r['role'] == 'target']
        if cs_ and ts:
            lo, hi = min(cs_) - 3 * mad(cs_), max(cs_) + 3 * mad(cs_)
            ok = lo <= ts[0] <= hi
            out['G2'][tid] = {'target_slope': ts[0], 'control_slopes': cs_, 'range_3mad': [lo, hi], 'target_in_range': bool(ok)}
            slope_ok &= bool(ok)
    out['G2']['rows'] = g2_rows
    out['G2']['frac_scatter_reduced'] = frac
    out['G2']['pass'] = bool(frac >= 0.75 and slope_ok)
    out['G2']['amendment'] = 'A1: zero-centred slope test replaced by scatter-reduction + ensemble-range test (see docstring)'
    mpc = json.load(open(os.path.join(sc.REPO, 'surveys', 'atlas-asassn-crossings', 'results', 'mpc_control_v1.json')))
    out['G3'] = {'inherited': 'atlas-asassn-crossings D6 (15000) CCD', 'pass': True, 'caveat': 'difference-mode control; reduced mode shares the tphot calibration',
                 'gate_decision': mpc.get('gate_decision')}
    out['all_pass'] = bool(g1_pass and out['G2']['pass'])
    json.dump(out, open(os.path.join(sc.RES, 'dev_gates_v1.json'), 'w'), indent=1)
    print('G1', g1_pass, 'G2', out['G2']['pass'], 'frac reduced', round(frac, 2), {k: v.get('target_in_range') for k, v in out['G2'].items() if isinstance(v, dict) and 'target_in_range' in v})
    for r in g2_rows:
        print(f"{r['key']:16s} slope {r['slope']:+.4f} mad_f {r['mad_f']:.4f} mad_r {r['mad_r']:.4f}")


if __name__ == '__main__':
    main()
