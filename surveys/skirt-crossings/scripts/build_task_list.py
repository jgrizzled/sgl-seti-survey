"""Task list for the outer-skirt survey (hypotheses D3, D6, D7, D8).

One full-history reduced-mode ATLAS task per star with server-side
proper motion: the 18 frozen targets (registry astrometry, epoch
2016.0, Gaia DR3 PM) and, per target, 8 controls drawn with the
frozen seed (20260907) from the recon pool of 16
(results/control_pool_v1.json; Gaia DR3 positions epoch 2016.0 and PM).
Output: results/task_list_v1.ecsv (+ _summary.json).
"""
import json
import os

import numpy as np
from astropy.table import Table

D = os.path.join(os.path.dirname(__file__), '..', 'results')
SEED = 20260907
N_CONTROLS = 8
TARGETS = ['gj-1111', 'wolf-1069', 'teegarden', 'gj-915', 'gj-3512', 'gj-3306', 'gj-3112', 'gj-293',
           'gj-1087', 'gj-9193', 'gj-2012', 'gj-13157', 'gj-518', 'gj-12724', 'gj-11547', 'gj-1276',
           'gj-11068', 'eps-ind-b']
DEV = ['teegarden', 'gj-2012']
LEDGER_ONLY = {'wolf-1069': ['0.9', '0.95'], 'gj-293': ['0.9', '0.95'], 'gj-13157': ['0.9'], 'gj-3112': ['0.9']}
MJD_MIN = 55000.0


def main():
    scope = json.load(open(os.path.join(D, 'target_scope_v1.json')))['targets']
    pools = json.load(open(os.path.join(D, 'control_pool_v1.json')))['targets']
    rng = np.random.default_rng(SEED)
    rows = []
    for tid in TARGETS:
        t = scope[tid]
        rows.append({'task_key': f'T-{tid}', 'target_id': tid, 'role': 'target', 'control_index': -1,
                     'gaia_source': str(t['gaia']['Source']), 'ra_deg': t['ra_deg'], 'dec_deg': t['dec_deg'],
                     'epoch_year': float(t['ref_epoch']), 'pm_ra_cosdec_mas_yr': float(t['pm_ra_cosdec_mas_yr']),
                     'pm_dec_mas_yr': float(t['pm_dec_mas_yr']), 'G': float(t['gaia']['Gmag']),
                     'bp_rp': float(t['gaia']['BP-RP']), 'split': 'dev' if tid in DEV else 'confirmatory',
                     'mjd_min': MJD_MIN, 'use_reduced': True})
        pool = pools[tid]['pool']
        idx = rng.choice(len(pool), size=N_CONTROLS, replace=False)
        for j, i in enumerate(sorted(idx)):
            c = pool[i]
            rows.append({'task_key': f'C-{tid}-{j}', 'target_id': tid, 'role': 'control', 'control_index': j,
                         'gaia_source': str(c['Source']), 'ra_deg': float(c['RA_ICRS']), 'dec_deg': float(c['DE_ICRS']),
                         'epoch_year': 2016.0, 'pm_ra_cosdec_mas_yr': float(c['pmRA'] or 0.0),
                         'pm_dec_mas_yr': float(c['pmDE'] or 0.0), 'G': float(c['Gmag']), 'bp_rp': float(c['BP-RP']),
                         'split': 'dev' if tid in DEV else 'confirmatory', 'mjd_min': MJD_MIN, 'use_reduced': True})
    tab = Table(rows=rows)
    tab.meta = {'seed': SEED, 'n_controls': N_CONTROLS, 'ledger_only': LEDGER_ONLY, 'dev': DEV,
                'hypotheses': 'surveys/skirt-crossings/hypotheses.md v1.0 (frozen 2026-09-07)'}
    tab.write(os.path.join(D, 'task_list_v1.ecsv'), overwrite=True)
    summ = {'n_tasks': len(tab), 'n_targets': len(TARGETS), 'n_controls': N_CONTROLS * len(TARGETS),
            'seed': SEED, 'dev': DEV, 'ledger_only': LEDGER_ONLY}
    json.dump(summ, open(os.path.join(D, 'task_list_v1_summary.json'), 'w'), indent=1)
    print(summ)


if __name__ == '__main__':
    main()
