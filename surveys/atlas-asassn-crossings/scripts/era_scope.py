"""Pre-freeze era scoping: intersect the frozen universal crossing list
with the ATLAS and ASAS-SN archive eras (geometry only — no archive
data touched).

Eras are the recon-measured values (notes/atlas_asassn_recon_2026-08-25.md):
ATLAS MJD 57227 (unit 02a first light) -> 61275 (probe currency);
southern units (03a/04a, needed for delta < -50) from MJD 59578.
ASAS-SN Sky Patrol v2 servable era JD 2456595 -> 2460841 (DB ceiling
at recon; re-measure at freeze).

Output: results/era_scope_v0.json (freeze-time eras) or, with
--atlas-end / --asassn-end / --tag, a re-scoped copy (v1 = coverage
start 2026-09-03: ATLAS end = MJD 61286, ASAS-SN v2 ceiling re-measured
by asassn_currency.py).
"""
import argparse
import collections
import json
import os

import numpy as np
from astropy.table import Table

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
EVENTS = os.path.join(REPO, 'crossings', 'universal_v1', 'events.ecsv')
OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'era_scope_v0.json')

ERAS = {
    'atlas': (57227.0, 61275.0),
    'asassn_v2': (56595.0, 60841.0),
}
ATLAS_SOUTH_START = 59578.0  # delta < -50 has no ATLAS epochs before this
RUNGS = [
    ('B', '1.2Rsun', 'b_rsun', 1.2),
    ('B', '2.5Rsun', 'b_rsun', 2.5),
    ('B', '0.1AU', 'b_au', 0.1),
    ('A', '0.1AU', 'b_au', 0.1),
    ('A', '1.0AU', 'b_au', 1.0),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--atlas-end', type=float, default=ERAS['atlas'][1])
    ap.add_argument('--asassn-end', type=float, default=ERAS['asassn_v2'][1])
    ap.add_argument('--tag', default='v0')
    args = ap.parse_args()
    ERAS['atlas'] = (ERAS['atlas'][0], args.atlas_end)
    ERAS['asassn_v2'] = (ERAS['asassn_v2'][0], args.asassn_end)
    out_path = OUT.replace('_v0.json', f'_{args.tag}.json')
    t = Table.read(EVENTS)
    mjd = np.array(t['t_ca_tdb_jd']) - 2400000.5
    is_a = np.array(t['link_direction']) == 'inbound'
    ax = np.array(t['axis_distance_au'])
    # amendment A1 (2026-09-03): the sunward axis crossing of each year
    # (A with ax < 0, B with ax > 0: channel position at elongation ~0)
    # is out of scope; 'events' counts searchable ones, the
    # '_incl_sunward' fields keep the v0 (link-direction-only) count.
    not_sunward = np.where(is_a, ax > 0, ax < 0)
    dec_pos = np.where(is_a, t['star_icrs_dec_deg'], t['relay_icrs_dec_deg'])
    b = {'b_rsun': np.array(t['b_min_solar_radii']), 'b_au': np.array(t['b_min_au'])}

    out = {'events_input': EVENTS.replace(REPO + '/', ''),
           'n_events_total': len(t), 'eras_mjd': ERAS,
           'atlas_south_start_mjd': ATLAS_SOUTH_START, 'channels': {}}
    for era_name, (lo, hi) in ERAS.items():
        inera = (mjd >= lo) & (mjd <= hi)
        rows = {}
        for ch, rung, key, lim in RUNGS:
            m_all = inera & (is_a if ch == 'A' else ~is_a) & (b[key] <= lim)
            m = m_all & not_sunward
            south = m & (dec_pos < -50)
            per_target = collections.Counter(str(x) for x in t['target_id'][m])
            rows[f'{ch}_{rung}'] = {
                'events': int(m.sum()),
                'events_incl_sunward': int(m_all.sum()),
                'targets': len(per_target),
                'south_pos_events': int(south.sum()),
                'south_pos_pre2022_events': int((south & (mjd < ATLAS_SOUTH_START)).sum()),
                'per_target': ({k: {'events': v,
                                    'b_min_rsun': round(float(b['b_rsun'][m & (t['target_id'] == k)].min()), 3),
                                    'v_perp_med_km_s': round(float(np.median(t['v_perp_km_s'][m & (t['target_id'] == k)])), 1)}
                                for k, v in sorted(per_target.items())}
                               if lim <= 0.1 or key == 'b_rsun' else 'omitted (wide rung)'),
            }
        out['channels'][era_name] = rows

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(out, f, indent=2)
    print(f'wrote {os.path.normpath(out_path)}')
    for era_name, rows in out['channels'].items():
        print(f'-- {era_name}:', {k: v['events'] for k, v in rows.items()})


if __name__ == '__main__':
    main()
