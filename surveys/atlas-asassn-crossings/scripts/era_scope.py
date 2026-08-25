"""Pre-freeze era scoping: intersect the frozen universal crossing list
with the ATLAS and ASAS-SN archive eras (geometry only — no archive
data touched).

Eras are the recon-measured values (notes/atlas_asassn_recon_2026-08-25.md):
ATLAS MJD 57227 (unit 02a first light) -> 61275 (probe currency);
southern units (03a/04a, needed for delta < -50) from MJD 59578.
ASAS-SN Sky Patrol v2 servable era JD 2456595 -> 2460841 (DB ceiling
at recon; re-measure at freeze).

Output: results/era_scope_v0.json
"""
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
    t = Table.read(EVENTS)
    mjd = np.array(t['t_ca_tdb_jd']) - 2400000.5
    is_a = np.array(t['link_direction']) == 'inbound'
    dec_pos = np.where(is_a, t['star_icrs_dec_deg'], t['relay_icrs_dec_deg'])
    b = {'b_rsun': np.array(t['b_min_solar_radii']), 'b_au': np.array(t['b_min_au'])}

    out = {'events_input': EVENTS.replace(REPO + '/', ''),
           'n_events_total': len(t), 'eras_mjd': ERAS,
           'atlas_south_start_mjd': ATLAS_SOUTH_START, 'channels': {}}
    for era_name, (lo, hi) in ERAS.items():
        inera = (mjd >= lo) & (mjd <= hi)
        rows = {}
        for ch, rung, key, lim in RUNGS:
            m = inera & (is_a if ch == 'A' else ~is_a) & (b[key] <= lim)
            south = m & (dec_pos < -50)
            per_target = collections.Counter(str(x) for x in t['target_id'][m])
            rows[f'{ch}_{rung}'] = {
                'events': int(m.sum()),
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

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, 'w') as f:
        json.dump(out, f, indent=2)
    print(f'wrote {os.path.normpath(OUT)}')
    for era_name, rows in out['channels'].items():
        print(f'-- {era_name}:', {k: v['events'] for k, v in rows.items()})


if __name__ == '__main__':
    main()
