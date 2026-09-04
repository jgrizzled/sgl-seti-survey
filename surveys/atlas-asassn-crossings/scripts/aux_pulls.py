"""Auxiliary ATLAS pulls for the threshold freeze (not survey statistics).

1. D6 positive control: MPC forced photometry of (15000) CCD (H 14.4,
   a 2.77 AU, condition code 0; V ~ 17.8-19.5 over an apparition),
   difference mode, MJD 60310 -> 61286 (two apparitions), through the
   identical task -> mask -> detrend chain later.
2. D2 saturation table material: 20-day reduced-mode (target-image
   tphot) series at each channel-A star position, placed 45-65 d after
   a mid-era event's t_ca - outside every real window of every rung
   (nearest real window edge > 30 d) - so the star's own o/c
   brightness and tphot saturation flags are measured on the survey
   substrate without touching a window-locked quantity.
Also snapshots the Gaia DR3 rows (G/BP/RP, PM, parallax) queried at
freeze prep for the seven stars as the catalog anchor.
Outputs under runs/atlas-asassn-crossings/aux/ with atlas_api sidecars.
"""
import json
import os
import sys

import numpy as np
from astropy.table import Table

sys.path.insert(0, os.path.dirname(__file__))
import atlas_api  # noqa: E402

REPO = atlas_api.REPO
AUX = os.path.join(REPO, 'runs', 'atlas-asassn-crossings', 'aux')
TASKS = os.path.join(os.path.dirname(__file__), '..', 'results', 'task_list_v1.ecsv')
ASTEROID = {'mpc_name': '15000', 'mjd_min': 60310.0, 'mjd_max': 61286.0,
            'use_reduced': 'false', 'comment': 'aux-mpc-15000'}
GAIA = {  # Gaia DR3 sync query 2026-09-04 (gea.esac.esa.int TAP), 36" cone at the 2016 star position
    'van-maanen': (2552928187080872832, 12.299163, 12.531230, 11.933642, 231.78),
    'wolf-359': (3864972938605115520, 11.038391, 13.770287, 9.585450, 415.18),
    'gj-1276': (2611561706216413696, 15.421229, 16.081415, 14.661167, 117.14),
    'teegarden': (35227046884571776, 12.263103, 15.318630, 10.790128, 260.99),
    'ross-128': (3796072592206250624, 9.601000, 11.361026, 8.327767, 296.31),
    'ross-154': (4075141768785646848, 9.126414, 10.732183, 7.898147, 336.03),
    'gj-908': (2739689239311660672, 8.152917, 9.216571, 7.120805, 169.22),
}


def main():
    os.makedirs(AUX, exist_ok=True)
    with open(os.path.join(AUX, 'gaia_dr3_targets.json'), 'w') as f:
        json.dump({k: dict(zip(('source_id', 'G', 'BP', 'RP', 'parallax_mas'), v))
                   for k, v in GAIA.items()}, f, indent=1)
    hdrs = atlas_api.headers()
    tasks = Table.read(TASKS)
    jobs = []
    for tid in GAIA:
        ev = tasks[(tasks['channel'] == 'A') & (tasks['target_id'] == tid)]
        i = int(np.argmin(np.abs(np.asarray(ev['t_ca_mjd']) - 59500)))
        e = ev[i]
        lo = float(e['t_ca_mjd']) + 45.0
        jobs.append((f'aux-reduced-{tid}', {'ra': float(e['ra_deg']), 'dec': float(e['dec_deg']),
                                            'mjd_min': lo, 'mjd_max': lo + 20.0,
                                            'use_reduced': 'true', 'comment': f'aux-reduced-{tid}'}))
    jobs.append(('aux-mpc-15000', ASTEROID))
    for key, p in jobs:
        out = os.path.join(AUX, key + '.txt')
        if os.path.exists(out + '.json'):
            print('done', key)
            continue
        j = atlas_api.submit(p, hdrs)
        print('submitted', key, j['url'], flush=True)
        sc = atlas_api.wait_and_fetch(j['url'], hdrs, out, query_params=p, poll_s=20)
        print('fetched', key, sc['result_bytes'], 'bytes', sc['wall_s'], 's', flush=True)


if __name__ == '__main__':
    main()
