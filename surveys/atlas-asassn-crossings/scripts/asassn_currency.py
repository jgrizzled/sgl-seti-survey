"""ASAS-SN Sky Patrol v2 currency re-measure (hypotheses §2 / D7).

Coverage-start action (resume checklist step 1): measure the servable
DB ceiling (max JD) of catalogued sources adjacent to the 7 narrow-rung
survey positions per channel (star side and antipode side), anonymous
v2 client, plain HTTP port 9006. Coverage records only — no window-
locked quantity is formed. Raw light curves are snapshotted under
runs/atlas-asassn-crossings/asassn/ (parquet + sidecar).

Must run as a file (download=True spawns worker processes).
"""
import hashlib
import json
import os
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from astropy.table import Table
from astropy.time import Time

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
EVENTS = os.path.join(REPO, 'crossings', 'universal_v1', 'events.ecsv')
RUN = os.path.join(REPO, 'runs', 'atlas-asassn-crossings', 'asassn')
OUT = os.path.join(os.path.dirname(__file__), '..', 'results')
TARGETS = ['van-maanen', 'wolf-359', 'gj-1276', 'teegarden', 'ross-128',
           'ross-154', 'gj-908']
ERA = (57227.0, 61286.0)
RADIUS_ARCSEC = 180.0


def positions():
    t = Table.read(EVENTS)
    mjd = Time(list(t['t_ca_utc']), format='isot', scale='utc').mjd
    era = (mjd >= ERA[0]) & (mjd <= ERA[1])
    out = []
    for tid in TARGETS:
        m = era & (np.asarray(t['target_id']) == tid) & (np.asarray(t['b_min_au']) <= 0.1)
        for ch, ld, cols in (('A', 'inbound', ('star_icrs_ra_deg', 'star_icrs_dec_deg')),
                             ('B', 'outbound', ('relay_icrs_ra_deg', 'relay_icrs_dec_deg'))):
            s = t[m & (np.asarray(t['link_direction']) == ld)]
            out.append((ch, tid, float(np.median(s[cols[0]])), float(np.median(s[cols[1]]))))
    return out


def main():
    from pyasassn.client import SkyPatrolClient
    os.makedirs(RUN, exist_ok=True)
    client = SkyPatrolClient(verbose=False)
    date = datetime.now(timezone.utc).strftime('%Y-%m-%d')
    rows = []
    for ch, tid, ra, dec in positions():
        cat = client.cone_search(ra, dec, RADIUS_ARCSEC, units='arcsec',
                                 catalog='master_list', download=False)
        n_src = 0 if cat is None else len(cat)
        rec = {'channel': ch, 'target_id': tid, 'ra_deg': ra, 'dec_deg': dec,
               'radius_arcsec': RADIUS_ARCSEC, 'n_sources': int(n_src),
               'jd_max': None, 'jd_min': None, 'n_epochs': 0,
               'n_epochs_good': 0, 'cameras': [], 'filters': []}
        if n_src:
            lcs = client.cone_search(ra, dec, RADIUS_ARCSEC, units='arcsec',
                                     catalog='master_list', download=True, threads=1)
            df = lcs.data if hasattr(lcs, 'data') else lcs
            df = pd.DataFrame(df)
            path = os.path.join(RUN, f'currency_{date}_{ch}_{tid}.parquet')
            df.to_parquet(path)
            sha = hashlib.sha256(open(path, 'rb').read()).hexdigest()
            with open(path + '.json', 'w') as f:
                json.dump({'query': {'ra_deg': ra, 'dec_deg': dec,
                                     'radius_arcsec': RADIUS_ARCSEC,
                                     'catalog': 'master_list',
                                     'host': 'http://asassn-lb01.ifa.hawaii.edu:9006'},
                           'fetched_utc': datetime.now(timezone.utc).isoformat(),
                           'rows': int(len(df)), 'sha256': sha,
                           'client': 'skypatrol 0.6.21 (pyasassn)'}, f, indent=2)
            if len(df):
                rec.update({'jd_max': float(df['jd'].max()), 'jd_min': float(df['jd'].min()),
                            'n_epochs': int(len(df)),
                            'n_epochs_good': int((df['quality'] == 'G').sum()) if 'quality' in df else None,
                            'cameras': sorted(map(str, df['camera'].unique())) if 'camera' in df else [],
                            'filters': sorted(map(str, df['phot_filter'].unique())) if 'phot_filter' in df else [],
                            'n_sources_with_epochs': int(df['asas_sn_id'].nunique()) if 'asas_sn_id' in df else None})
        print(json.dumps(rec), flush=True)
        rows.append(rec)
    jd_max = max((r['jd_max'] for r in rows if r['jd_max']), default=None)
    summary = {'measured_utc': datetime.now(timezone.utc).isoformat(),
               'jd_max_all': jd_max, 'mjd_max_all': (jd_max - 2400000.5) if jd_max else None,
               'iso_max_all': Time(jd_max, format='jd').isot if jd_max else None,
               'recon_ceiling_jd': 2460841.0, 'positions': rows}
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, f'asassn_currency_{date}.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    print('DB ceiling JD', jd_max, summary['iso_max_all'])


if __name__ == '__main__':
    main()
