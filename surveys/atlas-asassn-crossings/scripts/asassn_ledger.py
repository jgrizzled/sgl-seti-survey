"""ASAS-SN Sky Patrol v2 window coverage-fraction ledger v1 (hypotheses
§1 role (a), D7; recon open item "estimator").

Estimator: field-level union. Every catalogued source within the 3'
cone around a survey position (the currency snapshots,
runs/atlas-asassn-crossings/asassn/currency_*.parquet) lies on the
same 4.5-deg camera image as the position (8" pixels), so the set of
distinct image_ids (quality 'G') of those sources is the epoch list of
the field at the position. Per event x rung it counts in-window
distinct images (V / g), and the 8 pseudo-windows with data. Coverage
records only - not a null measurement at the antipode (channel B has
no v2 photometry at the point). For channel A the nearest catalogued
source to the star (separation, mean mag) is recorded as D2
saturation-cut material; the catalog rows are snapshotted.

Must run as a file (client spawns workers). Output:
results/asassn_ledger_v1_events.ecsv + asassn_ledger_v1_summary.json.
"""
import glob
import hashlib
import json
import os
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from astropy.table import Table

sys.path.insert(0, os.path.dirname(__file__))
from build_task_list import half_window_d, sep_arcsec, LADDER_A, LADDER_B  # noqa: E402

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
RES = os.path.join(os.path.dirname(__file__), '..', 'results')
RUN = os.path.join(REPO, 'runs', 'atlas-asassn-crossings', 'asassn')
TASKS = os.path.join(RES, 'task_list_v1.ecsv')
ASASSN_ERA = (56595.0, 60842.0)
PSEUDO = (-97, -71, -47, -23, 23, 47, 71, 97)
RADIUS_ARCSEC = 180.0


def snapshot_catalog(client, ch, tid, ra, dec):
    path = os.path.join(RUN, f'catalog_{ch}_{tid}.parquet')
    if os.path.exists(path):
        return pd.read_parquet(path)
    cat = client.cone_search(ra, dec, RADIUS_ARCSEC, units='arcsec',
                             catalog='master_list', download=False)
    cat = pd.DataFrame(cat)
    cat.to_parquet(path)
    with open(path + '.json', 'w') as f:
        json.dump({'query': {'ra_deg': ra, 'dec_deg': dec, 'radius_arcsec': RADIUS_ARCSEC,
                             'catalog': 'master_list'},
                   'fetched_utc': datetime.now(timezone.utc).isoformat(),
                   'rows': int(len(cat)),
                   'sha256': hashlib.sha256(open(path, 'rb').read()).hexdigest()}, f, indent=2)
    return cat


def main():
    from pyasassn.client import SkyPatrolClient
    client = SkyPatrolClient(verbose=False)
    tasks = Table.read(TASKS)
    tca_tasks = tasks[tasks['pos_role'] == 'tca']
    rows, nearest = [], {}
    for ch in ('B', 'A'):
        for tid in np.unique(tca_tasks['target_id']):
            ev = tca_tasks[(tca_tasks['channel'] == ch) & (tca_tasks['target_id'] == tid)]
            lc_path = glob.glob(os.path.join(RUN, f'currency_*_{ch}_{tid}.parquet'))
            lc = pd.read_parquet(lc_path[0]) if lc_path else pd.DataFrame()
            ra0, de0 = float(np.median(ev['ra_deg'])), float(np.median(ev['dec_deg']))
            cat = snapshot_catalog(client, ch, tid, ra0, de0)
            if ch == 'A' and len(cat):
                rc = 'ra_deg' if 'ra_deg' in cat else 'ra'
                dc = 'dec_deg' if 'dec_deg' in cat else 'dec'
                seps = [sep_arcsec(ra0, de0, float(r[rc]), float(r[dc])) for _, r in cat.iterrows()]
                i = int(np.argmin(seps))
                nearest[str(tid)] = {'asas_sn_id': str(cat.iloc[i]['asas_sn_id']),
                                     'sep_arcsec': round(seps[i], 1),
                                     **{c: (float(cat.iloc[i][c]) if pd.notna(cat.iloc[i][c]) else None)
                                        for c in cat.columns if 'mag' in c.lower()}}
            good = lc[lc['quality'] == 'G'] if len(lc) else lc
            epochs = (good.drop_duplicates('image_id')[['jd', 'image_id', 'camera', 'phot_filter']]
                      if len(good) else pd.DataFrame(columns=['jd', 'image_id', 'camera', 'phot_filter']))
            mjd = np.asarray(epochs['jd'], float) - 2400000.5
            filt = np.asarray(epochs['phot_filter'], str)
            for e in ev:
                tca = float(e['t_ca_mjd'])
                in_era = ASASSN_ERA[0] <= tca <= ASASSN_ERA[1]
                for rung, r in (LADDER_B if ch == 'B' else LADDER_A).items():
                    half = half_window_d(float(e['b_min_au']), r, float(e['v_perp_km_s']))
                    if not np.isfinite(half):
                        continue
                    dt = mjd - tca
                    inw = np.abs(dt) <= half
                    row = {'channel': ch, 'target_id': str(tid), 'event_id': str(e['event_id']),
                           'rung': rung, 't_ca_mjd': tca, 'in_asassn_era': in_era,
                           'half_window_d': half, 'n_images': int(inw.sum()),
                           'n_V': int((inw & (filt == 'V')).sum()), 'n_g': int((inw & (filt == 'g')).sum()),
                           'n_sources_cone': int(len(cat)),
                           'pseudo_with_data': int(sum((np.abs(dt - o) <= half).any() for o in PSEUDO))}
                    rows.append(row)
    out = Table(rows=rows)
    out.write(os.path.join(RES, 'asassn_ledger_v1_events.ecsv'), format='ascii.ecsv', overwrite=True)
    summ = {'estimator': 'field-level union of quality-G image_ids of catalogued sources within 3 arcmin',
            'era_mjd': list(ASASSN_ERA), 'channel_A_nearest_source': nearest, 'channels': {}}
    for ch in ('B', 'A'):
        summ['channels'][ch] = {}
        for rung in np.unique(out['rung'][out['channel'] == ch]):
            s = out[(out['channel'] == ch) & (out['rung'] == rung)]
            se = s[s['in_asassn_era']]
            cov = np.asarray(se['n_images']) > 0
            summ['channels'][ch][str(rung)] = {
                'events_total': len(s), 'events_in_era': len(se),
                'events_covered': int(cov.sum()),
                'median_images_covered': float(np.median(se['n_images'][cov])) if cov.any() else 0.0,
                'per_target': {str(t): [int((np.asarray(se['n_images'][se['target_id'] == t]) > 0).sum()),
                                        int((se['target_id'] == t).sum())] for t in np.unique(se['target_id'])}}
    with open(os.path.join(RES, 'asassn_ledger_v1_summary.json'), 'w') as f:
        json.dump(summ, f, indent=2)
    print(json.dumps(summ['channel_A_nearest_source'], indent=1))
    for ch, d in summ['channels'].items():
        for rung, v in d.items():
            print(ch, rung, {k: v[k] for k in v})


if __name__ == '__main__':
    main()
