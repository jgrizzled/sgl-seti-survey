"""Null-ensemble control pool (catalogue-level prep, recon stage).

For every photometrically eligible target (target_scope_v1.json:
ok_est / marginal_est) select a pool of matched field stars from Gaia
DR3 (VizieR I/355) within 1.0 deg (same ATLAS exposures, hence the same
solar-elongation / airmass / twilight history): |dG| <= 0.3 (widened
to 0.6 if the pool is thin), |d(BP-RP)| <= 0.4 (widened stepwise to
1.5 for the very red targets), not flagged VARIABLE, RUWE < 1.4,
PM < 100 mas/yr, and isolated (no Gaia neighbour within 10" brighter
than G_target + 3 -- tphot fits a single PSF). Pool = best 16 by
|dG| + |d(BP-RP)|; the freeze draws the 8 controls from the pool with
the declared seed. The target's own isolation is recorded the same
way. Output: results/control_pool_v1.json.
"""
import json
import os
import time

import numpy as np
import requests

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
SCOPE = os.path.join(os.path.dirname(__file__), '..', 'results', 'target_scope_v1.json')
GEOM = os.path.join(os.path.dirname(__file__), '..', 'results', 'skirt_geometry_v1.json')
OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'control_pool_v1.json')
TAP = 'https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync'
RADIUS_DEG, POOL = 1.0, 16


def tap(q, tries=3):
    for i in range(tries):
        try:
            r = requests.get(TAP, params={'REQUEST': 'doQuery', 'LANG': 'ADQL', 'FORMAT': 'json', 'QUERY': q}, timeout=180)
            r.raise_for_status()
            j = r.json()
            cols = [c['name'] for c in j['metadata']]
            return [dict(zip(cols, row)) for row in j['data']]
        except Exception as e:  # noqa: BLE001
            if i == tries - 1:
                raise
            time.sleep(10)


def neighbours(ra, dec, gmax, rad_arcsec=10.0):
    q = (f'SELECT Source, Gmag, RA_ICRS, DE_ICRS FROM "I/355/gaiadr3" WHERE Gmag < {gmax:.2f} AND '
         f"1=CONTAINS(POINT('ICRS',RA_ICRS,DE_ICRS), CIRCLE('ICRS',{ra},{dec},{rad_arcsec / 3600:.6f}))")
    return tap(q)


def main():
    scope = json.load(open(SCOPE))['targets']
    geom = json.load(open(GEOM))['targets']
    out = {}
    for tid, t in scope.items():
        if tid not in geom or t['o_screen'] not in ('ok_est', 'marginal_est'):
            continue
        g = t['gaia']
        G, col = g['Gmag'], g['BP-RP']
        # epoch-2016 catalogue position of the target (registry) for the cone centre
        ra, dec = t['ra_deg'], t['dec_deg']
        # target isolation (exclude itself; neighbours within 10" brighter than G+3)
        nb = [n for n in neighbours(ra, dec, G + 3) if str(n['Source']) != str(g['Source'])]
        rec = {'G': G, 'bp_rp': col, 'target_neighbours_10as': nb, 'pool': [], 'widening': None}
        for dG, dcol in [(0.3, 0.4), (0.6, 0.6), (0.6, 1.0), (0.8, 1.5), (1.0, 2.5), (1.0, 3.5), (1.5, 4.5)]:
            q = (f'SELECT Source, RA_ICRS, DE_ICRS, pmRA, pmDE, Gmag, "BP-RP", RUWE, VarFlag FROM "I/355/gaiadr3" WHERE '
                 f"1=CONTAINS(POINT('ICRS',RA_ICRS,DE_ICRS), CIRCLE('ICRS',{ra},{dec},{RADIUS_DEG})) "
                 f'AND Gmag BETWEEN {G - dG:.3f} AND {G + dG:.3f} AND "BP-RP" BETWEEN {col - dcol:.3f} AND {col + dcol:.3f} '
                 f"AND RUWE < 1.4 AND ABS(pmRA) < 100 AND ABS(pmDE) < 100 AND VarFlag != 'VARIABLE' "
                 f"AND Source != {g['Source']}")
            rows = tap(q)
            rows = [r for r in rows if r['Gmag'] is not None and r['BP-RP'] is not None]
            for r in rows:
                r['score'] = abs(r['Gmag'] - G) + abs(r['BP-RP'] - col)
                r['sep_deg'] = float(np.degrees(np.arccos(np.clip(
                    np.sin(np.radians(dec)) * np.sin(np.radians(r['DE_ICRS'])) +
                    np.cos(np.radians(dec)) * np.cos(np.radians(r['DE_ICRS'])) * np.cos(np.radians(ra - r['RA_ICRS'])), -1, 1))))
            rows.sort(key=lambda r: r['score'])
            pool = []
            for r in rows:
                if len(pool) >= POOL:
                    break
                if r['sep_deg'] < 0.01:      # not the target's own neighbourhood
                    continue
                nbs = [n for n in neighbours(r['RA_ICRS'], r['DE_ICRS'], r['Gmag'] + 3) if str(n['Source']) != str(r['Source'])]
                if nbs:
                    continue
                r['neighbours_10as'] = 0
                pool.append(r)
            rec['pool'] = pool
            rec['widening'] = {'dG': dG, 'dcol': dcol, 'n_candidates': len(rows)}
            if len(pool) >= POOL:
                break
        out[tid] = rec
        print(f"{tid:14s} G {G:5.2f} BP-RP {col:4.2f} pool {len(rec['pool'])} widening {rec['widening']} target-nbrs {len(nb)}", flush=True)
    json.dump({'radius_deg': RADIUS_DEG, 'pool_size': POOL, 'targets': out}, open(OUT, 'w'), indent=1)


if __name__ == '__main__':
    main()
