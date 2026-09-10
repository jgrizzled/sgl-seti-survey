"""Target photometric scope for the outer-skirt blended search (recon).

For every registry target: Gaia DR3 G/BP/RP (VizieR I/355 by source
id where the registry carries one; otherwise a 30" cone at the
epoch-2016 registry position), proper motion, and an ATLAS o-band
estimate from a transform anchored on the seven targets whose o was
measured on the substrate (atlas-asassn-crossings
results/saturation_cut_v1.json). The estimate is a screening quantity
only; the ATLAS D2 rule (sat 12.5 + 0.5 margin) is applied to the
measured reduced-mode value at the cut stage, never to the estimate.
Output: results/target_scope_v1.json.
"""
import json
import os
import sys

import numpy as np
import requests
import yaml

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
REG = os.path.join(REPO, 'registries', 'pilot_wise_2026.yaml')
SAT = os.path.join(REPO, 'surveys', 'atlas-asassn-crossings', 'results', 'saturation_cut_v1.json')
OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'target_scope_v1.json')
TAP = 'https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync'
# cone matches that are NOT the target (Gaia has no entry for the star)
OVERRIDE = {
    'sirius-a': ('excluded_bright', 'V -1.5, no Gaia entry; cone match is a background star'),
    'sirius-b': ('excluded_bright', 'blended with Sirius A on the ATLAS PSF'),
    'alpha-cen-a': ('excluded_bright', 'V 0.0, no Gaia entry; cone match is a background star'),
    'alpha-cen-b': ('excluded_bright', 'V 1.3, no Gaia entry; cone match is a background star'),
    'procyon-a': ('excluded_bright', 'V 0.4, no Gaia entry'),
    'procyon-b': ('excluded_bright', 'blended with Procyon A on the ATLAS PSF'),
    'fomalhaut': ('excluded_bright', 'V 1.2, no Gaia entry; cone match is a background star'),
    'wise-0855': ('no_optical', 'Y dwarf, no optical counterpart; cone match is a background star'),
}


def tap(q):
    r = requests.get(TAP, params={'REQUEST': 'doQuery', 'LANG': 'ADQL', 'FORMAT': 'json', 'QUERY': q}, timeout=120)
    r.raise_for_status()
    j = r.json()
    cols = [c['name'] for c in j['metadata']]
    return [dict(zip(cols, row)) for row in j['data']]


def main():
    reg = yaml.safe_load(open(REG))['targets']
    sat = json.load(open(SAT))["targets"]
    out = {}
    for tid, v in reg.items():
        st = v['state']['astrometry']
        gid = next((i['id'] for i in v.get('identifiers', []) if i['catalog'].startswith('Gaia')), None)
        rec = {'display_name': v.get('display_name'), 'endpoint_kind': v.get('endpoint_kind'),
               'ra_deg': st['ra_deg'], 'dec_deg': st['dec_deg'],
               'pm_total_arcsec_yr': float(np.hypot(st['pm_ra_cosdec_mas_per_yr'], st['pm_dec_mas_per_yr']) / 1000),
               'pm_ra_cosdec_mas_yr': st['pm_ra_cosdec_mas_per_yr'], 'pm_dec_mas_yr': st['pm_dec_mas_per_yr'],
               'parallax_mas': st['parallax_mas'], 'ref_epoch': st['reference_epoch_jyear'], 'gaia_dr3': gid}
        if gid:
            rows = tap(f'SELECT Source, Gmag, BPmag, RPmag, "BP-RP", RUWE FROM "I/355/gaiadr3" WHERE Source = {gid}')
        else:
            rows = tap(f'SELECT TOP 3 Source, Gmag, BPmag, RPmag, "BP-RP", RUWE, RA_ICRS, DE_ICRS FROM "I/355/gaiadr3" '
                       f"WHERE 1=CONTAINS(POINT('ICRS',RA_ICRS,DE_ICRS), CIRCLE('ICRS',{st['ra_deg']},{st['dec_deg']},0.01)) ORDER BY Gmag")
        rec['gaia'] = rows[0] if rows else None
        rec['gaia_match'] = 'by_id' if gid else ('cone_brightest' if rows else 'none')
        out[tid] = rec
        print(tid, rec['gaia_match'], rows[0] if rows else None, flush=True)
    # transform o(G, BP-RP) anchored on the substrate-measured seven
    anchors = []
    for tid, s in sat.items():
        o = s.get('o', {}).get('median_mag') if isinstance(s.get('o'), dict) else None
        g = out.get(tid, {}).get('gaia')
        # only substrate-unsaturated anchors (measured o >= 11.4; the
        # brighter three are saturation-compressed and would bias the fit)
        if o is not None and g and g['Gmag'] >= 11.0 and g and g.get('BP-RP') is not None:
            anchors.append((tid, g['Gmag'], g['BP-RP'], o))
    A = np.array([[1.0, c] for _, G, c, o in anchors]); y = np.array([o - G for _, G, c, o in anchors])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ coef
    for tid, rec in out.items():
        g = rec['gaia']
        if g and g.get('BP-RP') is not None:
            rec['o_est'] = float(g['Gmag'] + coef[0] + coef[1] * g['BP-RP'])
            rec['o_screen'] = 'excluded_est' if rec['o_est'] < 12.5 else ('marginal_est' if rec['o_est'] < 13.0 else 'ok_est')
        else:
            rec['o_est'] = None; rec['o_screen'] = 'no_gaia'
        if tid in OVERRIDE:
            rec['o_screen'], rec['o_est'], rec['override_reason'] = OVERRIDE[tid][0], None, OVERRIDE[tid][1]
    meta = {'transform': {'form': 'o - G = a + b*(BP-RP)', 'a': float(coef[0]), 'b': float(coef[1]),
                          'anchors': [{'target': t, 'G': G, 'bp_rp': c, 'o_measured': o} for t, G, c, o in anchors],
                          'rms_resid_mag': float(np.sqrt(np.mean(resid ** 2)))}}
    json.dump({'meta': meta, 'targets': out}, open(OUT, 'w'), indent=1)
    print(json.dumps(meta, indent=1))


if __name__ == '__main__':
    main()
