"""Channel-A saturation table (hypotheses D2, ATLAS o/c).

Decision quantity: the star's own o/c brightness measured on the survey
substrate - the 20-day reduced-mode (target-image tphot) series pulled
by aux_pulls.py at each channel-A position, outside every window - plus
tphot's own saturation/fit flags. Rule (D2, ZTF shape, ATLAS limit
sat ~ 12.5 with a 0.5 mag margin), per band:
  excluded : median m < 12.5, or > 50 % of rows fail the FAQ mask
             (saturated PSF fits fail err/chi/apfit terms); no rows -> no_data
  marginal : 12.5 <= median m < 13.0
  ok       : median m >= 13.0
Gaia DR3 G/BP/RP (aux/gaia_dr3_targets.json) recorded as the catalog
anchor; the transform-free measured value decides. Output:
results/saturation_cut_v1.json.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import atlas_api  # noqa: E402

REPO = atlas_api.REPO
AUX = os.path.join(REPO, 'runs', 'atlas-asassn-crossings', 'aux')
OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'saturation_cut_v1.json')
SAT, MARGIN = 12.5, 0.5
BANDS = ('o', 'c')


def status(m, fail_frac, n):
    if n == 0:
        return 'no_data'
    if fail_frac > 0.5 or (m is not None and m < SAT):
        return 'excluded'
    if m < SAT + MARGIN:
        return 'marginal'
    return 'ok'


def main():
    gaia = json.load(open(os.path.join(AUX, 'gaia_dr3_targets.json')))
    table = {}
    for tid, g in gaia.items():
        path = os.path.join(AUX, f'aux-reduced-{tid}.txt')
        sc = json.load(open(path + '.json'))
        df = atlas_api.read_result(path)
        rec = {'gaia': g, 'n_rows': int(len(df)), 'mjd_range': sc['query_params']['mjd_min'],
               'result_sha256': sc['result_sha256']}
        for b in BANDS:
            d = df[df.F == b]
            n = int(len(d))
            if n:
                ok = atlas_api.faq_quality_mask(d)
                fail = float(1 - ok.mean())
                # magnitude from positive fluxes on FAQ-passing rows; fall
                # back to all positive-flux rows if the mask empties it
                src = d[ok] if ok.any() else d
                pos = src[src.uJy > 0]
                m = float(np.median(pos.m)) if len(pos) else None
                err_frac = float((d.err != 0).mean())
                chin = float(d['chi/N'].median())
            else:
                fail, m, err_frac, chin = 1.0, None, 1.0, None
            rec[b] = {'n': n, 'median_mag': None if m is None else round(m, 2),
                      'faq_fail_frac': round(fail, 3), 'err_flag_frac': round(err_frac, 3),
                      'chi_per_n_median': None if chin is None else round(chin, 1),
                      'status': status(m, fail, n)}
        table[tid] = rec
        print(tid, {b: (rec[b]['median_mag'], rec[b]['faq_fail_frac'], rec[b]['status']) for b in BANDS},
              'G', g['G'], 'RP', g['RP'])
    with open(OUT, 'w') as f:
        json.dump({'rule': f'excluded m<{SAT} or FAQ-fail>50% or no rows; marginal {SAT}<=m<{SAT+MARGIN}; ok m>={SAT+MARGIN}',
                   'substrate': 'ATLAS reduced-mode tphot 20-d series outside all windows (aux_pulls.py)',
                   'targets': table}, f, indent=1)
    print('wrote', os.path.normpath(OUT))


if __name__ == '__main__':
    main()
