"""Empirical fixed-position PSF-fit response R(d) (amendment A3 support):
forced photometry of two quiet dev controls at deliberate Dec offsets
0.5", 1.0", 2.0" (mjd_min 59000), matched by Obs id to the on-position
series -> flux ratio vs offset in units of the per-exposure PSF sigma.
Fits R = exp(-d^2 / (k sigma^2)). Results: results/response_calib_v1.json.
  python response_calib.py --submit   (queue + fetch, ~40 min)
  python response_calib.py            (analyse)
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import search_core as sc  # noqa: E402

STARS = ['C-gj-2012-4', 'C-teegarden-4']
OFFSETS = [0.5, 1.0, 2.0]
OUT = os.path.join(sc.RUN, 'response')


def key(star, off):
    return f'RC-{star}-{off:.1f}'


def submit():
    hdrs = sc.atlas_api.headers()
    os.makedirs(OUT, exist_ok=True)
    jobs = []
    for star in STARS:
        row = sc.task_row(star)
        for off in OFFSETS:
            p = {'ra': float(row['ra_deg']), 'dec': float(row['dec_deg']) + off / 3600.0, 'mjd_min': 59000.0, 'use_reduced': 'true',
                 'radec_epoch_year': '2016.0', 'propermotion_ra': float(row['pm_ra_cosdec_mas_yr']), 'propermotion_dec': float(row['pm_dec_mas_yr']),
                 'comment': key(star, off)}
            j = sc.atlas_api.submit(p, hdrs)
            jobs.append((key(star, off), j, p))
            print('queued', key(star, off), j['url'], flush=True)
    for k, j, p in jobs:
        s = sc.atlas_api.wait_and_fetch(j['url'], hdrs, os.path.join(OUT, f'{k}.txt'), query_params=p)
        print('done', k, s['wall_s'], flush=True)


def analyse():
    res = {'stars': {}, 'fit': {}}
    x2, y = [], []
    for star in STARS:
        on = sc.atlas_api.read_result(os.path.join(sc.LC, f'{star}.txt'))
        on = on[sc.atlas_api.faq_quality_mask(on) & (on.F == 'o')].set_index('Obs')
        res['stars'][star] = {}
        for off in OFFSETS:
            p = os.path.join(OUT, f'{key(star, off)}.txt')
            if not os.path.exists(p):
                continue
            df = sc.atlas_api.read_result(p)
            df = df[(df.err == 0) & (df.F == 'o')].set_index('Obs')
            j = on.join(df, how='inner', lsuffix='_on', rsuffix='_off')
            j = j[(j.uJy_on > 0) & (j.uJy_off > -1e9)]
            sig = (1.86 * np.sqrt(j.maj_on.astype(float) * j.min_on.astype(float)) / 2.355)
            ratio = (j.uJy_off.astype(float) / j.uJy_on.astype(float)).values.astype(float)
            u2 = (off / sig.values) ** 2
            # bin in u2
            rec = {'n': int(len(j)), 'ratio_median': float(np.median(ratio)), 'u2_median': float(np.median(u2)),
                   'sigma_median_arcsec': float(np.median(sig))}
            res['stars'][star][str(off)] = rec
            good = np.isfinite(ratio) & (ratio > 0.05) & (ratio < 1.5)
            x2 += list(u2[good]); y += list(np.log(ratio[good]))
            print(star, off, rec, flush=True)
    x2, y = np.array(x2), np.array(y)
    x2_all, y_all = x2.copy(), y.copy()
    # small-offset regime (0.5" and 1.0" offsets, u2 < 0.6): the parallax
    # displacements are <= 0.3", where the PSF core sets the response
    x2, y = x2[x2 < 0.6], y[x2 < 0.6]
    # ln R = -u2/k  -> k = -sum(x2^2)/sum(x2 y)  (least squares through origin), robust by clipping
    for it in range(3):
        k = -np.sum(x2 * x2) / np.sum(x2 * y)
        resid = y + x2 / k
        m = np.abs(resid - np.median(resid)) < 3 * 1.4826 * np.median(np.abs(resid - np.median(resid)))
        x2, y = x2[m], y[m]
    k_all = -np.sum(x2_all * x2_all) / np.sum(x2_all * y_all)
    res['fit'] = {'k': float(k), 'n': int(len(x2)), 'k_all_offsets': float(k_all), 'regime': 'u2 < 0.6 (0.5\" and 1.0\" offsets)', 'model': 'R = exp(-d^2/(k sigma^2)), sigma = FWHM/2.355, FWHM = 1.86" sqrt(maj min)',
                  'R_at_0.26_arcsec_sigma_1.57': float(np.exp(-0.26 ** 2 / (k * 1.57 ** 2)))}
    print(res['fit'])
    json.dump(res, open(os.path.join(sc.RES, 'response_calib_v1.json'), 'w'), indent=1)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--submit', action='store_true'); a = ap.parse_args()
    submit() if a.submit else analyse()
