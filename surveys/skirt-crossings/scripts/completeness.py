"""Completeness by response-model step injection (threshold freeze
'completeness' block). For every searched unit and searched statistic:
a top-hat step of flux f (o-band AB grid) x R_i is added to every
IN exposure (both sides of the rung) of the target's real series
before normalisation; the frozen chain runs with a 20 % random
nightly-epoch dropout per draw; recovery = S > max(T, 0) against the
unit's real frozen T (from the search output). m90 = faintest grid
magnitude with >= 90 % recovery. Power P = F_line pi r_b^2.
Output: results/completeness_v1.json.
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import search_core as sc  # noqa: E402

LAMBDA_C, DLAMBDA = 690e-9, 260e-9
C = 2.99792458e8
DNU = C * DLAMBDA / LAMBDA_C ** 2
AU = 1.495978707e11


def f_line(m):
    return 3631e-26 * 10 ** (-0.4 * m) * DNU


def power_W(m, rb):
    return f_line(m) * np.pi * (float(rb) * AU) ** 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n-draw', type=int, default=None)
    ap.add_argument('--only', nargs='*')
    args = ap.parse_args()
    fz = json.load(open(os.path.join(sc.D, 'configs', 'threshold_freeze_v1.json')))
    cfg = fz['completeness']
    n_draw = args.n_draw or cfg['n_draw']
    mags = sorted(cfg['mag_grid'])   # ascending magnitude (bright to faint); stop once recovery reaches zero
    searches = {}
    for split in ('dev', 'confirmatory'):
        p = os.path.join(sc.RES, f'{split}_search_v1.json')
        if os.path.exists(p):
            for u in json.load(open(p))['units']:
                searches[(u['target_id'], u['rung_au'])] = u
    rng = np.random.default_rng(fz['seed'])
    out = {'seed': fz['seed'], 'n_draw': n_draw, 'dropout': cfg['dropout'], 'units': []}
    for u in fz['units']:
        if u['class'] != 'searched' or (args.only and u['target_id'] not in args.only):
            continue
        tid, rb = u['target_id'], u['rung_au']
        srch = searches.get((tid, rb))
        if srch is None:
            continue
        key = f'T-{tid}'
        mask = sc.inject_mask(tid, rb)
        rec = {'target_id': tid, 'rung_au': rb, 'split': u['split'], 'statistics': {}}
        Ts = {stat: max(srch['statistics'][stat]['T'], 0.0) for stat in u['statistics'] if srch['statistics'][stat]['searched']}
        curves = {stat: {} for stat in Ts}
        for m in mags:
            f_uJy = 10 ** (0.4 * (23.9 - m))
            hits = {stat: 0 for stat in Ts}
            for _ in range(n_draw):
                res = sc.star_stats(key, tid, 'o', rb, inject=(mask, f_uJy), dropout=(rng, cfg['dropout']))
                for stat in Ts:
                    r = res[stat]
                    if r['valid'] and r['S'] is not None and r['S'] > Ts[stat]:
                        hits[stat] += 1
            for stat in Ts:
                curves[stat][str(m)] = hits[stat] / n_draw
            if all(hits[s] == 0 for s in Ts):
                break
        for stat, T in Ts.items():
            ms = sorted((float(k), v) for k, v in curves[stat].items())
            m90 = None
            for m, rec_ in ms:
                if rec_ >= 0.9:
                    m90 = m
            rec['statistics'][stat] = {'T': T, 'recovery': curves[stat], 'm90': m90,
                                       'power_W_m90': power_W(m90, rb) if m90 is not None else None}
            print(f"{tid:12s} {rb} {stat:8s} T={T:.3f} m90={m90} P={None if m90 is None else f'{power_W(m90, rb) / 1e6:.2f} MW'}", flush=True)
        out['units'].append(rec)
        json.dump(out, open(os.path.join(sc.RES, 'completeness_v1.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
