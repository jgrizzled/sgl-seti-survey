"""D6 completeness by response-model injection (threshold freeze v1.0).

For every searched unit and searched statistic: inject a flat-chord
source of flux f (o-band AB magnitude grid) into the unit's real
series at a random off-window centre per covered event (a pseudo-
window that avoids the real 0.1 AU window and the pull edges, and
mirrors the coverage gate), with the per-epoch response R_i(offset,
FWHM) of the z = 550 AU track geometry shifted to the injection centre
(each mini-track position receives f x R_i for its own offset), then
evaluate the frozen statistic exactly as the search does (best-R
position per epoch, max over the nested z family) and test S > max(T,
0) against the unit's frozen threshold. 100 draws per magnitude, seed
20260904; m90 = faintest magnitude with >= 90 % recovery (linear
interpolation on the grid). Baseline medians and k are unaffected by
an in-window injection by construction (the evaluated window is
excluded from the baseline), so the injected residual is r + f R.

Power: monochromatic line through the o band (560-820 nm, effective
width 260 nm, lambda_c 690 nm): F_line = 3631 Jy 10^(-0.4 m) dnu;
channel B cone power P = F_line pi r^2 (d/z)^2 (relay at z beaming
through the rung's cone of radius r at the Sun, Earth at d = z - 1 AU);
channel A 10-m diffraction-limited uplink P = F_line pi (1.22 lambda/D
D_star)^2. The D6 gate (MPC control |residual| <= 0.2 mag) is read from
results/mpc_control_v1.json and stamped on every depth.
Output: results/completeness_v1.json (+ power_limits_v1.json).
"""
import argparse
import json
import os
import sys
from datetime import datetime, timezone

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import search_core as sc  # noqa: E402

RES = os.path.join(sc.D, 'results')
SEED = 20260904
N_DRAW = 100
MAGS = np.arange(16.0, 23.01, 0.5)
Z_INJ = 550.0
LAMBDA_C, DLAMBDA = 690e-9, 260e-9
C = 2.99792458e8
DNU = C * DLAMBDA / LAMBDA_C ** 2
AU = 1.495978707e11
RSUN = 6.957e8
RUNG_R = {'1.2Rsun': 1.2 * RSUN, '2.5Rsun': 2.5 * RSUN, '0.1AU': 0.1 * AU}
GAIA = json.load(open(os.path.join(sc.REPO, 'runs', 'atlas-asassn-crossings', 'aux', 'gaia_dr3_targets.json')))


def f_line(m):
    return 3631e-26 * 10 ** (-0.4 * m) * DNU


def power_W(ch, rung, tid, m, z=Z_INJ):
    F = f_line(m)
    if ch == 'B':
        d = (z - 1.0) * AU
        return F * np.pi * RUNG_R[rung] ** 2 * (d / (z * AU)) ** 2
    dstar = 1000.0 / GAIA[tid]['parallax_mas'] * 3.0857e16
    return F * np.pi * (1.22 * LAMBDA_C / 10.0 * dstar) ** 2


def draw_offset(rng, es, h):
    lo = es.half_widest + h + 1.0
    hi = 110.0 - h - 1.0
    for _ in range(30):
        o = rng.uniform(lo, hi) * rng.choice([-1, 1])
        tc = es.t_ca + o
        if any(((np.abs(df.MJD.values - tc) <= h).any()) for _, _, df in es.positions.values() if len(df)):
            return o
    return None


def inject_unit(es_list, rung, f_uJy, rng, stat_needed):
    """One draw: per event a random centre, injection, statistic.
    Returns (S_event, S_pulse, S_stack) maxed over z."""
    per_ev = []
    for es in es_list:
        h = es.half[rung]
        o = draw_offset(rng, es, h)
        if o is None:
            continue
        tc = es.t_ca + o
        inj = {}
        for pos, (ra, de, df) in es.positions.items():
            if len(df) == 0:
                continue
            inw = np.abs(df.MJD.values - tc) <= h
            if not inw.any():
                continue
            R = es.response(pos, df.MJD.values[inw], df.fwhm_arcsec.values[inw], tc, Z_INJ)
            inj[pos] = dict(zip(df.Obs.values[inw], f_uJy * R))
        es.inject = inj
        per_z = []
        for z in sc.Z_GRID:
            ws = es.window_stats(rung, tc, z)
            per_z.append(ws)
        es.inject = None
        per_ev.append(per_z)
    if not per_ev:
        return None
    best = {'S_event': -np.inf, 'S_pulse': -np.inf, 'S_stack': -np.inf}
    for iz in range(len(sc.Z_GRID)):
        wins = [pe[iz] for pe in per_ev if pe[iz] is not None]
        if not wins:
            continue
        best['S_event'] = max(best['S_event'], max(w['S_w'] for w in wins))
        best['S_pulse'] = max(best['S_pulse'], max(w['S_pulse_w'] for w in wins))
        best['S_stack'] = max(best['S_stack'], sum(w['num'] for w in wins) / np.sqrt(sum(w['den'] for w in wins)))
    return best


def m90_from_curve(mags, rec):
    mags, rec = np.asarray(mags), np.asarray(rec)
    ok = np.where(rec >= 0.9)[0]
    if len(ok) == 0:
        return None
    i = ok.max()
    if i + 1 < len(mags) and rec[i + 1] < 0.9:
        return float(mags[i] + (rec[i] - 0.9) / (rec[i] - rec[i + 1]) * (mags[i + 1] - mags[i]))
    return float(mags[i])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--only', nargs='*')
    ap.add_argument('--draws', type=int, default=N_DRAW)
    args = ap.parse_args()
    fz = sc.FREEZE
    gate = None
    try:
        mpc = json.load(open(os.path.join(RES, 'mpc_control_v1.json')))
        gate = mpc['bands'].get('o', {})
    except FileNotFoundError:
        pass
    conf = json.load(open(os.path.join(RES, 'confirmatory_search_v1.json')))
    dev = json.load(open(os.path.join(RES, 'dev_search_v1.json')))
    T_by_unit = {}
    for run in (conf, dev):
        for u in run['units']:
            if u['computed']:
                T_by_unit[(u['channel'], u['rung'], u['target_id'])] = {s: r['T'] for s, r in u['thresholds'].items()}
    out, power = {}, {}
    for ui, u in enumerate(fz['search_units']):
        if u['class'] != 'searched' or (args.only and u['target_id'] not in args.only):
            continue
        key = f"{u['channel']}|{u['rung']}|{u['target_id']}"
        stats = list(u['statistics'])
        T = T_by_unit[(u['channel'], u['rung'], u['target_id'])]
        es_list = [sc.EventSeries(u['channel'], u['target_id'], eid, u['band']) for eid in u['covered_event_ids']]
        for es in es_list:
            es.fast = True
        rng = np.random.default_rng(SEED + ui)
        curves = {s: [] for s in stats}
        mags_done = []
        zero_run = 0
        for m in MAGS:
            f_uJy = 10 ** (-0.4 * (m - 23.9))
            hits = {s: 0 for s in stats}
            n_ok = 0
            for _ in range(args.draws):
                b = inject_unit(es_list, u['rung'], f_uJy, rng, stats)
                if b is None:
                    continue
                n_ok += 1
                for s in stats:
                    thr = max(T[s], 0.0) if T[s] is not None else 0.0
                    hits[s] += int(b[s] > thr)
            mags_done.append(float(m))
            for s in stats:
                curves[s].append(hits[s] / max(n_ok, 1))
            print(key, m, {s: round(curves[s][-1], 2) for s in stats}, flush=True)
            if all(curves[s][-1] == 0.0 for s in stats):
                zero_run += 1
                if zero_run >= 2:
                    break
            else:
                zero_run = 0
        rec = {'unit': key, 'split': u['split'], 'statistics': {}, 'mags': mags_done, 'draws': args.draws, 'z_inj_au': Z_INJ}
        power[key] = {}
        for s in stats:
            m90 = m90_from_curve(mags_done, curves[s])
            rec['statistics'][s] = {'T': T[s], 'recovery': [round(x, 3) for x in curves[s]], 'm90_o_AB': m90}
            if m90 is not None:
                P = power_W(u['channel'], u['rung'], u['target_id'], m90)
                power[key][s] = {'m90_o_AB': round(m90, 2), 'F_line_W_m2': f_line(m90),
                                 ('P_cone_W' if u['channel'] == 'B' else 'P_tx_10m_W'): P}
        out[key] = rec
        with open(os.path.join(RES, 'completeness_v1.json'), 'w') as f:
            json.dump({'run_utc': datetime.now(timezone.utc).isoformat(), 'seed': SEED, 'freeze_hash': fz['freeze_content_hash'],
                       'mpc_gate_o': gate, 'band': {'lambda_c_nm': 690, 'dlambda_nm': 260, 'dnu_Hz': DNU},
                       'units': out}, f, indent=1, default=float)
        with open(os.path.join(RES, 'power_limits_v1.json'), 'w') as f:
            json.dump(power, f, indent=1, default=float)
    print('done')


if __name__ == '__main__':
    main()
