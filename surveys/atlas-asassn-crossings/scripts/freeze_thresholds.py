"""ATLAS + ASAS-SN crossings threshold freeze v1.0 (hypotheses v1.0 +
§12 amendments; user-approved recommendations 2026-09-04).

Declares, before any in-window statistic is formed, the statistics,
detrending, control constructions, per-statistic gates, threshold
rule, unit population, trials accounting, split and seed, bound to
content hashes of the frozen inputs. Conventions inherited from the
ZTF/PS1/PTF crossings freezes with the ATLAS substitutions.

Gates (approved 2026-09-04):
  S_pulse, S_event : every unit with >= 1 covered window (o band)
  S_stack          : units with >= 3 covered windows
  a statistic is a searched trial only if >= 4 of its controls are
  valid (mirror-gated pseudo-windows / pseudo-stacks); otherwise it is
  constraint-only for that unit. c band = chromatic annotation only.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.table import Table

REPO = Path(__file__).resolve().parents[3]
D = REPO / 'surveys' / 'atlas-asassn-crossings'
INPUTS = {
    'hypotheses': D / 'hypotheses.md',
    'task_list': D / 'results' / 'task_list_v1.ecsv',
    'coverage': D / 'results' / 'coverage_v1_events.ecsv',
    'asassn_ledger': D / 'results' / 'asassn_ledger_v1_events.ecsv',
    'saturation': D / 'results' / 'saturation_cut_v1.json',
}
SEED = 20260904
DEV_TARGETS = ('wolf-359', 'gj-908')
OFFSETS = (-97, -71, -47, -23, 23, 47, 71, 97)
MIN_STACK_EVENTS = 3
MIN_VALID_CONTROLS = 4
SEARCH_BAND = 'o'
RUNGS = {'B': ('1.2Rsun', '2.5Rsun', '0.1AU'), 'A': ('0.1AU',)}


def sha(p: Path) -> str:
    return 'sha256:' + hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    cov = Table.read(INPUTS['coverage'])
    sat = json.loads(INPUTS['saturation'].read_text())['targets']
    cov = cov[cov['pos_role'] == 'tca']
    units = []
    for ch in ('B', 'A'):
        for rung in RUNGS[ch]:
            s = cov[(cov['channel'] == ch) & (cov['rung'] == rung)]
            for tid in sorted(set(map(str, s['target_id']))):
                st = s[s['target_id'] == tid]
                # covered = >= 1 FAQ-passing o exposure in the real window
                # (B 0.1 AU: the mini-track union carries the same set at
                # the tca position plus the ends; the tca-position count
                # gates, the union is the searched epoch set)
                cov_ev = st[np.asarray(st[f'n_{SEARCH_BAND}']) > 0]
                n_cov = len(cov_ev)
                pw = np.array([[int(r[f'pw{i}_{SEARCH_BAND}']) for i in range(8)] for r in cov_ev]) \
                    if n_cov else np.zeros((0, 8), int)
                # control validity, mirror-gated:
                #  S_pulse / S_event: offset k valid iff >= 1 covered event
                #    has data at offset k (max over events, like the real
                #    statistic's max over covered windows)
                #  S_stack: pseudo-stack k valid iff >= MIN_STACK_EVENTS
                #    events have data at offset k
                per_offset_events = pw.sum(axis=0) if n_cov else np.zeros(8, int)
                valid_event = int((per_offset_events >= 1).sum())
                valid_stack = int((per_offset_events >= MIN_STACK_EVENTS).sum())
                stats, notes = {}, []
                sat_status = None
                if ch == 'A':
                    sat_status = sat[tid][SEARCH_BAND]['status']
                if n_cov >= 1 and sat_status != 'excluded':
                    for name in ('S_pulse', 'S_event'):
                        if valid_event >= MIN_VALID_CONTROLS:
                            stats[name] = {'n_valid_controls': valid_event,
                                           'p_crossing': round(1.0 / (valid_event + 1), 4)}
                        else:
                            notes.append(f'{name}: {valid_event} valid controls < {MIN_VALID_CONTROLS} -> constraint-only')
                    if n_cov >= MIN_STACK_EVENTS:
                        if valid_stack >= MIN_VALID_CONTROLS:
                            stats['S_stack'] = {'n_valid_controls': valid_stack,
                                                'p_crossing': round(1.0 / (valid_stack + 1), 4)}
                        else:
                            notes.append(f'S_stack: {valid_stack} valid pseudo-stacks < {MIN_VALID_CONTROLS} -> constraint-only')
                    else:
                        notes.append(f'S_stack: {n_cov} covered windows < {MIN_STACK_EVENTS} -> not formed')
                if sat_status == 'excluded':
                    cls = 'saturation_excluded'
                elif n_cov == 0:
                    cls = 'uncovered'
                elif stats:
                    cls = 'searched'
                else:
                    cls = 'constraint_only'
                units.append({
                    'target_id': tid, 'channel': ch, 'rung': rung, 'band': SEARCH_BAND,
                    'split': 'dev' if tid in DEV_TARGETS else 'confirmatory',
                    'class': cls, 'saturation_status': sat_status,
                    'n_events': len(st), 'n_covered_events': n_cov,
                    'in_window_exposures_o': int(np.sum(st[f'n_{SEARCH_BAND}'])),
                    'in_window_exposures_c': int(np.sum(st['n_c'])),
                    'covered_event_ids': [str(e) for e in cov_ev['event_id']],
                    'per_offset_covered_events': per_offset_events.tolist(),
                    'statistics': stats, 'n_trials': len(stats),
                    'expected_control_crossings': round(sum(v['p_crossing'] for v in stats.values()), 4),
                    'notes': notes,
                })
    searched = [u for u in units if u['class'] == 'searched']
    n_trials = sum(u['n_trials'] for u in searched)
    exp_cross = round(sum(u['expected_control_crossings'] for u in searched), 3)
    dev_trials = sum(u['n_trials'] for u in searched if u['split'] == 'dev')

    freeze = {
        'freeze_version': 'atlas-asassn-crossings-thresholds-v1.0',
        'frozen_at': '2026-09-04',
        'hypothesis_version': 'atlas-asassn-crossings-hypotheses-v1.0 + §12 amendments A1-A3',
        'seed': SEED,
        'input_hashes': {k: sha(p) for k, p in INPUTS.items()},
        'substrate': {
            'ATLAS': 'server-calibrated tphot difference-flux forced photometry (uJy, AB) at fixed '
                     'positions, windowed +/-110 d per event (task_list_v1); o band searched, c band '
                     'chromatic annotation only (3/44 and 10/55 grazing windows covered in c)',
            'ASAS-SN_v2': 'coverage-fraction ledger only (asassn_ledger_v1); no v2 photometry enters any '
                          'statistic at this freeze; channel-A on-star v2 series deferred to a PM-matched '
                          'follow-up (nearest-source matches are epoch-offset for the high-PM stars)',
            'observer': 'Earth-center universal list; site-vs-geocenter <= R_earth = 0.8 % of the tightest '
                        'rung, declared budget term',
        },
        'quality': {
            'primary_mask': 'FAQ recipe (atlas_api.faq_quality_mask): duJy<10000, err==0, 100<x,y<10460, '
                            '1.6<maj,min<5, -1<apfit<-0.1, mag5sig>17, Sky>17; H-filter rows dropped',
            'strict_mask': 'primary + mag5sig>18.5 + |MJD - template step (58417, 58882)| > 30 d; '
                           'strict re-run reported for every exceedance',
        },
        'detrending': {
            'rule': 'D5: per position and band, residual r_i = f_i - running median of off-window '
                    'primary-mask epochs within +/-15 d of t_i (30 d window; min 5 epochs, else the '
                    'off-window global median); in-window epochs never enter the median',
            'variance': 'v_i = k sigma_i^2 (duJy), k = median(r^2/sigma^2)/0.4549 over off-window '
                        'epochs after 3x3-sigma clipping, floor 1, per (position, band) - the '
                        'established recipe',
        },
        'statistics': {
            'unit': '(target, channel, rung, o band); channel B 0.1 AU unit epochs = union over the '
                    'D3 mini-track positions with per-epoch position selection (below)',
            'window': 'flat chord t_ca +/- sqrt(r^2 - b^2)/v_perp (task_list_v1 half-windows)',
            'weights': 'inverse variance w_i = R_i / v_i; R_i = response factor exp(-d_i^2 / 2 s_i^2), '
                       'd_i = offset between the task position and the predicted apparent source '
                       'position at epoch i (B: relay at z on the anti-star axis; A: star), '
                       's_i = FWHM_i/2.355 from the epoch maj/min (arcsec, 1.86"/px); epochs with '
                       'R_i < 0.2 dropped; B evaluated over the nested z family '
                       '(550, 1000, 2500, 5500, 10000 AU) choosing per epoch the mini-track position '
                       'of largest R_i, statistic = max over z; controls take the same max',
            'S_pulse': 'tertiary (pulse cell): max over in-window epochs of r_i sqrt(R_i) / sqrt(v_i); '
                       'one-sided positive',
            'S_event': 'secondary: max over covered windows of S_w = sum(w_i r_i) / sqrt(sum(w_i^2 v_i)); '
                       'one-sided positive',
            'S_stack': 'primary (the recurrence cell): the same weighted sum over all covered windows '
                       'of the unit pooled, S = sum_w sum(w r) / sqrt(sum_w sum(w^2 v)); requires '
                       f'>= {MIN_STACK_EVENTS} covered windows',
            'gates': {'S_pulse_S_event': '>= 1 covered window (>= 1 primary-mask o exposure)',
                      'S_stack': f'>= {MIN_STACK_EVENTS} covered windows',
                      'controls': f'>= {MIN_VALID_CONTROLS} valid controls per statistic, else '
                                  'constraint-only for that statistic'},
            'A_1.0AU': 'constraint-only, no ATLAS pull (hypotheses §3, standing window~season theorem); '
                       'coverage claim deferred to a ledger addendum',
            'A_saturation': 'D2 table results/saturation_cut_v1.json on the reduced-mode 20-d series: '
                            'excluded units carry no statistic; marginal units are searched and flagged',
        },
        'controls': {
            'kind': 'temporal pseudo-windows at the same fixed position (both channels)',
            'designated_offsets_days': list(OFFSETS),
            'redraws': 'none available inside the +/-110 d pull; overlap with a real window of any '
                       'rung is geometrically impossible (|offset| - half >= 17 d)',
            'validity': 'mirror-gated per statistic: for S_pulse/S_event an offset is a valid control '
                        'iff >= 1 covered event of the unit has primary-mask o data in that '
                        'pseudo-window (the control statistic takes the same max over events); for '
                        f'S_stack a pseudo-stack k is valid iff >= {MIN_STACK_EVENTS} events have data '
                        'at offset k; the same-offset pseudo-windows of all covered events form one '
                        'pseudo-stack, chord-weighted identically',
            'threshold_rule': 'T = max one-sided control statistic over the valid controls, per '
                              'statistic; exceedance S > max(T, 0); margin S - T reported',
            'per_trial_crossing_probability': '1/(n_valid + 1), tallied from measured validity - '
                                              'not the nominal 1/9',
        },
        'trials': {
            'n_searched_units': len(searched), 'n_trials': n_trials,
            'n_trials_dev': dev_trials, 'n_trials_confirmatory': n_trials - dev_trials,
            'expected_control_crossings': exp_cross,
            'disposition': 'every exceedance individually adjudicated by the veto ladder; no silent drops',
        },
        'veto_ladder': [
            '1. SkyBoT known-object census at every exceedance epoch (antipode fields sit in the '
            'opposition asteroid stream)',
            '2. intra-quad consistency: a window-locked signal must be present across the night\'s '
            'quad (a main-belt mover crosses the PSF in <~ 15 min)',
            '3. template-step annotation: exceedances within 30 d of MJD 58417 / 58882 re-run on the '
            'reduced-mode series (strict mask)',
            '4. recurrence: any S_event / S_pulse exceedance tested on the unit\'s other covered windows',
            '5. rate test vs the z-track prediction where the mini-track resolves it (B 0.1 AU)',
            'annotations (never sole grounds): c-band chromatic consistency where c epochs exist',
        ],
        'split': {'rule': 'D8 as approved 2026-09-04: dev = wolf-359 (D1 forced) + gj-908 (widest b, '
                          'thinnest A coverage); all other targets confirmatory; every rung and channel '
                          'of a dev target is dev',
                  'dev_targets': list(DEV_TARGETS), 'seed': SEED},
        'completeness': {
            'method': 'D6 response-model injections into the real off-window residual series of each '
                      'searched unit: signal uJy x R_i(offset, FWHM) with the chord temporal profile, '
                      '>= 100 draws per unit and statistic, seed 20260904; recovery = exceedance of the '
                      'frozen T',
            'positive_control': '(15000) CCD, MPC mode, MJD 60310-61286, through the identical mask -> '
                                'detrend -> statistic chain vs JPL Horizons predicted V (o via the '
                                'asteroid solar-colour offset); gate |residual| <= 0.2 mag, else a '
                                'small image-request subsample (D6)',
            'depth_honesty': 'grazing-rung stacks hold ~18-30 exposures (van-maanen 1.2 R_sun 18 over 4 '
                             'windows), ~20.6 stacked; 0.1 AU stacks ~21.5; per-window kW-class limits '
                             'through the grazing cones, not the hypotheses §9 "few hundred W" - stated '
                             'before search',
        },
        'alpha': {'fwer_alpha': 0.05, 'ks_alpha': 0.01},
        'search_units': units,
    }
    body = json.dumps(freeze, indent=1, sort_keys=True)
    freeze['freeze_content_hash'] = 'sha256:' + hashlib.sha256(body.encode()).hexdigest()
    out = D / 'configs' / 'threshold_freeze_v1.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(freeze, indent=1, sort_keys=True) + '\n')
    for u in units:
        print(f"{u['channel']} {u['rung']:7s} {u['target_id']:11s} [{u['split'][:4]}] {u['class']:20s} "
              f"cov {u['n_covered_events']:2d}/11 exp_o {u['in_window_exposures_o']:3d} "
              f"offsets {u['per_offset_covered_events']} "
              f"{'+'.join(f'{k}({v['n_valid_controls']})' for k, v in u['statistics'].items())} "
              f"{'; '.join(u['notes'])}")
    print(json.dumps(freeze['trials'], indent=1), freeze['freeze_content_hash'])


if __name__ == '__main__':
    main()
