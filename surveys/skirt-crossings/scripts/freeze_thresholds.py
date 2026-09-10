"""Outer-skirt threshold freeze v1.0 (hypotheses v1.0, user-approved
D1-D8 2026-09-07). Declares, before any in-beam statistic is formed,
the unit population, per-statistic searched/constraint-only class
from the coverage counts, the trials accounting, the split and seed,
bound to content hashes of the frozen inputs. Reads only
coverage_v1.json (counts) - no Delta is formed here.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))
import search_core as sc  # noqa: E402

D = Path(sc.D)
INPUTS = {'hypotheses': D / 'hypotheses.md', 'task_list': D / 'results' / 'task_list_v1.ecsv',
          'coverage': D / 'results' / 'coverage_v1.json', 'control_pool': D / 'results' / 'control_pool_v1.json',
          'geometry': D / 'results' / 'skirt_geometry_v1.json', 'search_core': Path(__file__).with_name('search_core.py')}
SEED = 20260907
DEV = ('teegarden', 'gj-2012')
MAGS = [float(m) for m in range(15, 24)] + [x + 0.5 for x in range(15, 23)]


def sha(p):
    return 'sha256:' + hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument('--tag', default='v1'); args = ap.parse_args()
    cov = json.load(open(INPUTS['coverage']))
    units = []
    for tid, rec in cov['targets'].items():
        if rec.get('status') == 'not_pulled':
            continue
        sat = rec['saturation']['o']['status']
        for rb in sc.RUNGS:
            r = rec['rungs'][rb]
            stats, notes = {}, []
            for stat in sc.STATS:
                q = r[stat]
                if q['searched']:
                    stats[stat] = {'n_valid_controls': q['n_valid_controls'], 'p_crossing': round(1 / (q['n_valid_controls'] + 1), 4),
                                   'n_valid_cycles_target': r['target_validity'][stat]['n_valid_cycles']}
                else:
                    why = ('ledger_only' if r['ledger_only'] else 'saturation_excluded' if sat == 'excluded'
                           else 'target_coverage' if not q['target_valid'] else f"{q['n_valid_controls']} valid controls < {sc.MIN_CONTROLS}")
                    notes.append(f'{stat}: {why}')
            cls = 'ledger_only' if r['ledger_only'] else 'saturation_excluded' if sat == 'excluded' else ('searched' if stats else 'constraint_only')
            units.append({'target_id': tid, 'rung_au': rb, 'band': 'o', 'split': 'dev' if tid in DEV else 'confirmatory',
                          'class': cls, 'saturation_status': sat, 'o_median': rec['saturation']['o']['median_mag'],
                          'n_nights_o': rec['n_nights_o'], 'statistics': stats, 'n_trials': len(stats),
                          'expected_control_crossings': round(sum(v['p_crossing'] for v in stats.values()), 4), 'notes': notes})
    searched = [u for u in units if u['class'] == 'searched']
    freeze = {
        'freeze_version': f'skirt-crossings-thresholds-{args.tag}',
        'frozen_at': datetime.now(timezone.utc).isoformat(),
        'hypothesis_version': 'skirt-crossings-hypotheses-v1.0 (D1-D8 approved 2026-09-07)',
        'seed': SEED, 'input_hashes': {k: sha(p) for k, p in INPUTS.items()},
        'substrate': 'ATLAS reduced-mode (target-image tphot) full-history forced photometry with server-side proper motion; o searched, c annotation',
        'sites': {k: v[0] for k, v in sc.SITES.items()},
        'geometry': {'floor_deg': cov['floor_deg'], 'rungs_au': list(sc.RUNGS), 'in_beam': 'r_E sin(eps) < r_b (both sides)',
                     'cycle': 'conjunction +/- 182.625 d', 'elongation_source': 'target curve for target and controls'},
        'quality': {'primary_mask': 'FAQ recipe + site prefix in {01,02,03,04}', 'response_gate': f'R >= {sc.R_GATE}',
                    'epoch': 'nightly-unit median of >= 2 exposures', 'airmass_layer': 'f = a + b (X-1) over all nightly epochs, per star and band'},
        'statistics': {'S_sym': 'sum(Delta_y/v_y)/sqrt(sum 1/v_y), IN both sides', 'S_skirt': 'same, S2 side only',
                       'S_year': 'max_y Delta_y/sqrt(v_y)', 'S_opp': 'annotation, opposition side only'},
        'gates': {'cycle_valid': f'n_in >= {sc.MIN_IN} (skirt {sc.MIN_SKIRT}) and n_out >= {sc.MIN_OUT}', 'min_cycles': sc.MIN_CYCLES,
                  'min_valid_controls': sc.MIN_CONTROLS, 'control_validity': 'same gate on the control series (mirror-gated)'},
        'threshold_rule': 'T = max over valid controls; exceedance S > max(T, 0); p_crossing = 1/(n_valid+1)',
        'controls': 'null ensemble: 8 field stars per target drawn with seed 20260907 from control_pool_v1 (task_list_v1)',
        'split': {'dev': list(DEV), 'confirmatory': sorted({u['target_id'] for u in units if u['split'] == 'confirmatory'})},
        'completeness': {'method': 'top-hat step of flux f x R_i added to every IN exposure (both sides) before normalisation; '
                                   'frozen chain; 100 draws per magnitude with 20 % random nightly-epoch dropout (pre-data '
                                   'clarification of hypotheses §8: adding photon noise to a real series would double-count it); '
                                   'm90 = faintest grid magnitude with >= 90 % of draws exceeding max(T, 0)',
                         'mag_grid': sorted(MAGS), 'dropout': 0.2, 'n_draw': 100,
                         'power': 'P = F_line pi r_b^2, o-band line 690 nm / 260 nm width, +/-0.1 mag scale (inherited ATLAS D6)'},
        'units': units,
        'population': {'n_units': len(units), 'n_searched_units': len(searched),
                       'n_trials': sum(u['n_trials'] for u in searched),
                       'n_trials_dev': sum(u['n_trials'] for u in searched if u['split'] == 'dev'),
                       'n_trials_confirmatory': sum(u['n_trials'] for u in searched if u['split'] == 'confirmatory'),
                       'expected_control_crossings': round(sum(u['expected_control_crossings'] for u in searched), 3),
                       'classes': {c: sum(1 for u in units if u['class'] == c) for c in ('searched', 'constraint_only', 'ledger_only', 'saturation_excluded')}},
    }
    body = json.dumps({k: v for k, v in freeze.items() if k != 'freeze_content_hash'}, sort_keys=True).encode()
    freeze['freeze_content_hash'] = 'sha256:' + hashlib.sha256(body).hexdigest()
    out = D / 'configs' / f'threshold_freeze_{args.tag}.json'
    out.write_text(json.dumps(freeze, indent=1))
    print(json.dumps(freeze['population'], indent=1))
    for u in units:
        print(u['target_id'], u['rung_au'], u['class'], u['n_trials'], u['notes'])


if __name__ == '__main__':
    main()
