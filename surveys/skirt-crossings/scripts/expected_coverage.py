"""Expected in-beam / out-of-beam nights per target, rung and side from
the geometry curves x the probe-measured ATLAS night efficiency vs
solar elongation (nights observed per available day, 5-deg bins,
both probes pooled). Screening only: the coverage stage measures the
real counts. Output: results/expected_coverage_v1.json."""
import json
import os

import numpy as np

D = os.path.join(os.path.dirname(__file__), '..', 'results')
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
GEOM = os.path.join(REPO, 'runs', 'skirt-crossings', 'geometry')
RUNGS = [0.85, 0.9, 0.95]


def main():
    probes = json.load(open(os.path.join(D, 'probe_elongation_v1.json')))
    scope = json.load(open(os.path.join(D, 'target_scope_v1.json')))['targets']
    geom = json.load(open(os.path.join(D, 'skirt_geometry_v1.json')))['targets']
    # efficiency(eps bin) = ok nights / available days, pooled over probes
    edges = np.arange(0, 181, 5)
    nights = np.zeros(len(edges) - 1); days = np.zeros(len(edges) - 1)
    for name, p in probes.items():
        h = p['hist_5deg_ok']
        # nights ~ ok rows / 4 (quad) is crude; use the per-bin night counts if present
        npy = p.get('nights_per_eps_bin_per_year', {})
        nyr = (p['mjd_range'][1] - p['mjd_range'][0]) / 365.25
        g = np.load(os.path.join(GEOM, 'wolf-359.npz'))   # beta ~ 0 like the probes
        for i, a in enumerate(edges[:-1]):
            key = f'{a}-{a + 5}'
            if key in npy:
                nights[i] += npy[key] * nyr
                days[i] += ((g['eps'] >= a) & (g['eps'] < a + 5)).sum() * (nyr / ((g['mjd'][-1] - g['mjd'][0]) / 365.25))
    eff = np.where(days > 0, nights / np.maximum(days, 1), np.nan)
    eff_tab = {f'{a}-{a + 5}': (None if np.isnan(e) else round(float(e), 3)) for a, e in zip(edges[:-1], eff)}
    out = {'efficiency_by_eps_bin': eff_tab, 'targets': {}}
    for tid, t in scope.items():
        if tid not in geom or t['o_screen'] not in ('ok_est', 'marginal_est'):
            continue
        g = np.load(os.path.join(GEOM, f'{tid}.npz'))
        eps, b_e, mjd = g['eps'], g['b_e'], g['mjd']
        nyr = (mjd[-1] - mjd[0]) / 365.25
        idx = np.clip(((eps) // 5).astype(int), 0, len(eff) - 1)
        e = np.nan_to_num(eff[idx])
        rec = {'ecl_lat_deg': geom[tid]['ecl_lat_deg'], 'o_est': t['o_est'], 'rungs': {}}
        for rb in RUNGS:
            inb = b_e < rb
            s2, a = eps < 90, eps > 90
            rec['rungs'][str(rb)] = {
                'skirt_nights_per_yr': round(float(e[inb & s2].sum() / nyr), 1),
                'opp_nights_per_yr': round(float(e[inb & a].sum() / nyr), 1),
                'out_nights_per_yr': round(float(e[~inb].sum() / nyr), 1),
                'skirt_days_per_yr': round(float((inb & s2).sum() / nyr), 1),
            }
        out['targets'][tid] = rec
    json.dump(out, open(os.path.join(D, 'expected_coverage_v1.json'), 'w'), indent=1)
    print('efficiency:', eff_tab)
    print(f"{'target':12s} {'beta':>6s} {'o_est':>5s} | 0.85 skirt/opp/out | 0.90 skirt/opp/out | 0.95 skirt/opp/out  (nights/yr)")
    for tid, r in sorted(out['targets'].items(), key=lambda kv: abs(kv[1]['ecl_lat_deg'])):
        cells = ' | '.join(f"{r['rungs'][str(rb)]['skirt_nights_per_yr']:5.1f}/{r['rungs'][str(rb)]['opp_nights_per_yr']:5.1f}/{r['rungs'][str(rb)]['out_nights_per_yr']:5.1f}" for rb in RUNGS)
        print(f"{tid:12s} {r['ecl_lat_deg']:6.1f} {r['o_est'] or 0:5.2f} | {cells}")


if __name__ == '__main__':
    main()
