"""Coverage / cut stage (hypotheses §4 D2-D3, §7 gate; no in-beam
statistic is formed here - counts only).

Per target: measured reduced-mode o/c brightness and the ATLAS D2
saturation rule; G1 proper-motion application check (forced position
vs the Gaia-propagated track at the first/last 30 epochs, both
mu_alpha readings); nightly-epoch counts per conjunction cycle and
rung for the target and its controls; validity per statistic; the
1 AU-rung coverage ledger (nights per S2 cycle at eps >= 50 deg, and
the A-side nights). Output: results/coverage_v1.json.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import search_core as sc  # noqa: E402

SAT, MARGIN = 12.5, 0.5
FLOOR = 50.0
LEDGER_ONLY = {'wolf-1069': ['0.9', '0.95'], 'gj-293': ['0.9', '0.95'], 'gj-13157': ['0.9'], 'gj-3112': ['0.9']}


def saturation(df):
    if not len(df):
        return {'n': 0, 'median_mag': None, 'status': 'no_data'}
    pos = df[df.uJy > 0]
    m = float(np.median(23.9 - 2.5 * np.log10(pos.uJy.values.astype(float)))) if len(pos) else None
    fail = 1 - len(df) / max(1, df.attrs.get('n_raw', len(df)))
    st = 'no_data' if m is None else ('excluded' if m < SAT else ('marginal' if m < SAT + MARGIN else 'ok'))
    return {'n': int(len(df)), 'median_mag': m, 'status': st, 'faq_fail_frac': float(fail)}


def pm_check(key):
    """Offsets between forced and propagated positions, both mu_alpha readings."""
    row = sc.task_row(key)
    df = sc.atlas_api.read_result(os.path.join(sc.LC, f'{key}.txt'))
    df = df[sc.atlas_api.faq_quality_mask(df)]
    if len(df) < 10:
        return None
    from astropy.coordinates import SkyCoord
    from astropy.time import Time
    import astropy.units as u
    dt = (df.MJD.values - Time(float(row['epoch_year']), format='jyear').mjd) / 365.25
    ra0, de0 = float(row['ra_deg']), float(row['dec_deg'])
    out = {}
    for reading, fac in (('pmra_cosdec', 1 / np.cos(np.radians(de0))), ('pmra_raw', 1.0)):
        ra_p = ra0 + float(row['pm_ra_cosdec_mas_yr']) / 3.6e6 * dt * fac
        de_p = de0 + float(row['pm_dec_mas_yr']) / 3.6e6 * dt
        d = SkyCoord(df.RA.values * u.deg, df.Dec.values * u.deg).separation(SkyCoord(ra_p * u.deg, de_p * u.deg)).arcsec
        o = np.argsort(df.MJD.values)
        out[reading] = {'first30_median_arcsec': float(np.median(d[o[:30]])), 'last30_median_arcsec': float(np.median(d[o[-30:]])),
                        'all_median_arcsec': float(np.median(d))}
    # does the forced position move at all?
    o = np.argsort(df.MJD.values)
    span = SkyCoord(df.RA.values[o[:30]].mean() * u.deg, df.Dec.values[o[:30]].mean() * u.deg).separation(
        SkyCoord(df.RA.values[o[-30:]].mean() * u.deg, df.Dec.values[o[-30:]].mean() * u.deg)).arcsec
    out['forced_position_span_arcsec'] = float(span)
    out['expected_span_arcsec'] = float(np.hypot(row['pm_ra_cosdec_mas_yr'], row['pm_dec_mas_yr']) / 1000 * (df.MJD.max() - df.MJD.min()) / 365.25)
    return out


def main():
    t = sc.tasks()
    targets = [str(x) for x in t['target_id'][t['role'] == 'target']]
    out = {'floor_deg': FLOOR, 'targets': {}}
    for tid in targets:
        key = f'T-{tid}'
        if not os.path.exists(os.path.join(sc.LC, f'{key}.txt')):
            out['targets'][tid] = {'status': 'not_pulled'}
            continue
        rec = {'split': str(sc.task_row(key)['split']), 'saturation': {}, 'pm_check': pm_check(key), 'rungs': {}, 'ledger_1au': {}}
        geom = sc.Geometry(tid)
        for band in ('o', 'c'):
            raw = sc.atlas_api.read_result(os.path.join(sc.LC, f'{key}.txt'))
            raw = raw[raw.F == band]
            df = sc.load_series(key, band)
            df.attrs['n_raw'] = len(raw)
            rec['saturation'][band] = saturation(df)
        df = sc.load_series(key, 'o')
        nt = sc.nightly(df)
        rec['n_nights_o'] = int(len(nt))
        # 1 AU ledger: nights per conjunction cycle on the S2 side at eps >= floor, and A-side nights
        if len(nt):
            eps = geom.eps_at(nt.mjd.values); cyc = geom.cycle_of(nt.mjd.values)
            for i in sorted(set(cyc[cyc >= 0])):
                m = cyc == i
                rec['ledger_1au'][int(i)] = {'conj_mjd': float(geom.conj[i]), 'conj_utc': sc.Time(geom.conj[i], format='mjd').iso[:10],
                                            'nights_s2_ge_floor': int((m & (eps < 90) & (eps >= FLOOR)).sum()),
                                            'nights_s2_lt_floor': int((m & (eps < 90) & (eps < FLOOR)).sum()),
                                            'min_eps_observed': float(eps[m].min()), 'nights_a_side': int((m & (eps > 90)).sum())}
        ckeys = [str(k) for k in t['task_key'][(t['target_id'] == tid) & (t['role'] == 'control')]]
        for rb in sc.RUNGS:
            cc = sc.cycle_counts(nt, geom, float(rb))
            r = {'target_counts': cc, 'target_validity': sc.validity_from_counts(cc), 'controls': {}, 'ledger_only': rb in LEDGER_ONLY.get(tid, [])}
            for k in ckeys:
                if os.path.exists(os.path.join(sc.LC, f'{k}.txt')):
                    cnt = sc.nightly(sc.load_series(k, 'o'))
                    r['controls'][k] = {'n_nights': int(len(cnt)), 'validity': sc.validity_from_counts(sc.cycle_counts(cnt, geom, float(rb)))}
                else:
                    r['controls'][k] = {'n_nights': 0, 'validity': None, 'status': 'not_pulled'}
            for stat in sc.STATS:
                nvc = sum(1 for c in r['controls'].values() if c['validity'] and c['validity'][stat]['valid'])
                tv = r['target_validity'][stat]['valid']
                r[stat] = {'n_valid_controls': nvc, 'target_valid': tv,
                           'searched': bool(tv and nvc >= sc.MIN_CONTROLS and not r['ledger_only'] and rec['saturation']['o']['status'] != 'excluded')}
            rec['rungs'][rb] = r
        out['targets'][tid] = rec
        print(tid, rec['saturation']['o']['status'], round(rec['saturation']['o']['median_mag'] or 0, 2), 'nights', rec['n_nights_o'],
              {rb: {s: (rec['rungs'][rb][s]['searched'], rec['rungs'][rb][s]['n_valid_controls']) for s in sc.STATS} for rb in sc.RUNGS},
              'pm', {k: round(v['all_median_arcsec'], 2) for k, v in rec['pm_check'].items() if isinstance(v, dict)} if rec['pm_check'] else None, flush=True)
    json.dump(out, open(os.path.join(sc.RES, 'coverage_v1.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
