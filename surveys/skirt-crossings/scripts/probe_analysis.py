"""Measure the ATLAS solar-elongation floor and low-elongation cadence
from the recon probes (runs/skirt-crossings/probe/*.txt) ->
results/probe_elongation_v1.json."""
import glob
import json
import os
import sys

import numpy as np
from astropy.coordinates import SkyCoord, get_sun, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'atlas-asassn-crossings', 'scripts'))
import atlas_api  # noqa: E402

PROBE = os.path.join(atlas_api.REPO, 'runs', 'skirt-crossings', 'probe')
OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'probe_elongation_v1.json')
# Obs-prefix -> site, VERIFIED from the probe (star/Sun altitudes at each
# epoch): 03 = Sutherland (South Africa), 04 = El Sauce (Chile) -- the
# 2026-08-25 ATLAS recon note had these two swapped.
SITES = {'01': ('Haleakala', 20.7076, -156.2569, 3055), '02': ('Mauna Loa', 19.5362, -155.5763, 3397),
         '03': ('Sutherland', -32.3796, 20.8107, 1798), '04': ('El Sauce', -30.4708, -70.7647, 1600)}


def main():
    out = {}
    for f in sorted(glob.glob(os.path.join(PROBE, '*.txt'))):
        name = os.path.basename(f)[:-4]
        sc = json.load(open(f + '.json'))
        q = sc['query_params']
        df = atlas_api.read_result(f)
        star = SkyCoord(q['ra'] * u.deg, q['dec'] * u.deg)
        t = Time(df.MJD.values, format='mjd')
        df['eps'] = get_sun(t).separation(star).deg
        df['unit'] = df.Obs.str[:2]
        ok = atlas_api.faq_quality_mask(df)
        df['ok'] = ok.values
        # airmass at each site
        am = np.full(len(df), np.nan)
        for code, (nm, lat, lon, h) in SITES.items():
            m = (df.unit == code).values
            if m.any():
                loc = EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=h * u.m)
                aa = star.transform_to(AltAz(obstime=t[m], location=loc))
                am[m] = aa.secz.value
        df['airmass'] = am
        g = df[df.ok]
        rec = {'query': q, 'wall_s': sc['wall_s'], 'n_rows': int(len(df)), 'n_ok': int(len(g)),
               'mjd_range': [float(df.MJD.min()), float(df.MJD.max())],
               'eps_min_all': float(df.eps.min()), 'eps_min_ok': float(g.eps.min()),
               'eps_p01_ok': float(np.percentile(g.eps, 1)), 'eps_p05_ok': float(np.percentile(g.eps, 5)),
               'per_unit': {}, 'hist_5deg_ok': {}, 'nights_per_eps_bin_per_year': {}, 'precision': {}}
        for code in sorted(set(df.unit)):
            d = g[g.unit == code]
            if len(d):
                rec['per_unit'][code] = {'name': SITES.get(code, ('?',))[0], 'n_ok': int(len(d)),
                                        'eps_min': float(d.eps.min()), 'eps_p05': float(np.percentile(d.eps, 5)),
                                        'airmass_median_eps_lt_70': float(np.nanmedian(d.airmass[d.eps < 70])) if (d.eps < 70).any() else None}
        edges = np.arange(0, 181, 5)
        h, _ = np.histogram(g.eps, bins=edges)
        rec['hist_5deg_ok'] = {f'{a}-{a + 5}': int(n) for a, n in zip(edges[:-1], h)}
        nyr = (df.MJD.max() - df.MJD.min()) / 365.25
        for a in range(0, 180, 5):
            d = g[(g.eps >= a) & (g.eps < a + 5)]
            nights = len(set(np.floor(d.MJD + 0.5)))
            rec['nights_per_eps_bin_per_year'][f'{a}-{a + 5}'] = round(nights / nyr, 1)
        for band in ('o', 'c'):
            d = g[g.F == band]
            if len(d) > 20:
                med = np.median(d.uJy)
                rec['precision'][band] = {'n': int(len(d)), 'median_mag': float(np.median(d.m)),
                                         'frac_scatter_mad': float(1.4826 * np.median(np.abs(d.uJy / med - 1))),
                                         'median_duJy_over_flux': float(np.median(d.duJy / med)),
                                         'frac_scatter_eps_lt_70': float(1.4826 * np.median(np.abs(d.uJy[d.eps < 70] / med - 1))) if (d.eps < 70).sum() > 10 else None,
                                         'frac_scatter_eps_gt_120': float(1.4826 * np.median(np.abs(d.uJy[d.eps > 120] / med - 1))) if (d.eps > 120).sum() > 10 else None}
        out[name] = rec
        print(name, json.dumps({k: rec[k] for k in ('n_rows', 'n_ok', 'eps_min_ok', 'eps_p01_ok', 'eps_p05_ok', 'wall_s')}))
        print(' nights/yr per eps bin:', rec['nights_per_eps_bin_per_year'])
        print(' per unit:', rec['per_unit'])
        print(' precision:', rec['precision'])
    json.dump(out, open(OUT, 'w'), indent=1)


if __name__ == '__main__':
    main()
