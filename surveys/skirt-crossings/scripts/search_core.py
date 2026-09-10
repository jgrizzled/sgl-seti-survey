"""Frozen statistics engine for the outer-skirt search (hypotheses v1.0
§6-§8; threshold freeze configs/threshold_freeze_v1.json). Every
construction here is the freeze text made executable.

Per star (target or one of its 8 controls) and band:
  series    : full-history reduced-mode task, FAQ primary mask,
              response R_i = exp(-d^2/4 sigma^2) from the forced
              position vs the apparent position (PM + PARALLAX,
              amendment A3); flux corrected by 1/R_i; gate R_i >= 0.9, site from the Obs prefix
              (01 Haleakala, 02 Mauna Loa, 03 Sutherland, 04 El Sauce
              -- verified in the recon), airmass per exposure
  nightly   : one epoch per unit-night (local-midnight night id),
              median of >= 2 exposures, fractional flux
              f = F / median(F) over all good exposures
  airmass   : r = f - (a + b (X - 1)) fitted over all nightly epochs
  geometry  : the TARGET's solar elongation eps(t) defines every set
              (controls share the exposures); b_e = r_E sin eps;
              IN(r_b) = b_e < r_b (skirt eps < 90 and opposition
              eps > 90), OUT = the quadrature band; cycles centred on
              each conjunction (+/- 182.6 d)
  per cycle : Delta_y = mean r(IN_y) - mean r(OUT_y);
              v_y = s^2(IN_y)/n_IN + s^2(OUT_y)/n_OUT
  statistics: S_sym  = sum(Delta_y / v_y) / sqrt(sum 1/v_y)
              S_skirt = same with IN, OUT restricted to eps < 90
              S_year = max_y Delta_y / sqrt(v_y)
              S_opp  = same as S_skirt on the eps > 90 side (annotation)
  validity  : cycle valid for S_sym/S_year when n_IN >= 5 and
              n_OUT >= 5; for S_skirt when n_skirt >= 3 and
              n_OUT(S2 side) >= 5; a series is valid for a statistic
              with >= 3 valid cycles (S_year: >= 1)
  controls  : T = max over valid controls; exceedance S > max(T, 0);
              searched trial needs >= 4 valid controls
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord, get_sun, EarthLocation, AltAz
from astropy.table import Table
from astropy.time import Time  # noqa: F401 (re-exported for cut_stage)
import astropy.units as u

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'atlas-asassn-crossings', 'scripts'))
import atlas_api  # noqa: E402

REPO = atlas_api.REPO
D = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
RES = os.path.join(D, 'results')
RUN = os.path.join(REPO, 'runs', 'skirt-crossings')
LC = os.path.join(RUN, 'lc')
GEOM = os.path.join(RUN, 'geometry')
SITES = {'01': ('Haleakala', 20.7076, -156.2569, 3055), '02': ('Mauna Loa', 19.5362, -155.5763, 3397),
         '03': ('Sutherland', -32.3796, 20.8107, 1798), '04': ('El Sauce', -30.4708, -70.7647, 1600)}
RUNGS = ('0.9', '0.95')
STATS = ('S_sym', 'S_skirt', 'S_year')
R_GATE = 0.9
_RC = os.path.join(D, 'results', 'response_calib_v1.json')
K_RESP = json.load(open(_RC))['fit']['k'] if os.path.exists(_RC) else 4.0   # A3: measured response exponent
MIN_IN, MIN_OUT, MIN_SKIRT = 5, 5, 3
MIN_CYCLES = {'S_sym': 3, 'S_skirt': 3, 'S_year': 1}
MIN_CONTROLS = 4
HALF_CYCLE = 182.625

_TASKS = None


def tasks():
    global _TASKS
    if _TASKS is None:
        _TASKS = Table.read(os.path.join(RES, 'task_list_v1.ecsv'))
    return _TASKS


def task_row(key):
    t = tasks()
    return t[np.where(t['task_key'] == key)[0][0]]


def night_id(mjd, unit):
    lon = np.array([SITES.get(u_, (None, 0, 0, 0))[2] for u_ in unit])
    return np.floor(mjd + lon / 360.0 + 0.5).astype(int)


_SERIES_CACHE = {}


def load_series(key, band, inject=None):
    """Masked exposure-level series for one task; inject = (mask_fn, f_uJy)
    adds f to the response-corrected flux of exposures selected by
    mask_fn(df). The astrometric/airmass part is cached per (key, band);
    injection acts on a copy."""
    ck = (key, band)
    if ck not in _SERIES_CACHE:
        _SERIES_CACHE[ck] = _load_series_uncached(key, band)
    df = _SERIES_CACHE[ck].copy()
    if inject is not None and len(df):
        mask_fn, f_uJy = inject
        sel = mask_fn(df)
        df.loc[sel, 'flux'] = df.loc[sel, 'flux'] + f_uJy       # transmitter at the star: same R, same correction
    return df


def _load_series_uncached(key, band):
    path = os.path.join(LC, f'{key}.txt')
    df = atlas_api.read_result(path)
    df = df[df.F == band].copy()
    if not len(df):
        return df
    row = task_row(key)
    df['unit'] = df.Obs.str[:2]
    ok = atlas_api.faq_quality_mask(df).values & df.unit.isin(SITES).values
    df = df[ok].copy()
    if not len(df):
        return df
    # response (amendment A3): forced position vs the star's apparent
    # position = Gaia track (PM) + parallactic displacement (registry
    # parallax, Earth's barycentric position); a fixed-position PSF fit
    # recovers exp(-d^2 / 4 sigma^2) of the flux (cross-correlation of
    # two Gaussians), so the measured flux is corrected by 1/R_i
    dt = (df.MJD.values - Time(float(row['epoch_year']), format='jyear').mjd) / 365.25
    ra0, de0 = float(row['ra_deg']), float(row['dec_deg'])
    ra_p, de_p = apparent_position(ra0, de0, float(row['pm_ra_cosdec_mas_yr']), float(row['pm_dec_mas_yr']),
                                   float(row['parallax_mas']), dt, df.MJD.values)
    d = SkyCoord(df.RA.values * u.deg, df.Dec.values * u.deg).separation(SkyCoord(ra_p * u.deg, de_p * u.deg)).arcsec
    fwhm = 1.86 * np.sqrt(df['maj'].values * df['min'].values)
    df['offset_arcsec'] = d
    df['R'] = np.exp(-d ** 2 / (K_RESP * (fwhm / 2.355) ** 2))
    df = df[df.R >= R_GATE].copy()
    if not len(df):
        return df
    # airmass per exposure at the star's own position
    t = Time(df.MJD.values, format='mjd')
    star = SkyCoord(ra0 * u.deg, de0 * u.deg)
    am = np.full(len(df), np.nan)
    for code, (nm, lat, lon, h) in SITES.items():
        m = (df.unit == code).values
        if m.any():
            loc = EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=h * u.m)
            am[m] = star.transform_to(AltAz(obstime=t[m], location=loc)).secz.value
    df['airmass'] = am
    df = df[(df.airmass > 0.99) & (df.airmass < 6)].copy()
    df['night'] = night_id(df.MJD.values, df.unit.values)
    df['flux'] = df.uJy.values.astype(float) / df.R.values   # A3 response correction
    return df


def apparent_position(ra0, de0, pmra_cosdec, pmdec, plx_mas, dt_yr, mjd):
    """ICRS apparent direction: PM-propagated + parallactic displacement
    (no aberration: the image WCS is ICRS). Returns (ra, dec) in deg."""
    from astropy.coordinates import get_body_barycentric
    ra = np.radians(ra0 + pmra_cosdec / 3.6e6 * dt_yr / np.cos(np.radians(de0)))
    de = np.radians(de0 + pmdec / 3.6e6 * dt_yr)
    s = np.vstack([np.cos(de) * np.cos(ra), np.cos(de) * np.sin(ra), np.sin(de)]).T
    r = get_body_barycentric('earth', Time(mjd, format='mjd')).xyz.to(u.AU).value.T   # ICRS AU
    plx = np.radians(plx_mas / 3.6e6)
    s2 = s - plx * r
    s2 /= np.linalg.norm(s2, axis=1)[:, None]
    return np.degrees(np.arctan2(s2[:, 1], s2[:, 0])) % 360.0, np.degrees(np.arcsin(s2[:, 2]))


def nightly(df):
    """Nightly-unit epochs: median fractional flux of >= 2 exposures."""
    if not len(df):
        return pd.DataFrame(columns=['night', 'unit', 'mjd', 'f', 'n', 'airmass'])
    med = np.median(df.flux.values)
    g = df.assign(f=df.flux / med).groupby(['unit', 'night'])
    out = g.agg(mjd=('MJD', 'mean'), f=('f', 'median'), n=('f', 'size'), airmass=('airmass', 'mean')).reset_index()
    return out[out.n >= 2].reset_index(drop=True)


def airmass_layer(nt):
    """r = f - (a + b (X-1)); returns nt with r, and (a, b)."""
    if len(nt) < 5:
        return nt.assign(r=nt.f - np.median(nt.f) if len(nt) else []), (np.nan, np.nan)
    A = np.vstack([np.ones(len(nt)), nt.airmass.values - 1.0]).T
    coef, *_ = np.linalg.lstsq(A, nt.f.values, rcond=None)
    return nt.assign(r=nt.f.values - A @ coef), (float(coef[0]), float(coef[1]))


class Geometry:
    """The target's elongation curve, conjunction cycles and rung sets."""

    def __init__(self, tid):
        g = np.load(os.path.join(GEOM, f'{tid}.npz'))
        self.mjd, self.eps, self.b_e, self.r_e = g['mjd'], g['eps'], g['b_e'], g['r_e']
        # conjunctions: local minima of eps
        e = self.eps
        idx = np.where((e[1:-1] < e[:-2]) & (e[1:-1] <= e[2:]))[0] + 1
        self.conj = self.mjd[idx]

    def eps_at(self, mjd):
        return np.interp(mjd, self.mjd, self.eps)

    def b_at(self, mjd):
        return np.interp(mjd, self.mjd, self.b_e)

    def cycle_of(self, mjd):
        """Index of the conjunction cycle containing each epoch (-1 if none)."""
        out = np.full(len(mjd), -1)
        for i, c in enumerate(self.conj):
            m = (mjd >= c - HALF_CYCLE) & (mjd < c + HALF_CYCLE)
            out[m] = i
        return out


def cycle_stats(nt, geom, rb):
    """Per-cycle Delta and variance for the three set definitions."""
    if not len(nt):
        return {}
    eps = geom.eps_at(nt.mjd.values)
    b = geom.b_at(nt.mjd.values)
    cyc = geom.cycle_of(nt.mjd.values)
    inb = b < rb
    s2 = eps < 90.0
    out = {}
    for i in sorted(set(cyc[cyc >= 0])):
        m = cyc == i
        rec = {}
        for name, inm, outm in (('sym', inb, ~inb), ('skirt', inb & s2, ~inb & s2), ('opp', inb & ~s2, ~inb & ~s2)):
            a, o = nt.r.values[m & inm], nt.r.values[m & outm]
            rec[name] = {'n_in': int(len(a)), 'n_out': int(len(o)),
                         'delta': float(a.mean() - o.mean()) if len(a) and len(o) else None,
                         'v': float((a.var(ddof=1) / len(a) if len(a) > 1 else np.nan) + (o.var(ddof=1) / len(o) if len(o) > 1 else np.nan)) if len(a) > 1 and len(o) > 1 else None}
        out[int(i)] = rec
    return out


def combine(cs):
    """Statistics + validity from per-cycle records."""
    res = {}
    for stat, name, min_in in (('S_sym', 'sym', MIN_IN), ('S_skirt', 'skirt', MIN_SKIRT), ('S_year', 'sym', MIN_IN), ('S_opp', 'opp', MIN_IN)):
        ds, vs = [], []
        for i, rec in cs.items():
            r = rec[name]
            if r['n_in'] >= min_in and r['n_out'] >= MIN_OUT and r['v'] and r['v'] > 0:
                ds.append(r['delta']); vs.append(r['v'])
        ds, vs = np.array(ds), np.array(vs)
        n = len(ds)
        if stat == 'S_year':
            S = float(np.max(ds / np.sqrt(vs))) if n else None
        else:
            S = float(np.sum(ds / vs) / np.sqrt(np.sum(1 / vs))) if n else None
        need = MIN_CYCLES.get(stat, 3)
        res[stat] = {'S': S, 'n_cycles': int(n), 'valid': bool(n >= need),
                     'delta_mean': float(np.sum(ds / vs) / np.sum(1 / vs)) if n else None}
    return res


def star_stats(key, tid, band, rb, inject=None, dropout=None):
    """Full chain for one star at one rung. dropout = (rng, frac) drops
    a random fraction of nightly epochs (completeness draws)."""
    df = load_series(key, band, inject=inject)
    nt = nightly(df)
    if dropout is not None and len(nt):
        rng, frac = dropout
        keep = rng.random(len(nt)) >= frac
        nt = nt[keep].reset_index(drop=True)
    nt, ab = airmass_layer(nt)
    geom = Geometry(tid)
    cs = cycle_stats(nt, geom, float(rb))
    res = combine(cs)
    res['n_nights'] = int(len(nt))
    res['airmass_fit'] = ab
    res['cycles'] = cs
    return res


def unit_stats(tid, band, rb, controls=None, **kw):
    """Target + control statistics for one unit."""
    t = tasks()
    ckeys = controls if controls is not None else [str(k) for k in t['task_key'][(t['target_id'] == tid) & (t['role'] == 'control')]]
    real = star_stats(f'T-{tid}', tid, band, rb, **kw)
    ctrl = {}
    for k in ckeys:
        if os.path.exists(os.path.join(LC, f'{k}.txt')):
            ctrl[k] = star_stats(k, tid, band, rb)
    return {'real': real, 'controls': ctrl}


def thresholds(us):
    """Max rule over valid controls, per statistic."""
    out = {}
    for stat in STATS:
        real = us['real'][stat]
        vals = {k: c[stat]['S'] for k, c in us['controls'].items() if c[stat]['valid'] and c[stat]['S'] is not None}
        T = max(vals.values()) if vals else None
        S = real['S'] if real['valid'] else None
        searched = real['valid'] and len(vals) >= MIN_CONTROLS
        out[stat] = {'S': S, 'T': T, 'n_valid_controls': len(vals), 'controls': {k: round(v, 3) for k, v in vals.items()},
                     'searched': bool(searched),
                     'exceedance': bool(searched and S is not None and S > max(T, 0.0)),
                     'margin': (None if S is None or T is None else round(S - T, 3)),
                     'S_opp': us['real']['S_opp']['S'], 'delta_mean': real.get('delta_mean')}
    return out


def inject_mask(tid, rb):
    geom = Geometry(tid)

    def fn(df):
        return geom.b_at(df.MJD.values) < float(rb)
    return fn


def cycle_counts(nt, geom, rb):
    """Coverage-stage counts only (no in-beam statistic is formed)."""
    if not len(nt):
        return {}
    eps = geom.eps_at(nt.mjd.values)
    b = geom.b_at(nt.mjd.values)
    cyc = geom.cycle_of(nt.mjd.values)
    inb = b < rb
    s2 = eps < 90.0
    out = {}
    for i in sorted(set(cyc[cyc >= 0])):
        m = cyc == i
        out[int(i)] = {'conj_mjd': float(geom.conj[i]),
                       'n_in_sym': int((m & inb).sum()), 'n_out_sym': int((m & ~inb).sum()),
                       'n_skirt': int((m & inb & s2).sum()), 'n_out_s2': int((m & ~inb & s2).sum()),
                       'n_opp': int((m & inb & ~s2).sum()), 'n_out_a': int((m & ~inb & ~s2).sum())}
    return out


def validity_from_counts(cc):
    v = {}
    for stat, key_in, key_out, min_in in (('S_sym', 'n_in_sym', 'n_out_sym', MIN_IN), ('S_skirt', 'n_skirt', 'n_out_s2', MIN_SKIRT), ('S_year', 'n_in_sym', 'n_out_sym', MIN_IN)):
        n = sum(1 for r in cc.values() if r[key_in] >= min_in and r[key_out] >= MIN_OUT)
        v[stat] = {'n_valid_cycles': n, 'valid': n >= MIN_CYCLES[stat]}
    return v
