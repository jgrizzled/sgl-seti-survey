"""Frozen statistics engine for the ATLAS crossings search
(configs/threshold_freeze_v1.json). Every construction here is the
freeze text made executable; the dev stage validates it, the
confirmatory stage runs it blind.

Per unit (target, channel, rung, band):
  series      : primary-mask epochs of the unit's task positions
                (B 0.1 AU: tca + ingress + egress merged by Obs id)
  baseline    : epochs outside the event's widest (0.1 AU) window and
                outside the evaluated (real or pseudo) window
  detrend     : r_i = f_i - running median of baseline epochs within
                +/-15 d (>= 5, else global baseline median)
  variance    : v_i = k sigma_i^2, k = median(r^2/sigma^2)/0.4549 over
                3x3-sigma-clipped baseline epochs, floor 1
  response    : R_i = exp(-d_i^2 / 2 s_i^2), d_i = task-position offset
                from the predicted apparent source at epoch i (B: relay
                at z on the anti-star axis, evaluated as if the event
                were centred on the evaluated window; A: the star),
                s_i = FWHM_i / 2.355, FWHM_i = 1.86" sqrt(maj min);
                epochs with R_i < 0.2 dropped; B takes per epoch the
                mini-track position of largest R_i and the max over z
  statistics  : S_w = sum(w r)/sqrt(sum(w^2 v)), w = R/v (window)
                S_pulse_w = max r sqrt(R)/sqrt(v)
                S_event = max_w S_w ; S_stack = pooled over windows ;
                S_pulse = max_w S_pulse_w
  controls    : the 8 designated offsets, mirror-gated validity,
                T = max valid control, exceedance S > max(T, 0)
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from astropy.table import Table

sys.path.insert(0, os.path.dirname(__file__))
import atlas_api  # noqa: E402
from astropy.coordinates import get_body_barycentric  # noqa: E402
from astropy.time import Time  # noqa: E402
from build_task_list import sep_arcsec  # noqa: E402


def relay_apparent_vec(star_ra, star_dec, z_au, t_mjd):
    """Vectorised build_task_list.relay_apparent (identical math)."""
    t = Time(np.atleast_1d(t_mjd), format='mjd', scale='utc')
    sun = get_body_barycentric('sun', t).xyz.to_value('AU')
    earth = get_body_barycentric('earth', t).xyz.to_value('AU')
    ra, de = np.radians(star_ra), np.radians(star_dec)
    u = np.array([np.cos(de) * np.cos(ra), np.cos(de) * np.sin(ra), np.sin(de)])[:, None]
    v = (sun - z_au * u) - earth
    v = v / np.linalg.norm(v, axis=0)
    return (np.degrees(np.arctan2(v[1], v[0])) % 360.0,
            np.degrees(np.arcsin(np.clip(v[2], -1, 1))))


def sep_arcsec_vec(ra1, de1, ra2, de2):
    d1, e1, d2, e2 = map(np.radians, (ra1, de1, ra2, de2))
    c = np.sin(e1) * np.sin(e2) + np.cos(e1) * np.cos(e2) * np.cos(d1 - d2)
    return np.degrees(np.arccos(np.clip(c, -1, 1))) * 3600

REPO = atlas_api.REPO
D = os.path.join(REPO, 'surveys', 'atlas-asassn-crossings')
LC = os.path.join(REPO, 'runs', 'atlas-asassn-crossings', 'lc')
FREEZE = json.load(open(os.path.join(D, 'configs', 'threshold_freeze_v1.json')))
TASKS = Table.read(os.path.join(D, 'results', 'task_list_v1.ecsv'))
EVENTS = Table.read(os.path.join(REPO, 'crossings', 'universal_v1', 'events.ecsv'))
OFFSETS = tuple(FREEZE['controls']['designated_offsets_days'])
Z_GRID = (550.0, 1000.0, 2500.0, 5500.0, 10000.0)
R_MIN, MED_HALF_D, MED_MIN_N, K_MIN_N = 0.2, 15.0, 5, 6
PIX = 1.86
HALF_COL = {'1.2Rsun': 'half_1p2Rsun_d', '2.5Rsun': 'half_2p5Rsun_d', '0.1AU': 'half_0p1AU_d'}
_STAR = {str(r['event_id']): (float(r['star_icrs_ra_deg']), float(r['star_icrs_dec_deg'])) for r in EVENTS}


def task_rows(ch, tid, eid):
    m = (TASKS['channel'] == ch) & (TASKS['target_id'] == tid) & (TASKS['event_id'] == eid)
    return TASKS[m]


def load_masked(key, band):
    df = atlas_api.read_result(os.path.join(LC, f'{key}.txt'))
    if len(df) == 0:
        return df
    df = df[atlas_api.faq_quality_mask(df) & (df.F == band)].copy()
    df['fwhm_arcsec'] = PIX * np.sqrt(df['maj'] * df['min'])
    return df.reset_index(drop=True)


class EventSeries:
    """All task positions of one (channel, event) merged by Obs id."""

    def __init__(self, ch, tid, eid, band):
        rows = task_rows(ch, tid, eid)
        self.ch, self.tid, self.eid = ch, tid, eid
        self.t_ca = float(rows[0]['t_ca_mjd'])
        self.half = {r: float(rows[0][c]) for r, c in HALF_COL.items() if np.isfinite(rows[0][c])}
        self.half_widest = max(self.half.values())
        self.positions = {}   # pos_role -> (ra, dec, dataframe)
        for r in rows:
            df = load_masked(str(r['task_key']), band)
            self.positions[str(r['pos_role'])] = (float(r['ra_deg']), float(r['dec_deg']), df)
        self.star = _STAR[eid]
        # union epoch table: one row per Obs id with per-position flux
        frames = []
        for pos, (ra, de, df) in self.positions.items():
            if len(df):
                frames.append(df[['Obs', 'MJD', 'uJy', 'duJy', 'fwhm_arcsec']].assign(pos=pos))
        self.long = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(
            columns=['Obs', 'MJD', 'uJy', 'duJy', 'fwhm_arcsec', 'pos'])
        self.mjd_by_obs = self.long.groupby('Obs')['MJD'].first()

    def baseline_mask(self, mjd, tc, h):
        return (np.abs(mjd - self.t_ca) > self.half_widest) & (np.abs(mjd - tc) > h)

    _inject = None  # optional {pos: {Obs: added uJy}} applied before detrending
    fast = False    # completeness: interpolate the relay track on a 0.02-d grid
    _grid = None
    _cache = None   # (pos, tc, h) -> detrended result; z-independent, cleared on inject change

    @property
    def inject(self):
        return self._inject

    @inject.setter
    def inject(self, value):
        self._inject = value
        self._cache = {}

    def track(self, z, t):
        """Apparent relay (ra, dec) at times t; exact, or gridded when fast."""
        if not self.fast:
            return relay_apparent_vec(self.star[0], self.star[1], z, t)
        if self._grid is None:
            self._grid = {}
        if z not in self._grid:
            tg = np.arange(self.t_ca - 240.0, self.t_ca + 240.0, 0.02)
            self._grid[z] = (tg,) + relay_apparent_vec(self.star[0], self.star[1], z, tg)
        tg, rg, dg = self._grid[z]
        return np.interp(t, tg, rg), np.interp(t, tg, dg)

    def detrended(self, pos, tc, h):
        """Residuals + rescaled variances for one position's series
        given the evaluated window (tc, h). Returns df with r, v, k."""
        if self._cache is None:
            self._cache = {}
        ck = (pos, round(tc, 6), round(h, 6))
        if ck in self._cache:
            return self._cache[ck]
        res = self._detrend(pos, tc, h)
        self._cache[ck] = res
        return res

    def _detrend(self, pos, tc, h):
        ra, de, df = self.positions[pos]
        if len(df) == 0:
            return df.assign(r=[], v=[]), 1.0, 0
        mjd = df.MJD.values
        f, s = df.uJy.values.astype(float), df.duJy.values.astype(float)
        if self.inject and pos in self.inject:
            f = f + np.array([self.inject[pos].get(o, 0.0) for o in df.Obs.values])
        base = self.baseline_mask(mjd, tc, h)
        gmed = np.median(f[base]) if base.sum() else 0.0
        med = np.empty_like(f)
        for i in range(len(f)):
            near = base & (np.abs(mjd - mjd[i]) <= MED_HALF_D)
            med[i] = np.median(f[near]) if near.sum() >= MED_MIN_N else gmed
        r = f - med
        rb, sb = r[base], s[base]
        k = 1.0
        if len(rb) >= K_MIN_N:
            keep = np.ones(len(rb), bool)
            for _ in range(3):
                z = rb / sb
                sd = 1.4826 * np.median(np.abs(z[keep] - np.median(z[keep])))
                keep = np.abs(z - np.median(z[keep])) <= 3 * max(sd, 1e-9)
            k = max(1.0, float(np.median((rb[keep] ** 2) / (sb[keep] ** 2)) / 0.4549))
        out = df.assign(r=r, v=k * s ** 2)
        return out, k, int(base.sum())

    def response(self, pos, mjd, fwhm, tc, z):
        """R_i for a position's epochs, event centred at tc (real: t_ca)."""
        ra, de, _ = self.positions[pos]
        if self.ch == 'A':
            return np.ones(len(mjd))
        shift = tc - self.t_ca
        pra, pde = self.track(z, np.asarray(mjd, float) - shift)
        d = sep_arcsec_vec(ra, de, pra, pde)
        s = fwhm / 2.355
        return np.exp(-d ** 2 / (2 * s ** 2))

    def window_stats(self, rung, tc, z):
        """(S_w, S_pulse_w, n_epochs, k_by_pos) for the window centred at
        tc (t_ca for real, t_ca+offset for controls) at relay distance z."""
        h = self.half[rung]
        best = {}   # Obs -> (R, r, v)
        ks = {}
        for pos in self.positions:
            det, k, nb = self.detrended(pos, tc, h)
            ks[pos] = (round(k, 3), nb)
            if len(det) == 0:
                continue
            inw = np.abs(det.MJD.values - tc) <= h
            if not inw.any():
                continue
            sub = det[inw]
            R = self.response(pos, sub.MJD.values, sub.fwhm_arcsec.values, tc, z)
            for obs, Ri, ri, vi in zip(sub.Obs.values, R, sub.r.values, sub.v.values):
                if Ri < R_MIN:
                    continue
                if obs not in best or Ri > best[obs][0]:
                    best[obs] = (Ri, ri, vi)
        if not best:
            return None
        R = np.array([b[0] for b in best.values()])
        r = np.array([b[1] for b in best.values()])
        v = np.array([b[2] for b in best.values()])
        w = R / v
        num, den = float(np.sum(w * r)), float(np.sum(w ** 2 * v))
        return {'num': num, 'den': den, 'S_w': num / np.sqrt(den),
                'S_pulse_w': float(np.max(r * np.sqrt(R) / np.sqrt(v))),
                'n': int(len(best)), 'k': ks}


def unit_stats(ch, tid, rung, band, event_ids, offsets=OFFSETS):
    """Real and control statistics for one unit. Returns dict."""
    zs = Z_GRID if ch == 'B' else (None,)
    series = {eid: EventSeries(ch, tid, eid, band) for eid in event_ids}
    centres = {'real': 0.0, **{f'off{o:+d}': float(o) for o in offsets}}
    out = {}
    for name, off in centres.items():
        per_z = []
        for z in zs:
            wins = {}
            for eid, es in series.items():
                ws = es.window_stats(rung, es.t_ca + off, z)
                if ws is not None:
                    wins[eid] = ws
            if not wins:
                per_z.append(None)
                continue
            num = sum(w['num'] for w in wins.values())
            den = sum(w['den'] for w in wins.values())
            per_z.append({'z': z, 'n_windows': len(wins),
                          'S_event': max(w['S_w'] for w in wins.values()),
                          'S_pulse': max(w['S_pulse_w'] for w in wins.values()),
                          'S_stack': num / np.sqrt(den),
                          'windows': {e: {'S_w': round(w['S_w'], 3), 'S_pulse_w': round(w['S_pulse_w'], 3),
                                          'n': w['n'], 'k': w['k']} for e, w in wins.items()}})
        valid = [p for p in per_z if p is not None]
        if not valid:
            out[name] = None
            continue
        best = {}
        for stat in ('S_event', 'S_pulse', 'S_stack'):
            i = int(np.argmax([p[stat] for p in valid]))
            best[stat] = {'S': float(valid[i][stat]), 'z': valid[i]['z'], 'n_windows': valid[i]['n_windows']}
        best['n_windows'] = max(p['n_windows'] for p in valid)
        best['detail_z550'] = valid[0]
        out[name] = best
    return out


def thresholds(stats, min_stack_events):
    """Apply the freeze's mirror-gated validity + max rule."""
    real = stats['real']
    res = {}
    for stat in ('S_event', 'S_pulse', 'S_stack'):
        vals, names = [], []
        for name, c in stats.items():
            if name == 'real' or c is None:
                continue
            need = min_stack_events if stat == 'S_stack' else 1
            if c['n_windows'] >= need:
                vals.append(c[stat]['S'])
                names.append(name)
        T = max(vals) if vals else None
        S = real[stat]['S'] if real else None
        res[stat] = {'S': S, 'T': T, 'n_valid_controls': len(vals),
                     'controls': dict(zip(names, [round(v, 3) for v in vals])),
                     'exceedance': (S is not None and T is not None and S > max(T, 0.0)),
                     'margin': (None if S is None or T is None else round(S - T, 3))}
    return res
