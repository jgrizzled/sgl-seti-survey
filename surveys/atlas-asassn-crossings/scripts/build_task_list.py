"""Build the ATLAS coverage-stage task list (hypotheses §6, decision D3).

Per in-era event (ATLAS era MJD 57227 -> 61286 = coverage start
2026-09-03), b <= 0.1 AU, seven narrow-rung targets, sunward
combinations excluded (hypotheses amendment A1: the universal list
holds two axis crossings per year; the one that puts the channel's
sky position at solar elongation ~0 is out of scope, leaving 11
annual events per target per channel - verified here by elongation):
  channel B (antipode, outbound & axis_distance < 0): three fixed
    positions - the tabulated apparent axis position at t_ca (serves
    the grazing rungs and the mini-track centre) plus the apparent
    z = 550 AU relay position at t_ca -/+ the 0.1 AU half-window
    (ingress / egress) - the D3 mini-track;
  channel A (star, inbound & axis_distance > 0): one fixed position,
    the per-event PM-propagated star position. The D3 "3-position
    mini-track" for A collapses by construction: the star has no
    z-track, and its motion over the 0.1 AU half-window (<= 5"/yr x
    ~6 d < 0.1") is far inside any resolvable offset, so the three
    positions coincide.
Every task is a windowed difference-flux forced-photometry request,
mjd in [t_ca - 110, t_ca + 110] d (window + baseline + the 8
pseudo-window offsets to +/-97 d). The A 1.0 AU constraint rung is not
in this list (constraint-only, frozen at draft).

Output: results/task_list_v1.ecsv + task_list_v1_summary.json.
"""
import json
import os

import numpy as np
from astropy.coordinates import SkyCoord, get_body_barycentric, get_sun
from astropy.table import Table
from astropy.time import Time

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
EVENTS = os.path.join(REPO, 'crossings', 'universal_v1', 'events.ecsv')
OUT = os.path.join(os.path.dirname(__file__), '..', 'results')
ERA = (57227.0, 61286.0)
TARGETS = ['van-maanen', 'wolf-359', 'gj-1276', 'teegarden', 'ross-128',
           'ross-154', 'gj-908']
KM_PER_AU = 1.495978707e8
RSUN_AU = 0.00465047
LADDER_B = {'1.2Rsun': 1.2 * RSUN_AU, '2.5Rsun': 2.5 * RSUN_AU, '0.1AU': 0.1}
LADDER_A = {'0.1AU': 0.1}
HALF_RANGE_D = 110.0
Z_TRACK_AU = 550.0


def relay_apparent(star_ra, star_dec, z_au, t_mjd):
    t = Time(t_mjd, format='mjd', scale='utc')
    sun = get_body_barycentric('sun', t).xyz.to_value('AU')
    earth = get_body_barycentric('earth', t).xyz.to_value('AU')
    ra, de = np.radians(star_ra), np.radians(star_dec)
    u = np.array([np.cos(de) * np.cos(ra), np.cos(de) * np.sin(ra), np.sin(de)])
    v = (sun - z_au * u) - earth
    v /= np.linalg.norm(v)
    return (float(np.degrees(np.arctan2(v[1], v[0])) % 360.0),
            float(np.degrees(np.arcsin(np.clip(v[2], -1, 1)))))


def half_window_d(b_au, r_au, v_perp):
    if b_au >= r_au:
        return np.nan
    return float(np.sqrt(r_au ** 2 - b_au ** 2) * KM_PER_AU / v_perp / 86400.0)


def sep_arcsec(ra1, de1, ra2, de2):
    d = np.radians([ra1, de1, ra2, de2])
    c = (np.sin(d[1]) * np.sin(d[3])
         + np.cos(d[1]) * np.cos(d[3]) * np.cos(d[0] - d[2]))
    return float(np.degrees(np.arccos(np.clip(c, -1, 1))) * 3600)


def main():
    t = Table.read(EVENTS)
    mjd = Time(list(t['t_ca_utc']), format='isot', scale='utc').mjd
    t['t_ca_mjd'] = mjd
    era = (mjd >= ERA[0]) & (mjd <= ERA[1])
    ax = np.asarray(t['axis_distance_au'])
    ld = np.asarray(t['link_direction'])
    b_au = np.asarray(t['b_min_au'])
    rows = []
    sunward_dropped = {}
    for tid in TARGETS:
        base = era & (np.asarray(t['target_id']) == tid) & (b_au <= 0.1)
        selA = t[base & (ld == 'inbound') & (ax > 0)]
        selB = t[base & (ld == 'outbound') & (ax < 0)]
        nsunA = int((base & (ld == 'inbound') & (ax < 0)).sum())
        nsunB = int((base & (ld == 'outbound') & (ax > 0)).sum())
        assert len(selA) == 11 and len(selB) == 11, (tid, len(selA), len(selB))
        assert nsunA == 11 and nsunB == 11, (tid, nsunA, nsunB)
        sunward_dropped[tid] = {'A': nsunA, 'B': nsunB}
        for ch, sel in (('B', selB), ('A', selA)):
            for ev in sel:
                b = float(ev['b_min_au'])
                vp = float(ev['v_perp_km_s'])
                tca = float(ev['t_ca_mjd'])
                halves = {k: half_window_d(b, r, vp) for k, r in
                          (LADDER_B if ch == 'B' else LADDER_A).items()}
                if ch == 'A':
                    poss = [('tca', float(ev['star_icrs_ra_deg']),
                             float(ev['star_icrs_dec_deg']), 0.0)]
                else:
                    ra0, de0 = float(ev['relay_icrs_ra_deg']), float(ev['relay_icrs_dec_deg'])
                    h = halves['0.1AU']
                    ra_i, de_i = relay_apparent(float(ev['star_icrs_ra_deg']),
                                                float(ev['star_icrs_dec_deg']),
                                                Z_TRACK_AU, tca - h)
                    ra_e, de_e = relay_apparent(float(ev['star_icrs_ra_deg']),
                                                float(ev['star_icrs_dec_deg']),
                                                Z_TRACK_AU, tca + h)
                    # sanity: the tabulated axis position must match the
                    # recomputed z=550 apparent position at t_ca
                    ra_c, de_c = relay_apparent(float(ev['star_icrs_ra_deg']),
                                                float(ev['star_icrs_dec_deg']),
                                                Z_TRACK_AU, tca)
                    chk = sep_arcsec(ra0, de0, ra_c, de_c)
                    assert chk < 2.0, (tid, ev['event_id'], chk)
                    poss = [('tca', ra0, de0, 0.0),
                            ('ingress', ra_i, de_i, -h), ('egress', ra_e, de_e, +h)]
                for pos, ra, de, dt in poss:
                    key = f"{ch}-{tid}-{str(ev['event_id'])[4:10]}-{pos}"
                    elong = float(get_sun(Time(tca + dt, format='mjd')).separation(
                        SkyCoord(ra, de, unit='deg')).deg)
                    assert elong > 150.0, (key, elong)
                    rows.append({
                        'task_key': key, 'channel': ch, 'target_id': tid,
                        'event_id': str(ev['event_id']), 'pos_role': pos,
                        'pos_dt_days': dt, 'z_track_au': Z_TRACK_AU if ch == 'B' else np.nan,
                        'ra_deg': ra, 'dec_deg': de, 'sun_elong_deg': elong,
                        'offset_from_tca_arcsec': sep_arcsec(ra, de, poss[0][1], poss[0][2]),
                        't_ca_mjd': tca, 'b_min_au': b,
                        'b_min_rsun': float(ev['b_min_solar_radii']),
                        'v_perp_km_s': vp,
                        'half_1p2Rsun_d': halves.get('1.2Rsun', np.nan),
                        'half_2p5Rsun_d': halves.get('2.5Rsun', np.nan),
                        'half_0p1AU_d': halves['0.1AU'],
                        'mjd_min': tca - HALF_RANGE_D, 'mjd_max': tca + HALF_RANGE_D,
                        'use_reduced': False,
                        'validity': str(ev['validity']),
                    })
    out = Table(rows=rows)
    os.makedirs(OUT, exist_ok=True)
    out.write(os.path.join(OUT, 'task_list_v1.ecsv'), format='ascii.ecsv', overwrite=True)
    summ = {'era_mjd': list(ERA), 'half_range_days': HALF_RANGE_D,
            'z_track_au': Z_TRACK_AU, 'n_tasks': len(out),
            'by_channel': {ch: int((out['channel'] == ch).sum()) for ch in 'AB'},
            'by_pos_role': {p: int((out['pos_role'] == p).sum()) for p in ('tca', 'ingress', 'egress')},
            'events_per_target_per_channel': 11,
            'sunward_events_dropped': sunward_dropped,
            'sun_elong_deg_range': [float(np.min(out['sun_elong_deg'])), float(np.max(out['sun_elong_deg']))],
            'mini_track_offset_arcsec': {
                'B_max': float(np.max(out['offset_from_tca_arcsec'][out['channel'] == 'B'])),
                'B_median_nonzero': float(np.median(out['offset_from_tca_arcsec'][(out['channel'] == 'B') & (out['pos_role'] != 'tca')]))},
            'half_windows_days': {
                k: [float(np.nanmin(out[c])), float(np.nanmax(out[c]))]
                for k, c in (('1.2Rsun', 'half_1p2Rsun_d'), ('2.5Rsun', 'half_2p5Rsun_d'), ('0.1AU', 'half_0p1AU_d'))},
            'channel_A_mini_track': 'collapsed to 1 position (geometric identity, see docstring)'}
    with open(os.path.join(OUT, 'task_list_v1_summary.json'), 'w') as f:
        json.dump(summ, f, indent=2)
    print(json.dumps(summ, indent=2))


if __name__ == '__main__':
    main()
