"""S2 1 AU outer-skirt geometry pass (plan §5.25; recon stage).

The S2 combination (uplink `inbound`, Earth on the anti-target side)
has the target star itself as the apparent source at solar elongation
eps, with Earth's distance from the Sun-star axis b_e = r_E sin(eps)
exactly (Earth on the anti-target side <=> eps < 90 deg; the channel-A
side is eps > 90 deg with the same b_e). A top-hat uplink beam of
radius r_b is therefore "on" whenever r_E sin(eps) < r_b, on BOTH
sides of the Sun: near conjunction (S2, the skirt) and near opposition
(A). This script tabulates, per target and ATLAS-era year, the
elongation curve and the in-beam night counts per candidate rung and
per candidate ATLAS solar-elongation floor, so the freeze can set the
rung ladder and the searchable population from measured numbers.

Inputs: crossings/universal_v1/events.ecsv (S2 rows give per-target
star ICRS positions and the yearly conjunction b_min for validation).
Output: results/skirt_geometry_v1.json (+ per-target daily curves in
runs/skirt-crossings/geometry/).
"""
import json
import os
import sys

import numpy as np
from astropy.coordinates import SkyCoord, get_sun, GeocentricTrueEcliptic
from astropy.table import Table
from astropy.time import Time
import astropy.units as u

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
EVENTS = os.path.join(REPO, 'crossings', 'universal_v1', 'events.ecsv')
OUT = os.path.join(os.path.dirname(__file__), '..', 'results', 'skirt_geometry_v1.json')
RUNDIR = os.path.join(REPO, 'runs', 'skirt-crossings', 'geometry')

ERA = (57227.0, 61286.0)          # ATLAS archive start -> ATLAS-survey era end (A2)
RUNGS_AU = [0.5, 0.7, 0.8, 0.85, 0.9, 0.95, 1.0]
FLOORS_DEG = [40, 45, 50, 55, 60, 70]   # candidate ATLAS solar-elongation floors


def main():
    t = Table.read(EVENTS)
    s2 = t[(t['link_direction'] == 'inbound') & (t['side'] == 'anti_target')]
    a = t[(t['link_direction'] == 'inbound') & (t['side'] == 'target')]
    targets = sorted(set(s2['target_id']))
    mjd = np.arange(ERA[0], ERA[1] + 1, 1.0)
    times = Time(mjd, format='mjd', scale='utc')
    sun = get_sun(times)
    r_e = sun.distance.to(u.AU).value           # Earth-Sun distance
    os.makedirs(RUNDIR, exist_ok=True)
    out = {'era_mjd': ERA, 'rungs_au': RUNGS_AU, 'floors_deg': FLOORS_DEG,
           'n_targets': len(targets), 'targets': {}}
    for tid in targets:
        rows = s2[s2['target_id'] == tid]
        # era-mean star position (PM over 11 yr is arcsec: irrelevant for eps)
        ra = float(np.median(rows['star_icrs_ra_deg']))
        dec = float(np.median(rows['star_icrs_dec_deg']))
        star = SkyCoord(ra * u.deg, dec * u.deg, frame='icrs')
        beta = float(star.transform_to(GeocentricTrueEcliptic(equinox='J2000')).lat.deg)
        eps = sun.separation(star).deg
        b_e = r_e * np.sin(np.radians(eps))
        # validation: yearly conjunction minima vs the events table b_min
        era_rows = rows[(rows['t_ca_tdb_jd'] - 2400000.5 > ERA[0]) & (rows['t_ca_tdb_jd'] - 2400000.5 < ERA[1])]
        val = []
        for r in era_rows:
            tca = r['t_ca_tdb_jd'] - 2400000.5
            i = int(np.argmin(np.abs(mjd - tca)))
            lo, hi = max(0, i - 200), min(len(mjd), i + 200)
            j = lo + int(np.argmin(b_e[lo:hi]))
            val.append({'t_ca_mjd': tca, 'b_min_table_au': float(r['b_min_au']),
                        'b_min_curve_au': float(b_e[j]), 'mjd_curve_min': float(mjd[j]),
                        'eps_min_deg': float(eps[j]), 'event_id': str(r['event_id'])})
        # per-rung in-beam day counts (era total), split by side and floor
        rung_stats = {}
        for rb in RUNGS_AU:
            inbeam = b_e < rb
            rec = {'days_inbeam_total': int(inbeam.sum()),
                   'days_inbeam_s2': int((inbeam & (eps < 90)).sum()),
                   'days_inbeam_a': int((inbeam & (eps > 90)).sum()),
                   'eps_edge_deg': float(np.degrees(np.arcsin(min(rb / r_e.mean(), 1.0)))),
                   'per_floor': {}}
            for fl in FLOORS_DEG:
                vis = eps >= fl
                rec['per_floor'][str(fl)] = {
                    'days_inbeam_s2_visible': int((inbeam & vis & (eps < 90)).sum()),
                    'days_outbeam_s2_visible': int((~inbeam & vis & (eps < 90)).sum()),
                    'days_inbeam_a': int((inbeam & (eps > 90)).sum()),
                    'days_outbeam_a': int((~inbeam & (eps > 90)).sum()),
                }
            rung_stats[str(rb)] = rec
        # skirt extent per year: fraction of the S2 half-year at eps >= floor
        yrs = {}
        for y in range(2015, 2027):
            m = (times.decimalyear >= y) & (times.decimalyear < y + 1) & (eps < 90)
            if m.sum() == 0:
                continue
            yrs[str(y)] = {fl: int((m & (eps >= fl)).sum()) for fl in FLOORS_DEG}
        out['targets'][tid] = {
            'ra_deg': ra, 'dec_deg': dec, 'ecl_lat_deg': beta,
            'eps_min_deg': float(eps.min()), 'eps_max_deg': float(eps.max()),
            'n_s2_events_era': len(era_rows), 'validation': val,
            'rungs': rung_stats, 's2_days_visible_per_year': yrs,
        }
        np.savez_compressed(os.path.join(RUNDIR, f'{tid}.npz'), mjd=mjd, eps=eps, b_e=b_e, r_e=r_e)
    # global validation summary
    resid = [abs(v['b_min_table_au'] - v['b_min_curve_au']) for tt in out['targets'].values() for v in tt['validation']]
    out['validation_summary'] = {'n': len(resid), 'max_abs_resid_au': float(max(resid)), 'median_abs_resid_au': float(np.median(resid))}
    with open(OUT, 'w') as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out['validation_summary']))
    print(f'{len(targets)} targets written to {OUT}')


if __name__ == '__main__':
    main()
