"""Recon probe: ATLAS solar-elongation floor and low-elongation cadence.

Two full-history reduced-mode tasks at NON-TARGET Gaia DR3 field stars
on the ecliptic (beta ~ 0, high Galactic latitude, G ~ 14, quiet,
low PM), so the epoch list samples the full 0-180 deg elongation
range under every ATLAS unit. Declared pre-freeze data contact: none
at any survey position; only epoch/elongation/quality counts are
formed. Raw results + sidecars under runs/skirt-crossings/probe/.
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'atlas-asassn-crossings', 'scripts'))
import atlas_api  # noqa: E402

OUTDIR = os.path.join(atlas_api.REPO, 'runs', 'skirt-crossings', 'probe')
PROBES = {  # Gaia DR3 (VizieR I/355) 2026-09-07
    'probe-ecl180': {'source_id': 3698956337898276352, 'ra': 180.05255145676, 'dec': 0.00422589136, 'G': 14.0855, 'bp_rp': 0.845},
    'probe-ecl000': {'source_id': 2546036791796380672, 'ra': 0.22557877329, 'dec': 0.09924125356, 'G': 14.0353, 'bp_rp': 0.755},
}


def main():
    hdrs = atlas_api.headers()
    os.makedirs(OUTDIR, exist_ok=True)
    tasks = {}
    for name, p in PROBES.items():
        params = {'ra': p['ra'], 'dec': p['dec'], 'mjd_min': 55000.0, 'use_reduced': True, 'comment': f'skirt-{name}'}
        j = atlas_api.submit(params, hdrs)
        tasks[name] = (j, params)
        print(name, 'queued', j['url'], flush=True)
    for name, (j, params) in tasks.items():
        sc = atlas_api.wait_and_fetch(j['url'], hdrs, os.path.join(OUTDIR, f'{name}.txt'), query_params={**params, **PROBES[name]})
        print(name, 'done', sc['wall_s'], 's', sc['result_bytes'], 'bytes', flush=True)


if __name__ == '__main__':
    main()
