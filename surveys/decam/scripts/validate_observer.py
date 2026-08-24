"""Validate the CTIO terrestrial-site observer against an independent
astropy computation (draft hypotheses §2; the ZTF pilot's Palomar check
found ≤ 0.016" agreement).

Two comparisons per epoch x endpoint, at fixed relay distance z:

1. *site effect*: angular shift of the sglseti locus point between the
   CTIO observer and the Earth-centre observer. Expected magnitude
   <= R_earth / z (16 mas at 550 AU), varying with hour angle.
2. *cross-check*: the same shift computed directly with astropy
   (EarthLocation -> GCRS site vector, relay at the Sun + z along the
   anti-star ray). |sglseti - astropy| is the validation residual.

Usage: uv run python surveys/decam/scripts/validate_observer.py
Writes surveys/decam/results/observer_validation.json.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import (GCRS, EarthLocation, SkyCoord,
                                 get_body_barycentric)
from astropy.time import Time

from sglseti import Role, adaptive_locus, load_target_registry

from sglsurvey.geometry import GeometryContext

sys.path.insert(0, str(Path(__file__).resolve().parent))
from decam_corridors import CTIO, PILOT_ENDPOINTS  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
REGISTRY_PATH = REPO / "registries" / "pilot_wise_2026.yaml"
OUT = REPO / "surveys" / "decam" / "results" / "observer_validation.json"

Z_TEST_AU = 550.0
#: Epochs spanning the archive era and the sidereal day (site vector
#: direction changes with LST).
EPOCHS_MJD = [56300.1, 56300.35, 58000.2, 60000.05, 60700.3]


def locus_point_at_z(ctx, target, t, z_au):
    """Interpolated (ra, dec) of the Rx locus at z_au (deg)."""
    al = adaptive_locus(target=target, role=Role.RX, observation_time=t,
                        observer=ctx.observer, relay_range=ctx.relay_range,
                        tolerance_arcsec=0.5, ephemeris=ctx.ephemeris,
                        model=ctx.model)
    zs = np.array([p.z_au for p in al.points])
    ra = np.array([p.icrs_ra_deg for p in al.points])
    dec = np.array([p.icrs_dec_deg for p in al.points])
    q = 1.0 / zs
    o = np.argsort(q)
    return (float(np.interp(1.0 / z_au, q[o], ra[o])),
            float(np.interp(1.0 / z_au, q[o], dec[o])))


def astropy_site_shift(dir_ra, dir_dec, t, z_au, site):
    """Expected apparent shift (site vs geocentre, arcsec) of a relay at
    the Sun + z_au along (dir_ra, dir_dec), computed with astropy only
    (no sglseti). Using the geocentric locus direction for the relay
    ray misplaces the relay by <= 1 AU transverse, which changes the
    16 mas site shift by < 0.1 uas — negligible for this check."""
    sun = get_body_barycentric("sun", t)
    earth = get_body_barycentric("earth", t)
    anti = SkyCoord(ra=dir_ra * u.deg, dec=dir_dec * u.deg).cartesian
    relay = sun + anti * (z_au * u.au)
    site_gcrs = site.get_gcrs_posvel(t)[0]
    d_geo = relay - earth
    d_site = d_geo - site_gcrs.without_differentials()
    v1 = d_geo.xyz / np.linalg.norm(d_geo.xyz)
    v2 = d_site.xyz / np.linalg.norm(d_site.xyz)
    return float(np.rad2deg(np.arccos(
        np.clip(float(np.dot(v1.value, v2.value)), -1, 1))) * 3600.0)


def ang_sep_arcsec(p1, p2):
    cosd = np.cos(np.deg2rad(0.5 * (p1[1] + p2[1])))
    return float(np.hypot((p1[0] - p2[0]) * cosd, p1[1] - p2[1]) * 3600.0)


def main() -> None:
    registry = load_target_registry(REGISTRY_PATH)
    ctx_site = GeometryContext.decam_default()
    ctx_geo = GeometryContext.wise_coarse_default()  # Earth centre
    site = EarthLocation.from_geodetic(CTIO["longitude_deg"] * u.deg,
                                       CTIO["latitude_deg"] * u.deg,
                                       CTIO["height_m"] * u.m)
    rows = []
    worst = 0.0
    for eid in PILOT_ENDPOINTS:
        target = registry[eid]
        for mjd in EPOCHS_MJD:
            t = Time(mjd, format="mjd")
            p_site = locus_point_at_z(ctx_site, target, t, Z_TEST_AU)
            p_geo = locus_point_at_z(ctx_geo, target, t, Z_TEST_AU)
            measured = ang_sep_arcsec(p_site, p_geo)
            expected = astropy_site_shift(p_geo[0], p_geo[1], t,
                                          Z_TEST_AU, site)
            resid = abs(measured - expected)
            worst = max(worst, resid)
            rows.append({"endpoint": eid, "mjd": mjd,
                         "measured_shift_arcsec": round(measured, 4),
                         "astropy_expected_arcsec": round(expected, 4),
                         "residual_arcsec": round(resid, 4)})
            print(f"{eid:16s} mjd {mjd:8.2f}: site shift "
                  f"{measured:.4f}\" vs astropy {expected:.4f}\" "
                  f"(residual {resid:.4f}\")", flush=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({
        "z_test_au": Z_TEST_AU, "site": CTIO,
        "worst_residual_arcsec": round(worst, 4), "rows": rows}, indent=2))
    print(f"\nworst residual {worst:.4f}\" -> "
          f"{'PASS' if worst < 0.02 else 'REVIEW'} (target <= 0.02\")")
    print(f"wrote {OUT.relative_to(REPO)}")


if __name__ == "__main__":
    main()
