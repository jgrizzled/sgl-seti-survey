"""STEREO-A sunward-channel geometry for the HI-1 survey.

Both sunward apparent sources are fixed ICRS directions (derived in the
LASCO tracks module and re-derived here for a general observer):

* S2 (uplink past the Sun; inbound, anti-target side): the target star
  itself, PM-propagated (parallax << the 72" pixel).
* S1 (downlink post-lens; outbound, target side) at the 0.1 AU rung:
  the relay direction = the star's antipode (the beam passes the Sun
  unlensed at b >> R_sun; the graze-point construction of the LASCO
  freeze reduces to the antipode direction exactly: src - obs = -s a_hat).

The in-beam condition is b_e(t) <= rung radius with b_e the transverse
offset of the spacecraft from the Sun-star axis; the HI-1 visibility
condition is the source's helioprojective position inside the camera
footprint (measured from real headers: HPLN in [-24.15, -3.90] deg,
HPLT in [-10.29, 9.97] deg for the 2010 test frame; the primary WCS is
axis-aligned to HPLN/HPLT, PC = identity).

Helioprojective-cartesian (Thompson 2006): with the observer->Sun unit
vector u, the solar rotation axis n (ICRS RA 286.13, Dec 63.87),
solar-east e = n x u (normalized), west w = -e, north p = u x e;
for a source direction s: HPLN = atan2(s.w, s.u), HPLT = asin(s.p).
Validated against the test-frame header WCS (see hi_census.py).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.coordinates import get_body_barycentric
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
OBS_TAB = REPO / "crossings" / "observers" / "stereoa_sc_ephemeris.npz"

RSUN_KM = 695_700.0
AU_KM = 1.495978707e8
RSUN_AU = RSUN_KM / AU_KM

# solar rotation axis, ICRS (IAU 2009 / Seidelmann et al.)
_SUN_POLE = np.array([np.cos(np.radians(63.87)) * np.cos(np.radians(286.13)),
                      np.cos(np.radians(63.87)) * np.sin(np.radians(286.13)),
                      np.sin(np.radians(63.87))])

# HI-1A footprint in helioprojective degrees (2010-06-15 test frame; an
# 8 px margin inside the CCD edge is applied by `in_fov`)
HI1_HPLN = (-24.152, -3.896)
HI1_HPLT = (-10.288, 9.968)
HI1_EDGE_MARGIN_DEG = 8 * 0.01998
HI1_HALF_W = ((HI1_HPLN[1] - HI1_HPLN[0]) / 2, (HI1_HPLT[1] - HI1_HPLT[0]) / 2)
HI1_CENTRE_EAST = (-14.0, -0.16)

# Pointing eras measured from the per-day headers (coverage stage,
# `results/coverage_v1.json`): the camera centre sits at HPLN -14 deg
# (east of the Sun: fixed directions approach conjunction) except for
# the rolled post-conjunction era at +14 deg (west: directions recede),
# with the 2014-08-19 -> 2015-11-16 safe-mode/conjunction gap between.
# MJD boundaries; HPLT drifts seasonally by +-1.9 deg (per-frame WCS).
POINTING_ERAS = [
    (54075.0, 56888.5, -14.0),   # 2006-12-06 .. 2014-08-19  east
    (56888.5, 57343.5, None),    # gap (no synoptic frames)
    (57343.5, 60172.5, +14.0),   # 2015-11-17 .. 2023-08-15  west (rolled 180)
    (60172.5, 62000.0, -14.0),   # 2023-08-16 ..             east
]

_obs = {k: np.array(v) for k, v in np.load(OBS_TAB).items()}   # eager (fork-safe)


def nominal_centre(mjd) -> tuple[float | None, float]:
    """(HPLN, HPLT) of the CCD centre from the era model; HPLN None in
    the gap. The coverage stage replaces this with per-day headers."""
    for lo, hi, c in POINTING_ERAS:
        if lo <= mjd < hi:
            return (c, HI1_CENTRE_EAST[1] if c is not None else 0.0)
    return (None, 0.0)



def stereoa_xyz_au(mjd) -> np.ndarray:
    mjd = np.atleast_1d(np.asarray(mjd, float))
    return np.stack([np.interp(mjd, _obs["mjd_utc"], _obs["xyz_au"][:, i])
                     for i in range(3)], axis=-1)


def sun_xyz_au(mjd) -> np.ndarray:
    t = Time(np.atleast_1d(mjd), format="mjd", scale="utc")
    return get_body_barycentric("sun", t).xyz.to_value("AU").T


def radec_to_vec(ra_deg, dec_deg) -> np.ndarray:
    ra, dec = np.radians(ra_deg), np.radians(dec_deg)
    return np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra),
                     np.sin(dec)], axis=-1)


def vec_to_radec(v) -> tuple[np.ndarray, np.ndarray]:
    v = np.asarray(v, float)
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    v = v / n
    return (np.degrees(np.arctan2(v[..., 1], v[..., 0])) % 360.0,
            np.degrees(np.arcsin(np.clip(v[..., 2], -1, 1))))


def source_direction(ev: dict, channel: str) -> np.ndarray:
    """Unit vector of the apparent source: the star (S2) or its antipode
    (S1). `ev` carries star_icrs_ra_deg/star_icrs_dec_deg."""
    a = radec_to_vec(ev["star_icrs_ra_deg"], ev["star_icrs_dec_deg"])
    return a if channel == "S2" else -a


def impact_parameter_au(ev: dict, mjd) -> np.ndarray:
    """b_e(t): transverse offset of STEREO-A from the Sun-star axis."""
    a = radec_to_vec(ev["star_icrs_ra_deg"], ev["star_icrs_dec_deg"])
    d = stereoa_xyz_au(mjd) - sun_xyz_au(mjd)
    perp = d - np.outer(d @ a, a)
    return np.linalg.norm(perp, axis=-1)


def helioprojective(s_hat, mjd) -> tuple[np.ndarray, np.ndarray]:
    """(HPLN, HPLT) in degrees of a fixed direction s_hat at epochs mjd,
    as seen from STEREO-A."""
    mjd = np.atleast_1d(np.asarray(mjd, float))
    u = sun_xyz_au(mjd) - stereoa_xyz_au(mjd)
    u /= np.linalg.norm(u, axis=-1, keepdims=True)
    e = np.cross(np.broadcast_to(_SUN_POLE, u.shape), u)
    e /= np.linalg.norm(e, axis=-1, keepdims=True)
    p = np.cross(u, e)
    s = np.broadcast_to(np.asarray(s_hat, float), u.shape)
    sx = -np.sum(s * e, axis=-1)      # west component
    sy = np.sum(s * p, axis=-1)       # north component
    sz = np.sum(s * u, axis=-1)
    return np.degrees(np.arctan2(sx, sz)), np.degrees(np.arcsin(np.clip(sy, -1, 1)))


def elongation_deg(s_hat, mjd) -> np.ndarray:
    mjd = np.atleast_1d(np.asarray(mjd, float))
    u = sun_xyz_au(mjd) - stereoa_xyz_au(mjd)
    u /= np.linalg.norm(u, axis=-1, keepdims=True)
    s = np.broadcast_to(np.asarray(s_hat, float), u.shape)
    return np.degrees(np.arccos(np.clip(np.sum(s * u, axis=-1), -1, 1)))


def in_fov(hpln, hplt, centre=HI1_CENTRE_EAST,
           margin_deg: float = HI1_EDGE_MARGIN_DEG) -> np.ndarray:
    """HI-1A footprint test about a CCD centre (HPLN, HPLT): the census
    uses the era model, the coverage stage each day's header centre."""
    if centre[0] is None:
        return np.zeros(np.shape(hpln), bool)
    return ((np.abs(hpln - centre[0]) <= HI1_HALF_W[0] - margin_deg)
            & (np.abs(hplt - centre[1]) <= HI1_HALF_W[1] - margin_deg))


def in_fov_era(hpln, hplt, mjd, margin_deg: float = HI1_EDGE_MARGIN_DEG) -> np.ndarray:
    """Footprint test with the era-model centre evaluated per epoch."""
    mjd = np.atleast_1d(mjd)
    out = np.zeros(len(mjd), bool)
    for lo, hi, c in POINTING_ERAS:
        m = (mjd >= lo) & (mjd < hi)
        if c is not None and m.any():
            out[m] = in_fov(np.asarray(hpln)[m], np.asarray(hplt)[m],
                            (c, HI1_CENTRE_EAST[1]), margin_deg)
    return out


def window_half_days(b_min_au: float, r_au: float, v_kms: float) -> float:
    if r_au <= b_min_au:
        return 0.0
    return float(np.sqrt(r_au**2 - b_min_au**2) * AU_KM / (v_kms * 86400.0))


def from_helioprojective(hpln_deg, hplt_deg, mjd) -> np.ndarray:
    """Inverse of `helioprojective`: the ICRS unit vector of the direction
    at (HPLN, HPLT) as seen from STEREO-A at one epoch."""
    u = sun_xyz_au(mjd)[0] - stereoa_xyz_au(mjd)[0]
    u /= np.linalg.norm(u)
    e = np.cross(_SUN_POLE, u)
    e /= np.linalg.norm(e)
    p = np.cross(u, e)
    w = -e
    ln, lt = np.radians(hpln_deg), np.radians(hplt_deg)
    return np.sin(ln) * np.cos(lt) * w + np.sin(lt) * p + np.cos(ln) * np.cos(lt) * u


def planet_directions(mjd, bodies=("mercury", "venus", "earth", "moon", "mars",
                                   "jupiter", "saturn")) -> dict:
    """ICRS unit vectors of solar-system bodies as seen from STEREO-A,
    vectorized over epochs: {body: (n, 3)}."""
    t = Time(np.atleast_1d(mjd), format="mjd", scale="utc")
    obs = stereoa_xyz_au(mjd)
    out = {}
    for b in bodies:
        v = get_body_barycentric(b, t).xyz.to_value("AU").T - obs
        out[b] = v / np.linalg.norm(v, axis=-1, keepdims=True)
    return out
