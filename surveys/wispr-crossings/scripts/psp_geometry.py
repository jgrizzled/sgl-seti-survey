"""PSP sunward-channel geometry for the WISPR observer-geometry pass
(plan §5.15 O5).

Same construction as the LASCO tracks module and STEREO `hi_geometry`,
for a general observer: with the star's unit vector a_hat, the observer
offset d = PSP - Sun, transverse part perp = d - (d.a_hat)a_hat, the
in-beam condition is b_e = |perp| <= rung radius and the apparent
source is the fixed ICRS direction

* S2 (uplink past the Sun; inbound, anti-target side): the star itself;
* S1 (downlink post-lens; outbound, target side): the star's antipode
  (src - obs = -(d.a_hat) a_hat exactly for the source at Sun + b_e p_hat,
  so the grazing-rung graze point and the 0.1 AU unlensed relay
  direction coincide in direction).

Either way the source's elongation from Sun center is
eps(t) = atan(b_e(t) / |r_along(t)|), r_along = d.a_hat — the inner-
heliosphere observer sees the same transverse offset at a far larger
angle than a 1 AU observer (b_e/1 AU -> b_e/0.05 AU).

WISPR field-of-view model: the WISPR Data Users Guide v5 (NRL, Sep
2025; `runs/wispr-crossings/recon/`) gives rectangular fields 40 deg
and 58 deg on a side, the inner telescope's sunward edge at 13.5 deg
from Sun center (from the spacecraft pointing vector, which is Sun
center whenever r < 0.25 AU) and the outer edge at 108.5 deg with a
3 deg overlap — WISPR-I 13.5-53.5, WISPR-O 50.5-108.5 deg. The fields
sit on the spacecraft's ram side, centred on the orbital plane
(Vourlidas et al. 2016; half-heights 20 / 26.5 deg here are the
literature values — the side and the vertical extent are recon
verification items from real L2 headers, which carry the S/C HCI
velocity and a ZPN WCS). Frame used here: u = Sun direction, t = prograde
tangential (velocity component perpendicular to u), n = t x u
(orbit normal, along the angular momentum); "ram longitude" lam = atan2(s.t, s.u) (positive on the
ram side), "orbit latitude" bet = asin(s.n). A source is in WISPR-I
when 13.5 <= lam <= 53.5 and |bet| <= 20, in WISPR-O when
50.5 <= lam <= 108.5 and |bet| <= 26.5. Encounters are the r < 0.25 AU arcs (the mission's
definition of the WISPR observing periods).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.coordinates import get_body_barycentric
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
OBS_TAB = REPO / "crossings" / "observers" / "psp_sc_ephemeris.npz"

RSUN_KM = 695_700.0
AU_KM = 1.495978707e8
RSUN_AU = RSUN_KM / AU_KM

ENCOUNTER_R_AU = 0.25
WISPR_I = {"lam_deg": (13.5, 53.5), "half_height_deg": 20.0}
WISPR_O = {"lam_deg": (50.5, 108.5), "half_height_deg": 26.5}

_obs = {k: np.array(v) for k, v in np.load(OBS_TAB).items()}   # eager
_MJD, _XYZ = _obs["mjd_utc"], _obs["xyz_au"]
# velocity by central differences on the 10-min table (AU/d)
_VEL = np.gradient(_XYZ, _MJD, axis=0)


def psp_xyz_au(mjd) -> np.ndarray:
    mjd = np.atleast_1d(np.asarray(mjd, float))
    return np.stack([np.interp(mjd, _MJD, _XYZ[:, i]) for i in range(3)], axis=-1)


def psp_vel_au_d(mjd) -> np.ndarray:
    mjd = np.atleast_1d(np.asarray(mjd, float))
    return np.stack([np.interp(mjd, _MJD, _VEL[:, i]) for i in range(3)], axis=-1)


def sun_xyz_au(mjd) -> np.ndarray:
    t = Time(np.atleast_1d(mjd), format="mjd", scale="utc")
    return get_body_barycentric("sun", t).xyz.to_value("AU").T


def radec_to_vec(ra_deg, dec_deg) -> np.ndarray:
    ra, dec = np.radians(ra_deg), np.radians(dec_deg)
    return np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra),
                     np.sin(dec)], axis=-1)


def vec_to_radec(v) -> tuple[np.ndarray, np.ndarray]:
    v = np.asarray(v, float)
    v = v / np.linalg.norm(v, axis=-1, keepdims=True)
    return (np.degrees(np.arctan2(v[..., 1], v[..., 0])) % 360.0,
            np.degrees(np.arcsin(np.clip(v[..., 2], -1, 1))))


def source_direction(ev: dict, channel: str) -> np.ndarray:
    a = radec_to_vec(ev["star_icrs_ra_deg"], ev["star_icrs_dec_deg"])
    return a if channel == "S2" else -a


def axis_geometry(ev: dict, mjd):
    """(b_e_au, r_along_au signed, r_helio_au) of PSP vs the Sun-star axis."""
    a = radec_to_vec(ev["star_icrs_ra_deg"], ev["star_icrs_dec_deg"])
    d = psp_xyz_au(mjd) - sun_xyz_au(mjd)
    along = d @ a
    perp = d - np.outer(along, a)
    return np.linalg.norm(perp, axis=-1), along, np.linalg.norm(d, axis=-1)


def heliocentric_r_au(mjd) -> np.ndarray:
    return np.linalg.norm(psp_xyz_au(mjd) - sun_xyz_au(mjd), axis=-1)


def ram_frame(mjd):
    """(u, t, n) unit vectors per epoch: Sun direction, prograde
    tangential, orbit normal."""
    mjd = np.atleast_1d(np.asarray(mjd, float))
    u = sun_xyz_au(mjd) - psp_xyz_au(mjd)
    u /= np.linalg.norm(u, axis=-1, keepdims=True)
    v = psp_vel_au_d(mjd)
    t = v - np.sum(v * u, axis=-1, keepdims=True) * u
    t /= np.linalg.norm(t, axis=-1, keepdims=True)
    n = np.cross(t, u)          # orbital angular-momentum direction (north)
    return u, t, n


def wispr_coords(s_hat, mjd):
    """(elongation, ram longitude, orbit latitude) in degrees of a fixed
    direction s_hat as seen from PSP at epochs mjd."""
    u, t, n = ram_frame(mjd)
    s = np.broadcast_to(np.asarray(s_hat, float), u.shape)
    su, st, sn = (np.sum(s * x, axis=-1) for x in (u, t, n))
    eps = np.degrees(np.arccos(np.clip(su, -1, 1)))
    lam = np.degrees(np.arctan2(st, su))
    bet = np.degrees(np.arcsin(np.clip(sn, -1, 1)))
    return eps, lam, bet


def in_wispr(lam, bet, cam: dict) -> np.ndarray:
    lo, hi = cam["lam_deg"]
    return (lam >= lo) & (lam <= hi) & (np.abs(bet) <= cam["half_height_deg"])


def encounters(r_max_au: float = ENCOUNTER_R_AU) -> list[dict]:
    """Perihelion passes from the observer table: number, perihelion
    epoch/distance, and the r < r_max arc (WISPR observing period)."""
    sub = slice(None, None, 6)          # hourly Sun samples, interpolated
    ts = Time(_MJD[sub], format="mjd", scale="utc")
    s = get_body_barycentric("sun", ts).xyz.to_value("AU").T
    sun = np.stack([np.interp(_MJD, _MJD[sub], s[:, i]) for i in range(3)], -1)
    r = np.linalg.norm(_XYZ - sun, axis=1)
    loc = np.nonzero((r[1:-1] < r[:-2]) & (r[1:-1] <= r[2:]))[0] + 1
    out = []
    for i in loc:
        if r[i] >= r_max_au:
            continue
        j0 = i
        while j0 > 0 and r[j0 - 1] < r_max_au:
            j0 -= 1
        j1 = i
        while j1 < len(r) - 1 and r[j1 + 1] < r_max_au:
            j1 += 1
        out.append({"encounter": f"E{len(out) + 1:02d}",
                    "mjd_perihelion": float(_MJD[i]),
                    "utc_perihelion": Time(_MJD[i], format="mjd").iso[:16],
                    "q_au": float(r[i]), "q_rsun": float(r[i] / RSUN_AU),
                    "mjd_start": float(_MJD[j0]), "mjd_end": float(_MJD[j1]),
                    "utc_start": Time(_MJD[j0], format="mjd").iso[:10],
                    "utc_end": Time(_MJD[j1], format="mjd").iso[:10],
                    "duration_days": float(_MJD[j1] - _MJD[j0])})
    return out


def window_half_days(b_min_au: float, r_au: float, v_kms: float) -> float:
    if r_au <= b_min_au:
        return 0.0
    return float(np.sqrt(r_au**2 - b_min_au**2) * AU_KM / (v_kms * 86400.0))


def planet_directions(mjd, bodies=("mercury", "venus", "earth", "mars",
                                   "jupiter", "saturn")) -> dict:
    t = Time(np.atleast_1d(mjd), format="mjd", scale="utc")
    obs = psp_xyz_au(mjd)
    out = {}
    for b in bodies:
        v = get_body_barycentric(b, t).xyz.to_value("AU").T - obs
        out[b] = v / np.linalg.norm(v, axis=-1, keepdims=True)
    return out
