"""Solar Orbiter sunward-channel geometry for the SoloHI observer-geometry
pass (plan §5.23).

Same construction as `psp_geometry` (WISPR) for a general observer:
with the star's unit vector a_hat, the observer offset d = SolO - Sun,
transverse part perp = d - (d.a_hat) a_hat, the in-beam condition is
b_e = |perp| <= rung radius and the apparent source is the fixed ICRS
direction

* S2 (uplink past the Sun; inbound, anti-target side): the star itself;
* S1 (downlink post-lens; outbound, target side): the star's antipode.

The source's elongation from Sun centre is eps = atan(b_e / |r_along|).

SoloHI field-of-view model, **measured from L2 headers** (this pass,
`runs/solohi-crossings/recon/`, 48 frames at 12 epochs 2021-12 ->
2026-04, r 0.29-1.01 AU, spacecraft roll -8 to +13 deg): the four
detector tiles mapped through the celestial ('A') WCS into the frame
u = Sun direction, t = prograde tangential (Horizons velocity), n = t x u
(orbit normal). Ram longitude lam = atan2(s.t, s.u), orbit latitude
bet = asin(s.n). The tile edges are constant-lam / constant-bet lines
to ~0.1 deg and are **stable in this orbit-plane frame** at every epoch
(the spacecraft rolls to keep the mosaic on the orbital plane; in the
solar-north / helioprojective frame the same edges wander by +-5 deg),
with ~1-3 deg pointing offsets on a few frames (2024-04-02, 2026-02-01,
2026-04-05) — a per-frame WCS is used at the survey stage; the census
uses the nominal tiles. **The field is on the anti-ram side** (lam < 0;
helioprojective longitude -5 to -45 deg, solar east), the opposite
side from WISPR, so the pre-t_ca half of every crossing window is the
visible one. Tiles: 1 (inner, south), 2 (inner, north), 3 (outer,
north), 4 (outer, south); inner edge 5.2 deg from Sun centre on the
orbital plane, outer edge 45 deg; 960 x 1024 2x2-binned pixels of
0.0206 deg = 74 arcsec; 5-exposure sums, XPOSURE ~29 s (r 0.38 AU).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.coordinates import get_body_barycentric
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
OBS_TAB = REPO / "crossings" / "observers" / "solo_sc_ephemeris.npz"

RSUN_KM = 695_700.0
AU_KM = 1.495978707e8
RSUN_AU = RSUN_KM / AU_KM

#: nominal tiles in (ram longitude, orbit latitude), degrees — header-
#: measured medians of the 2021-12 .. 2025-10 frames
TILES = {
    "1": {"lam_deg": (-25.6, -5.2), "bet_deg": (-19.3, -0.2), "descriptors": ("solohi-1ft",)},
    "2": {"lam_deg": (-24.3, -5.1), "bet_deg": (0.2, 20.5), "descriptors": ("solohi-2ft",)},
    "3": {"lam_deg": (-45.0, -24.7), "bet_deg": (0.3, 19.3), "descriptors": ("solohi-3ft", "solohi-3fg")},
    "4": {"lam_deg": (-45.1, -26.0), "bet_deg": (-20.5, -0.1), "descriptors": ("solohi-4ft", "solohi-4fg")},
}
INNER_EDGE_DEG = 5.2
OUTER_EDGE_DEG = 45.0

_obs = {k: np.array(v) for k, v in np.load(OBS_TAB).items()}
_MJD, _XYZ = _obs["mjd_utc"], _obs["xyz_au"]
_VEL = np.gradient(_XYZ, _MJD, axis=0)          # AU/d on the 10-min table


def solo_xyz_au(mjd) -> np.ndarray:
    mjd = np.atleast_1d(np.asarray(mjd, float))
    return np.stack([np.interp(mjd, _MJD, _XYZ[:, i]) for i in range(3)], axis=-1)


def solo_vel_au_d(mjd) -> np.ndarray:
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
    """(b_e_au, r_along_au signed, r_helio_au) of SolO vs the Sun-star axis."""
    a = radec_to_vec(ev["star_icrs_ra_deg"], ev["star_icrs_dec_deg"])
    d = solo_xyz_au(mjd) - sun_xyz_au(mjd)
    along = d @ a
    perp = d - np.outer(along, a)
    return np.linalg.norm(perp, axis=-1), along, np.linalg.norm(d, axis=-1)


def heliocentric_r_au(mjd) -> np.ndarray:
    return np.linalg.norm(solo_xyz_au(mjd) - sun_xyz_au(mjd), axis=-1)


def ram_frame(mjd):
    """(u, t, n) unit vectors per epoch: Sun direction, prograde
    tangential, orbit normal."""
    mjd = np.atleast_1d(np.asarray(mjd, float))
    u = sun_xyz_au(mjd) - solo_xyz_au(mjd)
    u /= np.linalg.norm(u, axis=-1, keepdims=True)
    v = solo_vel_au_d(mjd)
    t = v - np.sum(v * u, axis=-1, keepdims=True) * u
    t /= np.linalg.norm(t, axis=-1, keepdims=True)
    n = np.cross(t, u)
    return u, t, n


def solohi_coords(s_hat, mjd):
    """(elongation, ram longitude, orbit latitude) in degrees of a fixed
    direction s_hat as seen from Solar Orbiter at epochs mjd."""
    u, t, n = ram_frame(mjd)
    s = np.broadcast_to(np.asarray(s_hat, float), u.shape)
    su, st, sn = (np.sum(s * x, axis=-1) for x in (u, t, n))
    eps = np.degrees(np.arccos(np.clip(su, -1, 1)))
    lam = np.degrees(np.arctan2(st, su))
    bet = np.degrees(np.arcsin(np.clip(sn, -1, 1)))
    return eps, lam, bet


def in_tile(lam, bet, tile: dict) -> np.ndarray:
    lo, hi = tile["lam_deg"]
    b0, b1 = tile["bet_deg"]
    return (lam >= lo) & (lam <= hi) & (bet >= b0) & (bet <= b1)


def in_solohi(lam, bet) -> np.ndarray:
    out = np.zeros(np.shape(lam), bool)
    for tile in TILES.values():
        out |= in_tile(lam, bet, tile)
    return out


def tile_of(lam, bet) -> np.ndarray:
    """Tile id per epoch ('' when outside the mosaic)."""
    out = np.full(np.shape(lam), "", dtype="<U1")
    for k, tile in TILES.items():
        out[in_tile(lam, bet, tile) & (out == "")] = k
    return out


def perihelia() -> list[dict]:
    """Every perihelion pass of the observer table: number, epoch,
    distance, and the bracketing aphelia (orbit = aphelion-to-aphelion)."""
    sub = slice(None, None, 6)
    ts = Time(_MJD[sub], format="mjd", scale="utc")
    s = get_body_barycentric("sun", ts).xyz.to_value("AU").T
    sun = np.stack([np.interp(_MJD, _MJD[sub], s[:, i]) for i in range(3)], -1)
    r = np.linalg.norm(_XYZ - sun, axis=1)
    lo = np.nonzero((r[1:-1] < r[:-2]) & (r[1:-1] <= r[2:]))[0] + 1
    hi = np.nonzero((r[1:-1] > r[:-2]) & (r[1:-1] >= r[2:]))[0] + 1
    out = []
    for k, i in enumerate(lo):
        prev = hi[hi < i]
        nxt = hi[hi > i]
        j0 = int(prev[-1]) if len(prev) else 0
        j1 = int(nxt[0]) if len(nxt) else len(r) - 1
        out.append({"orbit": f"P{k + 1:02d}",
                    "mjd_perihelion": float(_MJD[i]),
                    "utc_perihelion": Time(_MJD[i], format="mjd").iso[:16],
                    "q_au": float(r[i]), "q_rsun": float(r[i] / RSUN_AU),
                    "Q_before_au": float(r[j0]), "Q_after_au": float(r[j1]),
                    "mjd_start": float(_MJD[j0]), "mjd_end": float(_MJD[j1]),
                    "utc_start": Time(_MJD[j0], format="mjd").iso[:10],
                    "utc_end": Time(_MJD[j1], format="mjd").iso[:10],
                    "period_days": float(_MJD[j1] - _MJD[j0])})
    return out


def planet_directions(mjd, bodies=("mercury", "venus", "earth", "mars",
                                   "jupiter", "saturn")) -> dict:
    t = Time(np.atleast_1d(mjd), format="mjd", scale="utc")
    obs = solo_xyz_au(mjd)
    out = {}
    for b in bodies:
        v = get_body_barycentric(b, t).xyz.to_value("AU").T - obs
        out[b] = v / np.linalg.norm(v, axis=-1, keepdims=True)
    return out
