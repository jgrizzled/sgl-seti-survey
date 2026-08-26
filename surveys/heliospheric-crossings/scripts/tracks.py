"""Per-event source tracks for the sunward channels (dev driver core).

For an event and a frame time t:
- axis direction a_hat = unit vector to the star (ICRS, from the event
  row; the axis passes through the Sun toward the star);
- observer offset d = SOHO - Sun (barycentric AU); its transverse part
  perp = d - (d.a_hat)a_hat gives b_e = |perp| (impact parameter, AU)
  and p_hat = perp/b_e (the azimuth the geometry study showed is
  preserved through the lens);
- the apparent source sits at 3D point Sun + b_e * p_hat: radius
  r_rsun = b_e/R_sun from Sun center, at the position angle of p_hat —
  for BOTH channels (S1: the limb graze point; S2: the star itself,
  whose apparent offset from the Sun is exactly the same transverse
  geometry).

The unit's photometric series lives at these solar-frame (r, PA)
positions: in-window the source occupies them; baseline frames sample
the same corona-frame track (the S2 star leaves the C3 FOV days after
t_ca, so a star-fixed baseline cannot exist); PA-ring controls rotate
p_hat about a_hat by the frozen offsets.

Validation: min over the window of b_e(t) must reproduce the events
table's b_min_au (same observer table, same axis model to first
order); run this module directly to check.
"""

from __future__ import annotations

import numpy as np
from astropy.coordinates import get_body_barycentric
from astropy.time import Time

from lasco_lib import RSUN_KM, AU_KM, soho_xyz_au

RSUN_AU = RSUN_KM / AU_KM


def _unit(v: np.ndarray) -> np.ndarray:
    return v / np.linalg.norm(v)


def radec_to_vec(ra_deg: float, dec_deg: float) -> np.ndarray:
    ra, dec = np.radians(ra_deg), np.radians(dec_deg)
    return np.array([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def vec_to_radec(v: np.ndarray) -> tuple[float, float]:
    v = _unit(v)
    return (float(np.degrees(np.arctan2(v[1], v[0])) % 360.0),
            float(np.degrees(np.arcsin(v[2]))))


_geom_cache: dict[float, tuple[np.ndarray, np.ndarray]] = {}


def warm_geometry_cache(mjds) -> None:
    """Vectorized precompute of (sun, soho) barycentric positions for a
    set of epochs — source_track per-call astropy cost (~10 ms) is the
    reduce bottleneck otherwise."""
    mjds = np.asarray(sorted(set(float(m) for m in mjds)))
    t = Time(mjds, format="mjd", scale="utc")
    sun = get_body_barycentric("sun", t).xyz.to_value("AU").T
    soho = np.stack([np.interp(mjds, __import__("lasco_lib")._soho["mjd_utc"],
                               __import__("lasco_lib")._soho["xyz_au"][:, i])
                     for i in range(3)], axis=1)
    for i, m in enumerate(mjds):
        _geom_cache[round(m, 6)] = (sun[i], soho[i])


def source_track(event: dict, mjd: float, dpa_deg: float = 0.0) -> dict:
    """Apparent source for one event at one time (optionally PA-rotated
    by dpa_deg about the axis — the ring-control construction).

    event: dict with star_icrs_ra_deg/star_icrs_dec_deg (axis target
    direction for the event's link geometry).
    Returns {ra, dec, b_e_au, r_rsun}.
    """
    hit = _geom_cache.get(round(float(mjd), 6))
    if hit is not None:
        sun, soho = hit
    else:
        t = Time(mjd, format="mjd", scale="utc")
        sun = get_body_barycentric("sun", t).xyz.to_value("AU")
        soho = soho_xyz_au(mjd)
    a_hat = radec_to_vec(event["star_icrs_ra_deg"], event["star_icrs_dec_deg"])
    d = soho - sun
    perp = d - np.dot(d, a_hat) * a_hat
    b_e = float(np.linalg.norm(perp))
    p_hat = perp / b_e
    if dpa_deg:
        th = np.radians(dpa_deg)
        # rotate p_hat about a_hat (Rodrigues)
        p_hat = (p_hat * np.cos(th) + np.cross(a_hat, p_hat) * np.sin(th)
                 + a_hat * np.dot(a_hat, p_hat) * (1 - np.cos(th)))
    src = sun + b_e * p_hat
    ra, dec = vec_to_radec(src - soho)
    return {"ra": ra, "dec": dec, "b_e_au": b_e, "r_rsun": b_e / RSUN_AU}


def window_epoch_gate(b_e_au: float, camera: str,
                      ann: dict[str, list[float]]) -> bool:
    lo, hi = ann[camera]
    return lo <= b_e_au / RSUN_AU <= hi


if __name__ == "__main__":
    # validation: reproduce b_min for a few grazing events
    import json
    from pathlib import Path

    from astropy.table import Table

    repo = Path(__file__).resolve().parents[3]
    t = Table.read(repo / "crossings" / "soho_v1" / "events.ecsv")
    sun_side = (np.asarray(t["link_direction"]) == "outbound") & (np.asarray(t["side"]) == "target")
    graze = t[sun_side & (np.asarray(t["b_min_solar_radii"]) < 2.5)]
    rng = np.random.default_rng(20260825)
    rows = graze[rng.choice(len(graze), 6, replace=False)]
    print("target        t_ca              b_min(table)  b_min(track)  dPA-90 sep check")
    for ev in rows:
        e = {k: float(ev[k]) for k in ("star_icrs_ra_deg", "star_icrs_dec_deg")}
        mjd_ca = Time(float(ev["t_ca_tdb_jd"]), format="jd", scale="tdb").utc.mjd
        mjds = mjd_ca + np.linspace(-0.6, 0.6, 49)
        bs = [source_track(e, m)["b_e_au"] for m in mjds]
        bmin = min(bs) / RSUN_AU
        # ring-control sanity: dpa=90 keeps the same radius
        s0 = source_track(e, mjd_ca)
        s9 = source_track(e, mjd_ca, dpa_deg=90.0)
        print(f"{ev['target_id']:12s} {ev['t_ca_utc'][:16]}  "
              f"{float(ev['b_min_solar_radii']):9.3f}  {bmin:11.3f}  "
              f"r0={s0['r_rsun']:.3f} r90={s9['r_rsun']:.3f}")
