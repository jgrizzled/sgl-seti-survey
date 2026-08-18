"""sglseti glue: discovery envelopes and coarse locus evaluation.

This layer owns every sglseti call the pipeline makes, so that the
geometry identities (model/ephemeris/observer versions, tolerances)
are recorded once and consistently on every record.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time

from sglseti import (AstropyEphemeris, Observer, RelayRange, Role,
                     adaptive_locus)

from sglsurvey.adapters.base import ConeRegion


@dataclass(frozen=True)
class GeometryContext:
    """Pinned geometry configuration for one analysis pass."""

    model: Any
    ephemeris: Any
    observer: Any
    relay_range: RelayRange
    tolerance_arcsec: float
    padding_arcsec: float
    confidence_level: float

    @classmethod
    def wise_coarse_default(cls) -> "GeometryContext":
        """The frozen WISE shakedown configuration (hypotheses v1.0):
        550-10,000 AU, Earth-center observer, +10" padding. At the coarse
        stage the propagated 99% locus uncertainty (<= 0.1" for every
        pilot endpoint) is folded into the padding rather than sampled.
        """
        from sglseti import Tusay2022Eq57V1

        return cls(model=Tusay2022Eq57V1(), ephemeris=AstropyEphemeris(),
                   observer=Observer.earth_center(),
                   relay_range=RelayRange(550.0, 10000.0),
                   tolerance_arcsec=5.0, padding_arcsec=10.0,
                   confidence_level=0.99)

    def identities(self) -> dict[str, Any]:
        return {
            "model_id": self.model.model_id,
            "model_version": self.model.model_version,
            "ephemeris_id": getattr(self.ephemeris, "ephemeris_id",
                                    "astropy_builtin"),
            "observer_id": self.observer.observer_id,
            "relay_range_au": [self.relay_range.z_min_au,
                               self.relay_range.z_max_au],
            "tolerance_arcsec": self.tolerance_arcsec,
            "padding_arcsec": self.padding_arcsec,
            "confidence_level": self.confidence_level,
        }


def locus_radec(ctx: GeometryContext, target: Any, role: Role, t: Time,
                *, tolerance_arcsec: float | None = None,
                ) -> tuple[np.ndarray, tuple[str, ...]]:
    """Adaptive locus polyline at one epoch: (N,2) ICRS deg + warnings."""
    al = adaptive_locus(
        target=target, role=role, observation_time=t, observer=ctx.observer,
        relay_range=ctx.relay_range,
        tolerance_arcsec=tolerance_arcsec or ctx.tolerance_arcsec,
        ephemeris=ctx.ephemeris, model=ctx.model,
    )
    pts = np.array([[p.icrs_ra_deg, p.icrs_dec_deg] for p in al.points])
    return pts, al.warnings


def enclosing_cone(points_deg: np.ndarray, extra_radius_deg: float,
                   ) -> ConeRegion:
    """Smallest practical cone around a point cloud: spherical centroid
    plus max angular distance, inflated by ``extra_radius_deg``."""
    sc = SkyCoord(ra=points_deg[:, 0], dec=points_deg[:, 1], unit="deg")
    xyz = np.stack(
        [np.asarray(sc.cartesian.x), np.asarray(sc.cartesian.y),
         np.asarray(sc.cartesian.z)], axis=1)
    center = xyz.mean(axis=0)
    center /= np.linalg.norm(center)
    ra_c = float(np.rad2deg(np.arctan2(center[1], center[0]))) % 360.0
    dec_c = float(np.rad2deg(np.arcsin(np.clip(center[2], -1.0, 1.0))))
    csc = SkyCoord(ra=ra_c, dec=dec_c, unit="deg")
    radius = float(csc.separation(sc).deg.max())
    return ConeRegion(ra_deg=float(csc.ra.deg), dec_deg=float(csc.dec.deg),
                      radius_deg=radius + extra_radius_deg)


def discovery_cone(ctx: GeometryContext, target: Any, roles: Sequence[Role],
                   mjd_start: float, mjd_stop: float,
                   *, sample_days: float = 60.0,
                   envelope_tolerance_arcsec: float = 30.0) -> ConeRegion:
    """Conservative discovery cone: union of coarse locus polylines
    sampled across the time range for every requested role, inflated by
    the sampling tolerance, the inter-sample corridor drift bound, and
    the search padding."""
    mjds = np.arange(mjd_start, mjd_stop, sample_days)
    mjds = np.append(mjds, mjd_stop)
    pts = []
    for role in roles:
        for mjd in mjds:
            p, _ = locus_radec(ctx, target, role, Time(mjd, format="mjd"),
                               tolerance_arcsec=envelope_tolerance_arcsec)
            pts.append(p)
    all_pts = np.concatenate(pts, axis=0)
    # Corridor drift between samples is bounded by the annual-parallax
    # rate at z_min (~2*pi*parallax/yr) over half a sampling interval.
    parallax_deg = np.rad2deg(1.0 / ctx.relay_range.z_min_au)
    drift_deg = np.pi * parallax_deg * (sample_days / 365.25)
    inflate = (envelope_tolerance_arcsec + ctx.padding_arcsec) / 3600.0
    return enclosing_cone(all_pts, inflate + drift_deg)
