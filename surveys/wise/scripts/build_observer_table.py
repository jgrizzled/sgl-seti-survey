"""Step C (observer): tabulate the WISE spacecraft barycentric position
per frame from the L1b cutout headers (hypotheses v2.0 §1.3).

SUN2SCX/Y/Z is the Sun-to-spacecraft vector in AU (J2000 equatorial,
verified below against astropy's Earth position and the header's
Earth-to-spacecraft vector); adding the Sun's barycentric position
from the pinned astropy ephemeris gives the SSB position that sglseti
observers return. Output: runs/wise/v2/observer/wise_sc_ephemeris.npz
(mjd_utc, xyz_au) with its content hash, used by every v2 geometry
call through ``sglsurvey.geometry.register_wise_spacecraft_observer``.

Usage: uv run python surveys/wise/scripts/build_observer_table.py
"""

from __future__ import annotations

import hashlib
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np
from astropy.coordinates import get_body_barycentric
from astropy.io import fits
from astropy.time import Time
from astropy.utils.exceptions import AstropyWarning

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v2common as C  # noqa: E402

OUT_DIR = C.RUN_DIR / "observer"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    t0 = time.monotonic()
    files = sorted(C.CUT_DIR.glob("cut-*-int-1b.fits"))
    print(f"{len(files)} int cutouts", flush=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AstropyWarning)
        for i, f in enumerate(files):
            h = fits.getheader(f)
            try:
                rows.append((float(h["MJD_OBS"]), float(h["SUN2SCX"]),
                             float(h["SUN2SCY"]), float(h["SUN2SCZ"]),
                             float(h["ERTH2SCX"]), float(h["ERTH2SCY"]),
                             float(h["ERTH2SCZ"])))
            except KeyError as exc:
                print(f"  {f.name}: missing {exc}", flush=True)
            if (i + 1) % 10000 == 0:
                print(f"  {i + 1}/{len(files)} ({(i + 1) / (time.monotonic() - t0):.0f}/s)",
                      flush=True)
    arr = np.array(rows)
    arr = arr[np.argsort(arr[:, 0])]
    # de-duplicate frames shared between bands (same MJD)
    _, keep = np.unique(np.round(arr[:, 0], 7), return_index=True)
    arr = arr[keep]
    mjd = arr[:, 0]
    t = Time(mjd, format="mjd", scale="utc")
    sun = get_body_barycentric("sun", t)
    earth = get_body_barycentric("earth", t)
    sun_xyz = np.stack([sun.x.to_value("AU"), sun.y.to_value("AU"), sun.z.to_value("AU")], 1)
    earth_xyz = np.stack([earth.x.to_value("AU"), earth.y.to_value("AU"), earth.z.to_value("AU")], 1)
    sc = sun_xyz + arr[:, 1:4]
    # frame check: Sun + SUN2SC - ERTH2SC should equal astropy's Earth
    resid_km = np.linalg.norm(sc - arr[:, 4:7] - earth_xyz, axis=1) * 1.495978707e8
    geo_km = np.linalg.norm(arr[:, 4:7], axis=1) * 1.495978707e8
    np.savez_compressed(OUT_DIR / "wise_sc_ephemeris.npz", mjd_utc=mjd, xyz_au=sc,
                        geocentric_km=geo_km)
    digest = "sha256:" + hashlib.sha256((OUT_DIR / "wise_sc_ephemeris.npz").read_bytes()).hexdigest()
    summary = {
        "n_frames": int(len(mjd)), "mjd_range": [float(mjd.min()), float(mjd.max())],
        "geocentric_radius_km": {"median": float(np.median(geo_km)),
                                 "min": float(geo_km.min()), "max": float(geo_km.max())},
        "frame_check_residual_km": {"median": float(np.median(resid_km)),
                                    "p99": float(np.percentile(resid_km, 99)),
                                    "max": float(resid_km.max())},
        "locus_shift_bound_mas_at_550au": float(np.max(geo_km) / 8.2281e10 * 206264806.0),
        "table_sha256": digest, "ephemeris": "astropy builtin (get_body_barycentric)",
    }
    (OUT_DIR / "summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
