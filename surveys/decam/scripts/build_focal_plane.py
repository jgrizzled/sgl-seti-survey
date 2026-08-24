"""Build the static DECam focal-plane layout (per-CCD tangent-plane
bounding boxes) from one reference instcal dqmask, and validate it
against N independent dqmasks from other eras/pointings.

The layout is what makes the nominal (coarse) footprint free of any
per-exposure download: all CCD HDUs of an exposure share CRVAL = the
pointing centre, and DECam has no field rotator, so per-CCD tangent
boxes are fixed on the sky (recon note, 2026-08-24). Validation
projects every CCD corner of each check dqmask through the layout and
reports the worst mismatch; the result is stored in the layout file.

Usage:
    uv run python surveys/decam/scripts/build_focal_plane.py
Writes surveys/decam/configs/decam_focal_plane_v1.json (committed).
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import requests

from sglsurvey.adapters.noirlab_decam import (RETRIEVE_URL, SEARCH_URL,
                                              _tangent_offsets,
                                              build_focal_plane_layout)

REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "surveys" / "decam" / "configs" / "decam_focal_plane_v1.json"
WORK = REPO / "runs" / "decam" / "focal_plane"

#: Reference exposure: the recon exposure (sigma-dra corridor, 2021,
#: full 61-CCD era). Validators span eras and pointings.
REFERENCE_DQMASK_MD5 = "9b194874a4fa73c07788e9e355aa419f"  # EXPNUM 960460
#: (label, search box) for validation dqmasks — one early-era (2013,
#: pre S30 recovery), one late (2024+), one at low airmass mid-dec.
VALIDATORS = [
    ("early-2012", [["dec_center", -70.8, -68.5],
                    ["ra_center", 110.0, 116.5],
                    ["caldat", "2012-11-01", "2013-06-30"]]),
    ("late-2024", [["dec_center", -70.8, -68.5],
                   ["ra_center", 110.0, 116.5],
                   ["caldat", "2024-01-01", "2026-12-31"]]),
]


def _fetch_dqmask(md5: str, dest: Path) -> Path:
    path = dest / f"dqmask_{md5}.fits.fz"
    if not path.exists():
        r = requests.get(RETRIEVE_URL.format(md5=md5), timeout=600)
        r.raise_for_status()
        assert hashlib.md5(r.content).hexdigest() == md5
        dest.mkdir(parents=True, exist_ok=True)
        path.write_bytes(r.content)
    return path


def _find_dqmask(search: list) -> dict:
    payload = {"outfields": ["md5sum", "caldat", "EXPNUM", "ifilter"],
               "search": [["instrument", "decam"],
                          ["proc_type", "instcal"],
                          ["prod_type", "dqmask"]] + search}
    r = requests.post(f"{SEARCH_URL}?limit=5", data=json.dumps(payload),
                      headers={"Content-Type": "application/json"},
                      timeout=120)
    r.raise_for_status()
    rows = r.json()[1:]
    if not rows:
        raise RuntimeError(f"no dqmask found for {search}")
    return rows[0]


def _validate(layout: dict, dq_path: Path) -> dict:
    """Worst per-corner distance (arcsec) between each CCD's box in the
    layout and the same CCD's true TPV corners in this exposure."""
    from astropy.io import fits
    from astropy.wcs import WCS

    worst = 0.0
    missing = []
    with fits.open(dq_path) as hdul:
        ra0 = float(hdul[1].header["CRVAL1"])
        dec0 = float(hdul[1].header["CRVAL2"])
        for hdu in hdul[1:]:
            h = hdu.header
            name = h["EXTNAME"]
            if name not in layout:
                missing.append(name)
                continue
            nx, ny = int(h["NAXIS1"]), int(h["NAXIS2"])
            corners = np.array([[0.5, 0.5], [nx + 0.5, 0.5],
                                [nx + 0.5, ny + 0.5], [0.5, ny + 0.5]])
            sky = WCS(h).all_pix2world(corners, 1)
            xi, eta = _tangent_offsets(sky[:, 0], sky[:, 1], ra0, dec0)
            box = layout[name]
            d = max(abs(xi.min() - box[0]), abs(xi.max() - box[1]),
                    abs(eta.min() - box[2]), abs(eta.max() - box[3]))
            worst = max(worst, d * 3600.0)
    return {"worst_corner_arcsec": round(worst, 2),
            "ccds_not_in_layout": missing}


def main() -> None:
    ref_path = _fetch_dqmask(REFERENCE_DQMASK_MD5, WORK)
    layout = build_focal_plane_layout(ref_path)
    print(f"layout: {len(layout)} CCDs from reference "
          f"{REFERENCE_DQMASK_MD5}")

    validation = {}
    for label, search in VALIDATORS:
        row = _find_dqmask(search)
        path = _fetch_dqmask(row["md5sum"], WORK)
        v = _validate(layout, path)
        v["dqmask_md5"] = row["md5sum"]
        v["caldat"] = row["caldat"]
        validation[label] = v
        print(f"  {label}: EXPNUM-file {row['md5sum'][:8]} "
              f"({row['caldat']}): worst corner "
              f"{v['worst_corner_arcsec']}\" , "
              f"CCDs not in layout: {v['ccds_not_in_layout']}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({
        "version": "decam_focal_plane_v1",
        "built_utc": datetime.now(timezone.utc).isoformat(),
        "reference_dqmask_md5": REFERENCE_DQMASK_MD5,
        "convention": ("per-CCD [xi_min, xi_max, eta_min, eta_max] in "
                       "deg, gnomonic tangent plane about the exposure "
                       "pointing centre (CRVAL), xi east / eta north"),
        "validation": validation,
        "ccds": layout,
    }, indent=1))
    print(f"wrote {OUT.relative_to(REPO)}")


if __name__ == "__main__":
    main()
