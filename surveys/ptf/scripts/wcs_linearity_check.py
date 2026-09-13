"""Recon check for hypotheses §7: how far the engine's linear-Jacobian
sampling (tangent-plane offsets from the stamp centre mapped through
the centre Jacobian) departs from the full TAN-SIP+PV inverse across
a 600-px PTF cutout. Samples every N-th mask cutout; reports the 99th
percentile and maximum of the displacement over a grid of offsets out
to +-300 px (the z = 550 AU locus half-span is ~190").
Writes surveys/ptf/results/wcs_linearity_check.json.

Usage: uv run python surveys/ptf/scripts/wcs_linearity_check.py [--every 20]
"""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

REPO = Path(__file__).resolve().parents[3]
MSK = REPO / "runs" / "ptf" / "products" / "msk"
OUT = REPO / "surveys" / "ptf" / "results" / "wcs_linearity_check.json"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--every", type=int, default=20)
    a = ap.parse_args()
    files = sorted(MSK.glob("cut600-*.fits"))[::a.every]
    rows = []
    offs = np.linspace(-300, 300, 13)
    gx, gy = np.meshgrid(offs, offs)
    for f in files:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with fits.open(f) as h:
                hdu = next(x for x in h if x.data is not None)
                w = WCS(hdu.header)
                ny, nx = hdu.data.shape
        x0, y0 = (nx - 1) / 2.0, (ny - 1) / 2.0
        ra0, dec0 = w.all_pix2world([[x0, y0]], 0)[0]
        cosd = np.cos(np.deg2rad(dec0))
        # Jacobian from 10" steps, as build.py does
        pp = w.all_world2pix([[ra0, dec0], [ra0 + 10 / 3600 / cosd, dec0], [ra0, dec0 + 10 / 3600]], 0)
        J = np.stack([(pp[1] - pp[0]) / 10.0, (pp[2] - pp[0]) / 10.0], axis=1)
        # true sky positions of grid pixels -> tangent offsets -> linear pixels
        sky = w.all_pix2world(np.stack([x0 + gx.ravel(), y0 + gy.ravel()], axis=1), 0)
        dra = (sky[:, 0] - ra0) * cosd * 3600.0
        ddec = (sky[:, 1] - dec0) * 3600.0
        lin = pp[0] + np.einsum("ij,nj->ni", J, np.stack([dra, ddec], axis=1))
        err = np.hypot(lin[:, 0] - (x0 + gx.ravel()), lin[:, 1] - (y0 + gy.ravel())) * 1.01
        r = np.hypot(gx.ravel(), gy.ravel())
        rows.append({"file": f.name, "max_arcsec": float(err.max()),
                     "p99_arcsec": float(np.percentile(err, 99)),
                     "max_within_200px": float(err[r <= 200].max()),
                     "sip": bool("SIP" in hdu.header.get("CTYPE1", "")),
                     "pv_terms": int(sum(1 for k in hdu.header if k.startswith("PV")))})
    mx = np.array([r["max_arcsec"] for r in rows])
    m200 = np.array([r["max_within_200px"] for r in rows])
    summary = {"n_cutouts": len(rows), "grid": "13x13 offsets to +-300 px",
               "max_arcsec": {"median": float(np.median(mx)), "p90": float(np.percentile(mx, 90)),
                              "max": float(mx.max())},
               "max_within_200px_arcsec": {"median": float(np.median(m200)),
                                           "p90": float(np.percentile(m200, 90)), "max": float(m200.max())},
               "ctype_sip_all": all(r["sip"] for r in rows),
               "rows": rows}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(summary, indent=1))
    print(json.dumps({k: v for k, v in summary.items() if k != "rows"}, indent=1))


if __name__ == "__main__":
    main()
