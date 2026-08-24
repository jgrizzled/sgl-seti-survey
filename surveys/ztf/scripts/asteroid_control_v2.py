"""ZTF v2 positive control: re-score the v1 asteroid (60000) control
cutouts (runs/ztf/control_v1) with the v2 rule (hypotheses v2.0 §1.5).
The predicted magnitude uses Horizons V with the v1 solar colours
(g − V = +0.25, V − r = +0.19). Output: runs/ztf/v2/control/60000/summary.json.

Usage: uv run python surveys/ztf/scripts/asteroid_control_v2.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.wcs import WCS

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from profile import PROFILE as P  # noqa: E402

from sglsurvey.adapters.irsa_ztf import ZtfExactFootprint  # noqa: E402
from sglsurvey.photometry import build_flux_map_ztf  # noqa: E402
from sglsurvey.control import rescore  # noqa: E402

V1 = P.v1_run_dir / "control_v1"
COLOUR = {"zg": 0.25, "zr": -0.19, "zi": -0.35}


def main():
    recs = [json.loads(l) for l in (V1 / "frames.jsonl").open()]
    ok = [r for r in recs if r["status"] == "ok"]
    frames, maps, positions, vpred = [], [], [], {b: [] for b in P.bands}
    for r in ok:
        nk = r["native_key"]
        stem = (f"ztf_{nk['filefracday']}_{nk['field']:06d}_{nk['filtercode']}_c{nk['ccdid']:02d}_o_q{nk['qid']}")
        sci = V1 / "cut" / f"cut120-{stem}_sciimg.fits"
        msk = V1 / "cut" / f"cut120-{stem}_mskimg.fits"
        diff = V1 / "cut" / f"cut120-{stem}_scimrefdiffimg.fits.fz"
        if not (sci.exists() and msk.exists()):
            continue
        fm = build_flux_map_ztf(sci, diff if diff.exists() else None, msk, ZtfExactFootprint.FATAL_MASK,
                                r["band"], r["mjd"], keep_inputs=False)
        ny, nx = fm.flux.shape
        sky = WCS(fm.header).wcs_pix2world([[(nx - 1) / 2.0, (ny - 1) / 2.0]], 0)[0]
        frames.append(r); maps.append(fm); positions.append((float(sky[0]), float(sky[1])))
        vpred[r["band"]].append(r["v_pred"])
    pred = {b: float(-2.5 * np.log10(np.mean(10 ** (-0.4 * np.array(v)))) + COLOUR[b]) for b, v in vpred.items() if v}
    rescore(P, frames, maps, positions, P.run_dir / "control" / "60000", "60000", predicted_mag=pred,
            extra={"source": "v1 control cutouts (runs/ztf/control_v1)", "colours": COLOUR})
    rescore(P, frames, maps, positions, P.run_dir / "control" / "60000_noclip", "60000 (no single-epoch clip)",
            predicted_mag=pred, extra={"source": "v1 control cutouts", "colours": COLOUR,
                                       "note": "the frozen clip removes the frames where a bright mover is individually detected"},
            clip_sigma=None)


if __name__ == "__main__":
    main()
