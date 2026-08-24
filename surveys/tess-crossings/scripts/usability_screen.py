"""Science-array usability screen (freeze amendment v1.1 input).

The coverage gate and refinement tested sector membership and cadence
counts; neither tested whether the cutout pixels are *science* pixels.
TESS CCDs expose collateral regions (leading virtual columns 1-44,
trailing 2093+, virtual rows 2049+) that TESScut serves happily with
aperture = 1 and ~zero flux. Three of the eight fetched cutouts turn
out to sit on or across the science edge — the PS1 nominal-vs-exact
footprint lesson in TESS form.

Rule (frozen here, pre-search — no signal statistic exists): the
source position must sit >= MARGIN_PX inside physical columns 45-2092
and rows 1-2048; for channel B additionally the +/-12 px ring extent
must stay inside. Writes results/usability_v1.json keyed by
(channel, target, event, sector).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from tesscut_lib import PRODUCTS  # noqa: E402

OUT = REPO / "surveys" / "tess-crossings" / "results"
SCI_X = (45, 2092)
SCI_Y = (1, 2048)
MARGIN_PX = 4          # kernel half-width + centroid slack
RING_PX = 12


def cube_for(label_prefix: str):
    for d in sorted(PRODUCTS.glob(f"{label_prefix}*")):
        p = sorted(d.glob("*.fits"))
        if p:
            return p[0]
    return None


def main():
    ref = Table.read(OUT / "coverage_refined_v1.ecsv")
    rows = []
    for r in ref:
        ch, tid = str(r["channel"]), str(r["target_id"])
        eid, sec = str(r["event_id"]), int(r["sector"])
        prefix = (f"B-{tid}-{eid[-6:]}-s{sec:04d}" if ch == "B"
                  else f"A-{tid}-s{sec:04d}")
        path = cube_for(prefix)
        if path is None:
            rows.append({"channel": ch, "target_id": tid,
                         "event_id": eid, "sector": sec,
                         "status": "no_cube"})
            continue
        with fits.open(path) as f:
            h = f[1].header
            x0 = int(h["1CRV4P"])
            y0 = int(h["2CRV4P"])
            ny, nx = f[2].data.shape
        # source sits at the cutout centre by construction
        cx, cy = x0 + nx // 2, y0 + ny // 2
        need = MARGIN_PX + (RING_PX if ch == "B" else 0)
        usable = (SCI_X[0] + need <= cx <= SCI_X[1] - need
                  and SCI_Y[0] + need <= cy <= SCI_Y[1] - need)
        src_on_sci = (SCI_X[0] <= cx <= SCI_X[1]
                      and SCI_Y[0] <= cy <= SCI_Y[1])
        rows.append({
            "channel": ch, "target_id": tid, "event_id": eid,
            "sector": sec, "cutout_x": [x0, x0 + nx - 1],
            "cutout_y": [y0, y0 + ny - 1],
            "source_ccd_xy": [cx, cy],
            "source_on_science": bool(src_on_sci),
            "status": "usable" if usable else "off_science_array",
        })
        print(rows[-1])
    (OUT / "usability_v1.json").write_text(
        json.dumps({"rule": {"science_cols": SCI_X,
                             "science_rows": SCI_Y,
                             "margin_px": MARGIN_PX,
                             "ring_extent_px_channel_B": RING_PX},
                    "rows": rows}, indent=1) + "\n")


if __name__ == "__main__":
    main()
