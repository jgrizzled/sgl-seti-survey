"""Shared corridor configuration for the ZTF scripts.

Corridor keys and endpoint→corridor mapping are the shared WISE ones
(sglsurvey/corridors.py); membership for ZTF is the set
of corridors graded visible by the ZTF overlay
(surveys/ztf/targets/overlay_v1.json, Dec > -28). Pilot v1.0 ran the
first three; the scale-up runs the overlay queue.
"""

import json
import sys
from collections import defaultdict
from pathlib import Path

_REPO = Path(__file__).resolve().parents[3]
from sglsurvey.corridors import CORRIDOR_OF as _WISE_CORRIDOR_OF  # noqa: E402

PILOT_CORRIDORS = ["ross128", "epsind", "proxima"]
_overlay_path = _REPO / "surveys" / "ztf" / "targets" / "overlay_v1.json"
if _overlay_path.exists():
    _ov = json.loads(_overlay_path.read_text())
    VISIBLE_CORRIDORS = [r["corridor"] for r in _ov["rows"] if r["ztf_visible"]]
    QUEUE = list(_ov["queue"])
    GRID_FLAG = {r["corridor"]: r.get("grid_flag") for r in _ov["rows"]}
else:  # pre-overlay fallback (pilot)
    VISIBLE_CORRIDORS, QUEUE, GRID_FLAG = list(PILOT_CORRIDORS), [], {}

CORRIDOR_OF = {e: c for e, c in _WISE_CORRIDOR_OF.items()
               if c in VISIBLE_CORRIDORS}

MEMBERS = defaultdict(list)
for _e, _c in CORRIDOR_OF.items():
    MEMBERS[_c].append(_e)
MEMBERS = dict(MEMBERS)

ALL_ENDPOINTS = list(CORRIDOR_OF)
ALL_CORRIDORS = list(MEMBERS)

#: Palomar P48 (hypotheses §4).
PALOMAR = {"name": "palomar-p48", "longitude_deg": -116.8650,
           "latitude_deg": 33.3563, "height_m": 1712.0}
