"""Shared corridor configuration for the Pan-STARRS1 scripts.

Corridor keys and endpoint->corridor mapping are the shared WISE ones
(surveys/wise/scripts/wise_corridors.py); membership for PS1 is the set
of corridors graded visible by the PS1 overlay
(surveys/panstarrs/targets/overlay_v1.json, Dec > -30). Pilot v1.0 ran
the first three; the scale-up runs the overlay queue.
"""

import sys
from collections import defaultdict
from pathlib import Path

_REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(_REPO / "surveys" / "wise" / "scripts"))
from wise_corridors import CORRIDOR_OF as _WISE_CORRIDOR_OF  # noqa: E402

PILOT_CORRIDORS = ["ross128", "epsind", "proxima"]
DEC_LIMIT = -30.0
_overlay_path = _REPO / "surveys" / "panstarrs" / "targets" / "overlay_v1.json"
if _overlay_path.exists():
    import json
    _ov = json.loads(_overlay_path.read_text())
    VISIBLE_CORRIDORS = [r["corridor"] for r in _ov["rows"]
                         if r["ps1_visible"] and r.get("warps", 0) > 0]
    QUEUE = list(_ov["queue"])
    OVERLAY = {r["corridor"]: r for r in _ov["rows"]}
else:  # pre-overlay fallback (pilot)
    VISIBLE_CORRIDORS, QUEUE, OVERLAY = list(PILOT_CORRIDORS), [], {}

CORRIDOR_OF = {e: c for e, c in _WISE_CORRIDOR_OF.items()
               if c in VISIBLE_CORRIDORS}

MEMBERS = defaultdict(list)
for _e, _c in CORRIDOR_OF.items():
    MEMBERS[_c].append(_e)
MEMBERS = dict(MEMBERS)

ALL_ENDPOINTS = list(CORRIDOR_OF)
ALL_CORRIDORS = list(MEMBERS)

#: PS1 Haleakala (hypotheses section 4); JPL Horizons site code F51.
HALEAKALA = {"name": "haleakala-ps1", "longitude_deg": -156.2559,
             "latitude_deg": 20.7071, "height_m": 3048.0, "mpc_code": "F51"}
