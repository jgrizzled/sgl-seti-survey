"""Shared corridor configuration for the PTF/iPTF Pipeline A scripts.

Corridor keys and endpoint->corridor mapping are the shared WISE ones
(sglsurvey/corridors.py). PTF observed from the same Palomar P48
telescope as ZTF, so membership is the ZTF-visible corridor set
(surveys/ztf/targets/overlay_v1.json, Dec > -28); the PTF overlay
(surveys/ptf/targets/overlay_v1.json, built from the coarse discovery
records) grades each corridor by its actual PTF coverage and marks
which are searchable under the frozen epoch floor.
"""

import json
from collections import defaultdict
from pathlib import Path

from sglsurvey.corridors import CORRIDOR_OF as _WISE_CORRIDOR_OF

_REPO = Path(__file__).resolve().parents[3]

_ztf_overlay = _REPO / "surveys" / "ztf" / "targets" / "overlay_v1.json"
_ov = json.loads(_ztf_overlay.read_text())
VISIBLE_CORRIDORS = [r["corridor"] for r in _ov["rows"] if r["ztf_visible"]]

_ptf_overlay = _REPO / "surveys" / "ptf" / "targets" / "overlay_v1.json"
if _ptf_overlay.exists():
    _pov = json.loads(_ptf_overlay.read_text())
    SEARCHABLE_CORRIDORS = [r["corridor"] for r in _pov["rows"] if r["searchable"]]
    GRADE = {r["corridor"]: r["grade"] for r in _pov["rows"]}
else:
    SEARCHABLE_CORRIDORS, GRADE = list(VISIBLE_CORRIDORS), {}

CORRIDOR_OF = {e: c for e, c in _WISE_CORRIDOR_OF.items()
               if c in VISIBLE_CORRIDORS}

MEMBERS = defaultdict(list)
for _e, _c in CORRIDOR_OF.items():
    MEMBERS[_c].append(_e)
MEMBERS = dict(MEMBERS)

ALL_ENDPOINTS = list(CORRIDOR_OF)
ALL_CORRIDORS = list(MEMBERS)
SEARCHABLE_ENDPOINTS = [e for e, c in CORRIDOR_OF.items() if c in SEARCHABLE_CORRIDORS]

#: Palomar P48 (same site as ZTF; GeometryContext.ztf_default()).
PALOMAR = {"name": "palomar-p48", "longitude_deg": -116.8650,
           "latitude_deg": 33.3563, "height_m": 1712.0}
