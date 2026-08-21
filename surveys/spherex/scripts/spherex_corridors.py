"""Shared corridor configuration for the SPHEREx scripts.

Corridor keys and endpoint→corridor mapping are the shared WISE ones
(surveys/wise/scripts/wise_corridors.py). SPHEREx is all-sky, so every
corridor is in-footprint; the *ordering* puts the ZTF-inaccessible
corridors (Dec < -28, from the ZTF overlay) first. Pilot v1.0 runs
three single-star corridors from that set.
"""

import json
import sys
from collections import defaultdict
from pathlib import Path

_REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(_REPO / "surveys" / "wise" / "scripts"))
from wise_corridors import CORRIDOR_OF as _WISE_CORRIDOR_OF  # noqa: E402

PILOT_CORRIDORS = ["lalande", "gj687", "sigmadra"]
PILOT_ENDPOINTS = ["lalande-21185", "gj-687", "sigma-dra"]

_ztf_overlay = _REPO / "surveys" / "ztf" / "targets" / "overlay_v1.json"
if _ztf_overlay.exists():
    _ov = json.loads(_ztf_overlay.read_text())
    ZTF_INACCESSIBLE = [r["corridor"] for r in _ov["rows"]
                        if not r["ztf_visible"]]
    ANTIPODE = {r["corridor"]: (r["antipode_ra_deg"], r["antipode_dec_deg"])
                for r in _ov["rows"]}
else:
    ZTF_INACCESSIBLE, ANTIPODE = list(PILOT_CORRIDORS), {}

CORRIDOR_OF = dict(_WISE_CORRIDOR_OF)
MEMBERS = defaultdict(list)
for _e, _c in CORRIDOR_OF.items():
    MEMBERS[_c].append(_e)
MEMBERS = dict(MEMBERS)

#: Scale-up order: ZTF-inaccessible corridors first, then the rest.
QUEUE = list(ZTF_INACCESSIBLE) + [c for c in MEMBERS if c not in ZTF_INACCESSIBLE]
ALL_ENDPOINTS = list(CORRIDOR_OF)
