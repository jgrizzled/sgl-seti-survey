"""Shared corridor configuration for the ZTF pilot scripts.

Registry ids are the shared survey-agnostic registry's. Corridor keys
match surveys/wise/scripts/wise_corridors.py so cross-archive joins
are by corridor name.
"""

from collections import defaultdict

CORRIDOR_OF = {
    # v1.0 pilot (hypotheses.md §1)
    "ross-128": "ross128",
    "eps-ind-a": "epsind",
    "proxima-cen": "proxima",
}

MEMBERS = defaultdict(list)
for _e, _c in CORRIDOR_OF.items():
    MEMBERS[_c].append(_e)
MEMBERS = dict(MEMBERS)

ALL_ENDPOINTS = list(CORRIDOR_OF)
ALL_CORRIDORS = list(MEMBERS)

#: Palomar P48 (hypotheses §4).
PALOMAR = {"name": "palomar-p48", "longitude_deg": -116.8650,
           "latitude_deg": 33.3563, "height_m": 1712.0}
