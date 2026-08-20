"""Shared corridor configuration for the WISE shakedown scripts.

One corridor per pilot system; multiple component endpoints share a
corridor (their arcs overlap at arcsecond scale). Registry v1.1.
"""

from collections import defaultdict

CORRIDOR_OF = {
    # v1.0 pilot
    "barnard-star": "barnard",
    "ross-154": "ross154",
    "lalande-21185": "lalande",
    "alpha-cen-a": "alphacen",
    "alpha-cen-b": "alphacen",
    "sirius-a": "sirius",
    "sirius-b": "sirius",
    # v1.1 expansion
    "proxima-cen": "proxima",
    "wolf-359": "wolf359",
    "ross-248": "ross248",
    "gj65-a": "gj65",
    "gj65-b": "gj65",
    # v1.2 batch 2 (universal list / WISE overlay queue)
    "ross-128": "ross128",
    "eps-ind-a": "epsind",
    "tau-cet": "taucet",
    "gj-54": "gj54",
    "teegarden": "teegarden",
    "lacaille-8760": "lacaille8760",
    "van-maanen": "vanmaanen",
    "gj-908": "gj908",
    "gj-784": "gj784",
    # v1.3 batch 3 (mid-confusion queue)
    "eps-eri": "epseri",
    "lacaille-9352": "lacaille9352",
    "gj-1061": "gj1061",
    "gj-12724": "gj12724",
    "wolf-1061": "wolf1061",
}

MEMBERS = defaultdict(list)
for _e, _c in CORRIDOR_OF.items():
    MEMBERS[_c].append(_e)
MEMBERS = dict(MEMBERS)

ALL_ENDPOINTS = list(CORRIDOR_OF)
ALL_CORRIDORS = list(MEMBERS)
