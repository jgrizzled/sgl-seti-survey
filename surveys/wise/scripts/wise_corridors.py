"""Shared corridor configuration for the WISE shakedown scripts.

One corridor per pilot system; multiple component endpoints share a
corridor (their arcs overlap at arcsecond scale). Registry v1.5.
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
    # v1.4 batch 4 + unblocked deferred
    "61-cyg-a": "cyg61", "61-cyg-b": "cyg61",
    "struve-2398-a": "struve2398", "struve-2398-b": "struve2398",
    "groombridge-34-a": "grb34", "groombridge-34-b": "grb34",
    "gj-1111": "gj1111",
    "luyten-star": "luyten",
    "kapteyn-star": "kapteyn",
    "lp-145-141": "lp145141",
    "gj-1221": "gj1221",
    "gj-9193": "gj9193",
    "gj-783": "gj783",
    "eps-ind-b": "epsind",
    "gj-11068": "gj11068",
    "wise-0855": "wise0855",
    "ez-aqr": "ezaqr",
    "luhman16-a": "luhman16", "luhman16-b": "luhman16",
    "procyon-a": "procyon", "procyon-b": "procyon",
    # v1.5 batch 5 (universal v2: 10 pc + picky-network baskets)
    "gj-876": "gj876",
    "gj-1002": "gj1002",
    "gj-832": "gj832",
    "gj-526": "gj526",
    "gj-581": "gj581",
    "gj-514": "gj514",
    "fomalhaut": "fomalhaut",
    "wolf-437": "wolf437",
    "gj-915": "gj915",
    "gj-518": "gj518",
    "gj-1276": "gj1276",
    "61-vir": "61vir",
    "gj-2012": "gj2012",
    "gj-11547": "gj11547",
    "lhs-1723": "lhs1723",
    "82-eri": "82eri",
    "gj-338-a": "gj338",
    "gj-338-b": "gj338",
    "gj-625": "gj625",
    "gj-293": "gj293",
    "gj-3306": "gj3306",
    "gj-3112": "gj3112",
    "gj-687": "gj687",
    "gj-674": "gj674",
    "gj-682": "gj682",
    "gj-251": "gj251",
    "sigma-dra": "sigmadra",
    "hd-219134": "hd219134",
    "gj-1087": "gj1087",
    "gj-318": "gj318",
    "wolf-1069": "wolf1069",
    "gj-588": "gj588",
    "gj-3512": "gj3512",
    "gj-13157": "gj13157",
    "gj-2066": "gj2066",
    "gj-367": "gj367",
    "gj-229-a": "gj229",
    "gj-667-c": "gj667",
    "ltt-1445-a": "ltt1445",
    "gj-66-a": "gj66",
    "gj-66-b": "gj66",
}

MEMBERS = defaultdict(list)
for _e, _c in CORRIDOR_OF.items():
    MEMBERS[_c].append(_e)
MEMBERS = dict(MEMBERS)

ALL_ENDPOINTS = list(CORRIDOR_OF)
ALL_CORRIDORS = list(MEMBERS)
