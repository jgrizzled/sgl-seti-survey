"""Shared corridor configuration for the DECam scripts.

Membership: the 15 southern corridors (antipode Dec < -30 deg) that
PS1 and ZTF cannot reach (plan §7 TODO) — the survey exists to close
the optical cell there. Corridor keys and endpoint->corridor mapping
come from the shared table (sglsurvey/corridors.py); a coverage overlay
(surveys/decam/targets/) will grade and order them after the pilot.

Pilot corridors are the plan-named three (§7): lalande (nearest
southern-corridor star), sigmadra and hd219134 (engineering-backbone
top picks). DECam instcal coverage at each antipode was verified in
notes/decam_recon_2026-08-24.md (57 / 185 / 81 exposures).
"""

from collections import defaultdict

from sglsurvey.corridors import CORRIDOR_OF as _ALL_CORRIDOR_OF

#: Corridors with antipodes below Dec -30 (registry v1.5 astrometry).
SOUTHERN_CORRIDORS = [
    "lalande", "ross248", "cyg61", "struve2398", "grb34", "gj1221",
    "gj338", "gj625", "gj687", "gj251", "sigmadra", "hd219134",
    "wolf1069", "gj3512", "gj13157",
]

PILOT_CORRIDORS = ["lalande", "sigmadra", "hd219134"]

CORRIDOR_OF = {e: c for e, c in _ALL_CORRIDOR_OF.items()
               if c in SOUTHERN_CORRIDORS}

MEMBERS = defaultdict(list)
for _e, _c in CORRIDOR_OF.items():
    MEMBERS[_c].append(_e)
MEMBERS = dict(MEMBERS)

ALL_ENDPOINTS = list(CORRIDOR_OF)
PILOT_ENDPOINTS = [e for e, c in CORRIDOR_OF.items()
                   if c in PILOT_CORRIDORS]

#: CTIO Blanco 4-m / DECam; JPL Horizons & MPC site code W84.
CTIO = {"name": "ctio-blanco-decam", "longitude_deg": -70.80655,
        "latitude_deg": -30.16928, "height_m": 2207.0, "mpc_code": "W84"}
