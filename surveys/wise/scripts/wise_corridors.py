"""Shared corridor configuration for the WISE shakedown scripts.

Since 2026-08-22 the table lives in ``sglsurvey.corridors`` (v2 plan
§8.5); this module re-exports it so the v1 scripts run unchanged.
"""

from sglsurvey.corridors import (ALL_CORRIDORS, ALL_ENDPOINTS,  # noqa: F401
                                 CORRIDOR_OF, MEMBERS)
