"""Shared DASCH/Starglass access + frozen-hypothesis constants
(hypotheses.md v1.0, FROZEN 2026-08-26).

All queries go through `post()` (curl-like UA — the WAF rejects
python-urllib's default; recon), are snapshotted verbatim to the
run directory, and sha256-recorded by the caller's manifest.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import time
import urllib.request
from pathlib import Path

BASE = "https://api.starglass.cfa.harvard.edu/public"
UA = "sgl-seti-survey/dasch-crossings (curl-compatible)"
REPO = Path(__file__).resolve().parents[3]
RUNS = REPO / "runs" / "dasch" / "v1"

# frozen scope (hypotheses §0) and rung ladder (§2)
TARGETS = ("van-maanen", "wolf-359", "teegarden", "gj-1276", "ross-128")
RSUN_AU = 6.957e5 / 1.496e8
RUNGS = {  # name -> radius in AU; channel applicability per §2
    "B_1.2Rs": 1.2 * RSUN_AU,
    "B_2.5Rs": 2.5 * RSUN_AU,
    "B_0.1AU": 0.1,
    "A_0.1AU": 0.1,
}
CHANNEL_LINK = {"A": "inbound", "B": "outbound"}
# §7 gates
TIMING_GATE_FRACTION = 0.5      # sigma_t <= 0.5 * half-width
DATE_ONLY_SIGMA_DAYS = 1.0      # 00:00:00-timestamp heuristic
EDGE_FLAG_DEG = 0.1             # near-edge exposures flagged for exact re-check
CLUSTER_TOL_DEG = 2.0 / 60.0    # locus clustering tolerance (coverage queries)


def sha256_file(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def post(endpoint: str, payload: dict, snapshot: Path,
         pause_s: float = 1.0, timeout: int = 300) -> list | dict:
    """POST to the public API; write the raw response to `snapshot`
    (reused if it exists — snapshots are immutable per run)."""
    if snapshot.exists():
        return json.loads(snapshot.read_bytes())
    req = urllib.request.Request(
        f"{BASE}/{endpoint.lstrip('/')}",
        data=json.dumps(payload).encode(),
        headers={"content-type": "application/json",
                 "accept": "application/json", "user-agent": UA})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        body = r.read()
    snapshot.parent.mkdir(parents=True, exist_ok=True)
    snapshot.write_bytes(body)
    time.sleep(pause_s)
    return json.loads(body)


def rows_of(response: list) -> list[dict]:
    """DASCH CSV-in-JSON -> list of dicts."""
    return list(csv.DictReader(io.StringIO("\n".join(response))))


def exposure_usable(r: dict) -> bool:
    """Frozen listing-level coverage gate (hypotheses §7): APASS
    photometric calibration present and non-logbook WCS."""
    return (r.get("limMagApass") not in (None, "", "0.000")
            and r.get("wcssource") in ("imwcs", "catalog"))


def exposure_interval_jd(r: dict):
    """[start, stop] JD-UTC from expdate + exptime (minutes; §7
    interval-based overlap). Returns None if no date."""
    from astropy.time import Time
    if not r.get("expdate"):
        return None
    t0 = Time(r["expdate"].replace("Z", ""), format="isot", scale="utc").jd
    dur_d = max(float(r.get("exptime") or 0.0), 0.1) / 60.0 / 24.0
    return t0, t0 + dur_d


def date_only(r: dict) -> bool:
    """00:00:00-timestamp heuristic => sigma_t = 1.0 d (§7)."""
    return r.get("expdate", "").endswith("T00:00:00.0Z")
