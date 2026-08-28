"""Rubin DP2 crossings blind confirmatory (threshold freeze v1.0+v1.1).

Runs the identical frozen chain (rubinlib — the dev-stage code path,
byte-for-byte) on the single confirmatory unit:
ross-128 B 0.1 AU r, visit 2025091300612. First data contact with the
unit visit in the programme. Snapshots under
runs/rubin-crossings/confirmatory_v1/.

Output: results/confirmatory_v1.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
import rubinlib as rl  # noqa: E402

from sglsurvey.adapters.rsp_tap import RspTapClient  # noqa: E402
from sglsurvey.snapshots import SnapshotStore  # noqa: E402

REPO = rl.REPO
SURVEY = REPO / "surveys" / "rubin-crossings"
OUT = SURVEY / "results" / "confirmatory_v1.json"
RUN = REPO / "runs" / "rubin-crossings" / "confirmatory_v1"


def main():
    freeze = json.loads(
        (SURVEY / "configs" / "threshold_freeze_v1.json").read_text())
    unit = freeze["search_family"]["units"][0]
    t = Table.read(REPO / "crossings" / "universal_v1" / "events.ecsv")
    ev = t[np.array(t["event_id"]) == unit["event_id"]][0]
    client = RspTapClient()
    store = SnapshotStore(RUN)

    # fresh visit-row pull (epoch + band from the archive, not the freeze)
    _, rows = client.query(
        "SELECT visit, band, expMidptMJD, expTime FROM dp2.Visit "
        f"WHERE visit = {unit['visit']}", store)
    assert len(rows) == 1 and str(rows[0][1]) == unit["band"]
    t_mjd = float(rows[0][2])

    rec = rl.run_chain(ev, unit["visit"], t_mjd, unit["band"],
                       client, store)
    rec["unit_id"] = unit["unit_id"]
    rec["band"] = unit["band"]
    rec["dt_days"] = t_mjd - unit["t_ca_mjd"]

    verdict = {
        "unit": unit["unit_id"],
        "status": rec["status"],
        "s_det": rec.get("s_det"),
        "threshold": rec.get("threshold"),
        "exceedance": rec.get("exceedance"),
        "n_trials": freeze["search_family"]["n_trials"],
        "expected_control_crossings":
            freeze["search_family"]["expected_control_crossings"],
    }
    out = {"freeze_sha_bound": "threshold_freeze_v1.json at run time",
           "unit_record": rec, "verdict": verdict}
    OUT.write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps(verdict, indent=1, default=float))
    if rec.get("association"):
        print("ASSOCIATION -> veto ladder required:",
              json.dumps(rec["association"], default=str))
    print("wrote", OUT.relative_to(REPO))


if __name__ == "__main__":
    main()
