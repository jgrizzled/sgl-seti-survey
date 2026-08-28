"""DASCH crossings blind confirmatory search.

Runs the 15 frozen confirmatory units (threshold_freeze_v1.json)
under the frozen chain v1.0 + v1.1 + v1.2 + v1.3. No thresholds or
constructions may change after this script touches in-window data.

Products: results/confirmatory_v1.json; snapshots under
runs/dasch/v1/confirmatory/.
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).parent))
import daschlib as dl        # noqa: E402
import search_lib as sl      # noqa: E402

D = Path(__file__).resolve().parents[1]
SNAP = dl.RUNS / "confirmatory"
CONF_UNITS = sl.FREEZE["units"]["confirmatory"]

REG = yaml.safe_load(open(dl.REPO / "registries"
                          / "pilot_wise_2026.yaml"))
TARGETS = REG.get("targets", REG)


def brief(res: dict) -> dict:
    return {k: res[k] for k in
            ("S_event", "T_event", "exceed_event", "S_stack",
             "T_stack", "exceed_stack") if k in res}


def main() -> int:
    out = {"generated_utc": datetime.now(timezone.utc).isoformat(
               timespec="seconds"),
           "chain": "hypotheses v1.0 + amendments v1.1/v1.2/v1.3; "
                    "threshold_freeze_v1",
           "freeze_sha256": dl.sha256_file(
               D / "configs" / "threshold_freeze_v1.json"),
           "hypotheses_sha256": dl.sha256_file(D / "hypotheses.md"),
           "units": {}}
    exceed = 0
    for unit in CONF_UNITS:
        target, rung = unit.split("/")
        chan = rung[0]
        print("==", unit, flush=True)
        if chan == "B":
            res = sl.run_hit_unit(unit, SNAP / "platephot", track=True)
        elif sl.A_ROUTING[target] is None:
            astro = TARGETS[target]["state"]["astrometry"]
            quies = sl.measure_quiescent(target, unit,
                                         SNAP / "offwindow_rate", astro)
            print("quiescent:", json.dumps(quies), flush=True)
            res = sl.run_hit_unit(unit, SNAP / "platephot", track=False,
                                  quiescent_mag=quies["quiescent_mag"])
            res["regime"] = "limits_only"
            res["quiescent"] = quies
        else:
            refcat, gsc, ref = sl.A_ROUTING[target]
            snap = SNAP / "lightcurve" / f"{target}_{refcat}_{ref}.json"
            resp = dl.post("dasch/dr7/lightcurve",
                           {"gsc_bin_index": gsc, "ref_number": ref,
                            "refcat": refcat}, snap)
            rows = [r for r in dl.rows_of(resp)
                    if (r.get("magcal_magdep") not in ("", "99.0")
                        and r.get("date_jd")
                        and r.get("reject_flag") in ("", "0")
                        and not (int(float(r["aflags"] or 0))
                                 & sl.TOO_BRIGHT_A)
                        and not (int(float(r["bflags"] or 0))
                                 & sl.SATURATED_B))]
            res = sl.run_lightcurve_unit(unit, rows)
            res["routing"] = {"refcat": refcat, "ref_number": ref}
            res["usable_rows"] = len(rows)
        out["units"][unit] = res
        exceed += int(bool(res.get("exceed_event"))) \
            + int(bool(res.get("exceed_stack")))
        print(json.dumps(brief(res)), flush=True)
    out["n_exceedances"] = exceed
    out["trials"] = 2 * len(CONF_UNITS)
    out["expected_control_crossings"] = round(
        2 * len(CONF_UNITS) / 9.0, 2)
    (D / "results" / "confirmatory_v1.json").write_text(
        json.dumps(out, indent=1) + "\n")
    print(f"DONE: {exceed} exceedances / {out['trials']} trials "
          f"(expected {out['expected_control_crossings']})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
