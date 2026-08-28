"""Rubin DP2 crossings threshold freeze v1.0.

Binds the frozen numerics to content hashes of hypotheses.md and the
coverage products, applies the section-6 saturation cut from the
frozen universal-list input (no new archive query), and fixes the
trial family. No DiaSource row is touched.

Output: configs/threshold_freeze_v1.json
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

SURVEY = Path(__file__).resolve().parents[1]
REPO = SURVEY.parents[1]


def sha(p: Path) -> str:
    return "sha256:" + hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> None:
    targets = json.loads((REPO / "targets" / "universal_v2.json")
                         .read_text())
    ross154 = next(s for s in targets["systems"]
                   if s["system"] == "Ross 154")
    gmag = float(ross154["components"][0]["gmag"])
    assert gmag < 16.5  # section-6 rule: excluded band

    payload = {
        "frozen_utc": "2026-08-26",
        "inputs": {
            "hypotheses_md": sha(SURVEY / "hypotheses.md"),
            "coverage_events": sha(SURVEY / "results"
                                   / "coverage_v1_events.ecsv"),
            "coverage_summary": sha(SURVEY / "results"
                                    / "coverage_v1_summary.json"),
            "coverage_md": sha(SURVEY / "results" / "coverage_v1.md"),
            "universal_targets": sha(REPO / "targets"
                                     / "universal_v2.json"),
            "universal_events": sha(REPO / "crossings" / "universal_v1"
                                    / "events.ecsv"),
        },
        "saturation_cut": {
            "rule": "hypotheses.md section 6 (G 16.5/17.5 bands)",
            "ross-154": {"gmag": gmag, "class": "excluded",
                         "source": "targets/universal_v2.json (frozen "
                                   "input; no new archive query)"},
            "consequence": "A 0.1 AU closes coverage-without-statistic "
                           "(saturation-limited, 16 on-detector visits)",
        },
        "search_family": {
            "units": [{
                "unit_id": "ross-128_B_0.1AU_r",
                "target_id": "ross-128",
                "event_id": "evt-5197f21cd67b",
                "channel": "B", "radius_au": 0.1, "band": "r",
                "t_ca_mjd": 60937.9414,
                "b_min_rsun": 1.855,
                "visit": 2025091300612,
                "visit_mjd": 60932.2843,
                "detector": 102,
                "maglim": 23.3661,
                "b_at_visit_rsun": 20.921,
                "statistics": ["S_det"],
            }],
            "n_trials": 1,
            "expected_control_crossings": 1.0 / 9.0,
            "fwer_alpha": 0.05,
        },
        "statistic": {
            "S_det": "max over gate-passing associated DiaSources of "
                     "psfFlux/psfFluxErr; 0 if none. Association: "
                     "within r_assoc of any deduplicated z-grid "
                     "position (positions closer than r_assoc merge).",
            "r_assoc_arcsec": 1.0,
            "z_grid_au": [550.0, 1000.0, 2500.0, 5500.0, 10000.0],
            "gates": "hypotheses.md section 7 (footprint; pixelFlags "
                     "exclusion set incl. archive injections; "
                     "psfFlux_flag; centroid_flag)",
            "static_exclusion_arcsec": 2.0,
            "static_exclusion_mag": "any-band Object psfMag < "
                                    "maglim(band) + 0.5",
        },
        "controls": {
            "n_per_unit": 8,
            "rule": "z-family offset pattern about the axis point, "
                    "rotated k*40 deg (k=1..8) about the axis point at "
                    "the visit epoch; a control failing the static "
                    "exclusion, the footprint gate, or the >= 2.5 "
                    "arcsec locus-avoidance (amendment v1.1) rotates "
                    "+5 deg until valid; threshold T = max control "
                    "S_det; exceedance S > max(T, 0)",
            "locus_avoid_arcsec": 2.5,
            "amendment": "v1.1 (hypotheses.md; was 10 arcsec in v1.0 — "
                         "geometrically impossible for inner-z "
                         "positions, found at dev)",
        },
        "dev_gates": [
            "association-rate null on off-window pseudo-units "
            "(all off-window visits at the ross-128 antipode field)",
            "control-rule census (rotation resolutions, static "
            "exclusions)",
            "positive control: SkyBoT-predicted known-object position "
            "on an off-window visit must associate a gate-passing "
            "DiaSource through the identical chain",
            "archive-injection census: pixelFlags_injected DiaSource "
            "population near the unit field (completeness route "
            "section 8; decides calibrated vs threshold-statement)",
        ],
    }
    out = SURVEY / "configs" / "threshold_freeze_v1.json"
    out.write_text(json.dumps(payload, indent=2) + "\n")
    print("wrote", out.relative_to(REPO))
    print("freeze sha:", sha(out))


if __name__ == "__main__":
    main()
