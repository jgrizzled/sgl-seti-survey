"""PTF corridor overlay v1: grade every ZTF-visible corridor by its
actual PTF level-1 coverage, from the coarse-discovery and precise-pass
records (no new archive queries). A corridor is *searchable* when at
least one (endpoint, role, band) cell holds >= MIN_EPOCHS usable or
partial precise-pass exposures (the engine's epoch floor); the grade
records the depth of that coverage and its calendar-month spread (the
parallax-phase proxy). Writes surveys/ptf/targets/overlay_v1.{json,md}.

Usage: uv run python surveys/ptf/targets/build_ptf_overlay.py
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "surveys" / "ptf" / "scripts"))
from ptf_corridors import ALL_CORRIDORS, CORRIDOR_OF, MEMBERS  # noqa: E402

from sglsurvey.records import read_records  # noqa: E402

RUNS = REPO / "runs" / "ptf"
MIN_EPOCHS = 5
OUT_JSON = REPO / "surveys" / "ptf" / "targets" / "overlay_v1.json"
OUT_MD = REPO / "surveys" / "ptf" / "targets" / "overlay_v1.md"


def main() -> None:
    obs = {r["observation_id"]: r for r in read_records(
        RUNS / "coarse_v1" / "records" / "observation.jsonl")}
    coarse = defaultdict(set)
    for r in read_records(RUNS / "coarse_v1" / "records" / "intersection_evaluation.jsonl"):
        coarse[CORRIDOR_OF[r["endpoint_id"]]].add(r["observation_id"])
    precise = defaultdict(lambda: defaultdict(set))   # corridor -> (endpoint, role, band) -> oids
    pp = RUNS / "precise_v1" / "records" / "intersection_evaluation.jsonl"
    if pp.exists():
        for r in read_records(pp):
            if r["stage"] == "precise" and r["usable"] in ("usable", "partial"):
                c = CORRIDOR_OF[r["endpoint_id"]]
                precise[c][(r["endpoint_id"], r["role"], obs[r["observation_id"]]["band"])].add(
                    r["observation_id"])
    cut_index = {}
    ci = RUNS / "precise_v1" / "records" / "cutout_index.jsonl"
    if ci.exists():
        for line in open(ci):
            rec = json.loads(line)
            cut_index[rec["observation_id"]] = rec
    rows = []
    for c in ALL_CORRIDORS:
        oids = coarse.get(c, set())
        cells = precise.get(c, {})
        by_band = defaultdict(set)
        for (e, role, band), s in cells.items():
            by_band[band] |= s
        usable_all = set().union(*by_band.values()) if by_band else set()
        mjds = np.array(sorted(obs[o]["t_mid_mjd_utc"] for o in usable_all))
        months = sorted({Time(m, format="mjd").datetime.strftime("%Y-%m") for m in mjds}) if len(mjds) else []
        doy = sorted({int(Time(m, format="mjd").datetime.strftime("%j")) // 30 for m in mjds}) if len(mjds) else []
        best = max((len(s) for s in cells.values()), default=0)
        n_cells_ok = sum(1 for s in cells.values() if len(s) >= MIN_EPOCHS)
        centre = None
        for o in usable_all:
            if o in cut_index:
                centre = (cut_index[o]["center_ra_deg"], cut_index[o]["center_dec_deg"])
                break
        if centre is None:
            for o in oids:
                w = obs[o]["wcs"]
                if w.get("crval1") is not None:
                    pass
        if best >= 20 and len(doy) >= 4:
            grade = "ok"
        elif best >= MIN_EPOCHS:
            grade = "thin" if len(doy) >= 2 else "single-phase"
        elif oids:
            grade = "below-floor"
        else:
            grade = "empty"
        rows.append({
            "corridor": c, "endpoints": MEMBERS[c],
            "antipode_ra_deg": centre[0] if centre else None,
            "antipode_dec_deg": centre[1] if centre else None,
            "n_exposures_discovered": len(oids),
            "n_usable_exposures": len(usable_all),
            "usable_by_band": {b: len(s) for b, s in sorted(by_band.items())},
            "best_cell_epochs": best, "n_cells_at_floor": n_cells_ok,
            "n_calendar_months": len(months), "n_doy_bins_30d": len(doy),
            "epoch_span_mjd": [float(mjds.min()), float(mjds.max())] if len(mjds) else None,
            "grade": grade, "searchable": best >= MIN_EPOCHS,
        })
    rows.sort(key=lambda r: (-r["best_cell_epochs"], r["corridor"]))
    out = {"built_utc": datetime.now(timezone.utc).isoformat(), "min_epochs": MIN_EPOCHS,
           "source": "runs/ptf/coarse_v1 + runs/ptf/precise_v1", "rows": rows}
    OUT_JSON.write_text(json.dumps(out, indent=1) + "\n")
    lines = ["# PTF corridor overlay v1", "",
             f"Built {out['built_utc'][:19]}Z from the coarse discovery and precise pass; "
             f"searchable = best (endpoint, role, band) cell >= {MIN_EPOCHS} usable exposures.", "",
             "| corridor | endpoints | discovered | usable | g | R | best cell | cells >= floor | months | 30-d bins | grade |",
             "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {r['corridor']} | {', '.join(r['endpoints'])} | {r['n_exposures_discovered']} | "
                     f"{r['n_usable_exposures']} | {r['usable_by_band'].get('g', 0)} | {r['usable_by_band'].get('R', 0)} | "
                     f"{r['best_cell_epochs']} | {r['n_cells_at_floor']} | {r['n_calendar_months']} | "
                     f"{r['n_doy_bins_30d']} | {r['grade']} |")
    n_s = sum(r["searchable"] for r in rows)
    lines += ["", f"Searchable corridors: {n_s} / {len(rows)}; grades: "
              + ", ".join(f"{g} {sum(1 for r in rows if r['grade'] == g)}"
                          for g in ("ok", "thin", "single-phase", "below-floor", "empty"))]
    OUT_MD.write_text("\n".join(lines) + "\n")
    print("\n".join(lines[-1:]))


if __name__ == "__main__":
    main()
