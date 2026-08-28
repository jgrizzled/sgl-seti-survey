"""DASCH crossings dev stage (thresholds.md v1.0).

Runs the three frozen dev units through the frozen machinery:
  - wolf-359/B_2.5Rs      B-chain hit statistic (track)
  - teegarden/A_0.1AU     limits-only regime decision + hit statistic
  - van-maanen/A_0.1AU    forced_dev; detected-regime lightcurve z

Products: results/dev_v1.json; snapshots under runs/dasch/v1/dev/.
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import daschlib as dl        # noqa: E402
import search_lib as sl      # noqa: E402

D = Path(__file__).resolve().parents[1]
SNAP = dl.RUNS / "dev"
STAR = {  # Gaia DR3 states (registry) propagated per pos_epoch at match
    "teegarden": (43.26964247679, 16.86437381898, 3429.083, -3805.541),
    "van-maanen": (12.29674025046, 5.37655660843, 1231.399, -2711.883),
}


def querycat_star(target: str) -> list[dict]:
    ra, dec, pmra, pmdec = STAR[target]
    snap = SNAP / "querycat" / f"{target}.json"
    resp = dl.post("dasch/dr7/querycat",
                   {"ra_deg": ra, "dec_deg": dec, "radius_arcsec": 60,
                    "refcat": "apass"}, snap)
    rows = dl.rows_of(resp)
    # merge rule: entries within R_MATCH of the PM-propagated star
    # position at each entry's pos_epoch
    out = []
    for r in rows:
        try:
            ep = float(r["pos_epoch"] or 2000.0)
        except ValueError:
            ep = 2000.0
        dt = ep - 2016.0
        pra = ra + pmra / 3.6e6 * dt / np.cos(np.radians(dec))
        pde = dec + pmdec / 3.6e6 * dt
        sep = 3600.0 * np.hypot(
            (float(r["ra_deg"]) - pra) * np.cos(np.radians(dec)),
            float(r["dec_deg"]) - pde)
        if sep <= sl.R_MATCH:
            r["_sep_arcsec"] = round(float(sep), 2)
            out.append(r)
    return out


def lightcurve_rows(target: str, entries: list[dict]) -> list[dict]:
    rows = []
    for e in entries:
        snap = SNAP / "lightcurve" / f"{target}_{e['ref_number']}.json"
        resp = dl.post("dasch/dr7/lightcurve",
                       {"gsc_bin_index": int(e["gsc_bin_index"]),
                        "ref_number": int(e["ref_number"]),
                        "refcat": "apass"}, snap)
        for r in dl.rows_of(resp):
            if (r.get("magcal_magdep") not in ("", "99.0")
                    and r.get("date_jd")
                    and r.get("reject_flag") in ("", "0")
                    and not (int(float(r["aflags"] or 0))
                             & sl.TOO_BRIGHT_A)
                    and not (int(float(r["bflags"] or 0))
                             & sl.SATURATED_B)):
                rows.append(r)
    return rows


def main() -> int:
    out = {"generated_utc": datetime.now(timezone.utc).isoformat(
        timespec="seconds"),
        "freeze_sha256": dl.sha256_file(
            D / "configs" / "threshold_freeze_v1.json"),
        "units": {}}

    # --- wolf-359 B 2.5 Rs: B-chain hits ---
    print("== wolf-359/B_2.5Rs", flush=True)
    res = sl.run_hit_unit("wolf-359/B_2.5Rs", SNAP / "platephot", True)
    out["units"]["wolf-359/B_2.5Rs"] = res
    print(json.dumps({k: res[k] for k in
                      ("S_event", "T_event", "exceed_event", "S_stack",
                       "T_stack", "exceed_stack",
                       "platephot_row_counts")}), flush=True)

    # --- teegarden A 0.1 AU: regime decision + limits-only hits ---
    print("== teegarden/A_0.1AU", flush=True)
    ent = querycat_star("teegarden")
    regime = "limits_only" if not ent else "decide_from_stdmag"
    if ent:
        stdmag = min(float(e["stdmag"]) for e in ent if e.get("stdmag"))
        # p90 in-window limiting mag measured from the search pulls
        regime = {"stdmag": stdmag}
    import yaml
    reg = yaml.safe_load(open(dl.REPO / "registries"
                              / "pilot_wise_2026.yaml"))
    astro = reg.get("targets", reg)["teegarden"]["state"]["astrometry"]
    quies = sl.measure_quiescent("teegarden", "teegarden/A_0.1AU",
                                 SNAP / "offwindow_rate", astro)
    print("quiescent:", json.dumps(quies), flush=True)
    res = sl.run_hit_unit("teegarden/A_0.1AU", SNAP / "platephot", False,
                          quiescent_mag=quies["quiescent_mag"])
    res["regime"] = regime
    res["quiescent"] = quies
    res["querycat_entries"] = len(ent)
    out["units"]["teegarden/A_0.1AU"] = res
    print(json.dumps({k: res[k] for k in
                      ("S_event", "T_event", "exceed_event", "S_stack",
                       "T_stack", "exceed_stack", "regime",
                       "querycat_entries")}), flush=True)

    # --- van-maanen A 0.1 AU (forced_dev): lightcurve z ---
    # Amendment v1.1: routed to the ATLAS refcat entry (the APASS
    # entry is a PM-less dummy; dev finding, hypotheses §12).
    print("== van-maanen/A_0.1AU (forced_dev, v1.1 routing)", flush=True)
    refcat, gsc, ref = sl.A_ROUTING["van-maanen"]
    snap = SNAP / "lightcurve" / f"van-maanen_{refcat}_{ref}.json"
    resp = dl.post("dasch/dr7/lightcurve",
                   {"gsc_bin_index": gsc, "ref_number": ref,
                    "refcat": refcat}, snap)
    rows = [r for r in dl.rows_of(resp)
            if (r.get("magcal_magdep") not in ("", "99.0")
                and r.get("date_jd")
                and r.get("reject_flag") in ("", "0")
                and not (int(float(r["aflags"] or 0)) & sl.TOO_BRIGHT_A)
                and not (int(float(r["bflags"] or 0))
                         & sl.SATURATED_B))]
    res = sl.run_lightcurve_unit("van-maanen/A_0.1AU", rows)
    res["routing"] = {"refcat": refcat, "ref_number": ref,
                      "amendment": "v1.1"}
    res["usable_rows"] = len(rows)
    out["units"]["van-maanen/A_0.1AU"] = res
    print(json.dumps({k: res[k] for k in
                      ("S_event", "T_event", "exceed_event", "S_stack",
                       "T_stack", "exceed_stack", "n_inwindow_rows",
                       "usable_rows")}), flush=True)

    (D / "results" / "dev_v1.json").write_text(
        json.dumps(out, indent=1) + "\n")
    print("wrote results/dev_v1.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
