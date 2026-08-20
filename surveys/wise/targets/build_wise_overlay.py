"""WISE overlay for the universal target list.

Consumes targets/universal_<version>.json (version from
targets/config.yaml) and produces the WISE-specific
work queue: per-corridor coverage and confusion grades, solution-gate
verdicts at WISE tolerance, and grandfathered/searched status.

WISE is all-sky, so no corridor is excluded; grades order the queue.
Confusion is measured, not proxied: one CatWISE2020 source count
(r = 0.2 deg) and one 2MASS brightest-Ks check (r = 0.3 deg) per
antipode, snapshotted to disk. Coverage uses |ecliptic beta| as the
proxy (exact epoch counts come free with each corridor's coarse
discovery run).

Usage: uv run python surveys/wise/targets/build_wise_overlay.py
"""

from __future__ import annotations

import csv
import io
import json
import time
from pathlib import Path

import requests
import yaml

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
VERSION = yaml.safe_load(open(REPO / "targets" / "config.yaml"))["version"]
UNIVERSAL = json.load(open(REPO / "targets" / f"universal_{VERSION}.json"))
TAP = "https://irsa.ipac.caltech.edu/TAP/sync"
SNAP = HERE / "corridor_measurements.json"


def tap_one(session, query):
    r = session.get(TAP, params={"QUERY": query, "FORMAT": "CSV"},
                    timeout=300)
    r.raise_for_status()
    rows = list(csv.DictReader(io.StringIO(r.text)))
    return rows[0] if rows else {}


def measure(session, ra, dec):
    n = tap_one(session, (
        "SELECT COUNT(*) AS n FROM catwise_2020 WHERE "
        f"CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.4f},"
        f"{dec:.4f},0.2))=1"))
    k = tap_one(session, (
        "SELECT MIN(k_m) AS kmin FROM fp_psc WHERE "
        f"CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.4f},"
        f"{dec:.4f},0.3))=1"))
    return {"catwise_n_r0p2": int(float(n.get("n") or 0)),
            "brightest_ks_r0p3": (float(k["kmin"])
                                  if k.get("kmin") else None)}


def main():
    session = requests.Session()
    meas = json.load(open(SNAP)) if SNAP.exists() else {}
    systems = UNIVERSAL["systems"]
    for i, s in enumerate(systems):
        if s["system"] in meas:
            continue
        t0 = time.monotonic()
        meas[s["system"]] = measure(session, s["anti_ra"], s["anti_dec"])
        SNAP.write_text(json.dumps(meas, indent=1))
        print(f"[{i + 1}/{len(systems)}] {s['system']}: "
              f"{meas[s['system']]} ({time.monotonic() - t0:.0f}s)",
              flush=True)

    # Grades. Confusion: density per deg^2 (area 0.1257 deg^2) and a
    # bright-star hazard when a Ks < 4 star sits in the corridor
    # (empirical: the Sirius corridor's contamination-limited depth).
    dens = {k: v["catwise_n_r0p2"] / 0.1257 for k, v in meas.items()}
    order = sorted(dens.values())
    lo, hi = order[len(order) // 3], order[2 * len(order) // 3]
    rows = []
    for s in systems:
        m = meas[s["system"]]
        d = dens[s["system"]]
        conf = "low" if d <= lo else ("mid" if d <= hi else "high")
        bright = (m["brightest_ks_r0p3"] is not None
                  and m["brightest_ks_r0p3"] < 4.0)
        cov = ("high" if s["abs_ecl_beta"] > 60
               else "mid" if s["abs_ecl_beta"] > 25 else "baseline")
        gate = s["solution_gate"]
        status = ("searched" if s["already_searched_wise"]
                  else "deferred" if gate in ("orbit_needed", "hard")
                  else "queued")
        rows.append({
            "system": s["system"], "dist_pc": s["dist_pc"],
            "n_tracks": s["n_tracks"], "baskets": s["baskets"],
            "universal_gate": gate,
            "catwise_density_deg2": round(d),
            "confusion": conf + ("+bright-star" if bright else ""),
            "coverage_proxy": cov, "status": status,
            "queue_key": (0 if status == "queued" else 1,
                          {"low": 0, "mid": 1, "high": 2}[conf]
                          + (2 if bright else 0),
                          s["dist_pc"]),
        })
    queue = sorted([r for r in rows if r["status"] == "queued"],
                   key=lambda r: r["queue_key"])
    for r in rows:
        r.pop("queue_key")
    out = {"universal_version": VERSION,
           "survey": "wise-neowise-merged-l1b",
           "exclusions": "none (all-sky survey)",
           "grandfathered": [r["system"] for r in rows
                             if r["status"] == "searched"],
           "systems": rows,
           "work_queue": [r["system"] for r in queue]}
    (HERE / f"overlay_{VERSION}.json").write_text(json.dumps(out, indent=1))
    write_md(rows, queue, HERE / f"overlay_{VERSION}.md")
    print(f"\n{len(rows)} systems: "
          f"{sum(1 for r in rows if r['status'] == 'searched')} searched, "
          f"{len(queue)} queued, "
          f"{sum(1 for r in rows if r['status'] == 'deferred')} deferred")
    print("\nWISE work queue (best first):")
    for r in queue:
        print(f"  {r['system']:24s} {r['dist_pc']:5.2f} pc "
              f"tracks={r['n_tracks']} confusion={r['confusion']:16s} "
              f"cov={r['coverage_proxy']}")
    print("\ndeferred:")
    for r in rows:
        if r["status"] == "deferred":
            print(f"  {r['system']:24s} gate={r['universal_gate']}")


def write_md(rows, queue, path):
    import datetime as dt
    L = ["---", f'title: "WISE overlay on universal {VERSION}"',
         f"date: {dt.date.today().isoformat()}",
         f'status: "{len(rows)} systems: '
         f'{sum(1 for r in rows if r["status"] == "searched")} searched, '
         f'{len(queue)} queued, '
         f'{sum(1 for r in rows if r["status"] == "deferred")} deferred"',
         "---", "", f"# WISE overlay — universal {VERSION}", "",
         f"Applies WISE-specific grades to `targets/universal_{VERSION}` "
         "without changing membership. All-sky: no exclusions; confusion "
         "(CatWISE2020 density, r = 0.2°) and bright-star hazard "
         "(2MASS Ks < 4 within 0.3°) order the queue; coverage proxy is "
         "|β_ecl|.", "",
         "## Work queue (best first)", "",
         "| # | System | d (pc) | Tracks | Confusion | Coverage | Baskets |",
         "| --- | --- | --- | --- | --- | --- | --- |"]
    for i, r in enumerate(queue, 1):
        L.append(f"| {i} | {r['system']} | {r['dist_pc']:.2f} | "
                 f"{r['n_tracks']} | {r['confusion']} | {r['coverage_proxy']} | "
                 f"{', '.join(r['baskets'])} |")
    L += ["", "## Deferred (solution gate)", ""]
    for r in rows:
        if r["status"] == "deferred":
            L.append(f"- {r['system']}: {r['universal_gate']}")
    L += ["", "## Searched (grandfathered from v1 registry)", "",
          ", ".join(r["system"] for r in rows if r["status"] == "searched")]
    path.write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
