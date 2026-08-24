"""Unified covered-window ledger (joint-crossings plan.md §1).

Aggregates the three finished crossings surveys into one table: every
(archive, channel, rung, event, band) with coverage, its statistical
disposition under that survey's frozen rules, and its depth where one
was calibrated. No statistic is recomputed; dispositions map 1:1 from
the frozen per-survey results (status vocabulary in plan.md).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

OUT = REPO / "surveys" / "joint-crossings" / "results"
ZTF = REPO / "surveys" / "ztf-crossings"
PS1 = REPO / "surveys" / "ps1-crossings"
WIS = REPO / "surveys" / "wise-crossings"

ROWS = []


def add(archive, channel, rung, tid, eid, t_ca, window_days, band,
        n_epochs, status, S=None, T=None, margin=None, m90=None,
        m90_cens=False, notes=""):
    ROWS.append({
        "archive": archive, "channel": channel, "rung_au": rung,
        "target_id": tid, "event_id": eid or "",
        "t_ca_mjd": np.nan if t_ca is None else float(t_ca),
        "window_days": np.nan if window_days is None else
        float(window_days),
        "band": band, "n_epochs": int(n_epochs or 0), "status": status,
        "S": np.nan if S is None else float(S),
        "T": np.nan if T is None else float(T),
        "margin": np.nan if margin is None else float(margin),
        "m90": np.nan if m90 is None else float(m90),
        "m90_censored": bool(m90_cens), "notes": notes})


def win_lookup(cov, channel_col=True):
    out = {}
    for r in cov:
        key = (str(r["event_id"]), float(r["radius_au"]))
        out[key] = (float(r["t_ca_mjd"]), float(r["window_days"]))
    return out


def ztf():
    cov = Table.read(ZTF / "results" / "coverage_v1_events.ecsv")
    wins = win_lookup(cov)
    res = json.loads((ZTF / "results" / "confirmatory_v1.json").read_text())
    comp = json.loads((ZTF / "results" / "completeness_v1.json").read_text())
    m90b = {(r["target_id"], r["band"], r.get("event_id", ""),
             float(r["radius_au"])): r for r in comp.get("B", [])}
    m90a = {(r["target_id"], r["band"], float(r["radius_au"])): r
            for r in comp.get("A", [])}
    for r in res.get("B", []):
        eid, rad = str(r.get("event_id", "")), float(r["radius_au"])
        t_ca, wd = wins.get((eid, rad), (r.get("t_ca_mjd"), None))
        n = int(r.get("n_epochs") or 0)
        if n == 0:
            st = "no_usable_data"
        elif r.get("exceedance"):
            st = "exceedance_adjudicated"
        elif r.get("kind") == "single_epoch" or r.get("status") == \
                "single_epoch":
            st = "single_epoch"
        else:
            st = "searched_null"
        c = m90b.get((r["target_id"], r["band"], eid, rad), {})
        add("ztf", "B", rad, r["target_id"], eid, t_ca, wd, r["band"],
            n, st, r.get("S"), r.get("T"),
            None if r.get("S") is None or r.get("T") is None
            else r["S"] - max(r["T"], 0.0),
            c.get("m90"), bool(c.get("grid_censored", False)))
    for r in res.get("A", []):
        rad = float(r.get("radius_au", 0.1))
        st = ("exceedance_adjudicated" if r.get("exceedance")
              else "searched_null" if r.get("status") == "searchable"
              else "constraint_only")
        c = m90a.get((r["target_id"], r["band"], rad), {})
        add("ztf", "A", rad, r["target_id"], "", None, None, r["band"],
            r.get("n_epochs") or r.get("n_windows_with_data") or 0, st,
            r.get("S"), r.get("T"), None, c.get("m90"),
            bool(c.get("grid_censored", False)),
            notes=str(r.get("status", "")))


def ps1():
    cov = Table.read(PS1 / "results" / "coverage_v1_events.ecsv")
    conf = json.loads((PS1 / "results" / "confirmatory_v1.json").read_text())
    dev = json.loads((PS1 / "results" / "dev_search_v1.json").read_text())
    comp = json.loads((PS1 / "results" / "completeness_v1.json").read_text())
    m90b = {(r["target_id"], r["band"], r["event_id"],
             float(r.get("radius_au", 0.1))): r for r in comp["B"]}
    m90a = {(r["target_id"], r["band"]): r for r in comp["A"]}
    B = cov[cov["channel"] == "B"]
    wins = {}
    for r in B:
        wins[(str(r["event_id"]), float(r["radius_au"]))] = (
            float(r["t_ca_mjd"]), float(r["window_days"]))
    for src, rows in (("conf", conf["B"]), ("dev", dev["B"])):
        for r in rows:
            eid, rad = r["event_id"], float(r.get("radius_au", 0.1))
            t_ca, wd = wins.get((eid, rad), (r.get("t_ca_mjd"), None))
            if r["status"] == "track_masked":
                st = "track_masked"
            elif r.get("exceedance"):
                st = ("retained_ambiguous"
                      if r["target_id"] == "gj-1276"
                      else "exceedance_vetoed")
            elif r.get("class") == "single_epoch":
                st = "single_epoch"
            else:
                st = "searched_null"
            c = m90b.get((r["target_id"], r["band"], eid, rad), {})
            add("ps1", "B", rad, r["target_id"], eid, t_ca, wd,
                r["band"], r.get("n_in_exposures_dedup",
                                 r.get("n_in_epochs_raw", 0)), st,
                r.get("S"), r.get("T"), r.get("margin"),
                c.get("m90"), bool(c.get("grid_censored", False)),
                notes=("dev" if src == "dev" else "") +
                      (";" + r["status"] if "partial" in r["status"]
                       else ""))
    for r in conf["A"]:
        c = m90a.get((r["target_id"], r["band"]), {})
        add("ps1", "A", 0.1, r["target_id"], "", None, None, r["band"],
            r["gates"].get("n_windows_with_data", 0), "constraint_only",
            r.get("S"), r.get("T"), None, c.get("m90"),
            bool(c.get("grid_censored", False)),
            notes=c.get("depth_kind", ""))
    A10 = cov[(cov["channel"] == "A") & (np.asarray(cov["radius_au"])
                                         == 1.0)]
    A10 = A10[np.asarray(A10["n_primary"]) > 0]
    for r in A10:
        add("ps1", "A", 1.0, str(r["target_id"]), str(r["event_id"]),
            float(r["t_ca_mjd"]), float(r["window_days"]), "any",
            int(r["n_primary"]), "coverage_only",
            notes="rung frozen constraint-only, zero trials")


def wise():
    cov = Table.read(WIS / "results" / "coverage_v1_events.ecsv")
    frz = json.loads((WIS / "configs"
                      / "threshold_freeze_v1.json").read_text())
    cls = {}
    for u in frz["search_units"]["single_window"]:
        cls[(u["target_id"], u["band"])] = "single_window"
    for u in frz["search_units"]["gate_blocked"]:
        cls[(u["target_id"], u["band"])] = "gate_blocked"
    covd = cov[np.asarray(cov["n_primary"]) > 0]
    for r in covd:
        tid = str(r["target_id"])
        for b in ("W1", "W2", "W3", "W4"):
            if int(r[f"n_{b}"]) == 0:
                continue
            note = cls.get((tid, b), "band_excluded_or_unclassed")
            add("wise", "A", 1.0, tid, str(r["event_id"]),
                float(r["t_ca_mjd"]), float(r["window_days"]), b,
                int(r[f"n_{b}"]), "coverage_only",
                notes=f"survey structural null; unit class: {note}")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ztf()
    ps1()
    wise()
    tab = Table(rows=ROWS)
    tab.write(OUT / "covered_window_ledger_v1.ecsv",
              format="ascii.ecsv", overwrite=True)
    from collections import Counter
    by_status = Counter(r["status"] for r in ROWS)
    by_archive = Counter(r["archive"] for r in ROWS)
    multi = {}
    for r in ROWS:
        multi.setdefault((r["target_id"], r["channel"]),
                         set()).add(r["archive"])
    shared = sorted(f"{t}/{ch}: {sorted(a)}"
                    for (t, ch), a in multi.items() if len(a) >= 2)
    summary = {
        "n_rows": len(ROWS), "by_archive": dict(by_archive),
        "by_status": dict(by_status),
        "n_target_channel_multi_archive": len(shared),
        "multi_archive": shared,
        "geometric_null_note": "WISE channels B and A-0.1 are "
            "elongation-null (0 observable events) and carry no rows",
    }
    (OUT / "ledger_v1_summary.json").write_text(
        json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: summary[k] for k in
                      ("n_rows", "by_archive", "by_status",
                       "n_target_channel_multi_archive")}, indent=1))


if __name__ == "__main__":
    main()
