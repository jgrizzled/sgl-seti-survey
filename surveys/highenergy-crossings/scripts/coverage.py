"""Coverage stage v1 (hypotheses.md §5 coverage gate, D2-D5): per unit and
window, the gated LAT livetime / per-lane exposure / Moon-excluded seconds
from the weekly spacecraft files, the searchable flag (lane-L livetime
>= 1 ks), and the kept pseudo-window count. No photon is touched.
Writes results/coverage_v1.json and results/coverage_v1.md."""
from __future__ import annotations

import json
from collections import Counter

import numpy as np
from astropy.time import Time

import hlib as H
import recon_scan as R


def main():
    events = R.load_events()
    units = H.all_units(events)
    out = {"utc": Time.now().isot, "freeze": H.clean(H.FREEZE), "units": []}
    for u in units:
        ws = H.unit_windows(u["channel"], u["target"], u["rung"], events)
        rows = []
        for w in ws:
            iv = H.intervals(u["ra"], u["dec"], w["mjd0"], w["mjd1"])
            live = float(iv["live"].sum()) if iv is not None else 0.0
            rec = {"event_id": w["event_id"], "t_ca_utc": w["t_ca_utc"], "b_rsun": w["b_rsun"], "dur_d": w["dur_d"],
                   "live_s": live, "live_raw_s": float(iv["live_raw"].sum()) if iv is not None else 0.0,
                   "moon_excluded_s": iv["moon_excluded_s"] if iv is not None else 0.0,
                   "min_moon_deg": float(iv["moon"].min()) if iv is not None else None,
                   "searchable": live >= H.FREEZE["window_livetime_min_s"],
                   "forced_dev": H.is_forced_dev(u["channel"], u["target"], w["t_ca_utc"])}
            for lane in H.FREEZE["lanes"]:
                rec["expo_" + lane] = float(iv["expo_" + lane].sum()) if iv is not None else 0.0
            n_pw = 0; n_pw_kept = 0
            if rec["searchable"]:
                for j in H.pseudo_offsets(u["rung"], w["dur_d"]):
                    off = j * max(w["dur_d"], 1.0)
                    iv2 = H.intervals(u["ra"], u["dec"], w["mjd0"] + off, w["mjd1"] + off)
                    n_pw += 1
                    if iv2 is not None and iv2["expo_L"].sum() >= H.FREEZE["pw_min_expo_frac"] * rec["expo_L"]:
                        n_pw_kept += 1
            rec["n_pw"] = n_pw; rec["n_pw_kept"] = n_pw_kept
            rows.append(rec)
        s = [r for r in rows if r["searchable"] and not r["forced_dev"]]
        u2 = dict(u, windows=rows, n_searchable=len(s),
                  live_searchable_s=float(sum(r["live_s"] for r in s)),
                  expo_L_searchable=float(sum(r["expo_L"] for r in s)),
                  expo_H_searchable=float(sum(r["expo_H"] for r in s)),
                  searchable_unit=len(s) > 0)
        out["units"].append(u2)
        print(f"[cov] {u['channel']} {u['target']:11s} {u['rung']:8s} windows {len(rows):3d} searchable {len(s):3d} "
              f"live {u2['live_searchable_s']/1e3:8.1f} ks  expoL {u2['expo_L_searchable']/1e4:8.1f} m2ks", flush=True)
    H.RES.mkdir(exist_ok=True)
    (H.RES / "coverage_v1.json").write_text(json.dumps(H.clean(out), indent=1))
    # ledger markdown
    md = ["# Coverage v1 (LAT arm) — 2026-09-07", "",
          "Gated livetime (θ ≤ 60°, DATA_QUAL > 0, LAT_CONFIG = 1, Moon/Sun > 8°) per window; searchable = lane-L livetime ≥ 1 ks.", "",
          "| Unit | Windows | Searchable | Σ live (ks) | Σ expo L (m² ks) | Σ expo H (m² ks) | Moon-excluded (ks) | dev |", "|---|---|---|---|---|---|---|---|"]
    for u in out["units"]:
        md.append(f"| {u['channel']} {u['target']} {u['rung']} | {u['n_windows']} | {u['n_searchable']} | {u['live_searchable_s']/1e3:.1f} | "
                  f"{u['expo_L_searchable']/1e4:.1f} | {u['expo_H_searchable']/1e4:.1f} | {sum(w['moon_excluded_s'] for w in u['windows'])/1e3:.1f} | {'dev' if u['dev'] else ''} |")
    c = Counter()
    for u in out["units"]:
        for w in u["windows"]:
            c[(u["rung"], "searchable" if w["searchable"] else "uncovered")] += 1
    md += ["", "Window census: " + ", ".join(f"{k[0]} {k[1]} {v}" for k, v in sorted(c.items())), ""]
    md += ["Coverage-without-statistic windows (< 1 ks):", ""]
    for u in out["units"]:
        for w in u["windows"]:
            if not w["searchable"]:
                md.append(f"- {u['channel']} {u['target']} {u['rung']} {w['t_ca_utc'][:10]}: {w['live_s']:.0f} s gated ({w['live_raw_s']:.0f} s raw, Moon-excluded {w['moon_excluded_s']:.0f} s)")
    (H.RES / "coverage_v1.md").write_text("\n".join(md))
    print("units", len(out["units"]), "searchable units", sum(u["searchable_unit"] for u in out["units"]))


if __name__ == "__main__":
    main()
