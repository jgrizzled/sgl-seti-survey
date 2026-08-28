"""DASCH crossings completeness (thresholds.md v1.0 §completeness).

Per searched unit and covered window: platephot-vs-querycat
field-star recovery in the shared subregion (r <= 600 arcsec) gives
the detection-efficiency-vs-magnitude curve (Wilson intervals) and
the per-window m90; the C1 analogue — `limMag*` columns alone never
qualify a constraint. Runs from the confirmatory/dev platephot
snapshots (no new photometric contact); querycat truth pulls are
catalog-metadata class.

Products: results/completeness_v1.json
"""
from __future__ import annotations

import json
import math
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import daschlib as dl        # noqa: E402
import search_lib as sl      # noqa: E402

D = Path(__file__).resolve().parents[1]
R_SUB = 600.0  # arcsec, shared truth/recovery aperture
BINS = np.arange(8.0, 19.0, 1.0)


def wilson(k: int, n: int, z: float = 1.6449):
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (max(0.0, c - h), min(1.0, c + h))


def truth_set(center, snapdir: Path, key: str):
    snap = snapdir / f"querycat_{key}.json"
    resp = dl.post("dasch/dr7/querycat",
                   {"ra_deg": center[0], "dec_deg": center[1],
                    "radius_arcsec": R_SUB, "refcat": "apass"}, snap)
    rows = sl.rows_of_safe(resp)
    out = {}
    cosd = np.cos(np.radians(center[1]))
    for r in rows:
        try:
            m = float(r["stdmag"])
        except (TypeError, ValueError):
            continue
        if m > 30:
            continue
        d = 3600 * np.hypot((float(r["ra_deg"]) - center[0]) * cosd,
                            float(r["dec_deg"]) - center[1])
        if d <= R_SUB:
            out[r["ref_number"]] = m
    return out


def main() -> int:
    out = {"generated_utc": datetime.now(timezone.utc).isoformat(
               timespec="seconds"), "units": {}}
    for stage, snapbase in (("confirmatory", dl.RUNS / "confirmatory"),
                            ("dev", dl.RUNS / "dev")):
        ppdir = snapbase / "platephot"
        if not ppdir.exists():
            continue
        units = (sl.FREEZE["units"]["confirmatory"]
                 if stage == "confirmatory" else
                 [u for u in sl.FREEZE["units"]["dev"]
                  if not u.startswith("van-maanen/A")])
        for unit in units:
            target, rung = unit.split("/")
            chan = rung[0]
            if chan == "A" and sl.A_ROUTING[target] is not None:
                continue  # lightcurve units: no platephot curve
            wins = sl.covered_windows(unit)
            per_window = []
            for w in wins:
                ev = sl.event_by_id(str(w["event_id"]))
                hw = float(w["half_width_days"])
                tca = float(ev["t_ca_tdb_jd"])
                exps = sl.inwindow_exposures(target, chan, tca - hw,
                                             tca + hw, hw)
                for e in exps:
                    t_mid = e["_t_mid_jd"]
                    if chan == "B":
                        center = sl.relay_apparent(ev, sl.Z_GRID[0],
                                                   t_mid)
                    else:
                        center = (float(ev["star_icrs_ra_deg"]),
                                  float(ev["star_icrs_dec_deg"]))
                    pid = (f"{e['series']}{int(e['platenum']):05d}")
                    snap = ppdir / (f"{pid}_s{e['solnum']}"
                                    f"_{center[0]:.3f}"
                                    f"_{center[1]:.3f}.json")
                    if not snap.exists():
                        continue
                    truth = truth_set(center, snapbase / "querycat_cx",
                                      f"{target}_{chan}"
                                      f"_{center[0]:.2f}_{center[1]:.2f}")
                    rows = sl.rows_of_safe(
                        json.loads(snap.read_bytes()))
                    cosd = np.cos(np.radians(center[1]))
                    rec = set()
                    for r in rows:
                        if not r.get("ref_number"):
                            continue
                        d = 3600 * np.hypot(
                            (float(r["ra_deg"]) - center[0]) * cosd,
                            float(r["dec_deg"]) - center[1])
                        if d <= R_SUB:
                            rec.add(r["ref_number"])
                    kbin = defaultdict(lambda: [0, 0])
                    for ref, m in truth.items():
                        b = int(np.digitize(m, BINS))
                        kbin[b][1] += 1
                        if ref in rec:
                            kbin[b][0] += 1
                    # m90: faintest bin center with eff>=0.9 and the
                    # next bin below 0.9 (monotone walk from bright end)
                    m90 = None
                    for b in sorted(kbin):
                        k, n = kbin[b]
                        if n < 5:
                            continue
                        eff = k / n
                        c = BINS[min(b, len(BINS) - 1)] - 0.5
                        if eff >= 0.9:
                            m90 = c
                        else:
                            break
                    per_window.append(dict(
                        event_id=str(w["event_id"]), plate=pid,
                        n_truth=len(truth), n_recovered=len(rec),
                        m90=m90,
                        lim_mag=e.get("limMagApass"),
                        bins={str(BINS[min(b, len(BINS)-1)] - 0.5):
                              dict(k=v[0], n=v[1],
                                   wilson90=wilson(v[0], v[1]))
                              for b, v in sorted(kbin.items())
                              if v[1] > 0}))
            m90s = [p["m90"] for p in per_window if p["m90"]]
            out["units"][unit] = dict(
                stage=stage, n_exposures=len(per_window),
                m90_median=(float(np.median(m90s)) if m90s else None),
                m90_best=(max(m90s) if m90s else None),
                windows=per_window)
            print(unit, "exp", len(per_window), "m90 median",
                  out["units"][unit]["m90_median"], "best",
                  out["units"][unit]["m90_best"], flush=True)
    (D / "results" / "completeness_v1.json").write_text(
        json.dumps(out, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
