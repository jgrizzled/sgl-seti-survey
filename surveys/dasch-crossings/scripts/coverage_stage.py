"""DASCH crossings coverage stage (hypotheses.md v1.0 §§7, 9).

Re-derives the frozen coverage record from fresh snapshot-disciplined
queryexps pulls: per searched-unit candidate, the covered windows
(usable exposures whose [start, stop] interval overlaps the
flat-chord window and pass the timing gate) and the
coverage-without-statistic ledger. Metadata only — no photometric
quantity is touched.

Products (surveys/dasch-crossings/results/):
  coverage_v1_windows.ecsv  one row per (unit, event) with coverage state
  coverage_v1_summary.json  per-unit tallies + manifest (input hashes,
                            snapshot hashes, package versions)
Raw queryexps snapshots: runs/dasch/v1/coverage/queryexps/ (ignored;
sha256 in the summary manifest).
"""
from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.time import Time

sys.path.insert(0, str(Path(__file__).parent))
import daschlib as dl  # noqa: E402

RESULTS = Path(__file__).resolve().parents[1] / "results"
SNAPDIR = dl.RUNS / "coverage" / "queryexps"
EVENTS = dl.REPO / "crossings" / "universal_1885_v1" / "events.ecsv"
ERA = (Time("1885-01-01", scale="utc").jd, Time("1993-01-01", scale="utc").jd)


def cluster_positions(events: Table, channel: str) -> list[dict]:
    """Greedy-cluster per-event query loci within CLUSTER_TOL_DEG so a
    century of PM drift maps to a handful of queryexps pulls.
    Channel A searches at the star (blend); B at the relay locus
    (antipode) — hypotheses §1."""
    col = "star" if channel == "A" else "relay"
    clusters: list[dict] = []
    for e in events:
        ra = float(e[f"{col}_icrs_ra_deg"])
        dec = float(e[f"{col}_icrs_dec_deg"])
        cosd = np.cos(np.radians(dec))
        hit = None
        for c in clusters:
            if (abs(dec - c["dec"]) < dl.CLUSTER_TOL_DEG
                    and abs((ra - c["ra"]) * cosd) < dl.CLUSTER_TOL_DEG):
                hit = c
                break
        if hit is None:
            hit = {"ra": ra, "dec": dec, "event_ids": []}
            clusters.append(hit)
        hit["event_ids"].append(str(e["event_id"]))
    return clusters


def main() -> int:
    ev = Table.read(EVENTS)
    era = (ev["t_ca_tdb_jd"] > ERA[0]) & (ev["t_ca_tdb_jd"] < ERA[1])
    ev = ev[era & np.isin(ev["target_id"], dl.TARGETS)]

    rows = []
    unit_summary: dict[str, dict] = {}
    snapshots: dict[str, str] = {}

    for target in dl.TARGETS:
        for chan, link in dl.CHANNEL_LINK.items():
            sub = ev[(ev["target_id"] == target)
                     & (ev["link_direction"] == link)]
            clusters = cluster_positions(sub, chan)
            # one queryexps pull per cluster
            exp_by_event: dict[str, np.ndarray] = {}
            for ci, c in enumerate(clusters):
                snap = SNAPDIR / f"{target}_{chan}_c{ci}.json"
                resp = dl.post("dasch/dr7/queryexps",
                               {"ra_deg": c["ra"], "dec_deg": c["dec"]}, snap)
                snapshots[snap.name] = dl.sha256_file(snap)
                recs = []
                for r in dl.rows_of(resp):
                    if not dl.exposure_usable(r):
                        continue
                    iv = dl.exposure_interval_jd(r)
                    if iv is None:
                        continue
                    sigma = (dl.DATE_ONLY_SIGMA_DAYS if dl.date_only(r)
                             else 0.001)
                    recs.append((iv[0], iv[1], sigma,
                                 float(r["limMagApass"])))
                arr = np.array(recs) if recs else np.empty((0, 4))
                for eid in c["event_ids"]:
                    exp_by_event[eid] = arr

            for rung, r_au in dl.RUNGS.items():
                if not rung.startswith(chan):
                    continue
                unit = f"{target}/{rung}"
                n_events = n_cov = n_exp = n_timing_lost = 0
                for e in sub:
                    b = float(e["b_min_au"])
                    if b >= r_au:
                        continue
                    n_events += 1
                    hw = (np.sqrt(r_au**2 - b**2) * 1.496e8
                          / float(e["v_perp_km_s"]) / 86400.0)
                    lo = float(e["t_ca_tdb_jd"]) - hw
                    hi = float(e["t_ca_tdb_jd"]) + hw
                    arr = exp_by_event.get(str(e["event_id"]))
                    if arr is None or not len(arr):
                        state, k = "no_usable_exposures", 0
                    else:
                        inw = (arr[:, 1] > lo) & (arr[:, 0] < hi)
                        gated = inw & (arr[:, 2]
                                       <= dl.TIMING_GATE_FRACTION * hw)
                        k = int(gated.sum())
                        if k:
                            state = "covered"
                            n_cov += 1
                            n_exp += k
                        elif inw.any():
                            state = "timing_gate_lost"
                            n_timing_lost += 1
                        else:
                            state = "uncovered"
                    rows.append(dict(
                        unit=unit, target=target, channel=chan, rung=rung,
                        event_id=str(e["event_id"]),
                        t_ca_utc=str(e["t_ca_utc"]),
                        b_min_rsun=float(e["b_min_solar_radii"]),
                        half_width_days=float(hw), state=state,
                        n_inwindow_exposures=k))
                unit_summary[unit] = dict(
                    events=n_events, covered_windows=n_cov,
                    inwindow_exposures=n_exp,
                    timing_gate_lost=n_timing_lost)

    RESULTS.mkdir(exist_ok=True)
    t = Table(rows=rows)
    t.write(RESULTS / "coverage_v1_windows.ecsv",
            format="ascii.ecsv", overwrite=True)
    import astropy
    summary = dict(
        generated_utc=datetime.now(timezone.utc).isoformat(
            timespec="seconds"),
        hypotheses="hypotheses.md v1.0 (FROZEN 2026-08-26)",
        input_events=str(EVENTS.relative_to(dl.REPO)),
        input_events_sha256=dl.sha256_file(EVENTS),
        era_jd=list(ERA), units=unit_summary,
        gates=dict(listing="limMagApass present + wcssource in "
                           "{imwcs,catalog}",
                   timing=f"sigma_t <= {dl.TIMING_GATE_FRACTION} * "
                          "half-width; date-only sigma "
                          f"{dl.DATE_ONLY_SIGMA_DAYS} d",
                   overlap="interval-based [start, start+exptime]"),
        snapshots=snapshots,
        versions=dict(astropy=astropy.__version__,
                      numpy=np.__version__),
    )
    (RESULTS / "coverage_v1_summary.json").write_text(
        json.dumps(summary, indent=1) + "\n")
    for u, s in unit_summary.items():
        print(f"{u:24s} events {s['events']:3d}  covered "
              f"{s['covered_windows']:3d}  exp {s['inwindow_exposures']:4d}"
              f"  timing-lost {s['timing_gate_lost']:3d}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
