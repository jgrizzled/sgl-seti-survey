"""GALEX crossings blind confirmatory run (frozen chain: hypotheses
v1.0 + v1.1 + v1.2, threshold freeze v1.0, dev assessment).

Phase `search` — the first and only sanctioned in-window contact:
for every confirmatory band-unit (5 units / 15 trials) fetch the
in-window photons at the frozen loci, compute the three frozen
statistics (S_period tick-jittered per v1.2), the 8-control
thresholds (B: same-visit pseudo-positions with the z-family pattern
mirrored across the unit's events; A: the designated off-window
same-duration segments), and record exceedances (S > max(T, 0)).
Statistics are written to results/confirmatory_v1.json BEFORE any
injection is drawn. Constraint-only lanes (gj-1276 A NUV/FUV,
wolf-359 A FUV; the D1a forced_dev S_rate on gj-1276 B 2010) are
computed and labeled, zero trials.

Phase `completeness` — photon-level injections onto the real unit
series through the identical chain, against the locked thresholds:
persistent rates, boxcar pulses (0.05 / 0.5 s), periodic trains with
orbital drift sampled U(+/-2.56e-5). -> results/completeness_v1.json
"""

from __future__ import annotations

import json
import sys
import time as _time
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from sglsurvey.adapters.mast_gphoton import (
    GphotonClient, galex_ms_to_iso, galex_ms_to_mjd, group_visits,
    mjd_to_galex_ms)
from sglsurvey.snapshots import SnapshotStore
import galexlib as gl
from dev_stage import Dev, FIELDS, UNIT_LIVE_S, boresight_mean

SURVEY = REPO / "surveys" / "galex-crossings"
OUT = SURVEY / "results"
RUN = REPO / "runs" / "galex-crossings"
SEED = 20260826
ZP_EFF = 19.592          # dev step ii, frozen input

#: (unit_key) -> channel, field, band, [(event t_ca date, trial?)],
#: status. D1a: S_rate uses only rate_events.
UNITS = {
    "gj-1276_B_NUV": dict(field="gj-1276_antipode", band="NUV",
                          events=["2007-03-03", "2010-03-03"],
                          rate_events=["2007-03-03"],
                          status="confirmatory"),
    "gj-1276_B_FUV": dict(field="gj-1276_antipode", band="FUV",
                          events=["2007-03-03"],
                          rate_events=["2007-03-03"],
                          status="confirmatory"),
    "wolf-359_A_NUV": dict(field="wolf-359_star", band="NUV",
                           events=["2007-03-03"],
                           rate_events=["2007-03-03"],
                           status="confirmatory"),
    "ross-128_A_NUV": dict(field="ross-128_star", band="NUV",
                           events=["2007-03-17"],
                           rate_events=["2007-03-17"],
                           status="confirmatory"),
    "ross-128_A_FUV": dict(field="ross-128_star", band="FUV",
                           events=["2007-03-17"],
                           rate_events=["2007-03-17"],
                           status="confirmatory"),
    # constraint-only lanes (dev segment gate)
    "gj-1276_A_NUV": dict(field="gj-1276_star", band="NUV",
                          events=["2007-09-05"],
                          rate_events=["2007-09-05"],
                          status="constraint-only"),
    "gj-1276_A_FUV": dict(field="gj-1276_star", band="FUV",
                          events=["2007-09-05"],
                          rate_events=["2007-09-05"],
                          status="constraint-only"),
    "wolf-359_A_FUV": dict(field="wolf-359_star", band="FUV",
                           events=["2007-03-03"],
                           rate_events=["2007-03-03"],
                           status="constraint-only"),
}
#: the D1a forced_dev rate lane, reported separately
D1A_EVENT = ("gj-1276_B_NUV", "2010-03-03")
#: deterministic per-unit jitter-seed offsets (str hash is
#: process-randomized and must not be used)
UNIT_OFF = {k: 37 * i for i, k in enumerate(UNITS)}


class Conf(Dev):
    """Dev machinery plus sanctioned in-window fetches."""

    def event_row(self, field, date):
        tid, ch, ld, _ = FIELDS[field]
        m = ((self.events["target_id"] == tid)
             & (self.events["link_direction"] == ld))
        for ev in self.events[m]:
            if str(ev["t_ca_utc"]).startswith(date):
                return ev
        raise KeyError((field, date))

    def window_of(self, ev):
        b = float(ev["b_min_au"])
        half_d = (np.sqrt(0.1 ** 2 - b ** 2) * gl.KM_PER_AU
                  / float(ev["v_perp_km_s"]) / 86400.0)
        t_ca = float(ev["t_ca_mjd"])
        return (mjd_to_galex_ms(t_ca - half_d),
                mjd_to_galex_ms(t_ca + half_d))

    def inwindow_visits(self, field, ev):
        _, _, _, (ra, dec) = FIELDS[field]
        w0, w1 = self.window_of(ev)
        rows = self.client.aspect_near(ra, dec, 40.0, w0, w1,
                                       self.store)
        times = sorted(set(t for t, _, _ in rows))
        return group_visits(times)

    def unit_positions(self, field, ev, t_ms):
        """Locus aperture positions: A -> [star]; B -> z-grid,
        deduplicated at half-aperture separation."""
        tid, ch, ld, _ = FIELDS[field]
        if ch == "A":
            return [(float(ev["star_icrs_ra_deg"]),
                     float(ev["star_icrs_dec_deg"]))]
        mjd = galex_ms_to_mjd(t_ms)
        pos = [gl.relay_apparent(float(ev["star_icrs_ra_deg"]),
                                 float(ev["star_icrs_dec_deg"]),
                                 z, mjd) for z in gl.Z_GRID_AU]
        keep = []
        for p in pos:
            if all(gl.ang_arcmin(p[0], p[1], q[0], q[1]) * 60.0
                   > gl.APERTURE_ARCSEC / 2.0 for q in keep):
                keep.append(p)
        return keep


def series_stats(conf, frags, band, positions, jitter_off,
                 min_live=30):
    """Fetch + compute the three statistics for one aperture set on
    one visit; statistic = max over the (deduplicated) apertures."""
    best = {"S_rate": 0.0, "S_burst": 0.0, "S_period": 0.0}
    detail = []
    for i, (ra, dec) in enumerate(positions):
        usable = conf.usable_at(frags, ra, dec)
        live = sum(1 for b in usable.values() if band in b)
        if live < min_live:
            detail.append({"live_s": live, "n": 0,
                           "note": "below min live"})
            continue
        ph = conf.photons_at(frags, band, ra, dec, usable)
        stamps = sorted(s for s, b in usable.items() if band in b)
        span = float(stamps[-1] - stamps[0] + 1)
        sp = gl.s_period(ph, span,
                         jitter_rng=np.random.default_rng(
                             SEED + jitter_off + i))
        st = {"S_rate": gl.s_rate(ph, live),
              "S_burst": gl.s_burst(ph, live),
              "S_period": sp if sp is not None else None}
        detail.append({"live_s": live, "n": int(len(ph)), **st})
        for k in best:
            if st[k] is not None:
                best[k] = max(best[k], st[k])
    return best, detail


def b_control_positions(conf, field, frags, locus_positions):
    """The frozen 8-pseudo-position rule; each control carries the
    z-family pattern rotated rigidly about the boresight."""
    from dev_stage import b_controls
    base = locus_positions[0]
    ctr_base, rotations = b_controls(conf, field, frags, base)
    rows = conf.aspect_frags(frags)
    bra, bde = boresight_mean(rows)
    out = []
    for (cra, cde) in ctr_base:
        # rigid translation of the z-family offsets to the control
        pat = []
        for (ra, dec) in locus_positions:
            dra = (ra - base[0])
            dde = (dec - base[1])
            pat.append((cra + dra, cde + dde))
        out.append(pat)
    return out, rotations


def a_control_segments(conf, field, band, need_s):
    """The designated 8 off-window segments (frozen selection)."""
    segs = []
    for v0, v1 in conf.visits(field):
        pos = conf.field_pos(field, (v0 + v1) / 2)
        usable = conf.usable_at([(v0, v1)], pos[0], pos[1])
        secs = sorted(s for s, b in usable.items() if band in b)
        for i in range(len(secs) // need_s):
            seg = secs[i * need_s:(i + 1) * need_s]
            segs.append(((v0, v1), pos, seg))
    if len(segs) < 8:
        return None
    idx = sorted(set(int(round(x)) for x in
                     np.linspace(0, len(segs) - 1, 8)))
    while len(idx) < 8:      # deterministic fill on collisions
        for j in range(len(segs)):
            if j not in idx:
                idx.append(j)
                break
        idx = sorted(idx)
    return [segs[j] for j in idx[:8]]


def a_segment_stats(conf, seg_entry, band, jitter_off):
    (v0, v1), pos, seg = seg_entry
    usable_all = conf.usable_at([(v0, v1)], pos[0], pos[1])
    ph_all = conf.photons_at([(v0, v1)], band, pos[0], pos[1],
                             usable_all)
    ids = set(int(s) for s in seg)
    ph = np.array([t for t in ph_all if int(gl.sec_id(t)) in ids],
                  dtype=np.int64)
    need = len(seg)
    sp = gl.s_period(ph, float(need),
                     jitter_rng=np.random.default_rng(
                         SEED + jitter_off))
    return {"S_rate": gl.s_rate(ph, need),
            "S_burst": gl.s_burst(ph, need),
            "S_period": sp if sp is not None else None,
            "n": int(len(ph))}


def phase_search(conf):
    results = {}
    for key, u in UNITS.items():
        field, band = u["field"], u["band"]
        tid, ch, ld, _ = FIELDS[field]
        print(f"== {key} ({u['status']}) ==", flush=True)
        ev_stats, ev_detail = {}, {}
        for date in u["events"]:
            ev = conf.event_row(field, date)
            visits = conf.inwindow_visits(field, ev)
            if not visits:
                ev_detail[date] = "no in-window visit"
                continue
            frags = visits
            pos = conf.unit_positions(field, ev,
                                      (frags[0][0] + frags[-1][1])
                                      // 2)
            best, detail = series_stats(conf, frags, band, pos,
                                        jitter_off=UNIT_OFF[key])
            ev_stats[date] = best
            ev_detail[date] = detail
            print(f"  event {date}: {best}", flush=True)
        # unit statistic = max over eligible events
        unit = {}
        for stat in ("S_rate", "S_burst", "S_period"):
            elig = (u["rate_events"] if stat == "S_rate"
                    else u["events"])
            vals = [ev_stats[d][stat] for d in elig if d in ev_stats]
            unit[stat] = max(vals) if vals else None
        # controls
        thr, ctrl_detail, rotations = {}, [], None
        if u["status"] == "confirmatory":
            if ch == "B":
                per_ctrl = {s: [] for s in unit}
                for date in u["events"]:
                    ev = conf.event_row(field, date)
                    frags = conf.inwindow_visits(field, ev)
                    if not frags:
                        continue
                    pos = conf.unit_positions(
                        field, ev, (frags[0][0] + frags[-1][1]) // 2)
                    patterns, rotations = b_control_positions(
                        conf, field, frags, pos)
                    for k, pat in enumerate(patterns):
                        best, _ = series_stats(
                            conf, frags, band, pat,
                            jitter_off=2000 + UNIT_OFF[key] + k)
                        for stat in per_ctrl:
                            elig = (u["rate_events"]
                                    if stat == "S_rate"
                                    else u["events"])
                            if date in elig:
                                if len(per_ctrl[stat]) <= k:
                                    per_ctrl[stat].append(best[stat])
                                else:
                                    per_ctrl[stat][k] = max(
                                        per_ctrl[stat][k],
                                        best[stat])
                ctrl_detail = per_ctrl
                thr = {s: (max(v) if v else None)
                       for s, v in per_ctrl.items()}
            else:
                need = UNIT_LIVE_S[field]
                segs = a_control_segments(conf, field, band, need)
                assert segs is not None, f"{key}: segment gate FAIL"
                sc = [a_segment_stats(conf, s, band,
                                      3000 + UNIT_OFF[key] + j)
                      for j, s in enumerate(segs)]
                ctrl_detail = sc
                thr = {stat: max((c[stat] for c in sc
                                  if c[stat] is not None),
                                 default=None)
                       for stat in ("S_rate", "S_burst", "S_period")}
        exceed = {}
        if u["status"] == "confirmatory":
            for stat in unit:
                su, T = unit[stat], thr.get(stat)
                if su is None:
                    exceed[stat] = "gate-not-met (constraint-only)"
                elif T is None:
                    exceed[stat] = "no-control (constraint-only)"
                else:
                    exceed[stat] = bool(su > max(T, 0))
        results[key] = {"status": u["status"], "band": band,
                        "events": ev_stats, "detail": ev_detail,
                        "unit": unit, "thresholds": thr,
                        "controls": ctrl_detail,
                        "control_rotations": rotations,
                        "exceedances": exceed}
        print(f"  unit {unit} vs T {thr} -> {exceed}", flush=True)
    # D1a lane annotation
    key, date = D1A_EVENT
    results["_d1a_forced_dev"] = {
        "lane": "S_rate on gj-1276 B 2010-03-03",
        "value": results[key]["events"].get(date, {}).get("S_rate"),
        "note": "forced_dev constraint-only (hypotheses D1a); "
                "excluded from the unit S_rate max and the trial "
                "family; pre-freeze contact must be cited at any "
                "adjudication"}
    (OUT / "confirmatory_v1.json").write_text(
        json.dumps({"chain": "hypotheses v1.0+v1.1+v1.2, threshold "
                             "freeze v1.0, dev assessment "
                             "(15 trials)",
                    "zp_eff": ZP_EFF, "results": results},
                   indent=2, default=float) + "\n")
    print("search phase locked -> confirmatory_v1.json")
    return results


def phase_completeness(conf, results):
    rng = np.random.default_rng(SEED + 77)
    out = {}
    for key, u in UNITS.items():
        if u["status"] != "confirmatory":
            continue
        field, band = u["field"], u["band"]
        # completeness on the unit's deepest event
        date = max(u["events"],
                   key=lambda d: (results[key]["events"]
                                  .get(d, {}).get("S_rate") or 0,
                                  d))
        # re-fetch the series (snapshotted; identical chain)
        ev = conf.event_row(field, date)
        frags = conf.inwindow_visits(field, ev)
        pos = conf.unit_positions(field, ev,
                                  (frags[0][0] + frags[-1][1]) // 2)
        ra, dec = pos[0]
        usable = conf.usable_at(frags, ra, dec)
        live = sum(1 for b in usable.values() if band in b)
        stamps = np.array(sorted(s * 1000 + 995
                                 for s, b in usable.items()
                                 if band in b), dtype=np.int64)
        if len(stamps) < 30:
            out[key] = {"note": "below min live"}
            continue
        span = (stamps[-1] + 1000 - stamps[0]) / 1000.0
        ph_u = conf.photons_at(frags, band, ra, dec, usable)
        T = results[key]["thresholds"]
        rows = {"event": date, "live_s": live,
                "persistent": [], "pulse": [], "train": []}

        def rec(times, stat, thr, joff):
            if thr is None:
                return None
            m = gl.merge_series(ph_u, times)
            if stat == "S_rate":
                return gl.s_rate(m, live) > max(thr, 0)
            if stat == "S_burst":
                return gl.s_burst(m, live) > max(thr, 0)
            sp = gl.s_period(m, span,
                             jitter_rng=np.random.default_rng(
                                 SEED + 5000 + joff))
            return (sp or 0.0) > max(thr, 0)

        for rate in (0.02, 0.04, 0.07, 0.12, 0.2, 0.35, 0.6, 1.0):
            n = sum(bool(rec(gl.inject_persistent(rng, rate, stamps),
                             "S_rate", T["S_rate"], 0))
                    for _ in range(25))
            rows["persistent"].append(
                {"rate_cts_s": rate, "recovered": n, "of": 25})
        for w in (0.05, 0.5):
            for n_ph in (2, 3, 4, 5, 6, 8, 12):
                n = sum(bool(rec(gl.inject_pulse(rng, n_ph, w,
                                                 stamps),
                                 "S_burst", T["S_burst"], 0))
                        for _ in range(25))
                rows["pulse"].append({"width_s": w,
                                      "n_photons": n_ph,
                                      "recovered": n, "of": 25})
        p_grid = [p for p in (0.05, 0.5, 5.0, 50.0)
                  if p <= span / 3.0]
        for p in p_grid:
            for n_tot in (15, 30, 60, 120):
                frac = n_tot * p / span
                n = 0
                for j in range(5):
                    drift = rng.uniform(-gl.ORBIT_DRIFT_RATE,
                                        gl.ORBIT_DRIFT_RATE)
                    n += bool(rec(
                        gl.inject_train(rng, p, frac, 0.1, stamps,
                                        drift_rate=drift),
                        "S_period", T["S_period"], j))
                rows["train"].append({"period_s": p,
                                      "n_injected": n_tot,
                                      "recovered": n, "of": 5})
            print(f"[c] {key} P={p} done", flush=True)
        out[key] = rows
        print(f"[c] {key} complete", flush=True)
    (OUT / "completeness_v1.json").write_text(
        json.dumps(out, indent=2, default=float) + "\n")
    print("completeness -> completeness_v1.json")


if __name__ == "__main__":
    conf = Conf()
    t0 = _time.monotonic()
    results = phase_search(conf)
    print(f"== search done in {_time.monotonic() - t0:.0f}s ==",
          flush=True)
    if "search-only" not in sys.argv:
        t0 = _time.monotonic()
        phase_completeness(conf, results)
        print(f"== completeness done in "
              f"{_time.monotonic() - t0:.0f}s ==", flush=True)
