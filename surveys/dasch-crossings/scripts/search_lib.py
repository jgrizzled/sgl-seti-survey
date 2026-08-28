"""Frozen search machinery shared by the dev and confirmatory stages
(thresholds.md v1.0). Everything here implements the frozen
constructions verbatim; stage scripts only choose units.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.time import Time

sys.path.insert(0, str(Path(__file__).parent))
import daschlib as dl  # noqa: E402

D = Path(__file__).resolve().parents[1]
FREEZE = json.loads((D / "configs" / "threshold_freeze_v1.json").read_text())
BC = FREEZE["b_construction"]
AC = FREEZE["a_construction"]
Z_GRID = BC["z_grid_au"]
R_MATCH = BC["match_radius_arcsec"]
# Amendment v1.2: SUSPECTED_DEFECT demoted from fatal to annotation
# (the B-chain positive control measured it killing a real transient).
FATAL = BC["fatal_aflags_mask"] & ~(1 << 25)
SUSPECTED_DEFECT = 1 << 25
RING = [tuple(v) for v in BC["ring_offsets_arcsec"]]
T_OFFS = AC["temporal_offsets_days"]
T_REDRAWS = AC["temporal_redraws_days"]
COVERAGE = Table.read(D / "results" / "coverage_v1_windows.ecsv")
EVENTS = Table.read(dl.REPO / "crossings" / "universal_1885_v1"
                    / "events.ecsv")
SNAP_COV = dl.RUNS / "coverage" / "queryexps"

TOO_BRIGHT_A = 1 << 29   # AFlags.TOO_BRIGHT
SATURATED_B = 1 << 2     # BFlags.SATURATED

#: Amendment v1.1 lightcurve routing (probed 2026-08-26): detected-
#: regime A units -> (refcat, gsc_bin_index, ref_number); None =
#: limits-only regime (no refcat entry tracks the star, or the star
#: sits below every plate limit).
A_ROUTING = {
    "van-maanen": ("atlas", 92397935, 9717808673),   # ATLAS carries PM
    "ross-128": ("apass", 85669484, 411474441004816),
    "gj-1276": None,      # APASS PM valid but stdmag 18.0 sub-limit
    "wolf-359": None,     # both refcat entries broken (v1.1)
    "teegarden": None,    # no refcat entry
}


def relay_apparent(ev, z_au: float, t_jd: float):
    """Apparent ICRS ra/dec of a relay at z_au on the anti-star axis
    (the PTF coverage_intersect construction, verbatim)."""
    from astropy.coordinates import get_body_barycentric
    t = Time(t_jd, format="jd", scale="utc")
    sun = get_body_barycentric("sun", t).xyz.to_value("AU")
    earth = get_body_barycentric("earth", t).xyz.to_value("AU")
    ra = np.radians(float(ev["star_icrs_ra_deg"]))
    de = np.radians(float(ev["star_icrs_dec_deg"]))
    u_star = np.array([np.cos(de) * np.cos(ra), np.cos(de) * np.sin(ra),
                       np.sin(de)])
    v = (sun - z_au * u_star) - earth
    v /= np.linalg.norm(v)
    return (float(np.degrees(np.arctan2(v[1], v[0])) % 360.0),
            float(np.degrees(np.arcsin(np.clip(v[2], -1, 1)))))


def event_by_id(eid: str):
    return EVENTS[EVENTS["event_id"] == eid][0]


def covered_windows(unit: str) -> Table:
    sub = COVERAGE[(COVERAGE["unit"] == unit)
                   & (COVERAGE["state"] == "covered")]
    return sub


def inwindow_exposures(target: str, chan: str, lo: float, hi: float,
                       hw: float) -> list[dict]:
    """Re-derive the usable in-window exposure rows (with plate ids)
    from the coverage-stage queryexps snapshots."""
    out, seen = [], set()
    for snap in sorted(SNAP_COV.glob(f"{target}_{chan}_c*.json")):
        for r in dl.rows_of(json.loads(snap.read_bytes())):
            if not dl.exposure_usable(r):
                continue
            iv = dl.exposure_interval_jd(r)
            if iv is None or not (iv[1] > lo and iv[0] < hi):
                continue
            sigma = dl.DATE_ONLY_SIGMA_DAYS if dl.date_only(r) else 0.001
            if sigma > dl.TIMING_GATE_FRACTION * hw:
                continue
            key = (r["series"], r["platenum"], r["solnum"], r["expnum"])
            if key in seen:
                continue
            seen.add(key)
            r["_t_mid_jd"] = 0.5 * (iv[0] + iv[1])
            out.append(r)
    return out


def _seg_dist_arcsec(pra, pdec, apts):
    """Min distance (arcsec) from points (pra, pdec) to the polyline
    through apts [(ra, dec), ...], small-angle tangent plane at
    apts[0]."""
    ra0, de0 = apts[0]
    cosd = np.cos(np.radians(de0))
    P = np.column_stack([(np.asarray(pra) - ra0) * cosd * 3600.0,
                         (np.asarray(pdec) - de0) * 3600.0])
    V = np.array([[(a - ra0) * cosd * 3600.0, (d - de0) * 3600.0]
                  for a, d in apts])
    best = np.full(len(P), np.inf)
    for i in range(len(V) - 1):
        a, b = V[i], V[i + 1]
        ab = b - a
        L2 = ab @ ab
        t = np.clip(((P - a) @ ab) / L2, 0, 1) if L2 > 0 else 0.0
        proj = a + np.outer(np.atleast_1d(t), ab)
        best = np.minimum(best, np.linalg.norm(P - proj, axis=1))
    return best


def _row_ok(r: dict, strict: bool = False) -> bool:
    af = int(float(r["aflags"] or 0))
    bf = int(float(r["bflags"] or 0))
    rej = r.get("reject_flag", "")
    mask = FATAL | (BC["strict_extra_mask"] if strict else 0)
    return (not (af & mask) and not (af & TOO_BRIGHT_A)
            and not (bf & SATURATED_B)
            and (rej in ("", "0")))


def platephot_hits(exposure: dict, ev, t_mid: float, snapdir: Path,
                   center: tuple[float, float], stats: dict,
                   track: bool) -> dict:
    """One platephot pull; returns hit counts at the locus and the 8
    ring offsets. `track`=True uses the z-family polyline (channel B);
    False a single point (limits-only A at the star)."""
    plate_id = f"{exposure['series']}{int(exposure['platenum']):05d}"
    snap = snapdir / (f"{plate_id}_s{exposure['solnum']}"
                      f"_{center[0]:.3f}_{center[1]:.3f}.json")
    resp = dl.post("dasch/dr7/platephot",
                   {"plate_id": plate_id,
                    "solution_number": int(exposure["solnum"]),
                    "center_ra_deg": center[0],
                    "center_dec_deg": center[1], "refcat": "apass"}, snap)
    stats["snapshots"][snap.name] = dl.sha256_file(snap)
    rows = dl.rows_of(resp) if isinstance(resp, list) else []
    stats["platephot_rows"].append(len(rows))
    ok = [r for r in rows if _row_ok(r)]
    cand = [r for r in ok if not r.get("ref_number")]
    # refcat-matched flux-anomalous class: needs the querycat join —
    # evaluated at adjudication; primary hit class = uncatalogued.
    if not cand:
        return {"locus": 0, "ring": [0] * len(RING), "hits_detail": []}
    pra = np.array([float(r["ra_deg"]) for r in cand])
    pdec = np.array([float(r["dec_deg"]) for r in cand])
    if track:
        apts = [relay_apparent(ev, z, t_mid) for z in Z_GRID]
    else:
        apts = [center, center]
    d0 = _seg_dist_arcsec(pra, pdec, apts)
    hits_detail = [
        dict(dist_arcsec=round(float(d0[i]), 1),
             mag=cand[i].get("magcal_magdep"),
             aflags=cand[i].get("aflags"),
             suspected_defect=bool(int(float(cand[i]["aflags"] or 0))
                                   & SUSPECTED_DEFECT),
             fwhm_pix=cand[i].get("fwhm_pix"),
             ellipticity=cand[i].get("ellipticity"),
             plate=plate_id)
        for i in np.where(d0 <= R_MATCH)[0]]
    res = {"locus": int((d0 <= R_MATCH).sum()), "ring": [],
           "hits_detail": hits_detail}
    cosd = np.cos(np.radians(center[1]))
    for dx, dy in RING:
        shifted = [(a + dx / 3600.0 / cosd, d + dy / 3600.0)
                   for a, d in apts]
        res["ring"].append(int((_seg_dist_arcsec(pra, pdec, shifted)
                                <= R_MATCH).sum()))
    return res


def run_hit_unit(unit: str, snapdir: Path, track: bool,
                 quiescent_mag: float | None = None) -> dict:
    """B-channel (track=True) or limits-only A (track=False) unit:
    S_event / S_stack hit statistics with ring controls.

    Amendment v1.3 (limits-only): `quiescent_mag` = the star's
    measured off-window quiescent level; a locus hit not >= 1.0 mag
    brighter is classified `quiescent_star` (annotated, excluded
    from S). Ring hits are never gated (conservative)."""
    target, rung = unit.split("/")
    chan = rung[0]
    wins = covered_windows(unit)
    stats = {"snapshots": {}, "platephot_rows": []}
    per_window = []
    for w in wins:
        ev = event_by_id(str(w["event_id"]))
        hw = float(w["half_width_days"])
        tca = float(ev["t_ca_tdb_jd"])
        exps = inwindow_exposures(target, chan, tca - hw, tca + hw, hw)
        wl, wr, wd = 0, [0] * len(RING), []
        used = []
        for e in exps:
            t_mid = e["_t_mid_jd"]
            if track:
                center = relay_apparent(ev, Z_GRID[0], t_mid)
            else:
                center = (float(ev["star_icrs_ra_deg"]),
                          float(ev["star_icrs_dec_deg"]))
            h = platephot_hits(e, ev, t_mid, snapdir, center, stats, track)
            if not track:
                # v1.3 quiescent-star gate; rule 3 fallback = own
                # plate limit when no off-window level exists
                ref = (quiescent_mag if quiescent_mag is not None
                       else float(e["limMagApass"]))
                kept = []
                for hd_ in h["hits_detail"]:
                    try:
                        bright = float(hd_["mag"]) <= ref - 1.0
                    except (TypeError, ValueError):
                        bright = True
                    hd_["quiescent_star"] = not bright
                    if bright:
                        kept.append(hd_)
                h["locus"] = len(kept)
            wl += h["locus"]
            wr = [a + b for a, b in zip(wr, h["ring"])]
            wd.extend(h["hits_detail"])
            used.append(f"{e['series']}{int(e['platenum']):05d}"
                        f"/{e['solnum']}")
        per_window.append(dict(event_id=str(w["event_id"]),
                               t_ca=str(w["t_ca_utc"]),
                               n_exposures=len(exps), locus_hits=wl,
                               ring_hits=wr, hits_detail=wd,
                               plates=used))
    S_event = max((w["locus_hits"] for w in per_window), default=0)
    S_stack = sum(w["locus_hits"] for w in per_window)
    T_event = max((max(w["ring_hits"][i] for w in per_window)
                   for i in range(len(RING))), default=0) \
        if per_window else 0
    ring_tot = [sum(w["ring_hits"][i] for w in per_window)
                for i in range(len(RING))]
    T_stack = max(ring_tot, default=0)
    return dict(unit=unit, statistic="hits",
                S_event=S_event, T_event=T_event,
                exceed_event=bool(S_event > max(T_event, 0)),
                S_stack=S_stack, T_stack=T_stack,
                exceed_stack=bool(S_stack > max(T_stack, 0)),
                ring_totals=ring_tot, windows=per_window,
                platephot_row_counts=sorted(set(
                    stats["platephot_rows"])),
                snapshots=stats["snapshots"])


def _robust_z(inw: np.ndarray, base: np.ndarray) -> float:
    mad = np.median(np.abs(base - np.median(base))) * 1.4826
    if mad <= 0 or not len(inw):
        return np.nan
    return float((np.median(base) - np.median(inw))
                 / (mad / np.sqrt(len(inw))))  # brightening positive


def run_lightcurve_unit(unit: str, lc_rows: list[dict]) -> dict:
    """Detected-regime channel A: robust-z window excess with the 8
    temporal pseudo-window controls. Amendment v1.1: era-local
    baseline (±5 yr per window, padded windows excluded) and
    Stouffer S_stack. `lc_rows` = usable rows of the routed refcat
    entry (v1.1 routing rule)."""
    wins = covered_windows(unit)
    tca = np.array([float(event_by_id(str(w["event_id"]))["t_ca_tdb_jd"])
                    for w in wins])
    hw = np.array([float(w["half_width_days"]) for w in wins])
    t = np.array([float(r["date_jd"]) for r in lc_rows])
    m = np.array([float(r["magcal_magdep"]) for r in lc_rows])
    LOCAL_D = 5 * 365.25

    def stat(shift: float):
        pad = np.zeros(len(t), bool)
        for c, h in zip(tca + shift, hw):
            pad |= np.abs(t - c) <= 2 * h
        zs, n_tot = [], 0
        for c, h in zip(tca + shift, hw):
            inw = np.abs(t - c) <= h
            if not inw.any():
                continue
            base = m[(np.abs(t - c) <= LOCAL_D) & ~pad]
            if len(base) < AC["min_offwindow_epochs"]:
                continue
            z = _robust_z(m[inw], base)
            if np.isfinite(z):
                zs.append(z)
                n_tot += int(inw.sum())
        if not zs:
            return np.nan, np.nan, 0
        return (float(np.max(zs)),
                float(np.sum(zs) / np.sqrt(len(zs))), n_tot)

    S_event, S_stack, n_in = stat(0.0)
    ctrl = []
    redraws = list(T_REDRAWS)
    for off in T_OFFS:
        se, ss, k = stat(float(off))
        while (not np.isfinite(se) or k == 0) and redraws:
            se, ss, k = stat(float(redraws.pop(0)))
        ctrl.append((se, ss))
    ev_ctrl = [c[0] for c in ctrl if np.isfinite(c[0] if c[0] is not
               None else np.nan)]
    st_ctrl = [c[1] for c in ctrl if np.isfinite(c[1] if c[1] is not
               None else np.nan)]
    T_event = max(ev_ctrl) if ev_ctrl else np.nan
    T_stack = max(st_ctrl) if st_ctrl else np.nan
    return dict(unit=unit, statistic="lightcurve_z",
                n_inwindow_rows=n_in,
                S_event=None if not np.isfinite(S_event) else S_event,
                T_event=float(T_event),
                exceed_event=bool(np.isfinite(S_event)
                                  and S_event > max(T_event, 0)),
                S_stack=None if not np.isfinite(S_stack) else S_stack,
                T_stack=float(T_stack),
                exceed_stack=bool(np.isfinite(S_stack)
                                  and S_stack > max(T_stack, 0)),
                controls=[[None if not np.isfinite(x) else float(x)
                           for x in c] for c in ctrl])


def measure_quiescent(target: str, unit: str, snapdir: Path,
                      registry_astrometry: dict,
                      max_plates: int = 12) -> dict:
    """Amendment v1.3 rule 1: the star's off-window quiescent level
    from uncatalogued detections within R_MATCH of the propagated
    position on the deepest off-window usable plates."""
    import csv as _csv
    import io as _io
    a = registry_astrometry
    chan = unit.split("/")[1][0]
    cov = covered_windows(unit)
    tcas = [float(event_by_id(str(w["event_id"]))["t_ca_tdb_jd"])
            for w in cov]
    pool, seen = [], set()
    for snap in sorted(SNAP_COV.glob(f"{target}_{chan}_c*.json")):
        for r in _csv.DictReader(_io.StringIO(
                "\n".join(json.loads(snap.read_bytes())))):
            if not dl.exposure_usable(r) or dl.date_only(r):
                continue
            key = (r["series"], r["platenum"], r["solnum"])
            if key in seen:
                continue
            seen.add(key)
            iv = dl.exposure_interval_jd(r)
            if iv is None or any(abs(iv[0] - c) < 12 for c in tcas):
                continue
            pool.append(r)
    pool.sort(key=lambda r: -float(r["limMagApass"]))
    mags, n_hit, n_test = [], 0, 0
    for r in pool[:max_plates]:
        iv = dl.exposure_interval_jd(r)
        t_mid = 0.5 * (iv[0] + iv[1])
        yr = (t_mid - 2451545.0) / 365.25 + 2000.0
        dt = yr - 2016.0
        pra = (a["ra_deg"] + a["pm_ra_cosdec_mas_per_yr"] / 3.6e6 * dt
               / np.cos(np.radians(a["dec_deg"])))
        pde = a["dec_deg"] + a["pm_dec_mas_per_yr"] / 3.6e6 * dt
        pid = f"{r['series']}{int(r['platenum']):05d}"
        snap = snapdir / f"quiescent_{pid}_s{r['solnum']}.json"
        resp = dl.post("dasch/dr7/platephot",
                       {"plate_id": pid,
                        "solution_number": int(r["solnum"]),
                        "center_ra_deg": pra, "center_dec_deg": pde,
                        "refcat": "apass"}, snap)
        rows = rows_of_safe(resp)
        n_test += 1
        ok = [x for x in rows if _row_ok(x) and not x.get("ref_number")]
        if not ok:
            continue
        d = np.array([3600 * np.hypot(
            (float(x["ra_deg"]) - pra) * np.cos(np.radians(pde)),
            float(x["dec_deg"]) - pde) for x in ok])
        j = int(np.argmin(d))
        if d[j] <= R_MATCH:
            n_hit += 1
            try:
                mags.append(float(ok[j]["magcal_magdep"]))
            except (TypeError, ValueError):
                pass
    return dict(n_tested=n_test, n_detections=n_hit,
                quiescent_mag=(float(np.median(mags)) if mags else None),
                mags=sorted(round(m, 2) for m in mags))


def rows_of_safe(resp):
    return dl.rows_of(resp) if isinstance(resp, list) else []
