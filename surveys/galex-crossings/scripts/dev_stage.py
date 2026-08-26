"""GALEX crossings dev stage (threshold freeze v1.0, D8: pseudo-units
only — every visit touched is guarded off-window against the union of
that field's in-era 0.1 AU flat-chord windows; the confirmatory
in-window visits are unreachable by construction).

Steps (results merged into results/dev_v1.json):
  i    flag-64 astrometric verification (amendment v1.1 item)
  ii   live-time/aperture calibration -> effective ZP + <= 0.2 mag gate
  iii  control-validity census (B rotation rule; A segment supply)
  iv   flare census at the three A stars (veto template input;
       conditional S_burst positive control)
  v    end-to-end pseudo-units: statistics + controls + injections

GALEX off-window data at these fields is largely petal-pattern
fragments (60-110 s dwells); dev works on observation *blocks*
(fragments merged within 30 min) for calibration/astrometry and on
the one contiguous 1,380 s flag-64 segment (2009-03-08) for the
B-long pseudo-unit.
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

SURVEY = REPO / "surveys" / "galex-crossings"
OUT = SURVEY / "results"
RUN = REPO / "runs" / "galex-crossings"
DEV_JSON = OUT / "dev_v1.json"

ERA_UTC = ("2003-06-07T05:02:17.995", "2013-05-01T18:02:22.995")
SEED = 20260826

FIELDS = {
    "gj-1276_antipode": ("gj-1276", "B", "outbound", (163.5170, 6.7936)),
    "gj-1276_star": ("gj-1276", "A", "inbound", (343.4781, -6.7833)),
    "wolf-359_star": ("wolf-359", "A", "inbound", (164.1111, 7.0082)),
    "ross-128_star": ("ross-128", "A", "inbound", (176.9365, 0.8016)),
}
UNIT_LIVE_S = {"gj-1276_star": 92, "wolf-359_star": 97,
               "ross-128_star": 110, "gj-1276_antipode": 109}


class Dev:
    def __init__(self):
        self.client = GphotonClient()
        self.store = SnapshotStore(RUN)
        t = Table.read(REPO / "crossings" / "universal_v1" /
                       "events.ecsv")
        mjd = Time(list(t["t_ca_utc"]), format="isot", scale="utc").mjd
        t["t_ca_mjd"] = mjd
        era = Time(list(ERA_UTC), format="isot", scale="utc").mjd
        side = np.asarray(t["axis_distance_au"])
        keep = ((mjd >= era[0]) & (mjd <= era[1])
                & (t["validity"] == "valid")
                & (((t["link_direction"] == "inbound") & (side > 0))
                   | ((t["link_direction"] == "outbound") & (side < 0))))
        self.events = t[keep]
        self.era_ms = [mjd_to_galex_ms(x) for x in era]
        self.windows = {
            name: gl.windows_ms(self.events, tid, ld, 0.1,
                                mjd_to_galex_ms)
            for name, (tid, ch, ld, _) in FIELDS.items()}
        self.tracks = {name: gl.star_track(self.events, tid, ld)
                       for name, (tid, ch, ld, _) in FIELDS.items()}
        self._visits = {}
        self._aspect = {}
        self._mcat = {}

    # -- off-window visit fragments and blocks -------------------------

    def visits(self, field):
        if field not in self._visits:
            _, _, _, (ra, dec) = FIELDS[field]
            rows = self.client.aspect_near(ra, dec, 40.0,
                                           self.era_ms[0],
                                           self.era_ms[1], self.store)
            times = sorted(set(t for t, _, _ in rows))
            out = []
            for v0, v1 in group_visits(times):
                try:
                    gl.assert_off_window(v0, v1 + 1000,
                                         self.windows[field])
                except RuntimeError:
                    continue
                out.append((v0, v1))
            self._visits[field] = out
        return self._visits[field]

    def blocks(self, field, min_total_ms=800_000,
               gap_ms=1_800_000):
        blocks, cur = [], []
        for v in self.visits(field):
            if cur and v[0] - cur[-1][1] > gap_ms:
                blocks.append(cur)
                cur = []
            cur.append(v)
        if cur:
            blocks.append(cur)
        return [b for b in blocks
                if sum(v1 - v0 for v0, v1 in b) >= min_total_ms]

    # -- cached fetches -------------------------------------------------

    def aspect_frags(self, frags):
        rows = []
        for v0, v1 in frags:
            key = (v0, v1)
            if key not in self._aspect:
                self._aspect[key] = self.client.aspect_range(
                    v0, v1 + 1000, self.store)
            rows.extend(self._aspect[key])
        return rows

    def usable_at(self, frags, ra, dec):
        """{second_id: band} across fragments, v1.1 gates."""
        out = {}
        for t_ms, ra0, de0, band, flag in self.aspect_frags(frags):
            s = int(gl.sec_id(t_ms))
            if s in out or flag % 2 != 0:
                continue
            if gl.ang_arcmin(ra0, de0, ra, dec) > gl.USABLE_ARCMIN:
                continue
            out[s] = band
        return out

    def photons_at(self, frags, band, ra, dec, usable):
        times = []
        for v0, v1 in frags:
            ph = gl.photons_aperture(self.client, self.store, band,
                                     ra, dec, v0, v1 + 1000, usable)
            times.extend(ph)
        return np.array(sorted(times), dtype=np.int64)

    def mcat(self, field, half=0.25):
        if field not in self._mcat:
            _, _, _, (ra, dec) = FIELDS[field]
            self._mcat[field] = self.client.mcat_box(ra, dec, half,
                                                     self.store,
                                                     limit=500)
        return self._mcat[field]

    def field_pos(self, field, t_ms):
        tid, ch, ld, _ = FIELDS[field]
        mjd = galex_ms_to_mjd(t_ms)
        ra_s, de_s = self.tracks[field](mjd)
        if ch == "A":
            return ra_s, de_s
        return gl.relay_apparent(ra_s, de_s, 550.0, mjd)


def merge(step, payload):
    data = json.loads(DEV_JSON.read_text()) if DEV_JSON.exists() else {}
    data[step] = payload
    DEV_JSON.write_text(json.dumps(data, indent=2, default=float)
                        + "\n")
    print(f"[{step}] written to {DEV_JSON.name}")


def boresight_mean(rows):
    return (float(np.mean([r for _, r, _, _, _ in rows])),
            float(np.mean([d for _, _, d, _, _ in rows])))


# ---------------------------------------------------------------- i --

def step_i(dev):
    """Flag-64 astrometry on off-window blocks at the antipode."""
    field = "gj-1276_antipode"
    cands = sorted((float(n), float(r), float(d))
                   for _, n, f, r, d in dev.mcat(field)
                   if n is not None and 16.5 <= float(n) <= 19.5)
    results = []
    for frags in dev.blocks(field):
        rows = dev.aspect_frags(frags)
        f64 = {int(gl.sec_id(t)): (r0, d0) for t, r0, d0, b, fl in rows
               if fl == 64 and "NUV" in b}
        if len(f64) < 200:
            continue
        bra, bde = boresight_mean(rows)
        for m, sra, sdec in cands:
            if gl.ang_arcmin(bra, bde, sra, sdec) > 25.0:
                continue
            usable = {s: "NUV" for s in f64}
            keep_ra, keep_de = [], []
            for v0, v1 in frags:
                raw = dev.client.photons_box(
                    "NUV", sra, sdec, 11.0 / 3600.0, v0, v1 + 1000,
                    dev.store)
                for t, r, d, fl in raw:
                    if fl == 0 and int(gl.sec_id(t)) in f64:
                        keep_ra.append(r)
                        keep_de.append(d)
            if len(keep_ra) < 30:
                continue
            pra, pde = np.array(keep_ra), np.array(keep_de)
            cosd = np.cos(np.radians(sdec))
            off = float(np.hypot((pra.mean() - sra) * cosd,
                                 pde.mean() - sdec) * 3600.0)
            rms = float(np.sqrt(np.mean(
                ((pra - pra.mean()) * cosd) ** 2
                + (pde - pde.mean()) ** 2)) * 3600.0)
            results.append({
                "block_utc": galex_ms_to_iso(frags[0][0]),
                "nuv_mag": m, "n_photons": len(keep_ra),
                "centroid_offset_arcsec": off, "rms_arcsec": rms,
                "n_flag64_s": len(f64)})
            print(f"[i] {galex_ms_to_iso(frags[0][0])[:16]} NUV "
                  f"{m:.2f}: {len(keep_ra)} ph, offset {off:.2f}\", "
                  f"rms {rms:.2f}\"", flush=True)
        if len(results) >= 10:
            break
    offs = [r["centroid_offset_arcsec"] for r in results]
    gate = bool(len(offs) >= 3 and np.median(offs) <= 2.0)
    merge("i_flag64_astrometry", {
        "sources": results,
        "median_offset_arcsec": float(np.median(offs)) if offs else None,
        "gate_median_le_2_arcsec": gate})
    print(f"[i] GATE {'PASS' if gate else 'FAIL'} over {len(offs)} "
          f"source-blocks")


# --------------------------------------------------------------- ii --

def step_ii(dev):
    """Live-time/aperture calibration + blank-sky background."""
    field = "gj-1276_antipode"
    rng = np.random.default_rng(SEED)
    mcat = dev.mcat(field)
    raw = [(float(n) if n is not None and float(n) > -99 else 99.0,
            float(r), float(d)) for _, n, f, r, d in mcat]
    # the MCAT carries near-duplicate per-visit rows of the same star:
    # cluster within 3 arcsec, keep the brightest row per star
    raw.sort()
    allpos = []
    for m, r, d in raw:
        cosd = np.cos(np.radians(d))
        if any(np.hypot((r - r2) * cosd, d - d2) * 3600.0 < 3.0
               for m2, r2, d2 in allpos):
            continue
        allpos.append((m, r, d))
    iso = []
    for m, r, d in allpos:
        if not (16.5 <= m <= 20.5):
            continue
        cosd = np.cos(np.radians(d))
        ok = all(np.hypot((r - r2) * cosd, d - d2) * 3600.0 >= 12.0
                 or m2 >= m + 1.0
                 or (abs(r - r2) + abs(d - d2)) < 1e-7
                 for m2, r2, d2 in allpos)
        if ok:
            iso.append((m, r, d))
    allpos = [(r, d, m) for m, r, d in allpos]
    print(f"[ii] {len(allpos)} deduped MCAT stars, {len(iso)} "
          f"isolated calibrator candidates", flush=True)
    cal_rows, bg_rates = [], []
    for frags in dev.blocks(field)[:4]:
        rows = dev.aspect_frags(frags)
        bra, bde = boresight_mean(rows)
        n_bg, tries = 0, 0
        while n_bg < 8 and tries < 200:
            tries += 1
            th = rng.uniform(0, 2 * np.pi)
            rr = rng.uniform(5, 20) / 60.0
            bra2 = bra + rr * np.cos(th) / np.cos(np.radians(bde))
            bde2 = bde + rr * np.sin(th)
            cosd = np.cos(np.radians(bde2))
            if any(np.hypot((bra2 - r2) * cosd, bde2 - d2) * 3600 < 20
                   for r2, d2, m2 in allpos if m2 < 22):
                continue
            usable = dev.usable_at(frags, bra2, bde2)
            live = sum(1 for b in usable.values() if "NUV" in b)
            if live < 200:
                continue
            ph = dev.photons_at(frags, "NUV", bra2, bde2, usable)
            bg_rates.append(len(ph) / live)
            n_bg += 1
        bg = float(np.median(bg_rates)) if bg_rates else np.nan
        for m, r, d in iso:
            if gl.ang_arcmin(bra, bde, r, d) > 25.0:
                continue
            usable = dev.usable_at(frags, r, d)
            live = sum(1 for b in usable.values() if "NUV" in b)
            if live < 200:
                continue
            ph = dev.photons_at(frags, "NUV", r, d, usable)
            net = len(ph) / live - bg
            if net <= 0:
                continue
            cal_rows.append({"block_utc": galex_ms_to_iso(frags[0][0]),
                             "nuv_mag": m, "live_s": live,
                             "rate": len(ph) / live, "net": net,
                             "zp_eff": float(m + 2.5 * np.log10(net))})
        print(f"[ii] {galex_ms_to_iso(frags[0][0])[:16]}: bg "
              f"{bg:.3f}, {len(cal_rows)} cal rows", flush=True)
    zps = np.array([c["zp_eff"] for c in cal_rows])
    med = float(np.median(zps)) if len(zps) else np.nan
    mad = (float(1.4826 * np.median(np.abs(zps - med)))
           if len(zps) else np.nan)
    gate = bool(len(zps) >= 8 and mad <= 0.2)
    merge("ii_calibration", {
        "n_calibrators": len(cal_rows), "zp_eff_median": med,
        "zp_eff_robust_scatter": mad,
        "bg_rate_cts_s_aperture_median":
            float(np.median(bg_rates)) if bg_rates else None,
        "nominal_zp": 20.08,
        "aperture_deadtime_term_mag":
            round(20.08 - med, 3) if len(zps) else None,
        "rows": cal_rows, "gate_scatter_le_0p2": gate})
    print(f"[ii] GATE {'PASS' if gate else 'FAIL'}: ZP_eff {med:.3f} "
          f"± {mad:.3f} ({len(zps)} stars)")


# -------------------------------------------------------------- iii --

def b_controls(dev, field, frags, locus):
    """The frozen 8-pseudo-position rule; returns positions and the
    rotation counts."""
    rows = dev.aspect_frags(frags)
    bra, bde = boresight_mean(rows)
    bright = [(float(r), float(d)) for _, n, f, r, d in dev.mcat(field)
              if n is not None and -99 < float(n) < 21.0]
    r_loc = gl.ang_arcmin(bra, bde, locus[0], locus[1]) / 60.0
    th0 = np.degrees(np.arctan2(
        locus[1] - bde, (locus[0] - bra) * np.cos(np.radians(bde))))
    ctrls, rotations = [], []
    for k in range(1, 9):
        th = th0 + 40.0 * k
        nrot = 0
        while nrot <= 24:
            a = np.radians(th)
            cra = bra + r_loc * np.cos(a) / np.cos(np.radians(bde))
            cde = bde + r_loc * np.sin(a)
            cosd = np.cos(np.radians(cde))
            bad = (any(np.hypot((cra - r2) * cosd, cde - d2) * 3600
                       < 30 for r2, d2 in bright)
                   or gl.ang_arcmin(cra, cde, locus[0], locus[1]) * 60
                   < 60)
            if not bad:
                break
            th += 5.0
            nrot += 1
        ctrls.append((cra, cde))
        rotations.append(nrot)
    return ctrls, rotations


def step_iii(dev):
    out = {}
    field = "gj-1276_antipode"
    b_rows = []
    for frags in dev.blocks(field)[:3]:
        locus = dev.field_pos(field, (frags[0][0] + frags[-1][1]) / 2)
        ctrls, rot = b_controls(dev, field, frags, locus)
        b_rows.append({"block_utc": galex_ms_to_iso(frags[0][0]),
                       "rotations": rot,
                       "n_valid": sum(1 for n in rot if n <= 24)})
        print(f"[iii] B {galex_ms_to_iso(frags[0][0])[:16]}: "
              f"rotations {rot}", flush=True)
    out["B_pseudo_positions"] = b_rows
    a_rows = {}
    for field in ("gj-1276_star", "wolf-359_star", "ross-128_star"):
        need = UNIT_LIVE_S[field]
        nseg, per_visit = 0, []
        for v0, v1 in dev.visits(field):
            pos = dev.field_pos(field, (v0 + v1) / 2)
            usable = dev.usable_at([(v0, v1)], pos[0], pos[1])
            live = sum(1 for b in usable.values() if "NUV" in b)
            k = live // need
            nseg += k
            if k:
                per_visit.append([galex_ms_to_iso(v0), live, int(k)])
        a_rows[field] = {"segment_s": need, "n_segments": int(nseg),
                         "visits": per_visit,
                         "gate_ge_8": bool(nseg >= 8)}
        print(f"[iii] A {field}: {nseg} x {need}s segments "
              f"({'PASS' if nseg >= 8 else 'FAIL'})", flush=True)
    out["A_segments"] = a_rows
    merge("iii_controls", out)


# --------------------------------------------------------------- iv --

def step_iv(dev):
    out = {}
    for field in ("gj-1276_star", "wolf-359_star", "ross-128_star"):
        rows, rates = [], []
        for v0, v1 in dev.visits(field):
            pos = dev.field_pos(field, (v0 + v1) / 2)
            usable = dev.usable_at([(v0, v1)], pos[0], pos[1])
            for band in ("NUV", "FUV"):
                live = sum(1 for b in usable.values() if band in b)
                if live < 60:
                    continue
                ph = dev.photons_at([(v0, v1)], band, pos[0], pos[1],
                                    usable)
                r = len(ph) / live
                sb = gl.s_burst(ph, live)
                rows.append({"visit_utc": galex_ms_to_iso(v0),
                             "band": band, "live_s": live,
                             "n_photons": len(ph), "rate": float(r),
                             "s_burst": float(sb)})
                if band == "NUV":
                    rates.append(r)
        med = float(np.median(rates)) if rates else np.nan
        flares = [r for r in rows
                  if r["s_burst"] > 6.0
                  or (r["band"] == "NUV" and med > 0
                      and r["rate"] > 3.0 * med)]
        out[field] = {"visits": rows, "nuv_rate_median": med,
                      "flare_candidates": flares}
        print(f"[iv] {field}: {len(rows)} visit-bands, median NUV "
              f"{med:.3f} cts/s, {len(flares)} flare candidates",
              flush=True)
    merge("iv_flare_census", out)


# ---------------------------------------------------------------- v --

def step_v(dev):
    rng = np.random.default_rng(SEED)
    out = {}
    field = "gj-1276_antipode"
    # the contiguous long flag-64 segment (2009-03-08)
    v0, v1 = max(dev.visits(field), key=lambda v: v[1] - v[0])
    frags = [(v0, v1)]
    locus = dev.field_pos(field, (v0 + v1) / 2)
    usable = dev.usable_at(frags, locus[0], locus[1])
    live = sum(1 for b in usable.values() if "NUV" in b)
    stamps = np.array(sorted(s * 1000 + 995 for s, b in usable.items()
                             if "NUV" in b), dtype=np.int64)
    span = (stamps[-1] + 1000 - stamps[0]) / 1000.0
    ph_u = dev.photons_at(frags, "NUV", locus[0], locus[1], usable)
    ctrls, rot = b_controls(dev, field, frags, locus)
    stats_c = {"S_rate": [], "S_burst": [], "S_period": []}
    for cra, cde in ctrls:
        cu = dev.usable_at(frags, cra, cde)
        cl = sum(1 for b in cu.values() if "NUV" in b)
        cp = dev.photons_at(frags, "NUV", cra, cde, cu)
        stats_c["S_rate"].append(gl.s_rate(cp, cl))
        stats_c["S_burst"].append(gl.s_burst(cp, cl))
        sp = gl.s_period(cp, span)
        stats_c["S_period"].append(sp if sp is not None else 0.0)
        print(f"[v] B control: {len(cp)} ph / {cl}s", flush=True)
    T = {k: max(v) for k, v in stats_c.items()}
    su = {"S_rate": gl.s_rate(ph_u, live),
          "S_burst": gl.s_burst(ph_u, live),
          "S_period": gl.s_period(ph_u, span) or 0.0}
    out["B_pseudo_unit"] = {
        "visit_utc": galex_ms_to_iso(v0), "live_s": live,
        "n_photons": len(ph_u), "unit": su, "thresholds": T,
        "control_rotations": rot, "controls": stats_c,
        "exceedances": {k: bool(su[k] > max(T[k], 0)) for k in su}}
    print(f"[v] B pseudo-unit: {su} vs T {T}", flush=True)

    def recovered(times, stat, thr):
        m = gl.merge_series(ph_u, times)
        if stat == "S_rate":
            return gl.s_rate(m, live) > max(thr, 0)
        if stat == "S_burst":
            return gl.s_burst(m, live) > max(thr, 0)
        return (gl.s_period(m, span) or 0.0) > max(thr, 0)

    inj = {"persistent": [], "pulse": [], "train": []}
    for rate in (0.02, 0.05, 0.1, 0.2, 0.5):
        n_rec = sum(recovered(gl.inject_persistent(rng, rate, stamps),
                              "S_rate", T["S_rate"])
                    for _ in range(40))
        inj["persistent"].append({"rate_cts_s": rate,
                                  "recovered": n_rec, "of": 40})
        print(f"[v] inj persistent {rate}: {n_rec}/40", flush=True)
    for n_ph in (3, 5, 8, 12, 20):
        n_rec = sum(recovered(gl.inject_pulse(rng, n_ph, 0.5, stamps),
                              "S_burst", T["S_burst"])
                    for _ in range(40))
        inj["pulse"].append({"n_photons": n_ph, "width_s": 0.5,
                             "recovered": n_rec, "of": 40})
        print(f"[v] inj pulse {n_ph}ph: {n_rec}/40", flush=True)
    for period in (0.1, 1.0, 10.0):
        for frac in (0.05, 0.15, 0.4):
            for drift in (0.0, gl.ORBIT_DRIFT_RATE):
                n_rec = sum(recovered(
                    gl.inject_train(rng, period, frac, 0.1, stamps,
                                    drift_rate=drift or None),
                    "S_period", T["S_period"]) for _ in range(6))
                inj["train"].append({"period_s": period,
                                     "photons_per_cycle": frac,
                                     "drift": drift,
                                     "recovered": n_rec, "of": 6})
                print(f"[v] inj train P={period} f={frac} "
                      f"d={drift:g}: {n_rec}/6", flush=True)
    out["B_injections"] = inj

    # A-short pseudo-unit: wolf-359 star segment machinery
    field = "wolf-359_star"
    need = UNIT_LIVE_S[field]
    v0, v1 = max(dev.visits(field), key=lambda v: v[1] - v[0])
    pos = dev.field_pos(field, (v0 + v1) / 2)
    usable = dev.usable_at([(v0, v1)], pos[0], pos[1])
    secs = np.array(sorted(s for s, b in usable.items()
                           if "NUV" in b), dtype=np.int64)
    nseg = len(secs) // need
    segs = [secs[i * need:(i + 1) * need] for i in range(nseg)]
    rng.shuffle(segs)
    unit_seg, ctrl_segs = segs[0], segs[1:9]
    ph_all = dev.photons_at([(v0, v1)], "NUV", pos[0], pos[1], usable)

    def seg_stats(seg):
        ids = set(int(s) for s in seg)
        ph = np.array([t for t in ph_all
                       if int(gl.sec_id(t)) in ids], dtype=np.int64)
        sp = gl.s_period(ph, float(need))
        return {"S_rate": gl.s_rate(ph, need),
                "S_burst": gl.s_burst(ph, need),
                "S_period": sp if sp is not None else 0.0,
                "n": len(ph)}

    su = seg_stats(unit_seg)
    sc = [seg_stats(s) for s in ctrl_segs]
    Ta = {k: max(c[k] for c in sc)
          for k in ("S_rate", "S_burst", "S_period")}
    out["A_pseudo_unit"] = {
        "visit_utc": galex_ms_to_iso(v0), "segment_s": need,
        "n_segments_available": nseg, "unit": su,
        "thresholds": Ta, "controls": sc,
        "exceedances": {k: bool(su[k] > max(Ta[k], 0)) for k in Ta}}
    print(f"[v] A pseudo-unit ({su['n']} ph/seg): {su} vs {Ta}",
          flush=True)
    merge("v_pseudo_units", out)


STEPS = {"i": step_i, "ii": step_ii, "iii": step_iii, "iv": step_iv,
         "v": step_v}

if __name__ == "__main__":
    dev = Dev()
    for arg in (sys.argv[1:] or ["i", "ii", "iii", "iv", "v"]):
        t0 = _time.monotonic()
        STEPS[arg](dev)
        print(f"== step {arg} done in {_time.monotonic() - t0:.0f}s ==",
              flush=True)
