"""Dev-stage search (threshold freeze v1.0, dev split only).

Port of surveys/ps1-crossings/scripts/dev_search.py to the PTF
scie-direct substrate. Stages: --plan (enumerate dev-unit epochs;
dev_plan_v1.json), --fetch (resumable scie/msk cutout downloads +
one MAST DR2 mean-star calibrator cone per dev position), --stats
(star-calibrated matched-filter photometry, frozen statistics,
dev_search_v1.json), --satgate (frozen saturation verification:
bright-star cutouts, dmask bit-8/bleed core test, incl. the ross-128
excluded-class exercise; saturation_verify_v1.json), --asteroid
(positive control: a SkyBoT-selected numbered asteroid through the
identical chain; asteroid_control_v1.json).

Dev scope (freeze D8): gj-1276 A 0.1 AU R (the one searched dev
unit, S_event, 1 trial), wolf-359 B 0.1 AU R (single-epoch class,
machinery only), ross-128 A R (saturation-excluded class, exercised
in --satgate).

Concrete realizations fixed at dev (the freeze left them
implementation-defined; validated here and reported):
  * calibrators: one MAST DR2 `mean` cone (r = 0.10 deg) per dev
    position, snapshotted; predicted PTF mags through the frozen
    transform (g = gMeanPSFMag; R = r - 0.153(r-i) - 0.117);
    calibrator quality nDetections >= 5, finite mags, err <= 0.1,
    predicted mag in [15.5, 19.5] (R) / [16.0, 20.0] (g) and PS1
    color inside the transform's dwarf locus (r-i in [0.0, 0.8] for
    R; g-r in [0.2, 1.2] for g) — the unrestricted set inflated the
    per-frame scatter to ~0.25 via off-locus red stars (measured at
    dev; restriction drops a test frame to 0.047), and all cutouts
    are 256 px so the restricted density still meets the >= 5 gate;
  * per-frame ZP = median over calibrators with good_frac >= 0.7 and
    S/N >= 5 of (m_pred + 2.5 log10 flux); gate >= 5 stars and robust
    scatter <= 0.2 mag (freeze D3), else the frame is unusable;
    fluxes rescaled to ZP_REF 25 AB;
  * channel-A star position per epoch: linear interpolation in time
    of the per-event star_icrs positions from the universal events
    table (PM propagation between tabulated events);
  * channel-A baseline for a given (real or pseudo) window = usable
    epochs outside the real window AND outside the evaluated window;
    k = median(r^2/v)/0.4549 (3x3sigma clip, floor 1, n >= 6) from
    the baseline residuals about their median;
  * exposure duplicates cannot occur (one row per CCD exposure,
    asserted at coverage); no dedup bin needed.
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import requests
from astropy.table import Table

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from sglsurvey.adapters.base import ConeRegion, CutoutSpec, MjdRange
from sglsurvey.adapters.irsa_ptf import (PtfLevel1Adapter,
                                         MASK_FATAL_TEMPLATE)
from sglsurvey.adapters.mast_ps1 import MEAN_COLS, catalog_cone
from sglsurvey.photometry import build_flux_map_ptf
from sglsurvey.snapshots import SnapshotStore

from coverage_intersect import load_events, relay_apparent, window  # noqa: E402

D = REPO / "surveys" / "ptf-crossings"
FREEZE = json.loads((D / "configs" / "threshold_freeze_v1.json").read_text())
RUN = REPO / "runs" / "ptf-crossings"
PRODUCTS = RUN / "products" / "dev"
ERA = MjdRange(54891.0, 57051.0)
Z_GRID = (550.0, 1000.0, 2500.0, 5500.0, 10000.0)
RING = [(20.0, 0.0), (-20.0, 0.0), (30.0, 0.0), (-30.0, 0.0),
        (40.0, 0.0), (-40.0, 0.0), (0.0, 30.0), (0.0, -30.0)]
OFFSETS_D = [-97, -71, -47, -23, 23, 47, 71, 97]
REDRAWS_D = [-127, -113, 113, 127]
CUT_B_IN, CUT_B_OFF, CUT_A, CUT_SAT = 256, 256, 256, 96
WEIGHT_CAP = 20.0
ZP_REF = 25.0
GOOD_FRAC_MIN = 0.7
CAL_CONE_DEG = 0.10
CAL_MIN, CAL_SCATTER_MAX = 5, 0.2
SAT_EST = {"R": 14.0, "g": 14.5}       # frozen point-source estimates

DEV = {
    "A_search": {"target_id": "gj-1276", "band": "R",
                 "event_id": "evt-4b29abd25d25"},
    "B_single": {"target_id": "wolf-359", "band": "R",
                 "event_id": "evt-86fd5e8fad87"},
    "A_excluded": {"target_id": "ross-128", "band": "R"},
}

_events = None


def events():
    global _events
    if _events is None:
        _events = load_events()
    return _events


def star_track(tid):
    """Linear (RA, Dec) vs MJD from the per-event star positions."""
    sub = events()["A"]
    sub = sub[np.asarray([str(x) == tid for x in sub["target_id"]])]
    t = np.asarray(sub["t_ca_mjd"], float)
    ra = np.asarray(sub["star_icrs_ra_deg"], float)
    de = np.asarray(sub["star_icrs_dec_deg"], float)
    o = np.argsort(t)
    t, ra, de = t[o], ra[o], de[o]

    def at(mjd):
        return (float(np.interp(mjd, t, ra)),
                float(np.interp(mjd, t, de)))
    return at


def event_row(ch, tid, eid):
    sub = events()[ch]
    m = (np.asarray([str(x) == tid for x in sub["target_id"]])
         & np.asarray([str(x) == eid for x in sub["event_id"]]))
    return sub[m][0]


def nkey(o):
    return json.dumps(o.native_key, sort_keys=True)


_obs_cache = {}


def discover_at(adapter, store, name, ra, dec, radius=0.03):
    if name not in _obs_cache:
        _obs_cache[name] = list(adapter.discover(
            ConeRegion(ra, dec, radius), ERA, store))
    return _obs_cache[name]


def _dest(group, key):
    kd = json.loads(key)
    return PRODUCTS / group / f"e{kd['expid']}_c{kd['ccdid']}"


def plan():
    adapter = PtfLevel1Adapter()
    store = SnapshotStore(RUN)
    out = {}

    # -- B single-epoch unit: wolf-359 antipode
    u = DEV["B_single"]
    ev = event_row("B", u["target_id"], u["event_id"])
    lo, hi = window(ev, 0.1)
    anti = (float(ev["relay_icrs_ra_deg"]), float(ev["relay_icrs_dec_deg"]))
    obs = discover_at(adapter, store, "wolf359-anti", *anti)
    ins, offs = [], []
    for o in obs:
        if o.band != u["band"]:
            continue
        fp = adapter.nominal_footprint(o)
        t = o.t_mid_mjd_utc
        if lo <= t <= hi:
            pos = {z: relay_apparent(ev, z, t) for z in Z_GRID}
            if any(fp.contains(r_, d_) for r_, d_ in pos.values()):
                ins.append({"key": nkey(o), "t": t,
                            "pos": {str(z): list(p)
                                    for z, p in pos.items()}})
        elif fp.contains(*anti):
            offs.append({"key": nkey(o), "t": t})
    out["B"] = {**u, "t_ca_mjd": float(ev["t_ca_mjd"]),
                "window": [lo, hi], "anti": list(anti),
                "n_in": len(ins), "n_off": len(offs),
                "in_epochs": ins, "off_epochs": offs}

    # -- A searched unit: gj-1276 star
    u = DEV["A_search"]
    ev = event_row("A", u["target_id"], u["event_id"])
    lo, hi = window(ev, 0.1)
    at = star_track(u["target_id"])
    ra0, de0 = at(float(ev["t_ca_mjd"]))
    obs = discover_at(adapter, store, "gj1276-star", ra0, de0)
    eps = []
    for o in obs:
        if o.band != u["band"]:
            continue
        t = o.t_mid_mjd_utc
        ra, dec = at(t)
        if adapter.nominal_footprint(o).contains(ra, dec):
            eps.append({"key": nkey(o), "t": t, "pos": [ra, dec],
                        "in_window": bool(lo <= t <= hi)})
    out["A"] = {**u, "t_ca_mjd": float(ev["t_ca_mjd"]),
                "window": [lo, hi],
                "n_epochs": len(eps),
                "n_in": sum(1 for e in eps if e["in_window"]),
                "epochs": eps}
    (D / "configs" / "dev_plan_v1.json").write_text(
        json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: {kk: v[kk] for kk in v
                          if kk not in ("in_epochs", "off_epochs",
                                        "epochs")}
                      for k, v in out.items()}, indent=2))


def _load_plan():
    return json.loads((D / "configs" / "dev_plan_v1.json").read_text())


def fetch():
    adapter = PtfLevel1Adapter()
    store = SnapshotStore(RUN)
    p = _load_plan()
    jobs = []
    omap_b = {nkey(o): o for o in discover_at(
        adapter, store, "wolf359-anti", *p["B"]["anti"])}
    for e in p["B"]["in_epochs"]:
        ra, dec = e["pos"]["2500.0"]
        jobs.append(("B", e["key"], omap_b[e["key"]], ra, dec, CUT_B_IN))
    for e in p["B"]["off_epochs"]:
        jobs.append(("B", e["key"], omap_b[e["key"]],
                     p["B"]["anti"][0], p["B"]["anti"][1], CUT_B_OFF))
    at = star_track(p["A"]["target_id"])
    ra0, de0 = at(p["A"]["t_ca_mjd"])
    omap_a = {nkey(o): o for o in discover_at(
        adapter, store, "gj1276-star", ra0, de0)}
    for e in p["A"]["epochs"]:
        jobs.append(("A", e["key"], omap_a[e["key"]],
                     e["pos"][0], e["pos"][1], CUT_A))
    print(f"{len(jobs)} epoch-cutouts to fetch", flush=True)
    ok = fail = 0
    for i, (grp, key, o, ra, dec, size) in enumerate(jobs):
        dest = _dest(grp, key)
        cut = CutoutSpec(ra_deg=ra, dec_deg=dec, size_pix=size)
        try:
            adapter.fetch(o, ("sci", "msk"), dest, cutout=cut)
            ok += 1
        except FileNotFoundError:
            print(f"  404 {dest.name}", flush=True)
            fail += 1
        except Exception as exc:
            print(f"  ! {dest.name}: {exc}", flush=True)
            fail += 1
        if (i + 1) % 25 == 0:
            print(f"  {i+1}/{len(jobs)} (ok={ok} fail={fail})", flush=True)
    print(f"cutouts done: {ok}/{len(jobs)} ok, {fail} failed", flush=True)

    # calibrator cones (one per dev position, incl. ross-128 star for
    # the satgate)
    r128 = star_track("ross-128")(56368.5)
    for name, (ra, dec) in {
            "wolf359-anti": tuple(p["B"]["anti"]),
            "gj1276-star": (ra0, de0),
            "ross128-star": r128}.items():
        dest = RUN / "calibrators" / f"{name}.json"
        if dest.exists():
            continue
        rows, snaps = catalog_cone("mean", ra, dec, CAL_CONE_DEG,
                                   MEAN_COLS, store)
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(json.dumps(
            {"center": [ra, dec], "radius_deg": CAL_CONE_DEG,
             "snapshots": snaps, "rows": rows}) + "\n")
        print(f"calibrators {name}: {len(rows)} mean objects", flush=True)


def mould_r(r, i):
    return r - 0.153 * (r - i) - 0.117


def load_calibrators(name, band):
    doc = json.loads((RUN / "calibrators" / f"{name}.json").read_text())
    return calibrators_from_rows(doc["rows"], band)


def calibrators_from_rows(rows, band):
    lo, hi = (15.5, 19.5) if band == "R" else (16.0, 20.0)
    out = []
    for r in rows:
        try:
            nd = int(float(r["nDetections"]))
            rm, im = float(r["rMeanPSFMag"]), float(r["iMeanPSFMag"])
            gm = float(r["gMeanPSFMag"])
            if band == "R":
                re_, ie = (float(r["rMeanPSFMagErr"]),
                           float(r["iMeanPSFMagErr"]))
                if min(rm, im) < -100 or max(re_, ie) > 0.1:
                    continue
                if not (0.0 <= rm - im <= 0.8):
                    continue
                m = mould_r(rm, im)
            else:
                if (min(gm, rm) < -100
                        or float(r["gMeanPSFMagErr"]) > 0.1):
                    continue
                if not (0.2 <= gm - rm <= 1.2):
                    continue
                m = gm
        except (KeyError, TypeError, ValueError):
            continue
        if nd >= 5 and lo <= m <= hi:
            out.append((float(r["raMean"]), float(r["decMean"]), m))
    return out


_fm_cache = {}


def fluxmap(grp, key, o, band, cal_name, size):
    ck = (grp, key, size)
    if ck in _fm_cache:
        return _fm_cache[ck]
    dest = _dest(grp, key)
    res = None
    if dest.exists():
        pref = f"cut{size}-"
        sci = next((f for f in dest.iterdir()
                    if f.name.startswith(pref) and "_scie_" in f.name
                    and f.name.endswith(".fits")), None)
        msk = next((f for f in dest.iterdir()
                    if f.name.startswith(pref) and "_mask_" in f.name),
                   None)
        if sci is not None and msk is not None:
            try:
                fm = build_flux_map_ptf(sci, msk, MASK_FATAL_TEMPLATE,
                                        band, o.t_mid_mjd_utc)
                zps = []
                for cra, cde, cm in load_calibrators(cal_name, band):
                    f, v, g = fm.sample(cra, cde)
                    if (np.isfinite(g[0]) and g[0] >= GOOD_FRAC_MIN
                            and np.isfinite(f[0]) and np.isfinite(v[0])
                            and v[0] > 0 and f[0] > 5 * np.sqrt(v[0])):
                        zps.append(cm + 2.5 * np.log10(f[0]))
                if len(zps) >= CAL_MIN:
                    zps = np.array(zps)
                    zp = float(np.median(zps))
                    sc = 1.4826 * float(np.median(np.abs(zps - zp)))
                    if sc <= CAL_SCATTER_MAX:
                        s = 10 ** (0.4 * (ZP_REF - zp))
                        fm.flux = fm.flux * s
                        fm.var = fm.var * s * s
                        fm.zp_star = zp
                        fm.zp_scatter = sc
                        fm.n_cal = len(zps)
                        res = fm
                    else:
                        res = ("gate_scatter", sc, len(zps))
                else:
                    res = ("gate_ncal", None, len(zps))
            except Exception as exc:
                print(f"  fluxmap fail {dest.name}: {exc}", flush=True)
                res = ("error", None, None)
    _fm_cache[ck] = res
    return res


def _wstack(res):
    """res: list of (residual, var). Weighted one-sided stack."""
    r = np.array([x[0] for x in res], float)
    v = np.array([x[1] for x in res], float)
    ok = np.isfinite(r) & np.isfinite(v) & (v > 0)
    if not ok.any():
        return np.nan
    r, v = r[ok], v[ok]
    w = 1.0 / v
    w = np.minimum(w, WEIGHT_CAP * np.median(w))
    return float(np.sum(w * r) / np.sqrt(np.sum(w)))


def stats():
    adapter = PtfLevel1Adapter()
    store = SnapshotStore(RUN)
    p = _load_plan()
    results = {}

    # ---------------- channel B single-epoch unit -------------------
    b = p["B"]
    ev = event_row("B", b["target_id"], b["event_id"])
    omap = {nkey(o): o for o in discover_at(
        adapter, store, "wolf359-anti", *b["anti"])}
    gate = defaultdict(int)
    off = []
    for e in b["off_epochs"]:
        fm = fluxmap("B", e["key"], omap[e["key"]], b["band"],
                     "wolf359-anti", CUT_B_OFF)
        if not hasattr(fm, "sample"):
            gate[fm[0] if fm else "missing"] += 1
            continue
        f, v, g = fm.sample(b["anti"][0], b["anti"][1])
        if np.isfinite(g[0]) and g[0] >= GOOD_FRAC_MIN:
            off.append((f[0], v[0]))
    k = 1.0
    k_status = f"undefined_n{len(off)}"
    if len(off) >= 6:
        fo = np.array([f for f, _ in off])
        vo = np.array([v for _, v in off])
        for _ in range(3):
            med = np.median(fo)
            sd = 1.4826 * np.median(np.abs(fo - med)) or 1.0
            m = np.abs(fo - med) <= 3 * sd
            fo, vo = fo[m], vo[m]
        if len(fo) >= 6:
            med = np.median(fo)
            k = max(1.0, float(np.median((fo - med) ** 2 / vo) / 0.4549))
            k_status = f"ok_n{len(fo)}"
    infms = []
    for e in b["in_epochs"]:
        fm = fluxmap("B", e["key"], omap[e["key"]], b["band"],
                     "wolf359-anti", CUT_B_IN)
        if hasattr(fm, "sample"):
            infms.append((omap[e["key"]], fm))
        else:
            gate[fm[0] if fm else "missing"] += 1

    def track_S(dx, dy):
        best = -np.inf
        for z in Z_GRID:
            sm = []
            for o, fm in infms:
                ra, dec = relay_apparent(ev, z, o.t_mid_mjd_utc)
                ra += dx / 3600.0 / max(np.cos(np.radians(dec)), .05)
                dec += dy / 3600.0
                f, v, g = fm.sample(ra, dec)
                if np.isfinite(g[0]) and g[0] >= GOOD_FRAC_MIN:
                    sm.append((f[0], v[0] * k))
            if sm:
                best = max(best, _wstack(sm))
        return best if np.isfinite(best) else np.nan

    S = track_S(0.0, 0.0)
    controls = [c for c in (track_S(dx, dy) for dx, dy in RING)
                if np.isfinite(c)]
    T = max(controls) if controls else np.nan

    # static annotation from ptf_objects (TAP snapshot)
    static_min = None
    try:
        url = "https://irsa.ipac.caltech.edu/TAP/sync"
        q = ("SELECT ra,dec,ngoodobs FROM ptf_objects WHERE "
             f"CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',"
             f"{b['anti'][0]:.5f},{b['anti'][1]:.5f},0.02))=1")
        from datetime import datetime, timezone
        resp = requests.get(url, params={"QUERY": q, "FORMAT": "csv"},
                            timeout=300)
        store.store(service_url=url, query=q,
                    request_utc=datetime.now(timezone.utc).isoformat(),
                    response_bytes=resp.content, row_count=None,
                    http_status=resp.status_code)
        import csv as _csv
        import io as _io
        rows = [r for r in _csv.DictReader(_io.StringIO(resp.text))
                if float(r["ngoodobs"] or 0) >= 3]
        dmin = np.inf
        for e in b["in_epochs"]:
            for z in Z_GRID:
                ra, dec = e["pos"][str(z)]
                for r in rows:
                    d = np.hypot((float(r["ra"]) - ra)
                                 * np.cos(np.radians(dec)),
                                 float(r["dec"]) - dec) * 3600.0
                    dmin = min(dmin, d)
        static_min = None if not np.isfinite(dmin) else round(dmin, 2)
    except Exception as exc:
        static_min = f"query_failed: {exc}"

    results["B"] = {
        **{k_: b[k_] for k_ in ("target_id", "band", "event_id",
                                "t_ca_mjd")},
        "class": "single_epoch", "n_trials": 0,
        "n_in_epochs": len(infms), "n_off_usable": len(off),
        "calibration_gate_failures": dict(gate),
        "k_rescale": round(k, 3), "k_status": k_status,
        "S": None if np.isnan(S) else round(S, 3),
        "n_controls": len(controls),
        "controls": [round(c, 3) for c in controls],
        "T": None if not controls else round(T, 3),
        "margin": (None if not controls or np.isnan(S)
                   else round(S - max(T, 0.0), 3)),
        "static_min_arcsec": static_min,
        "exceedance_reportable": bool(
            np.isfinite(S) and controls and S > max(T, 0.0)),
    }

    # ---------------- channel A searched unit -----------------------
    a = p["A"]
    at = star_track(a["target_id"])
    ra0, de0 = at(a["t_ca_mjd"])
    omap = {nkey(o): o for o in discover_at(
        adapter, store, "gj1276-star", ra0, de0)}
    gate = defaultdict(int)
    eps = []           # (t, flux, var) usable, ZP 25
    for e in a["epochs"]:
        fm = fluxmap("A", e["key"], omap[e["key"]], a["band"],
                     "gj1276-star", CUT_A)
        if not hasattr(fm, "sample"):
            gate[fm[0] if fm else "missing"] += 1
            continue
        f, v, g = fm.sample(e["pos"][0], e["pos"][1])
        if np.isfinite(g[0]) and g[0] >= GOOD_FRAC_MIN:
            eps.append((e["t"], f[0], v[0]))
    lo, hi = a["window"]
    half = (hi - lo) / 2.0

    def window_S(t_c):
        wlo, whi = t_c - half, t_c + half
        inw = [(f, v) for t, f, v in eps if wlo <= t <= whi]
        base = [(f, v) for t, f, v in eps
                if not (wlo <= t <= whi) and not (lo <= t <= hi)]
        if len(inw) < 2 or len(base) < 6:
            return None, len(inw), len(base)
        fb = np.array([f for f, _ in base])
        vb = np.array([v for _, v in base])
        for _ in range(3):
            med = np.median(fb)
            sd = 1.4826 * np.median(np.abs(fb - med)) or 1.0
            m = np.abs(fb - med) <= 3 * sd
            fb, vb = fb[m], vb[m]
        if len(fb) < 6:
            return None, len(inw), len(base)
        med = float(np.median(fb))
        k_ = max(1.0, float(np.median((fb - med) ** 2 / vb) / 0.4549))
        res = [(f - med, v * k_) for f, v in inw]
        return _wstack(res), len(inw), len(base)

    t_ca = a["t_ca_mjd"]
    S_A, n_in, n_base = window_S(t_ca)
    ctrl, invalid = [], []
    for off_d in OFFSETS_D + REDRAWS_D:
        if len(ctrl) >= 8:
            break
        s, n_i, _ = window_S(t_ca + off_d)
        if s is not None and n_i >= 2:
            ctrl.append((off_d, round(s, 3)))
        else:
            invalid.append(off_d)
    T_A = max(s for _, s in ctrl) if ctrl else np.nan
    a_status = ("ok" if S_A is not None and len(ctrl) >= 8
                else "constraint_only_gates")
    results["A"] = {
        **{k_: a[k_] for k_ in ("target_id", "band", "event_id",
                                "t_ca_mjd")},
        "class": "search_unit", "n_trials": 1,
        "n_epochs_usable": len(eps),
        "n_in_window_usable": n_in, "n_baseline": n_base,
        "calibration_gate_failures": dict(gate),
        "status": a_status,
        "S": None if S_A is None else round(S_A, 3),
        "n_valid_offsets": len(ctrl),
        "controls": ctrl, "invalid_offsets": invalid,
        "T": None if not ctrl else round(float(T_A), 3),
        "margin": (None if not ctrl or S_A is None
                   else round(S_A - max(T_A, 0.0), 3)),
        "exceedance": bool(S_A is not None and ctrl
                           and len(ctrl) >= 8 and S_A > max(T_A, 0.0)),
    }
    zp_stats = [(_fm_cache[k].zp_scatter, _fm_cache[k].n_cal)
                for k in _fm_cache
                if hasattr(_fm_cache[k], "zp_scatter")]
    results["calibration"] = {
        "n_frames_calibrated": len(zp_stats),
        "median_zp_scatter": (round(float(np.median(
            [s for s, _ in zp_stats])), 3) if zp_stats else None),
        "median_n_calibrators": (int(np.median([n for _, n in zp_stats]))
                                 if zp_stats else None),
    }
    (D / "results" / "dev_search_v1.json").write_text(
        json.dumps(results, indent=1) + "\n")
    print(json.dumps(results, indent=1))


def satgate():
    """Frozen saturation verification: per band, DR2 mean stars near
    the frozen boundary measured on in-band cutouts; dmask bit 8
    (saturated, 256) / bit 6 (bleed, 64) core test vs magnitude; plus
    the ross-128 excluded-class exercise (R ~ 9.9 star core)."""
    from astropy.io import fits
    from astropy.wcs import WCS as AWCS

    adapter = PtfLevel1Adapter()
    store = SnapshotStore(RUN)
    p = _load_plan()
    at = star_track("gj-1276")
    positions = {
        "wolf359-anti": tuple(p["B"]["anti"]),
        "gj1276-star": at(p["A"]["t_ca_mjd"]),
    }
    r128 = star_track("ross-128")(56368.5)
    out = {"frozen_sat_estimate": SAT_EST,
           "frozen_exclude": {"R": 14.5, "g": 15.0}, "bands": {}}
    SATBITS = 256 | 64
    for band in ("R", "g"):
        rows = []
        for name, (cra, cde) in positions.items():
            doc = json.loads((RUN / "calibrators" / f"{name}.json"
                              ).read_text())
            cand = []
            for r in doc["rows"]:
                try:
                    if band == "R":
                        m = mould_r(float(r["rMeanPSFMag"]),
                                    float(r["iMeanPSFMag"]))
                    else:
                        m = float(r["gMeanPSFMag"])
                    if m < 5 or int(float(r["nDetections"])) < 5:
                        continue
                except (KeyError, TypeError, ValueError):
                    continue
                cand.append((float(r["raMean"]), float(r["decMean"]), m))
            cand = np.array(cand) if cand else np.empty((0, 3))
            picks = []
            for m0 in (SAT_EST[band] - 0.8, SAT_EST[band],
                       SAT_EST[band] + 0.8):
                if len(cand):
                    i = int(np.argmin(np.abs(cand[:, 2] - m0)))
                    picks.append(tuple(cand[i]))
            obs = [o for o in discover_at(adapter, store, name, cra, cde)
                   if o.band == band]
            for ra, dec, m in dict.fromkeys(picks):
                o = next((o for o in obs
                          if adapter.nominal_footprint(o).contains(ra, dec)),
                         None)
                if o is None:
                    continue
                dest = PRODUCTS / "satgate" / f"{name}-{band}-{m:.2f}"
                try:
                    ps = adapter.fetch(o, ("msk",), dest,
                                       cutout=CutoutSpec(ra, dec, CUT_SAT))
                except Exception:
                    continue
                with fits.open(ps.products[0].path) as hd:
                    hdu = next(h for h in hd if h.data is not None)
                    mdat = (np.asarray(hdu.data).astype(np.int64) & 0xFFFF)
                    w = AWCS(hdu.header)
                x, y = w.all_world2pix([[ra, dec]], 0)[0]
                xi, yi = int(round(x)), int(round(y))
                r5 = 5
                core = mdat[max(0, yi - r5):yi + r5 + 1,
                            max(0, xi - r5):xi + r5 + 1]
                rows.append({"field": name, "mag_pred": round(m, 2),
                             "sat_bit8": bool((core & 256).any()),
                             "bleed_bit6": bool((core & 64).any()),
                             "sat_or_bleed": bool((core & SATBITS).any()),
                             "expid": o.native_key["expid"]})
        satm = [r["mag_pred"] for r in rows if r["sat_or_bleed"]]
        clean = [r["mag_pred"] for r in rows if not r["sat_or_bleed"]]
        out["bands"][band] = {
            "stars": rows,
            "faintest_with_sat_bits": max(satm) if satm else None,
            "brightest_clean": min(clean) if clean else None,
        }
    # ross-128 excluded-class exercise
    obs = [o for o in discover_at(adapter, store, "ross128-star", *r128)
           if o.band == "R"]
    exc = None
    for o in obs:
        if not adapter.nominal_footprint(o).contains(*r128):
            continue
        dest = PRODUCTS / "satgate" / "ross128-star-exercise"
        try:
            ps = adapter.fetch(o, ("msk",), dest,
                               cutout=CutoutSpec(r128[0], r128[1],
                                                 CUT_SAT))
        except Exception:
            continue
        with fits.open(ps.products[0].path) as hd:
            hdu = next(h for h in hd if h.data is not None)
            mdat = (np.asarray(hdu.data).astype(np.int64) & 0xFFFF)
            w = AWCS(hdu.header)
        x, y = w.all_world2pix([[r128[0], r128[1]]], 0)[0]
        xi, yi = int(round(x)), int(round(y))
        core = mdat[max(0, yi - 5):yi + 6, max(0, xi - 5):xi + 6]
        exc = {"expid": o.native_key["expid"],
               "R_est": 9.86,
               "sat_bit8": bool((core & 256).any()),
               "bleed_bit6": bool((core & 64).any()),
               "verdict": ("exclusion_confirmed"
                           if (core & SATBITS).any()
                           else "UNEXPECTED_clean_core")}
        break
    out["ross128_excluded_class_exercise"] = exc
    (D / "results" / "saturation_verify_v1.json").write_text(
        json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


def asteroid():
    """Positive control: a numbered main-belt asteroid crossing the
    wolf-359 antipode field, forced-photometered through the identical
    chain at its SkyBoT per-epoch positions, vs predicted V + (V-R) =
    0.4 +/- 0.1 (declared class-color term)."""
    from datetime import datetime, timezone

    adapter = PtfLevel1Adapter()
    store = SnapshotStore(RUN)
    p = _load_plan()
    anti = tuple(p["B"]["anti"])
    obs = [o for o in discover_at(adapter, store, "wolf359-anti", *anti)
           if o.band == "R"]
    by_night = defaultdict(list)
    for o in obs:
        by_night[int(o.t_mid_mjd_utc)].append(o)
    nights = sorted(by_night, key=lambda n: -len(by_night[n]))

    sess = requests.Session()
    SKYBOT = "https://vo.imcce.fr/webservices/skybot/skybotconesearch_query.php"

    def skybot(mjd):
        params = {"EPOCH": f"{mjd + 2400000.5:.6f}", "RA": f"{anti[0]:.5f}",
                  "DEC": f"{anti[1]:.5f}", "SR": "0.5", "-mime": "text",
                  "-loc": "675", "-filter": "120"}
        resp = sess.get(SKYBOT, params=params, timeout=120)
        store.store(service_url=SKYBOT,
                    query="&".join(f"{k}={v}" for k, v in params.items()),
                    request_utc=datetime.now(timezone.utc).isoformat(),
                    response_bytes=resp.content, row_count=None,
                    http_status=resp.status_code)
        rows = []
        for line in resp.text.splitlines():
            if line.startswith(("#", "-")) or "|" not in line:
                continue
            f = [x.strip() for x in line.split("|")]
            try:
                num = int(f[0])
                from astropy.coordinates import Angle
                ra = Angle(f[2] + " hours").degree
                dec = Angle(f[3] + " degrees").degree
                mv = float(f[5])
                err = float(f[6])
            except (ValueError, IndexError):
                continue
            if err > 1.0:
                continue
            rows.append({"num": num, "name": f[1], "ra": ra, "dec": dec,
                         "mv": mv})
        return rows

    # find one numbered asteroid V<=19 (ephemeris err <= 1 arcsec)
    # present in >= 4 exposures over the best-sampled nights
    seen = defaultdict(list)     # num -> [(o, ra, dec, mv)]
    for n in nights[:4]:
        for o in by_night[n]:
            for r in skybot(o.t_mid_mjd_utc):
                if r["mv"] > 19.0:
                    continue
                if adapter.nominal_footprint(o).contains(r["ra"], r["dec"]):
                    seen[(r["num"], r["name"])].append(
                        (o, r["ra"], r["dec"], r["mv"]))
    if not seen:
        print("no suitable asteroid found")
        return
    (num, name), hits = max(seen.items(), key=lambda kv: len(kv[1]))
    hits = hits[:8]
    # dedicated calibrator cone at the asteroid's median position (it
    # roams up to ~0.5 deg from the antipode, outside the field cone)
    mra = float(np.median([h[1] for h in hits]))
    mde = float(np.median([h[2] for h in hits]))
    cal_path = RUN / "calibrators" / f"asteroid-{num}.json"
    if not cal_path.exists():
        rows_c, snaps = catalog_cone("mean", mra, mde, 0.15,
                                     MEAN_COLS, store)
        cal_path.write_text(json.dumps(
            {"center": [mra, mde], "radius_deg": 0.15,
             "snapshots": snaps, "rows": rows_c}) + "\n")
    ast_cal = calibrators_from_rows(
        json.loads(cal_path.read_text())["rows"], "R")
    print(f"asteroid ({num}) {name}: {len(hits)} hits, "
          f"{len(ast_cal)} calibrators", flush=True)
    meas = []
    for j, (o, ra, dec, mv) in enumerate(hits):
        dest = PRODUCTS / "asteroid" / f"{num}-{j}"
        try:
            adapter.fetch(o, ("sci", "msk"), dest,
                          cutout=CutoutSpec(ra, dec, CUT_A))
        except Exception as exc:
            print(f"  fetch fail {j}: {exc}")
            continue
        pref = f"cut{CUT_A}-"
        sci = next(f for f in dest.iterdir()
                   if f.name.startswith(pref) and "_scie_" in f.name)
        msk = next(f for f in dest.iterdir()
                   if f.name.startswith(pref) and "_mask_" in f.name)
        fm = build_flux_map_ptf(sci, msk, MASK_FATAL_TEMPLATE, "R",
                                o.t_mid_mjd_utc)
        zps = []
        for cra, cde, cm in ast_cal:
            f, v, g = fm.sample(cra, cde)
            if (np.isfinite(g[0]) and g[0] >= GOOD_FRAC_MIN
                    and np.isfinite(f[0]) and v[0] > 0
                    and f[0] > 5 * np.sqrt(v[0])):
                zps.append(cm + 2.5 * np.log10(f[0]))
        if len(zps) < CAL_MIN:
            continue
        zp = float(np.median(np.array(zps)))
        f, v, g = fm.sample(ra, dec)
        if not (np.isfinite(g[0]) and g[0] >= GOOD_FRAC_MIN
                and f[0] > 0):
            continue
        m_meas = zp - 2.5 * np.log10(f[0])
        m_pred = mv - 0.4
        meas.append({"mjd": round(o.t_mid_mjd_utc, 5),
                     "expid": o.native_key["expid"],
                     "V_pred": mv, "R_pred": round(m_pred, 3),
                     "R_meas": round(float(m_meas), 3),
                     "resid": round(float(m_meas - m_pred), 3),
                     "n_cal": len(zps)})
    resid = np.array([m["resid"] for m in meas])
    out = {"asteroid": {"number": num, "name": name},
           "vr_assumed": 0.4, "vr_class_scatter": 0.1,
           "n_epochs": len(meas), "epochs": meas,
           "median_resid": (round(float(np.median(resid)), 3)
                            if len(meas) else None),
           "resid_rms": (round(float(np.std(resid)), 3)
                         if len(meas) else None)}
    (D / "results" / "asteroid_control_v1.json").write_text(
        json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    if "--plan" in sys.argv:
        plan()
    elif "--fetch" in sys.argv:
        fetch()
    elif "--stats" in sys.argv:
        stats()
    elif "--satgate" in sys.argv:
        satgate()
    elif "--asteroid" in sys.argv:
        asteroid()
    else:
        print("usage: dev_search.py --plan|--fetch|--stats|--satgate"
              "|--asteroid")
