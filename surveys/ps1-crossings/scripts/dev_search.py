"""Dev-target search (threshold freeze v1.0, dev split only).

Port of surveys/ztf-crossings/scripts/dev_search.py to the PS1
warp-direct substrate. Stages: --plan (enumerate dev channel-B units,
in-window and off-window epochs; write dev_plan_v1.json), --fetch
(resumable img/wt/msk fitscut cutout downloads), --stats
(star-calibrated matched-filter photometry, frozen statistics,
dev_search_v1.json), --satgate (the frozen saturation verification
gate: bright-star warp cutouts + CELL.SATURATION,
saturation_verify_v1.json).

Channel A has an empty dev split (all searchable targets are singleton
strata); its constructions stand as imported from ZTF v1.0+1.1+1.2 and
only the saturation gate runs at dev.

Concrete realizations fixed at dev (the freeze left them
implementation-defined; they are validated here and reported):
  * skycell duplicates (one exposure listed on two overlapping
    skycells) are collapsed at the sample level by best good_frac,
    exposures grouped at 0.0002 d (17 s; TTI spacing is >= ~600 s) --
    the ps1-v2 duplicate convention;
  * the channel-B off-window sample for the empirical variance rescale
    k = median(r^2/v)/0.4549 (3x3sigma clip, floor 1) is every same-band
    primary era epoch of the field outside the event window, sampled at
    the event's fixed tabulated antipode position (the track
    parameterization is only meaningful in-window; a fixed clean field
    point measures the same matched-filter scatter);
  * per-epoch flux scale: zp_star from the retained v1 star calibration
    (runs/panstarrs/calib_t0_59800/zeropoints.jsonl, deterministic
    observation ids), fallback FPA.ZP + 2.5 log10(EXPTIME) + per-band
    median star-vs-header offset; all fluxes rescaled to ZP 25 AB.
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.table import Table

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from sglsurvey.adapters.base import ConeRegion, CutoutSpec, MjdRange
from sglsurvey.adapters.mast_ps1 import Ps1WarpAdapter, Ps1ExactFootprint
from sglsurvey.photometry import build_flux_map_ps1
from sglsurvey.snapshots import SnapshotStore
from sglsurvey.vetting import load_ps1_mean

from coverage_intersect import load_events, relay_apparent, source_radec, window  # noqa: E402

D = REPO / "surveys" / "ps1-crossings"
FREEZE = json.loads((D / "configs" / "threshold_freeze_v1.json").read_text())
COV = Table.read(D / "results" / "coverage_v1_events.ecsv")
RUN = REPO / "runs" / "ps1-crossings"
PRODUCTS = RUN / "products" / "dev"
SKYCELL_CACHE = RUN / "skycell_wcs.json"
V1_CAL = REPO / "runs" / "panstarrs"
SCREEN = V1_CAL / "screen_v1"
ERA = MjdRange(54900.0, 57300.0)
Z_GRID = (550.0, 1000.0, 2500.0, 5500.0, 10000.0)
RING = [(20.0, 0.0), (-20.0, 0.0), (30.0, 0.0), (-30.0, 0.0),
        (40.0, 0.0), (-40.0, 0.0), (0.0, 30.0), (0.0, -30.0)]
CUT_IN, CUT_OFF, CUT_SAT = 640, 128, 96      # pix
DUP_BIN_DAYS = 0.0002
WEIGHT_CAP = 20.0
ZP_REF = 25.0
GOOD_FRAC_MIN = 0.7
STATIC_ARCSEC = 2.0
DEV_B = FREEZE["split"]["B_wide"]["dev"]
SAT_E = {"g": 14.0, "r": 14.0, "i": 14.0, "z": 13.0, "y": 12.0}

_events_cache = None


def events_b():
    global _events_cache
    if _events_cache is None:
        _events_cache = load_events()["B"]
    return _events_cache


def corridor_name(tid):
    return tid.replace("-", "")


def dev_units():
    us = [u for u in FREEZE["search_units"]["B"]["units"]
          + FREEZE["search_units"]["B"]["single_epoch"]
          if u["target_id"] in DEV_B and u["radius_au"] == 0.1]
    return us


def target_cone(tid):
    sub = events_b()[np.asarray([str(x) == tid
                                 for x in events_b()["target_id"]])]
    ras = np.array([source_radec("B", e)[0] for e in sub])
    des = np.array([source_radec("B", e)[1] for e in sub])
    ra0, de0 = float(np.median(ras)), float(np.median(des))
    cosd = max(np.cos(np.radians(de0)), 0.05)
    radius = max(np.ptp(ras) * cosd, np.ptp(des)) / 2.0 + 0.15
    return ConeRegion(ra0, de0, radius), sub


_obs_cache = {}


def discover(adapter, store, tid):
    if tid not in _obs_cache:
        cone, _ = target_cone(tid)
        _obs_cache[tid] = list(adapter.discover(cone, ERA, store))
    return _obs_cache[tid]


def okq(o):
    return o.quality_flags.get("badflag", 0) == 0


def event_row(sub, eid):
    m = np.asarray([str(x) == eid for x in sub["event_id"]])
    return sub[m][0]


def nkey(o):
    return json.dumps(o.native_key, sort_keys=True)


def plan():
    adapter = Ps1WarpAdapter(skycell_cache=SKYCELL_CACHE)
    store = SnapshotStore(RUN)
    out = {"A": [], "A_note": "dev split empty (freeze v1.0); "
                              "constructions imported from ZTF chain; "
                              "saturation gate runs via --satgate",
           "B": []}
    for u in dev_units():
        tid, band, eid = u["target_id"], u["band"], u["event_id"]
        cone, sub = target_cone(tid)
        ev = event_row(sub, eid)
        lo, hi = window(ev, 0.1)
        anti = (float(ev["relay_icrs_ra_deg"]),
                float(ev["relay_icrs_dec_deg"]))
        obs = discover(adapter, store, tid)
        ins, offs = [], []
        for o in obs:
            if o.band != band or not okq(o):
                continue
            if o.extra["skycell"] not in adapter._skycells:
                continue
            fp = adapter.nominal_footprint(o)
            t = o.t_mid_mjd_utc
            if lo <= t <= hi:
                pos = {z: relay_apparent(ev, z, t) for z in Z_GRID}
                if any(fp.contains(r_, d_) for r_, d_ in pos.values()):
                    ins.append({"key": nkey(o), "t": t,
                                "skycell": o.extra["skycell"],
                                "pos": {str(z): list(p)
                                        for z, p in pos.items()}})
            elif fp.contains(*anti):
                offs.append({"key": nkey(o), "t": t,
                             "skycell": o.extra["skycell"]})
        out["B"].append({
            "target_id": tid, "band": band, "event_id": eid,
            "t_ca_mjd": float(ev["t_ca_mjd"]),
            "window": [lo, hi], "anti": list(anti),
            "n_in": len(ins), "n_off": len(offs),
            "in_epochs": ins, "off_epochs": offs})
    (D / "configs" / "dev_plan_v1.json").write_text(
        json.dumps(out, indent=1) + "\n")
    print(json.dumps([{k: p[k] for k in
                       ("target_id", "band", "event_id", "n_in", "n_off")}
                      for p in out["B"]], indent=2))


def _load_plan():
    return json.loads((D / "configs" / "dev_plan_v1.json").read_text())


def _dest(tid, key):
    kd = json.loads(key)
    return (PRODUCTS / f"B-{tid}"
            / f"{kd['projcell']:04d}.{kd['subcell']:03d}.{kd['filter']}"
              f".{kd['mjd_tag']}")


def fetch():
    adapter = Ps1WarpAdapter(skycell_cache=SKYCELL_CACHE)
    store = SnapshotStore(RUN)
    plan_doc = _load_plan()
    jobs = {}
    for p in plan_doc["B"]:
        tid = p["target_id"]
        omap = {nkey(o): o for o in discover(adapter, store, tid)}
        for e in p["in_epochs"]:
            o = omap[e["key"]]
            ra, dec = e["pos"]["2500.0"]
            jobs.setdefault((tid, e["key"], CUT_IN), (o, ra, dec))
        for e in p["off_epochs"]:
            o = omap[e["key"]]
            ra, dec = p["anti"]
            jobs.setdefault((tid, e["key"], CUT_OFF), (o, ra, dec))
    print(f"{len(jobs)} epoch-cutouts to fetch", flush=True)
    ok = fail = 0
    for i, ((tid, key, size), (o, ra, dec)) in enumerate(sorted(jobs.items())):
        dest = _dest(tid, key)
        cut = CutoutSpec(ra_deg=ra, dec_deg=dec, size_pix=size)
        try:
            for kind in ("img", "wt", "msk"):
                adapter.fetch(o, (kind,), dest, cutout=cut)
            ok += 1
        except FileNotFoundError:
            print(f"  404 {dest.name}", flush=True)
            fail += 1
        except Exception as exc:
            print(f"  ! {dest.name}: {exc}", flush=True)
            fail += 1
        if (i + 1) % 25 == 0:
            print(f"  {i+1}/{len(jobs)} (ok={ok} fail={fail})", flush=True)
    print(f"done: {ok}/{len(jobs)} ok, {fail} failed", flush=True)


_zp_cache = None


def zeropoints():
    global _zp_cache
    if _zp_cache is None:
        zp = {}
        for pth in (V1_CAL / "calib_t0_59800" / "zeropoints.jsonl",
                    V1_CAL / "calib_v1" / "zeropoints.jsonl"):
            if pth.exists():
                for line in pth.read_text().splitlines():
                    if line.strip():
                        r = json.loads(line)
                        zp.setdefault(r["observation_id"], r)
        delta = defaultdict(list)
        for r in zp.values():
            if (r.get("zp_star") is not None and np.isfinite(r["zp_star"])
                    and r.get("zp_hdr") and r.get("exptime")):
                delta[r["band"]].append(
                    r["zp_star"] - r["zp_hdr"]
                    - 2.5 * np.log10(r["exptime"]))
        _zp_cache = (zp, {b: float(np.median(v))
                          for b, v in delta.items()})
    return _zp_cache


def _zp_eff(o, fm):
    zp, delta = zeropoints()
    r = zp.get(o.observation_id)
    if r and r.get("zp_star") is not None and np.isfinite(r["zp_star"]):
        return float(r["zp_star"]), "star"
    hdr = fm.header
    z = hdr.get("FPA.ZP", hdr.get("HIERARCH FPA.ZP"))
    ext = float(getattr(fm, "exptime", 0.0) or 0.0)
    if z is None or ext <= 0:
        return None, "none"
    return (float(z) + 2.5 * np.log10(ext)
            + delta.get(o.band, -0.5)), "fallback"


_fm_cache = {}


def _fluxmap(tid, o, size):
    ck = (tid, o.observation_id, size)
    if ck in _fm_cache:
        return _fm_cache[ck]
    dest = _dest(tid, nkey(o))
    res = None
    if dest.exists():
        pref = f"cut{size}-"
        img = next((p for p in dest.iterdir()
                    if p.name.startswith(pref)
                    and p.name.endswith(".fits")
                    and ".wt." not in p.name and ".mask." not in p.name),
                   None)
        wt = next((p for p in dest.iterdir()
                   if p.name.startswith(pref) and ".wt." in p.name), None)
        msk = next((p for p in dest.iterdir()
                    if p.name.startswith(pref) and ".mask." in p.name),
                   None)
        if img is not None and msk is not None:
            try:
                fm = build_flux_map_ps1(img, wt, msk,
                                        Ps1ExactFootprint.FATAL_MASK,
                                        o.band, o.t_mid_mjd_utc)
                zpv, src = _zp_eff(o, fm)
                if zpv is not None:
                    s = 10 ** (0.4 * (ZP_REF - zpv))
                    fm.flux = fm.flux * s
                    fm.var = fm.var * s * s
                    fm.zp_source = src
                    res = fm
            except Exception as exc:
                print(f"  fluxmap fail {dest.name}: {exc}", flush=True)
    _fm_cache[ck] = res
    return res


def _snr_stack(samples):
    f = np.array([s[0] for s in samples], float)
    v = np.array([s[1] for s in samples], float)
    okm = np.isfinite(f) & np.isfinite(v) & (v > 0)
    if not okm.any():
        return np.nan, 0.0
    f, v = f[okm], v[okm]
    w = 1.0 / v
    w = np.minimum(w, WEIGHT_CAP * np.median(w))
    S = float(np.sum(w * f) / np.sqrt(np.sum(w)))
    n_eff = float(np.sum(w) ** 2 / np.sum(w ** 2))
    return S, n_eff


def _dedup_best(samples):
    """samples: list of (tbin, good_frac, flux, var); collapse skycell
    duplicates keeping best good_frac per exposure-time bin."""
    best = {}
    for tb, g, f, v in samples:
        if tb not in best or g > best[tb][0]:
            best[tb] = (g, f, v)
    return [(f, v) for _, f, v in best.values()]


def stats():
    adapter = Ps1WarpAdapter(skycell_cache=SKYCELL_CACHE)
    store = SnapshotStore(RUN)
    plan_doc = _load_plan()
    results = {"A": [], "A_note": plan_doc["A_note"], "B": []}
    cat_cache = {}

    for p in plan_doc["B"]:
        tid, band, eid = p["target_id"], p["band"], p["event_id"]
        _, sub = target_cone(tid)
        ev = event_row(sub, eid)
        omap = {nkey(o): o for o in discover(adapter, store, tid)}

        # off-window samples at the fixed antipode point -> rescale k
        off = []
        zp_src = defaultdict(int)
        for e in p["off_epochs"]:
            o = omap.get(e["key"])
            if o is None:
                continue
            fm = _fluxmap(tid, o, CUT_OFF)
            if fm is None:
                continue
            zp_src[fm.zp_source] += 1
            f, v, g = fm.sample(p["anti"][0], p["anti"][1])
            if np.isfinite(g[0]) and g[0] >= GOOD_FRAC_MIN:
                off.append((round(o.t_mid_mjd_utc / DUP_BIN_DAYS),
                            g[0], f[0], v[0]))
        offd = _dedup_best(off)
        k = 1.0
        k_status = "undefined_n<6"
        if len(offd) >= 6:
            fo = np.array([f for f, _ in offd])
            vo = np.array([v for _, v in offd])
            keep = np.isfinite(fo) & np.isfinite(vo) & (vo > 0)
            fo, vo = fo[keep], vo[keep]
            for _ in range(3):
                med = np.median(fo)
                sd = 1.4826 * np.median(np.abs(fo - med)) or 1.0
                m = np.abs(fo - med) <= 3 * sd
                fo, vo = fo[m], vo[m]
            if len(fo) >= 6:
                med = np.median(fo)
                k = max(1.0, float(np.median((fo - med) ** 2 / vo)
                                   / 0.4549))
                k_status = f"ok_n{len(fo)}"

        # in-window flux maps
        infms = []
        for e in p["in_epochs"]:
            o = omap.get(e["key"])
            if o is None:
                continue
            fm = _fluxmap(tid, o, CUT_IN)
            if fm is not None:
                zp_src[fm.zp_source] += 1
                infms.append((o, fm))

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
                        sm.append((round(o.t_mid_mjd_utc / DUP_BIN_DAYS),
                                   g[0], f[0], v[0] * k))
                smd = _dedup_best(sm)
                if smd:
                    S, _ = _snr_stack(smd)
                    best = max(best, S)
            return best if np.isfinite(best) else np.nan

        S = track_S(0.0, 0.0)
        controls = [track_S(dx, dy) for dx, dy in RING]
        controls = [c for c in controls if np.isfinite(c)]
        T = max(controls) if controls else np.nan

        # static-source annotation: any in-window track node within 2"
        corr = corridor_name(tid)
        if corr not in cat_cache:
            try:
                cat_cache[corr] = load_ps1_mean(SCREEN, corr, band="r")
            except Exception:
                cat_cache[corr] = None
        cat = cat_cache[corr]
        static_min = None
        if cat is not None:
            dmin = np.inf
            for e in p["in_epochs"]:
                for z in Z_GRID:
                    ra, dec = e["pos"][str(z)]
                    d = np.hypot((cat.ra - ra)
                                 * np.cos(np.radians(dec)),
                                 cat.dec - dec) * 3600.0
                    m = cat.ndet >= 3
                    if m.any():
                        dmin = min(dmin, float(np.min(d[m])))
            static_min = None if not np.isfinite(dmin) else round(dmin, 2)

        n_dedup = len(_dedup_best(
            [(round(o.t_mid_mjd_utc / DUP_BIN_DAYS), 1.0, 0.0, 1.0)
             for o, _ in infms]))
        exceed = bool(np.isfinite(S) and S > max(T if np.isfinite(T)
                                                 else -np.inf, 0.0))
        status = ("track_masked" if np.isnan(S)
                  else "ok" if len(controls) == 8
                  else f"controls_partial_{len(controls)}of8")
        results["B"].append({
            "target_id": tid, "band": band, "event_id": eid,
            "t_ca_mjd": p["t_ca_mjd"], "status": status,
            "n_in_epochs_raw": len(p["in_epochs"]),
            "n_in_epochs_fetched": len(infms),
            "n_in_exposures_dedup": n_dedup,
            "n_off_samples": len(offd), "k_rescale": round(k, 3),
            "k_status": k_status,
            "zp_sources": dict(zp_src),
            "S": None if np.isnan(S) else round(S, 3),
            "n_controls": len(controls),
            "controls": [round(c, 3) for c in controls],
            "T": None if not controls else round(T, 3),
            "margin": (None if not controls or np.isnan(S)
                       else round(S - max(T, 0.0), 3)),
            "static_min_arcsec": static_min,
            "static_annotated": (static_min is not None
                                 and static_min < STATIC_ARCSEC),
            "class": ("search_unit" if n_dedup >= 2 else "single_epoch"),
            "exceedance": exceed})

    out = D / "results" / "dev_search_v1.json"
    out.write_text(json.dumps(results, indent=1) + "\n")
    print(json.dumps(results, indent=1))


def satgate():
    """Frozen saturation verification gate: per band, bright DR2 mean
    stars in the dev corridors measured on in-band warp cutouts;
    empirical SAT/STARCORE core bits vs the frozen exclusion levels,
    plus the header CELL.SATURATION-derived point-source m_sat."""
    from astropy.io import fits
    from astropy.wcs import WCS as AWCS

    adapter = Ps1WarpAdapter(skycell_cache=SKYCELL_CACHE)
    store = SnapshotStore(RUN)
    out = {"frozen_exclude": SAT_E, "bands": {}}
    SATBITS = 32 | 4096          # SAT | STARCORE
    for band in "grizy":
        rows = []
        for tid in DEV_B:
            corr = corridor_name(tid)
            try:
                cat = load_ps1_mean(SCREEN, corr, band=band)
            except Exception:
                continue
            m = np.isfinite(cat.mag) & (cat.mag > 5) & (cat.mag < 16.5) \
                & (cat.ndet >= 3)
            idx0 = np.where(m)[0]
            # two magnitude bins: near the empirical mask boundary and
            # inside the frozen marginal zone
            idx = np.concatenate([
                idx0[np.argsort(np.abs(cat.mag[idx0] - 13.5))][:4],
                idx0[np.argsort(np.abs(cat.mag[idx0] - 14.75))][:4]])
            _, uniq_i = np.unique(idx, return_index=True)
            idx = idx[np.sort(uniq_i)]
            obs = [o for o in discover(adapter, store, tid)
                   if o.band == band and okq(o)
                   and o.extra["skycell"] in adapter._skycells]
            for i in idx:
                ra, dec = float(cat.ra[i]), float(cat.dec[i])
                o = next((o for o in obs
                          if adapter.nominal_footprint(o).contains(ra, dec)),
                         None)
                if o is None:
                    continue
                dest = PRODUCTS / "satgate" / f"{corr}-{band}-{i}"
                try:
                    cut = CutoutSpec(ra_deg=ra, dec_deg=dec,
                                     size_pix=CUT_SAT)
                    ps_i = adapter.fetch(o, ("img",), dest, cutout=cut)
                    ps_m = adapter.fetch(o, ("msk",), dest, cutout=cut)
                except Exception:
                    continue
                with fits.open(ps_i.products[0].path) as hd:
                    hdu = next(h for h in hd if h.data is not None)
                    hdr = dict(hdu.header)
                    img_wcs = AWCS(hdu.header)
                with fits.open(ps_m.products[0].path) as hd:
                    mh = next(h for h in hd if h.data is not None)
                    mdat = np.where(np.isfinite(mh.data), mh.data, 0
                                    ).astype(np.int64)
                x, y = img_wcs.wcs_world2pix([[ra, dec]], 0)[0]
                xi, yi = int(round(x)), int(round(y))
                r = 5
                core = mdat[max(0, yi - r):yi + r + 1,
                            max(0, xi - r):xi + r + 1]
                satbits = bool((core & SATBITS).any()) if core.size else None
                sat_only = bool((core & 32).any()) if core.size else None
                starcore_only = (bool((core & 4096).any())
                                 if core.size else None)
                satlvl = hdr.get("CELL.SATURATION",
                                 hdr.get("HIERARCH CELL.SATURATION"))
                see = float(hdr.get("CHIP.SEEING",
                                    hdr.get("HIERARCH CHIP.SEEING", 5.0))
                            or 5.0)
                zp_h = hdr.get("FPA.ZP", hdr.get("HIERARCH FPA.ZP"))
                ext = float(hdr.get("EXPTIME", 0) or 0)
                m_sat = None
                if satlvl and zp_h and ext > 0:
                    _, delta = zeropoints()
                    zp_eff = (float(zp_h) + 2.5 * np.log10(ext)
                              + delta.get(band, -0.5))
                    sig = see / 2.355
                    peak_frac = 1.0 / (2 * np.pi * sig * sig)
                    m_sat = round(zp_eff - 2.5 * np.log10(
                        float(satlvl) / peak_frac), 2)
                rows.append({"corridor": corr, "mag": round(float(
                    cat.mag[i]), 2), "sat_core_bits": satbits,
                    "sat_bit": sat_only, "starcore_bit": starcore_only,
                    "cell_saturation": satlvl, "seeing_pix": see,
                    "m_sat_header": m_sat,
                    "warp": o.native_key["mjd_tag"]})
        sat_mags = [r["mag"] for r in rows if r["sat_core_bits"]]
        clean = [r["mag"] for r in rows if r["sat_core_bits"] is False]
        out["bands"][band] = {
            "stars": rows,
            "faintest_with_sat_bits": max(sat_mags) if sat_mags else None,
            "brightest_clean": min(clean) if clean else None,
            "median_m_sat_header": (float(np.median(
                [r["m_sat_header"] for r in rows
                 if r["m_sat_header"] is not None]))
                if any(r["m_sat_header"] is not None for r in rows)
                else None),
            "verdict": "unverified" if not rows else "see_stars"}
    (D / "results" / "saturation_verify_v1.json").write_text(
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
    else:
        print("usage: dev_search.py --plan | --fetch | --stats | --satgate")
