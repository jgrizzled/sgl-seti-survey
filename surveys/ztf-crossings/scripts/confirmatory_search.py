"""Confirmatory search: channel B under freeze v1.0, channel A under
v1.0 + amendments v1.1 + v1.2, on the confirmatory split only.

Implementation correction over the dev scripts (declared in the
report): channel-A forced photometry and cutout centers use per-epoch
PROPAGATED star positions (registry linear astrometry: PM + parallax
displacement), since PM up to ~8 arcsec/yr over the 8.5-yr era dwarfs
the PSF. The frozen statistic ("forced photometry at the star
position") is unchanged - the star position at epoch t is the
propagated one.

Stages: --plan | --fetch | --stats
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.time import Time

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from dev_search import (BANDCODE, COV, CUT_A, CUT_B, D, ERA, EVENTS, FREEZE,
                        RUN, WEIGHT_CAP, Z_GRID, _snr_stack, target_windows)
from dev_search import build_flux_map_ztf  # noqa
from dev_search_v11 import (FIT_MAX, FIT_MIN, REG, parallax_factors,
                            valid_offsets_same_rung)
from coverage_intersect import relay_apparent
from sglsurvey.adapters.base import ConeRegion, CutoutSpec
from sglsurvey.adapters.irsa_ztf import (MASK_FATAL_TEMPLATE,
                                         ZtfNominalFootprint, ZtfSciAdapter)
from sglsurvey.snapshots import SnapshotStore

PRODUCTS = RUN / "products" / "conf"
CONF_A = set(FREEZE["split"]["A"]["confirmatory"])
DEV_B = set(FREEZE["split"]["B_wide"]["dev"])
PLAN_PATH = D / "configs" / "conf_plan_v1.json"
EV_BY_ID = {str(e["event_id"]): e for e in EVENTS}


def star_pos(tid, mjd):
    """Propagated apparent ICRS position (deg) at epoch(s) mjd."""
    st = REG[tid]["state"]["astrometry"]
    mjd = np.atleast_1d(np.asarray(mjd, float))
    dt = (Time(mjd, format="mjd", scale="utc").jyear
          - st["reference_epoch_jyear"])
    p_ra, p_de = parallax_factors(tid, mjd)          # arcsec
    cosd = np.cos(np.radians(st["dec_deg"]))
    ra = st["ra_deg"] + (st["pm_ra_cosdec_mas_per_yr"] / 1000.0 * dt
                         + p_ra) / 3600.0 / max(cosd, 1e-6)
    de = st["dec_deg"] + (st["pm_dec_mas_per_yr"] / 1000.0 * dt
                          + p_de) / 3600.0
    return ra, de


def a_unit_list():
    units = [dict(u, kind="unit") for u in FREEZE["search_units"]["A"]["units"]
             if u["target_id"] in CONF_A]
    units += [dict(u, kind="single_window")
              for u in FREEZE["search_units"]["A"]["single_window"]
              if u["target_id"] in CONF_A]
    return units


def b_row_list():
    rows = [dict(u, kind="unit") for u in FREEZE["search_units"]["B"]["units"]
            if u["target_id"] not in DEV_B]
    rows += [dict(u, kind="single_epoch")
             for u in FREEZE["search_units"]["B"]["single_epoch"]
             if u["target_id"] not in DEV_B]
    return rows


def quality_ok(o):
    q = o.quality_flags
    return not q.get("bad_quality") and (q.get("seeing") or 99) <= 4.0


def discover_cones():
    adapter, store = ZtfSciAdapter(), SnapshotStore(RUN)
    cones = {}
    for u in a_unit_list():
        tid = u["target_id"]
        if ("A", tid) not in cones:
            ra, de = star_pos(tid, 0.5 * (ERA.start_mjd_utc + ERA.stop_mjd_utc))
            cones[("A", tid)] = (float(ra[0]), float(de[0]), 0.06)
    for u in b_row_list():
        tid = u["target_id"]
        if ("B", tid) not in cones:
            ev = EV_BY_ID[u["event_id"]]
            cones[("B", tid)] = (float(ev["relay_icrs_ra_deg"]),
                                 float(ev["relay_icrs_dec_deg"]), 0.15)
    omaps = {}
    for key, (ra, de, rad) in sorted(cones.items()):
        obs = list(adapter.discover(ConeRegion(ra, de, rad), ERA, store))
        omaps[key] = {json.dumps(o.native_key, sort_keys=True): o for o in obs}
        print(f"[{key[0]}] {key[1]}: {len(obs)} rows", file=sys.stderr,
              flush=True)
    return adapter, omaps


def a_epochs(tid, band, omap):
    eps = []
    for k, o in omap.items():
        if o.band != BANDCODE[band] or not quality_ok(o):
            continue
        ra, de = star_pos(tid, o.t_mid_mjd_utc)
        if ZtfNominalFootprint(o).contains(float(ra[0]), float(de[0])):
            eps.append((k, o))
    return sorted(eps, key=lambda x: x[1].t_mid_mjd_utc)


def build_plan(adapter=None, omaps=None):
    if omaps is None:
        adapter, omaps = discover_cones()
    plan = {"A": [], "B": []}
    for u in a_unit_list():
        tid, band, rad = u["target_id"], u["band"], u["radius_au"]
        offs = valid_offsets_same_rung(tid, rad)
        allw = [(lo, hi) for lo, hi, _, _ in target_windows(tid)]
        rung = [(lo, hi) for lo, hi, r_, _ in target_windows(tid) if r_ == rad]
        eps = a_epochs(tid, band, omaps[("A", tid)])
        offw = [k for k, o in eps
                if not any(lo <= o.t_mid_mjd_utc <= hi for lo, hi in allw)]
        step = max(1, len(offw) // FIT_MAX)
        fit_keys = offw[::step][:FIT_MAX]
        spans = list(rung) + [(lo + o, hi + o) for o in offs
                              for lo, hi in rung]
        stat_keys = [k for k, o in eps
                     if any(lo <= o.t_mid_mjd_utc <= hi for lo, hi in spans)]
        plan["A"].append({
            "target_id": tid, "band": band, "radius_au": rad,
            "kind": u["kind"], "valid_offsets": offs,
            "rung_windows": [[lo, hi] for lo, hi in rung],
            "all_windows": [[lo, hi] for lo, hi in allw],
            "epoch_keys": sorted(set(fit_keys + stat_keys)),
            "n_fit": len(fit_keys), "n_stat": len(stat_keys)})
    for u in b_row_list():
        tid, band, eid = u["target_id"], u["band"], u["event_id"]
        rad = u["radius_au"]
        m = ((COV["channel"] == "B")
             & np.asarray([str(x) == eid for x in COV["event_id"]])
             & (np.asarray(COV["radius_au"]) == rad))
        row = COV[m][0]
        lo = row["t_ca_mjd"] - row["window_days"] / 2
        hi = row["t_ca_mjd"] + row["window_days"] / 2
        ev = EV_BY_ID[eid]
        keys = []
        for k, o in omaps[("B", tid)].items():
            if o.band != BANDCODE[band] or not quality_ok(o):
                continue
            if not (lo <= o.t_mid_mjd_utc <= hi):
                continue
            pos = [relay_apparent(ev, z, o.t_mid_mjd_utc) for z in Z_GRID]
            if any(ZtfNominalFootprint(o).contains(r_, d_) for r_, d_ in pos):
                keys.append(k)
        plan["B"].append({
            "target_id": tid, "band": band, "event_id": eid,
            "radius_au": rad, "kind": u["kind"],
            "window": [float(lo), float(hi)], "epoch_keys": sorted(keys)})
    PLAN_PATH.write_text(json.dumps(plan, indent=1) + "\n")
    nA = len({(p["target_id"], k) for p in plan["A"] for k in p["epoch_keys"]})
    nB = len({(p["target_id"], k) for p in plan["B"] for k in p["epoch_keys"]})
    print(json.dumps({"A_rows": len(plan["A"]), "B_rows": len(plan["B"]),
                      "A_epochs": nA, "B_epochs": nB}, indent=1))
    return plan, adapter, omaps


def _dest(ch, tid, o):
    return PRODUCTS / f"{ch}-{tid}" / str(o.native_key["pid"])


def fetch():
    plan, adapter, omaps = build_plan()
    tasks = []
    for p in plan["B"]:
        tid = p["target_id"]
        for k in p["epoch_keys"]:
            o = omaps[("B", tid)][k]
            ev = EV_BY_ID[p["event_id"]]
            ra, de = relay_apparent(ev, 2500.0, o.t_mid_mjd_utc)
            tasks.append(("B", tid, o, ra, de, CUT_B))
    for p in plan["A"]:
        tid = p["target_id"]
        for k in p["epoch_keys"]:
            o = omaps[("A", tid)][k]
            ra, de = star_pos(tid, o.t_mid_mjd_utc)
            tasks.append(("A", tid, o, float(ra[0]), float(de[0]), CUT_A))
    seen, todo = set(), []
    for t in tasks:
        key = str(_dest(t[0], t[1], t[2]))
        if key not in seen:
            seen.add(key)
            todo.append(t)
    todo = [t for t in todo if not _dest(t[0], t[1], t[2]).exists()]
    print(f"{len(todo)} epochs to fetch", flush=True)
    ok = fail = 0
    for i, (ch, tid, o, ra, de, size) in enumerate(todo):
        dest = _dest(ch, tid, o)
        cut = CutoutSpec(ra_deg=ra, dec_deg=de, size_pix=size)
        good = True
        for kind in ("sci", "msk", "diff"):
            try:
                adapter.fetch(o, (kind,), dest, cutout=cut)
            except FileNotFoundError:
                if kind != "diff":
                    good = False
                    break
            except Exception as e:
                print(f"  ! {dest}: {e}", flush=True)
                good, fail = False, fail + 1
                break
        ok += bool(good)
        if (i + 1) % 250 == 0:
            print(f"  {i+1}/{len(todo)} ok={ok} fail={fail}", flush=True)
    print(f"done {ok}/{len(todo)} fail={fail}", flush=True)


def _fm(ch, tid, o):
    dest = _dest(ch, tid, o)
    if not dest.exists():
        return None
    sci = next((x for x in dest.iterdir()
                if x.name.endswith("sciimg.fits")), None)
    msk = next((x for x in dest.iterdir()
                if x.name.endswith("mskimg.fits")), None)
    diff = next((x for x in dest.iterdir()
                 if x.name.endswith("scimrefdiffimg.fits.fz")), None)
    if sci is None or msk is None:
        return None
    try:
        return build_flux_map_ztf(sci, diff, msk, MASK_FATAL_TEMPLATE,
                                  o.band, o.t_mid_mjd_utc)
    except Exception:
        return None


def stats():
    plan = json.loads(PLAN_PATH.read_text())
    adapter, omaps = discover_cones()
    out = {"A": [], "B": []}

    for p in plan["A"]:
        tid, band = p["target_id"], p["band"]
        omap = omaps[("A", tid)]
        samples = {}
        for k in p["epoch_keys"]:
            o = omap[k]
            fm = _fm("A", tid, o)
            if fm is None:
                continue
            ra, de = star_pos(tid, o.t_mid_mjd_utc)
            f, v, g = fm.sample(float(ra[0]), float(de[0]))
            if np.isfinite(g[0]) and g[0] >= 0.7 and np.isfinite(f[0]):
                samples[k] = (o.t_mid_mjd_utc, float(f[0]), float(v[0]),
                              float(o.quality_flags.get("seeing") or 2.0))
        allw = [tuple(w) for w in p["all_windows"]]
        offi = [(k, s) for k, s in samples.items()
                if not any(lo <= s[0] <= hi for lo, hi in allw)]
        gate = {"n_fit": len(offi)}
        coef = med_see = None
        k_scale = 1.0
        if len(offi) >= FIT_MIN:
            t_a = np.array([s[0] for _, s in offi])
            f_a = np.array([s[1] for _, s in offi])
            v_a = np.array([s[2] for _, s in offi])
            see = np.array([s[3] for _, s in offi])
            pra, pde = parallax_factors(tid, t_a)
            med_see = float(np.median(see))
            B = np.column_stack([np.ones_like(t_a), (t_a - 59800.0) / 365.25,
                                 pra, pde, see - med_see])
            keep = np.ones(len(t_a), bool)
            for _ in range(3):
                coef, *_ = np.linalg.lstsq(B[keep], f_a[keep], rcond=None)
                res = f_a - B @ coef
                sd = 1.4826 * np.median(np.abs(res[keep] - np.median(res[keep])))
                keep = np.abs(res - np.median(res[keep])) <= 3 * sd
            gate["n_fit_kept"] = int(keep.sum())
            r2v = ((f_a - B @ coef) ** 2 / v_a)[keep]
            k_scale = max(1.0, float(np.median(r2v) / 0.4549))
            gate["k_scale"] = round(k_scale, 1)
        else:
            gate["n_fit_kept"] = len(offi)

        def resid(key):
            t, f, v, s = samples[key]
            if coef is None:
                return f, v
            pra, pde = parallax_factors(tid, t)
            b = np.array([1.0, (t - 59800.0) / 365.25, float(pra[0]),
                          float(pde[0]), s - med_see])
            return f - float(b @ coef), v * k_scale

        def span_S(spans):
            per = []
            for lo, hi in spans:
                sm = [resid(k) for k, s in samples.items() if lo <= s[0] <= hi]
                if sm:
                    per.append(_snr_stack(sm))
            if not per:
                return np.nan, 0
            Sw = np.array([s for s, _ in per])
            return float(np.sum(Sw) / np.sqrt(len(Sw))), len(per)

        rung = [tuple(w) for w in p["rung_windows"]]
        S, nw = span_S(rung)
        ctrls = []
        for off in p["valid_offsets"]:
            Sc, _ = span_S([(lo + off, hi + off) for lo, hi in rung])
            if np.isfinite(Sc):
                ctrls.append(round(Sc, 3))
        T = max(ctrls) if ctrls else np.nan
        medc = float(np.median(ctrls)) if ctrls else np.nan
        searchable = (p["kind"] == "unit" and len(p["valid_offsets"]) == 8
                      and gate["n_fit_kept"] >= FIT_MIN
                      and np.isfinite(medc) and abs(medc) <= 1.0
                      and nw >= 2)
        exceed = bool(searchable and np.isfinite(S)
                      and S > max(T if np.isfinite(T) else 0.0, 0.0))
        out["A"].append({
            "target_id": tid, "band": band, "radius_au": p["radius_au"],
            "kind": p["kind"], "n_valid_offsets": len(p["valid_offsets"]),
            "gate": gate, "n_windows_with_data": nw,
            "S": None if np.isnan(S) else round(S, 3),
            "controls": ctrls,
            "T": None if np.isnan(T) else round(T, 3),
            "median_control": None if np.isnan(medc) else round(medc, 3),
            "margin": (None if np.isnan(S) or not ctrls
                       else round(S - max(T, 0.0), 3)),
            "status": "searchable" if searchable else "constraint_only",
            "exceedance": exceed})

    RING = [(20.0, 0.0), (-20.0, 0.0), (30.0, 0.0), (-30.0, 0.0),
            (40.0, 0.0), (-40.0, 0.0), (0.0, 30.0), (0.0, -30.0)]
    for p in plan["B"]:
        tid, band, eid = p["target_id"], p["band"], p["event_id"]
        ev = EV_BY_ID[eid]
        omap = omaps[("B", tid)]
        fms = []
        for k in p["epoch_keys"]:
            o = omap[k]
            fm = _fm("B", tid, o)
            if fm is not None:
                fms.append((o, fm))

        def track_S(dx, dy):
            best = -np.inf
            for z in Z_GRID:
                sm = []
                for o, fm in fms:
                    ra, de = relay_apparent(ev, z, o.t_mid_mjd_utc)
                    ra += dx / 3600.0 / max(np.cos(np.radians(de)), .05)
                    de += dy / 3600.0
                    f, v, g = fm.sample(ra, de)
                    if np.isfinite(g[0]) and g[0] >= 0.7:
                        sm.append((f[0], v[0]))
                if sm:
                    s_, _ = _snr_stack(sm)
                    best = max(best, s_)
            return best if np.isfinite(best) else np.nan

        S = track_S(0.0, 0.0)
        ctrls = [c for c in (track_S(dx, dy) for dx, dy in RING)
                 if np.isfinite(c)]
        T = max(ctrls) if ctrls else np.nan
        n_ep = len(fms)
        searchable = p["kind"] == "unit" and n_ep >= 2 and len(ctrls) == 8
        exceed = bool(searchable and np.isfinite(S) and np.isfinite(T)
                      and T > 0 and S / T > 1)
        out["B"].append({
            "target_id": tid, "band": band, "event_id": eid,
            "radius_au": p["radius_au"], "kind": p["kind"], "n_epochs": n_ep,
            "S": None if np.isnan(S) else round(S, 3),
            "n_controls": len(ctrls),
            "T": None if np.isnan(T) else round(T, 3),
            "R": (None if np.isnan(S) or np.isnan(T) or T <= 0
                  else round(S / T, 3)),
            "status": "searchable" if searchable else p["kind"],
            "exceedance": exceed})

    (D / "results" / "confirmatory_v1.json").write_text(
        json.dumps(out, indent=1) + "\n")
    exc = [r for r in out["A"] + out["B"] if r["exceedance"]]
    print(json.dumps({"A_rows": len(out["A"]), "B_rows": len(out["B"]),
                      "exceedances": exc}, indent=1))


if __name__ == "__main__":
    if "--plan" in sys.argv:
        build_plan()
    elif "--fetch" in sys.argv:
        fetch()
    elif "--stats" in sys.argv:
        stats()
