"""Injection completeness (threshold freeze v1.0 §completeness).

Exact stamp-response injections (sglsurvey.inject.stamp_response — the
matched filter is linear, so the response to an added Moffat beta=3
stamp at the warp's CHIP.SEEING FWHM is computed locally and exactly)
through the identical star-calibrated chain, for every searchable
channel-B unit (dev + confirmatory, all rungs incl. the van-maanen
grazing family and the single-epoch class) and channel-A unit
(searchable; constraint-only A units get a clearly-labeled reference
depth against their available controls).

Per unit: magnitude grid 14-22 (0.5 steps), 200 bootstrap draws per
grid point. A draw picks the injected z from the frozen grid (channel
B), draws per-epoch background fluxes with replacement from the unit's
clipped off-window sample pool, adds F x R (R = exact per-epoch
stamp response sampled on every z' track, so the nested-z max sees the
injected source exactly as the search would), and recovers iff
S_inj > max(T, 0) with the unit's actual frozen threshold. m90 =
interpolated 90% recovery; values at 22.0 are grid-censored lower
limits. Chord temporal model: d = 1 top-hat across the window
(hypotheses v1.0) — every in-window epoch carries the source.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from sglsurvey.adapters.mast_ps1 import Ps1ExactFootprint, Ps1WarpAdapter
from sglsurvey.inject import MoffatPSF, stamp_response
from sglsurvey.photometry import build_flux_map_ps1
from sglsurvey.snapshots import SnapshotStore

import dev_search as ds
import confirmatory_search as cs
from coverage_intersect import relay_apparent

D = ds.D
RUN = ds.RUN
MAGS = np.arange(14.0, 22.01, 0.5)
NDRAW = 200
RNG = np.random.default_rng(20260824)
Z_GRID = ds.Z_GRID


def build_keep(tid, o, size, base_dir):
    """FluxMap with stored kernel/denominator (keep_inputs) + ZP-25
    scale factor; returns (fm, s) or None."""
    ds.PRODUCTS = base_dir
    dest = ds._dest(tid, ds.nkey(o))
    if not dest.exists():
        return None
    pref = f"cut{size}-"
    img = next((p for p in dest.iterdir() if p.name.startswith(pref)
                and p.name.endswith(".fits") and ".wt." not in p.name
                and ".mask." not in p.name), None)
    wt = next((p for p in dest.iterdir() if p.name.startswith(pref)
               and ".wt." in p.name), None)
    msk = next((p for p in dest.iterdir() if p.name.startswith(pref)
                and ".mask." in p.name), None)
    if img is None or msk is None:
        return None
    try:
        fm = build_flux_map_ps1(img, wt, msk, Ps1ExactFootprint.FATAL_MASK,
                                o.band, o.t_mid_mjd_utc, keep_inputs=True)
    except Exception:
        return None
    zpv, _ = ds._zp_eff(o, fm)
    if zpv is None:
        return None
    return fm, 10 ** (0.4 * (ds.ZP_REF - zpv))


def snr_stack_arr(f, v):
    w = 1.0 / v
    w = np.minimum(w, ds.WEIGHT_CAP * np.median(w))
    return float(np.sum(w * f) / np.sqrt(np.sum(w)))


def m90_from(recovery):
    """recovery: array over MAGS (fraction). Returns (m90, censored)."""
    r = np.asarray(recovery)
    if r[0] < 0.9:
        return None, False           # never reaches 90% even at 14
    below = np.where(r < 0.9)[0]
    if not below.size:
        return float(MAGS[-1]), True
    i = below[0]                     # first grid point under 90%
    m_lo, m_hi = MAGS[i - 1], MAGS[i]
    r_lo, r_hi = r[i - 1], r[i]
    frac = (r_lo - 0.9) / max(r_lo - r_hi, 1e-9)
    return float(m_lo + frac * (m_hi - m_lo)), False


def complete_b(plan_row, res_row, base_dir, adapter, store):
    tid, band, eid = (plan_row["target_id"], plan_row["band"],
                      plan_row["event_id"])
    _, sub = ds.target_cone(tid)
    ev = ds.event_row(sub, eid)
    omap = {ds.nkey(o): o for o in ds.discover(adapter, store, tid)}
    k = res_row["k_rescale"]
    T = res_row.get("T")
    thresh = max(T if T is not None else -np.inf, 0.0)

    # off-window pool (ZP-25 scaled, clipped)
    ds.PRODUCTS = base_dir
    pool = []
    for e in plan_row["off_epochs"]:
        o = omap.get(e["key"])
        if o is None:
            continue
        fm = ds._fluxmap(tid, o, ds.CUT_OFF)
        if fm is None:
            continue
        f, v, g = fm.sample(plan_row["anti"][0], plan_row["anti"][1])
        if np.isfinite(g[0]) and g[0] >= ds.GOOD_FRAC_MIN:
            pool.append(f[0])
    pool = np.asarray(pool)
    for _ in range(3):
        med = np.median(pool)
        sd = 1.4826 * np.median(np.abs(pool - med)) or 1.0
        pool = pool[np.abs(pool - med) <= 3 * sd]
    if pool.size < 4:
        return {"status": "no_off_pool"}

    # per-epoch: base samples, variances, responses R[z_inj][z'] --
    # duplicates collapsed by best good_frac exactly as the search
    ep = []          # (tbin, g, v_scaled, R (nz, nz))
    for e in plan_row["in_epochs"]:
        o = omap.get(e["key"])
        if o is None:
            continue
        built = build_keep(tid, o, ds.CUT_IN, base_dir)
        if built is None:
            continue
        fm, s = built
        pos = np.array([relay_apparent(ev, z, o.t_mid_mjd_utc)
                        for z in Z_GRID])
        f, v, g = fm.sample(pos[:, 0], pos[:, 1])
        if not (np.isfinite(g) & (g >= ds.GOOD_FRAC_MIN)).any():
            continue
        pix = fm.world2pix(pos[:, 0], pos[:, 1])
        psf = MoffatPSF(fwhm_pix=fm.fwhm_pix, beta=3.0, band=band)
        R = np.zeros((len(Z_GRID), len(Z_GRID)))
        for zi in range(len(Z_GRID)):
            stamp, ox, oy = psf.render(pix[zi, 0], pix[zi, 1], half=16)
            rw = stamp_response(fm, stamp, ox, oy)
            R[zi] = rw.sample(pix[:, 0], pix[:, 1])
        gz = np.where(np.isfinite(g), g, -1.0)
        ep.append({"tbin": round(o.t_mid_mjd_utc / ds.DUP_BIN_DAYS),
                   "g": gz, "v": v * s * s * k, "R": R})
    if not ep:
        return {"status": "no_usable_epochs"}

    # dedup per z' by best good_frac (matches the search's selection)
    nz = len(Z_GRID)
    sel = {}         # z' -> list of epoch indices
    for zj in range(nz):
        best = {}
        for i, e in enumerate(ep):
            if e["g"][zj] < ds.GOOD_FRAC_MIN or not np.isfinite(e["v"][zj]):
                continue
            tb = e["tbin"]
            if tb not in best or e["g"][zj] > ep[best[tb]]["g"][zj]:
                best[tb] = i
        sel[zj] = sorted(best.values())
    if not any(sel.values()):
        return {"status": "no_usable_epochs"}

    rec = np.zeros(len(MAGS))
    for _ in range(NDRAW):
        zi = int(RNG.integers(nz))
        bg = {i: RNG.choice(pool) for i in range(len(ep))}
        for mi, m in enumerate(MAGS):
            F = 10 ** (0.4 * (ds.ZP_REF - m))
            best = -np.inf
            for zj in range(nz):
                idx = sel[zj]
                if not idx:
                    continue
                f = np.array([bg[i] + F * ep[i]["R"][zi, zj] for i in idx])
                v = np.array([ep[i]["v"][zj] for i in idx])
                best = max(best, snr_stack_arr(f, v))
            if best > thresh:
                rec[mi] += 1
    rec /= NDRAW
    m90, cens = m90_from(rec)
    return {"status": "ok", "m90": m90, "grid_censored": cens,
            "recovery": [round(float(r), 3) for r in rec],
            "threshold": round(thresh, 3),
            "n_epochs_used": len(ep), "n_off_pool": int(pool.size)}


def complete_a(plan_row, res_row, adapter, store):
    tid, band = plan_row["target_id"], plan_row["band"]
    omap = {ds.nkey(o): o for o in cs.discover_a(adapter, store, tid)}
    k = res_row["k_rescale"]
    T = res_row.get("T")
    thresh = max(T if T is not None else -np.inf, 0.0)
    if not np.isfinite(thresh):
        return {"status": "no_controls"}
    wins = plan_row["windows"]
    base = res_row.get("baseline_zp25")
    if base is None:
        return {"status": "no_baseline"}

    ds.PRODUCTS = cs.CONF_PRODUCTS
    resid_pool, inwin = [], []
    for e in plan_row["epochs"]:
        o = omap.get(e["key"])
        if o is None:
            continue
        is_in = any(w["lo"] <= e["t"] <= w["hi"] for w in wins)
        if is_in:
            built = build_keep(tid, o, cs.CUT_A, cs.CONF_PRODUCTS)
            if built is None:
                continue
            fm, s = built
            pos = np.array([e["pos"]])
            f, v, g = fm.sample(pos[:, 0], pos[:, 1])
            if not (np.isfinite(g[0]) and g[0] >= ds.GOOD_FRAC_MIN):
                continue
            pix = fm.world2pix(pos[:, 0], pos[:, 1])
            psf = MoffatPSF(fwhm_pix=fm.fwhm_pix, beta=3.0, band=band)
            stamp, ox, oy = psf.render(pix[0, 0], pix[0, 1], half=16)
            rw = stamp_response(fm, stamp, ox, oy)
            R = float(rw.sample(pix[:, 0], pix[:, 1])[0])
            wi = next(j for j, w in enumerate(wins)
                      if w["lo"] <= e["t"] <= w["hi"])
            inwin.append({"win": wi, "v": float(v[0]) * s * s * k, "R": R})
        else:
            fm = ds._fluxmap(tid, o, cs.CUT_A)
            if fm is None:
                continue
            f, v, g = fm.sample(e["pos"][0], e["pos"][1])
            if np.isfinite(g[0]) and g[0] >= ds.GOOD_FRAC_MIN:
                resid_pool.append(f[0] - base)
    pool = np.asarray(resid_pool)
    for _ in range(3):
        med = np.median(pool)
        sd = 1.4826 * np.median(np.abs(pool - med)) or 1.0
        pool = pool[np.abs(pool - med) <= 3 * sd]
    if pool.size < 4 or not inwin:
        return {"status": "insufficient_data"}

    rec = np.zeros(len(MAGS))
    for _ in range(NDRAW):
        bg = RNG.choice(pool, size=len(inwin))
        for mi, m in enumerate(MAGS):
            F = 10 ** (0.4 * (ds.ZP_REF - m))
            per_win = {}
            for i, e in enumerate(inwin):
                per_win.setdefault(e["win"], []).append(
                    (bg[i] + F * e["R"], e["v"]))
            Sws = []
            for sm in per_win.values():
                f = np.array([a for a, _ in sm])
                v = np.array([b for _, b in sm])
                Sws.append(snr_stack_arr(f, v))
            S = float(np.sum(Sws) / np.sqrt(len(Sws)))
            if S > thresh:
                rec[mi] += 1
    rec /= NDRAW
    m90, cens = m90_from(rec)
    return {"status": "ok", "m90": m90, "grid_censored": cens,
            "recovery": [round(float(r), 3) for r in rec],
            "threshold": round(thresh, 3),
            "n_in_epochs_used": len(inwin), "n_resid_pool": int(pool.size)}


def main():
    adapter = Ps1WarpAdapter(skycell_cache=ds.SKYCELL_CACHE)
    store = SnapshotStore(RUN)
    dev_plan = json.loads((D / "configs" / "dev_plan_v1.json").read_text())
    conf_plan = json.loads(cs.PLAN_PATH.read_text())
    dev_res = json.loads(
        (D / "results" / "dev_search_v1.json").read_text())
    conf_res = json.loads(
        (D / "results" / "confirmatory_v1.json").read_text())
    out = {"mag_grid": [float(m) for m in MAGS], "n_draws": NDRAW,
           "B": [], "A": []}

    for plan_doc, res_doc, base in (
            (dev_plan, dev_res, RUN / "products" / "dev"),
            (conf_plan, conf_res, cs.CONF_PRODUCTS)):
        for p in plan_doc["B"]:
            r = next((r for r in res_doc["B"]
                      if r["target_id"] == p["target_id"]
                      and r["band"] == p["band"]
                      and r["event_id"] == p["event_id"]
                      and r.get("radius_au", 0.1)
                      == p.get("radius_au", 0.1)), None)
            if r is None or r["status"] == "track_masked":
                continue
            print(f"B {p['target_id']} {p['band']} {p['event_id'][:14]} "
                  f"r={p.get('radius_au', 0.1)}", flush=True)
            c = complete_b(p, r, base, adapter, store)
            out["B"].append({"target_id": p["target_id"], "band": p["band"],
                             "event_id": p["event_id"],
                             "radius_au": p.get("radius_au", 0.1),
                             "class": r["class"], **c})

    for p in conf_plan["A"]:
        r = next((r for r in conf_res["A"]
                  if r["target_id"] == p["target_id"]
                  and r["band"] == p["band"]), None)
        if r is None:
            continue
        print(f"A {p['target_id']} {p['band']} [{r['status']}]", flush=True)
        c = complete_a(p, r, adapter, store)
        out["A"].append({"target_id": p["target_id"], "band": p["band"],
                         "status_search": r["status"],
                         "depth_kind": ("calibrated" if r["status"] ==
                                        "searchable" else
                                        "reference_constraint_only"), **c})

    (D / "results" / "completeness_v1.json").write_text(
        json.dumps(out, indent=1) + "\n")
    for ch in ("B", "A"):
        for r in out[ch]:
            print(ch, r["target_id"], r["band"],
                  r.get("event_id", "")[:14], r.get("radius_au", ""),
                  "m90:", r.get("m90"), "cens:", r.get("grid_censored"),
                  r["status"])


if __name__ == "__main__":
    main()
