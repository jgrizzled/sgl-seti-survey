"""Injection completeness (threshold freeze v1.0 §completeness).

Exact stamp-response injections (sglsurvey.inject.stamp_response — the
matched filter is linear, so the response to an added Moffat beta=3
stamp at the frame's SEEING FWHM is computed locally and exactly)
through the identical star-calibrated chain, for every searched
confirmatory channel-B unit. The temporal model is the frozen d = 1
persistent-recurrent chord: the injected relay is present in every
covered event of the unit (top-hat across each window), so S_event
and S_stack see it exactly as the search would.

Per unit: magnitude grid 15–22 AB (0.5 steps), 200 draws per grid
point. A draw picks the injected z from the recoverable subset of the
frozen grid (shared across the unit's events — one physical relay
distance), draws per-epoch background fluxes with replacement from
the unit's clipped off-window pool, adds F x R (R = exact per-epoch
stamp response sampled on every z' track node, so the nested-z max
sees the injected source exactly), and recovers iff any of the
unit's statistics exceeds its measured frozen threshold
(S > max(T, 0) with T from confirmatory_v1.json). m90 = interpolated
90 % recovery; 22.0 values are grid-censored lower limits. A z_inj
masked in every usable epoch is declared dead and excluded (reported).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from sglsurvey.adapters.irsa_ptf import (PtfLevel1Adapter,
                                         MASK_FATAL_TEMPLATE)
from sglsurvey.inject import MoffatPSF, stamp_response
from sglsurvey.photometry import build_flux_map_ptf
from sglsurvey.snapshots import SnapshotStore

import dev_search as ds
import confirmatory_search as cs
from coverage_intersect import relay_apparent

D = ds.D
RUN = ds.RUN
MAGS = np.arange(15.0, 22.01, 0.5)
NDRAW = 200
RNG = np.random.default_rng(20260826)
Z_GRID = ds.Z_GRID


def build_keep(grp, key, o, band, calname):
    """keep_inputs FluxMap + ZP-25 scale factor, same gate as search."""
    dest = ds._dest(grp, key)
    if not dest.exists():
        return None
    pref = f"cut{cs.CUT}-"
    sci = next((f for f in dest.iterdir()
                if f.name.startswith(pref) and "_scie_" in f.name
                and f.name.endswith(".fits")), None)
    msk = next((f for f in dest.iterdir()
                if f.name.startswith(pref) and "_mask_" in f.name), None)
    if sci is None or msk is None:
        return None
    try:
        fm = build_flux_map_ptf(sci, msk, MASK_FATAL_TEMPLATE, band,
                                o.t_mid_mjd_utc, keep_inputs=True)
    except Exception:
        return None
    zps = []
    for cra, cde, cm in ds.load_calibrators(calname, band):
        f, v, g = fm.sample(cra, cde)
        if (np.isfinite(g[0]) and g[0] >= ds.GOOD_FRAC_MIN
                and np.isfinite(f[0]) and v[0] > 0
                and f[0] > 5 * np.sqrt(v[0])):
            zps.append(cm + 2.5 * np.log10(f[0]))
    if len(zps) < ds.CAL_MIN:
        return None
    zps = np.array(zps)
    zp = float(np.median(zps))
    if 1.4826 * float(np.median(np.abs(zps - zp))) > ds.CAL_SCATTER_MAX:
        return None
    return fm, 10 ** (0.4 * (ds.ZP_REF - zp))


def snr_stack_arr(f, v):
    w = 1.0 / v
    w = np.minimum(w, ds.WEIGHT_CAP * np.median(w))
    return float(np.sum(w * f) / np.sqrt(np.sum(w)))


def m90_from(recovery):
    r = np.asarray(recovery)
    if r[0] < 0.9:
        return None, False
    below = np.where(r < 0.9)[0]
    if not below.size:
        return float(MAGS[-1]), True
    i = below[0]
    m_lo, m_hi = MAGS[i - 1], MAGS[i]
    r_lo, r_hi = r[i - 1], r[i]
    frac = (r_lo - 0.9) / max(r_lo - r_hi, 1e-9)
    return float(m_lo + frac * (m_hi - m_lo)), False


def complete_unit(u_plan, u_res, adapter, store):
    tid, band = u_plan["target_id"], u_plan["band"]
    grp = f"conf-{tid}"
    calname = cs.cal_name(tid)
    omap = {ds.nkey(o): o for o in ds.discover_at(
        adapter, store, f"conf-{tid}", *u_plan["anti"])}
    k = u_res["k_rescale"]
    stats = {n: s for n, s in u_res["statistics"].items() if "S" in s}
    if not stats:
        return {"status": "no_searched_statistic"}
    thresh = {n: max(s["T"] if s["T"] is not None else -np.inf, 0.0)
              for n, s in stats.items()}

    # off-window pool (ZP-25, clipped)
    pool = []
    for e in u_plan["off_epochs"]:
        fm = ds.fluxmap(grp, e["key"], omap[e["key"]], band, calname,
                        cs.CUT)
        if not hasattr(fm, "sample"):
            continue
        f, v, g = fm.sample(u_plan["anti"][0], u_plan["anti"][1])
        if np.isfinite(g[0]) and g[0] >= ds.GOOD_FRAC_MIN:
            pool.append(f[0])
    pool = np.asarray(pool)
    for _ in range(3):
        if not pool.size:
            break
        med = np.median(pool)
        sd = 1.4826 * np.median(np.abs(pool - med)) or 1.0
        pool = pool[np.abs(pool - med) <= 3 * sd]
    if pool.size < 4:
        return {"status": "no_background_pool"}

    tab = cs.target_events(tid)
    nz = len(Z_GRID)
    events = []       # per event: list of epoch dicts
    for e in u_plan["events"]:
        row = tab[np.asarray([str(x) == e["event_id"]
                              for x in tab["event_id"]])][0]
        eps = []
        for epch in e["in_epochs"]:
            o = omap.get(epch["key"])
            built = build_keep(grp, epch["key"], o, band, calname)
            if built is None:
                continue
            fm, s = built
            pos = np.array([relay_apparent(row, z, o.t_mid_mjd_utc)
                            for z in Z_GRID])
            f, v, g = fm.sample(pos[:, 0], pos[:, 1])
            if not (np.isfinite(g) & (g >= ds.GOOD_FRAC_MIN)).any():
                continue
            pix = fm.world2pix(pos[:, 0], pos[:, 1])
            psf = MoffatPSF(fwhm_pix=fm.fwhm_pix, beta=3.0, band=band)
            R = np.zeros((nz, nz))
            for zi in range(nz):
                stamp, ox, oy = psf.render(pix[zi, 0], pix[zi, 1],
                                           half=16)
                rw = stamp_response(fm, stamp, ox, oy)
                R[zi] = rw.sample(pix[:, 0], pix[:, 1])
            eps.append({"g": np.where(np.isfinite(g), g, -1.0),
                        "v": v * s * s * k, "R": R})
        events.append(eps)
    if not any(events):
        return {"status": "no_usable_epochs"}

    # usable (event, epoch) per z'
    def usable(eps, zj):
        return [i for i, e in enumerate(eps)
                if e["g"][zj] >= ds.GOOD_FRAC_MIN
                and np.isfinite(e["v"][zj])]

    alive = [zi for zi in range(nz)
             if any(e["R"][zi, zj] > 1e-6
                    for eps in events for zj in range(nz)
                    for i in usable(eps, zj) for e in [eps[i]])]
    if not alive:
        return {"status": "all_z_masked"}

    multi_idx = [j for j, eps in enumerate(events) if len(eps) >= 2]
    cov_idx = [j for j, eps in enumerate(events) if len(eps) >= 1]

    rec = np.zeros(len(MAGS))
    for _ in range(NDRAW):
        zi = alive[int(RNG.integers(len(alive)))]
        bgs = [[RNG.choice(pool) for _ in eps] for eps in events]
        for mi, m in enumerate(MAGS):
            F = 10 ** (0.4 * (ds.ZP_REF - m))
            got = False
            if "S_event" in stats:
                best = -np.inf
                for j in multi_idx:
                    eps = events[j]
                    for zj in range(nz):
                        idx = usable(eps, zj)
                        if len(idx) < 2:
                            continue
                        f = np.array([bgs[j][i]
                                      + F * eps[i]["R"][zi, zj]
                                      for i in idx])
                        v = np.array([eps[i]["v"][zj] for i in idx])
                        best = max(best, snr_stack_arr(f, v))
                got |= best > thresh["S_event"]
            if not got and "S_stack" in stats:
                best = -np.inf
                for zj in range(nz):
                    f, v = [], []
                    for j in cov_idx:
                        eps = events[j]
                        for i in usable(eps, zj):
                            f.append(bgs[j][i]
                                     + F * eps[i]["R"][zi, zj])
                            v.append(eps[i]["v"][zj])
                    if f:
                        best = max(best, snr_stack_arr(np.array(f),
                                                       np.array(v)))
                got |= best > thresh["S_stack"]
            if got:
                rec[mi] += 1
    rec /= NDRAW
    m90, cens = m90_from(rec)
    return {"status": "ok", "m90": (None if m90 is None
                                    else round(m90, 2)),
            "grid_censored": cens,
            "recovery": [round(float(r), 3) for r in rec],
            "thresholds": {n: round(t, 3) for n, t in thresh.items()},
            "z_alive": [float(Z_GRID[z]) for z in alive],
            "z_dead": [float(Z_GRID[z]) for z in range(nz)
                       if z not in alive],
            "n_events_multi": len(multi_idx),
            "n_events_covered": len(cov_idx),
            "n_epochs_used": int(sum(len(e) for e in events)),
            "n_off_pool": int(pool.size)}


def main():
    adapter = PtfLevel1Adapter()
    store = SnapshotStore(RUN)
    plan = cs._plan()
    res = json.loads((D / "results" / "confirmatory_v1.json").read_text())
    res_by = {(r["target_id"], r["band"], r["radius_au"]): r
              for r in res["units"]}
    out = {"mags": MAGS.tolist(), "ndraw": NDRAW, "units": []}
    for u_plan in plan["units"]:
        key = (u_plan["target_id"], u_plan["band"],
               float(u_plan["radius_au"]))
        u_res = res_by[key]
        print(f"injecting {key}", flush=True)
        c = complete_unit(u_plan, u_res, adapter, store)
        c.update({"target_id": key[0], "band": key[1],
                  "radius_au": key[2]})
        out["units"].append(c)
        print(f"  -> {c['status']} m90={c.get('m90')}", flush=True)
    (D / "results" / "completeness_v1.json").write_text(
        json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
