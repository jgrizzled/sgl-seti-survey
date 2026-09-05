"""SPHEREx template-absorption control (plan §5.15 item O3 (ii)).

Question. The v3 static-sky template (per corridor x detector, per 3"
node: f = a + b (lam - lam0) fitted across epochs, two-pass 3-sigma
clip) was fitted to data that contain any real source, so a slow real
source is partly absorbed into (a, b) at the nodes it visits and the
search — which subtracts the template at sampling — sees only the
remainder. The v2 injections add the source to the *template-
subtracted* samples and so do not model this loss.

Measurement. For each (endpoint, role, detector) and a fixed ladder of
sources on the real trajectory, the injected per-epoch flux is rendered
with the exposure's own PSF plane and the matched-filter response
window (the v2 injection machinery, unchanged) at every template node
within R_NODE of the source position, the template is REFITTED with the
injected fluxes added (the full two-pass clipped fit at the affected
nodes), and the absorbed fraction at the source's own track is

    f_abs = sum_e w_e T_inj(track_e, lam_e) / sum_e w_e delta_e(track_e),

with T_inj = (refit - original) template evaluated bilinearly at the
track and w_e the epoch's inverse-variance weight at the track. The
equivalent completeness shift is dm = -2.5 log10(1 - f_abs): an
injection-calibrated m90 overstates the depth for a slow persistent
source by dm, so the corrected limit is m90 - dm.

Ladder (per cell): the 8 z-interval centres of the completeness
tiling at mu = 0 (a pure parallax track — the source spends the most
time per node, the conservative case) at three magnitudes (the cell's
v1 m90, 1.5 mag brighter, 1 mag fainter — the clip non-linearity
check), plus mu = (+1, +1) and (-1, -1) "/yr at the two most distant z
at m90 (the residual-motion check); temporal model persistent (every
epoch on); cross-track 0. Persistent is the conservative case, but
visit- and block-scale sources are absorbed just like it: a node's
epochs come from 2–3 visits and the moving source is near a node during
one visit only, so what matters is the fraction of the node's weight in
that visit, not the source's duty cycle — only exposure-flicker (random
on/off within a visit) is absorbed proportionally less (~duty x).
Validated end to end (template_control.py `validate`: source added to
the images, templates regenerated, tensor rebuilt) to |Δf_abs| < 0.01
in all six detectors on van-maanen/rx at z = 1364 AU.

Products: runs/spherex/template_control/absorb/<corridor>__<band>.json
per group; `summarise` aggregates over cells into
surveys/spherex/results/template_absorption.{json,md}.
"""

from __future__ import annotations

import json
import sys
import time as _time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import template_control as TC  # noqa: E402
import profile as spx  # noqa: E402
from spherex_corridors import MEMBERS  # noqa: E402
from sglsurvey import completeness_stage as CS  # noqa: E402

P = spx.PROFILE
REPO = TC.REPO
OUT = TC.OUT / "absorb"
R_NODE_ARCSEC = 15.0          # nodes within this radius of the source receive injected flux (response < 2e-4 beyond)
MAG_OFFSETS = (-1.5, 0.0, 1.0)  # relative to the cell's v1 m90 (persistent)
MU_CHECK = ((1.0, 1.0), (-1.0, -1.0))
N_MU_CHECK_Z = 2              # the most distant z-interval centres get the mu check
MODEL = "persistent"


def ladder(endpoint, role, band):
    """The frozen source ladder of one cell: list of dicts (z_au, iz, mu, mag, model)."""
    edges = CS.interval_edges(P)
    zc = np.sqrt(edges[:-1] * edges[1:])                 # geometric interval centres, ascending z
    q = 1.0 / P.z_grid
    m90 = spx.v1_m90(endpoint, role, band)
    if m90 is None:
        m90 = P.default_m90[band]
    out = []
    for z in zc:
        iz = int(np.argmin(np.abs(q - 1.0 / z)))
        for dm in MAG_OFFSETS:
            out.append({"z_au": float(P.z_grid[iz]), "iz": iz, "mu": (0.0, 0.0), "mag": float(m90 + dm), "model": MODEL,
                        "mag_offset": dm})
    for z in zc[-N_MU_CHECK_Z:]:
        iz = int(np.argmin(np.abs(q - 1.0 / z)))
        for mu in MU_CHECK:
            out.append({"z_au": float(P.z_grid[iz]), "iz": iz, "mu": mu, "mag": float(m90), "model": MODEL, "mag_offset": 0.0})
    return out, m90


def gnomonic_xy(ec, ra, dec):
    """Sky -> fractional template-node coordinates (StaticTemplate.__call__ convention)."""
    r0, d0 = np.deg2rad(ec.ra_c), np.deg2rad(ec.dec_c)
    r, d = np.deg2rad(ra), np.deg2rad(dec)
    cosc = np.sin(d0) * np.sin(d) + np.cos(d0) * np.cos(d) * np.cos(r - r0)
    xi = np.rad2deg(np.cos(d) * np.sin(r - r0) / cosc) * 3600
    eta = np.rad2deg((np.cos(d0) * np.sin(d) - np.sin(d0) * np.cos(d) * np.cos(r - r0)) / cosc) * 3600
    return ec.node_xy(xi, eta)


def bilinear(vals_full, x, y, side):
    """Bilinear lookup on a (side*side,) node array at fractional node
    coords; NaN where any corner is undefined or outside."""
    if not (0 <= x <= side - 1.001 and 0 <= y <= side - 1.001):
        return np.nan
    x0, y0 = int(np.floor(x)), int(np.floor(y)); fx, fy = x - x0, y - y0
    c = [vals_full[y0 * side + x0], vals_full[y0 * side + x0 + 1], vals_full[(y0 + 1) * side + x0], vals_full[(y0 + 1) * side + x0 + 1]]
    if not np.all(np.isfinite(c)):
        return np.nan
    return c[0] * (1 - fx) * (1 - fy) + c[1] * fx * (1 - fy) + c[2] * (1 - fx) * fy + c[3] * fx * fy


def run_group(corridor, band, recs, obs_by_id, manifest):
    """All cells of one corridor x detector: build the epoch cache, fit
    the original template, render every ladder source's injected flux
    at the nearby nodes per epoch, refit per source, measure absorption."""
    t0 = _time.monotonic()
    ec = TC.EpochCache(corridor, band, recs, obs_by_id, load=False)
    oid_to_k = {r["observation_id"]: k for k, r in enumerate(recs)}   # provisional (before load drops failures)
    side = ec.side
    r_node = R_NODE_ARCSEC / TC.SPACING_ARCSEC
    cells = []
    for e in P.members(corridor):
        for role in ("rx", "tx"):
            tp = P.tensor_dir / f"{e}__{role}.npz"
            if not tp.exists():
                continue
            with np.load(tp) as d:
                bi = P.band_idx[band]
                eb = np.where(d["band_idx"] == bi)[0]
                if len(eb) < P.min_epochs:
                    continue
                centre = tuple(float(x) for x in d["centre"])   # the tensor's own tangent point (a detector template's centre)
                oids = [str(o) for o in d["oid"][eb]]
                d0 = d["d0"][eb].astype(float)               # (E_b, nz, 2) arcsec offsets (dRA cos d, dDec)
                mjd = d["mjd"][eb].astype(float)
                # the stack's own per-epoch weights at the real trajectory (every ladder source sits on a grid node)
                f_t = d["f"][0, eb].astype(float); v_t = d["v"][0, eb].astype(float); g_t = d["g"][0, eb].astype(float)
                cap_t = d["frame_cap"][eb].astype(float)
                with np.errstate(all="ignore"):
                    valid = np.isfinite(f_t) & np.isfinite(v_t) & (v_t > 0) & (g_t >= P.min_good_frac) & (np.abs(f_t) / np.sqrt(np.where(v_t > 0, v_t, np.inf)) <= P.clip_sigma)
                w_t = np.where(valid, np.minimum(1.0 / np.where(valid, v_t, 1.0), cap_t[:, None, None, None]), 0.0)
            ks = np.array([oid_to_k.get(o, -1) for o in oids])
            sel = ks >= 0
            if sel.sum() < P.min_epochs:
                continue
            eb_k, ks = np.arange(len(eb))[sel], ks[sel]
            srcs, m90 = ladder(e, role, band)
            cells.append({"endpoint": e, "role": role, "band": band, "m90_v1": m90, "n_epochs": int(sel.sum()), "_centre": centre,
                          "_oids": oids, "_w": w_t, "_eb": eb_k, "_ks": ks, "_d0": d0, "_mjd": mjd, "sources": srcs,
                          "_delta": [dict() for _ in srcs], "_track": [dict() for _ in srcs]})
    cells_ok = [c for c in cells if "sources" in c]
    if not cells_ok:
        return {"corridor": corridor, "band": band, "cells": [c for c in cells if "error" in c], "seconds": _time.monotonic() - t0}
    # ---- node subset: every node within R_NODE (+1 spacing) of any ladder source at any epoch ----
    sub = set()
    for c in cells_ok:
        cx, cy = c["_centre"]; cosd = np.cos(np.deg2rad(cy))
        for s in c["sources"]:
            dt_yr = (c["_mjd"] - P.t0_mjd) / 365.25
            off = c["_d0"][:, s["iz"], :] + np.asarray(s["mu"])[None, :] * dt_yr[:, None]
            sx, sy = gnomonic_xy(ec, cx + off[:, 0] / 3600.0 / cosd, cy + off[:, 1] / 3600.0)
            for x, y in zip(sx, sy):
                if not (np.isfinite(x) and np.isfinite(y)):
                    continue
                x0, x1 = int(np.floor(x - r_node - 1)), int(np.ceil(x + r_node + 1))
                y0, y1 = int(np.floor(y - r_node - 1)), int(np.ceil(y + r_node + 1))
                if x1 < 0 or y1 < 0 or x0 >= side or y0 >= side:
                    continue
                xs = np.arange(max(0, x0), min(side - 1, x1) + 1); ys = np.arange(max(0, y0), min(side - 1, y1) + 1)
                sub.update((ys[:, None] * side + xs[None, :]).ravel().tolist())
    ec.load(np.fromiter(sub, np.int64))
    if not len(ec.f):
        return None
    a0_l, b0_l, N0_l, rchi0_l = ec.fit()
    a0 = np.full(ec.nn, np.nan); b0 = np.full(ec.nn, np.nan); rchi0 = np.full(ec.nn, np.nan)
    a0[ec.nodes] = a0_l; b0[ec.nodes] = b0_l; rchi0[ec.nodes] = rchi0_l
    loaded = ec.local >= 0
    oid_to_k = {o: k for k, o in enumerate(ec.oid)}
    for c in cells_ok:                       # re-map tensor epochs to the loaded epoch order
        ks = np.array([oid_to_k.get(o, -1) for o in c["_oids"]])
        sel = ks >= 0
        c["_eb"], c["_ks"] = np.arange(len(c["_oids"]))[sel], ks[sel]
    # ---- render the injected flux per epoch (one flux map per epoch, all sources) ----
    n_maps = 0
    for k, oid in enumerate(ec.oid):
        needed = [(c, np.where(c["_ks"] == k)[0]) for c in cells_ok]
        needed = [(c, i[0]) for c, i in needed if len(i)]
        if not needed:
            continue
        fm = spx.build_map(oid, obs_by_id[oid], manifest.get(oid), True, None)
        if fm is None:
            continue
        n_maps += 1
        up = int(getattr(fm, "upsample", 1))
        ny, nx = fm.flux.shape
        psf = P.make_psf(band, fm)
        for c, ie in needed:
            eidx = c["_eb"][ie]
            dt_yr = (c["_mjd"][eidx] - P.t0_mjd) / 365.25
            for si, s in enumerate(c["sources"]):
                off = c["_d0"][eidx, s["iz"]] + np.array(s["mu"]) * dt_yr           # arcsec, about the tensor centre (build-stage convention)
                cx, cy = c["_centre"]
                ra_s = cx + off[0] / 3600.0 / np.cos(np.deg2rad(cy)); dec_s = cy + off[1] / 3600.0
                pix = fm.world2pix(np.array([ra_s]), np.array([dec_s]))[0]
                x, y = float(pix[0]), float(pix[1])
                if not (np.isfinite(x) and np.isfinite(y) and 0 <= x < nx and 0 <= y < ny):
                    continue
                amp = 10.0 ** (0.4 * (P.zp_ref - s["mag"]))
                stamp, ox, oy = psf.render(x / up, y / up, half=P.stamp_half)
                rw = P.stamp_response(fm, stamp, ox, oy)
                # nodes within R_NODE of the source
                sx, sy = gnomonic_xy(ec, np.array([ra_s]), np.array([dec_s]))
                sx, sy = float(sx[0]), float(sy[0])
                xs = np.arange(max(0, int(np.floor(sx - r_node))), min(side - 1, int(np.ceil(sx + r_node))) + 1)
                ys = np.arange(max(0, int(np.floor(sy - r_node))), min(side - 1, int(np.ceil(sy + r_node))) + 1)
                near = (ys[:, None] * side + xs[None, :]).ravel()
                near = near[loaded[near]]
                if len(near) == 0:
                    continue
                li = ec.local[near]
                pn = fm.world2pix(ec.RA[li], ec.DEC[li])
                resp = rw.sample(pn[:, 0], pn[:, 1])
                c["_delta"][si][k] = (near, (amp * resp).astype(np.float32))
                imu = (int(np.argmin(np.abs(P.mu_grid - s["mu"][0]))), int(np.argmin(np.abs(P.mu_grid - s["mu"][1]))))
                c["_track"][si][k] = (float(amp * rw.sample(np.array([x]), np.array([y]))[0]), sx, sy,
                                      float(c["_w"][eidx, s["iz"], imu[0], imu[1]]))
        del fm
    # ---- refit per source and measure ----------------------------------------------
    results = []
    for c in cells_ok:
        cres = {"endpoint": c["endpoint"], "role": c["role"], "band": band, "m90_v1": c["m90_v1"],
                "n_epochs_band": c["n_epochs"], "sources": []}
        for si, s in enumerate(c["sources"]):
            dl = c["_delta"][si]; tr = c["_track"][si]
            if len(dl) < P.min_epochs:
                cres["sources"].append({**s, "status": "too_few_epochs_in_map", "n_epochs_on": len(dl)})
                continue
            nodes = np.unique(np.concatenate([v[0] for v in dl.values()]))
            pos = {n: i for i, n in enumerate(nodes)}
            delta = []
            for k in range(len(ec.f)):
                arr = np.zeros(len(nodes), np.float32)
                if k in dl:
                    near, val = dl[k]
                    arr[[pos[n] for n in near]] = val
                delta.append(arr)
            a1, b1, N1, rchi1 = ec.fit(delta=delta, nodes=nodes)
            # scatter to full-grid arrays for the bilinear lookup
            da = np.full(ec.nn, np.nan); db = np.full(ec.nn, np.nan)
            da[nodes] = a1 - a0[nodes]; db[nodes] = b1 - b0[nodes]
            # where the original template is undefined the search subtracts nothing: absorption 0
            undefined = ~np.isfinite(a0[nodes])
            da[nodes[undefined]] = 0.0; db[nodes[undefined]] = 0.0
            # ... and where the refit newly defines a node (cannot happen: N unchanged) keep NaN -> skipped
            num = den = 0.0; per_epoch = []
            for k, (dtrack, sx, sy, w) in tr.items():
                if dtrack <= 0 or w <= 0:      # w: the v2 stack's capped, clipped weight of this epoch at the node
                    continue
                t_inj = bilinear(da, sx, sy, side) + bilinear(db, sx, sy, side) * (ec.lam[k] - ec.lam0)
                if not np.isfinite(t_inj):
                    continue
                num += w * t_inj; den += w * dtrack
                per_epoch.append(t_inj / dtrack)
            if den <= 0 or len(per_epoch) < P.min_epochs:
                cres["sources"].append({**s, "status": "no_weighted_epochs", "n_epochs_on": len(dl)})
                continue
            f_abs = num / den
            rr = rchi1 / rchi0[nodes]; rr = rr[np.isfinite(rr)]
            cres["sources"].append({**s, "status": "ok", "n_epochs_on": len(dl), "n_epochs_weighted": len(per_epoch),
                                    "f_abs": float(f_abs), "f_abs_epoch_median": float(np.median(per_epoch)),
                                    "f_abs_epoch_p90": float(np.quantile(per_epoch, 0.9)),
                                    "delta_mag": float(-2.5 * np.log10(max(1.0 - f_abs, 1e-6))),
                                    "rchi2_ratio_median": float(np.median(rr)) if rr.size else None,
                                    "n_nodes_refit": int(len(nodes))})
        results.append(cres)
    return {"corridor": corridor, "band": band, "n_epochs": len(ec.f), "n_maps_rendered": n_maps,
            "template": {"ra_c": ec.ra_c, "dec_c": ec.dec_c, "side": side, "lam0": ec.lam0},
            "cells": results + [c for c in cells if "error" in c], "seconds": _time.monotonic() - t0}


_G = {}


def _init(obs_by_id, groups, manifest):
    _G.update(obs_by_id=obs_by_id, groups=groups, manifest=manifest)


def _run(key):
    corridor, band = key
    out = OUT / f"{corridor}__{band}.json"
    try:
        res = run_group(corridor, band, _G["groups"][key], _G["obs_by_id"], _G["manifest"])
        if res is None:
            return f"[{corridor}/{band}] no epochs"
        out.write_text(json.dumps(res, indent=1, default=float))
        n_ok = sum(1 for c in res["cells"] for s in c.get("sources", []) if s.get("status") == "ok")
        fa = [s["f_abs"] for c in res["cells"] for s in c.get("sources", []) if s.get("status") == "ok"]
        return (f"[{corridor}/{band}] {len(res['cells'])} cells, {n_ok} sources ok, median f_abs "
                f"{np.median(fa) if fa else float('nan'):.3f} (max {max(fa) if fa else float('nan'):.3f}); {res['seconds']:.0f}s")
    except Exception as exc:
        import traceback
        return f"[{corridor}/{band}] FAILED {exc}\n{traceback.format_exc()}"


def absorb(set_name="dev", workers=8, only_missing=True):
    OUT.mkdir(parents=True, exist_ok=True)
    fz = P.load_freeze()
    if set_name == "all":
        corr = set(P.corridors())
    else:
        corr = set(fz["split"]["development" if set_name == "dev" else "confirmatory"]["corridors"])
    obs_by_id, groups = TC.load_groups()
    _, _, manifest = P.load_inputs()
    keys = sorted(k for k in groups if k[0] in corr)
    if only_missing:
        keys = [k for k in keys if not (OUT / f"{k[0]}__{k[1]}.json").exists()]
    keys.sort(key=lambda k: -len(groups[k]))
    print(f"[absorb/{set_name}] {len(keys)} corridor x detector groups", flush=True)
    log = OUT.parent / "absorb.log"
    if workers <= 1:
        _init(obs_by_id, groups, manifest)
        for k in keys:
            msg = _run(k); print(msg, flush=True); open(log, "a").write(msg + "\n")
        return
    import multiprocessing as mp
    with mp.get_context("fork").Pool(workers, initializer=_init, initargs=(obs_by_id, groups, manifest)) as pool:
        for msg in pool.imap_unordered(_run, keys):
            print(msg, flush=True); open(log, "a").write(msg + "\n")


# -- summary --------------------------------------------------------------------

def summarise():
    fz = P.load_freeze()
    dev = set(fz["split"]["development"]["endpoints"])
    edges = CS.interval_edges(P)
    rows = []
    for p in sorted(OUT.glob("*.json")):
        g = json.loads(p.read_text())
        for c in g["cells"]:
            for s in c.get("sources", []):
                if s.get("status") != "ok":
                    continue
                zi = int(np.clip(np.searchsorted(edges, s["z_au"], side="right") - 1, 0, len(edges) - 2))
                rows.append({"endpoint": c["endpoint"], "role": c["role"], "band": c["band"], "set": "dev" if c["endpoint"] in dev else "confirmatory",
                             "z_interval": zi, "z_au": s["z_au"], "mu": tuple(s["mu"]), "mag_offset": s["mag_offset"],
                             "f_abs": s["f_abs"], "delta_mag": s["delta_mag"], "rchi2_ratio": s.get("rchi2_ratio_median"),
                             "n_epochs_on": s["n_epochs_on"]})
    if not rows:
        raise SystemExit("no absorption products")
    R = rows
    def q(vals, p): return float(np.quantile(vals, p)) if len(vals) else None
    summ = {"n_measurements": len(R), "n_cells": len({(r["endpoint"], r["role"], r["band"]) for r in R}),
            "design": {"R_NODE_ARCSEC": R_NODE_ARCSEC, "MAG_OFFSETS": MAG_OFFSETS, "MU_CHECK": MU_CHECK, "model": MODEL},
            "by_band_zinterval": {}, "by_band": {}, "mag_linearity": {}, "mu_check": {}, "worst_cells": []}
    base = [r for r in R if r["mu"] == (0.0, 0.0) and r["mag_offset"] == 0.0]
    for b in P.bands:
        rb = [r for r in base if r["band"] == b]
        dm = [r["delta_mag"] for r in rb]
        summ["by_band"][b] = {"n": len(rb), "delta_mag_median": q(dm, 0.5), "delta_mag_p90": q(dm, 0.9), "delta_mag_max": max(dm) if dm else None,
                              "f_abs_median": q([r["f_abs"] for r in rb], 0.5)}
        for zi in range(len(edges) - 1):
            rz = [r for r in rb if r["z_interval"] == zi]
            dm = [r["delta_mag"] for r in rz]
            summ["by_band_zinterval"][f"{b}/{zi}"] = {"z_interval_au": [round(float(edges[zi]), 1), round(float(edges[zi + 1]), 1)], "n": len(rz),
                                                     "f_abs_median": q([r["f_abs"] for r in rz], 0.5), "delta_mag_median": q(dm, 0.5),
                                                     "delta_mag_p90": q(dm, 0.9), "delta_mag_max": max(dm) if dm else None}
    for dmo in MAG_OFFSETS:
        rm = [r for r in R if r["mu"] == (0.0, 0.0) and r["mag_offset"] == dmo]
        summ["mag_linearity"][str(dmo)] = {"n": len(rm), "f_abs_median": q([r["f_abs"] for r in rm], 0.5),
                                           "f_abs_p90": q([r["f_abs"] for r in rm], 0.9), "rchi2_ratio_median": q([r["rchi2_ratio"] for r in rm if r["rchi2_ratio"] is not None], 0.5)}
    for mu in MU_CHECK:
        rm = [r for r in R if r["mu"] == mu]
        ref = [r for r in base if r["z_interval"] >= len(edges) - 1 - N_MU_CHECK_Z]
        summ["mu_check"][str(mu)] = {"n": len(rm), "f_abs_median": q([r["f_abs"] for r in rm], 0.5),
                                     "f_abs_median_mu0_same_z": q([r["f_abs"] for r in ref], 0.5)}
    worst = sorted(base, key=lambda r: -r["delta_mag"])[:12]
    summ["worst_cells"] = [{k: r[k] for k in ("endpoint", "role", "band", "z_au", "f_abs", "delta_mag", "n_epochs_on", "set")} for r in worst]
    for s in ("dev", "confirmatory"):
        rs = [r for r in base if r["set"] == s]
        summ[f"{s}_delta_mag"] = {"n": len(rs), "median": q([r["delta_mag"] for r in rs], 0.5), "p90": q([r["delta_mag"] for r in rs], 0.9)}
    res_dir = P.survey_dir / "results"
    res_dir.mkdir(parents=True, exist_ok=True)
    (res_dir / "template_absorption.json").write_text(json.dumps(summ, indent=1, default=float))
    lines = ["---", 'title: "SPHEREx template-absorption control (generated by template_absorb.summarise)"', "---", "",
             f"{summ['n_measurements']} measurements on {summ['n_cells']} cells (persistent source on the real track, "
             f"template refitted with the source present). Δm = −2.5 log10(1 − f_abs) is the depth overstatement of the "
             f"injection-calibrated m90 for a slow persistent source; corrected limit = m90 − Δm.", "",
             "| detector | z interval (AU) | n | median f_abs | median Δm | p90 Δm | max Δm |", "| --- | --- | --- | --- | --- | --- | --- |"]
    for b in P.bands:
        for zi in range(len(edges) - 1):
            r = summ["by_band_zinterval"][f"{b}/{zi}"]
            f = lambda x, nd=3: "—" if x is None else f"{x:.{nd}f}"  # noqa: E731
            lines.append(f"| {b} | {r['z_interval_au'][0]}–{r['z_interval_au'][1]} | {r['n']} | {f(r['f_abs_median'])} | {f(r['delta_mag_median'])} | {f(r['delta_mag_p90'])} | {f(r['delta_mag_max'])} |")
    lines += ["", "Magnitude linearity (µ = 0, all z; f_abs should not depend on the injected magnitude below the clip):", "",
              "| mag − m90 | n | median f_abs | p90 f_abs | median rchi2 ratio at refit nodes |", "| --- | --- | --- | --- | --- |"]
    for k, v in summ["mag_linearity"].items():
        lines.append(f"| {k} | {v['n']} | {v['f_abs_median']:.4f} | {v['f_abs_p90']:.4f} | {v['rchi2_ratio_median']:.3f} |")
    lines += ["", "Residual-motion check (two most distant z intervals):", "", "| µ (\"/yr) | n | median f_abs | µ = 0 same z |", "| --- | --- | --- | --- |"]
    for k, v in summ["mu_check"].items():
        lines.append(f"| {k} | {v['n']} | {v['f_abs_median']:.4f} | {v['f_abs_median_mu0_same_z']:.4f} |")
    lines += ["", "Largest Δm (µ = 0, m90):", "", "| cell | set | z (AU) | f_abs | Δm | epochs on |", "| --- | --- | --- | --- | --- | --- |"]
    for w in summ["worst_cells"]:
        lines.append(f"| {w['endpoint']}/{w['role']}/{w['band']} | {w['set']} | {w['z_au']:.0f} | {w['f_abs']:.3f} | {w['delta_mag']:.3f} | {w['n_epochs_on']} |")
    (res_dir / "template_absorption.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


# -- correction of the frozen limits (report-level) ------------------------------

def correct():
    """Per cell x z interval: the injection-calibrated m90 of the frozen
    per-detector (runs/spherex/v2) and joint (runs/spherex/joint6)
    Constraint records, and the template-absorption-corrected limit
    m90 - dm, with dm from the ladder source at the cell's m90 and mu = 0
    in that z interval (joint cell: mean f_abs over the six detectors).
    Records are not modified; writes surveys/spherex/results/
    template_absorption_corrected.{json,md}."""
    from sglsurvey.records import read_records
    edges = CS.interval_edges(P)
    fabs = {}     # (endpoint, role, band, zi) -> f_abs
    for p in sorted(OUT.glob("*.json")):
        g = json.loads(p.read_text())
        for c in g["cells"]:
            for s in c.get("sources", []):
                if s.get("status") == "ok" and tuple(s["mu"]) == (0.0, 0.0) and s["mag_offset"] == 0.0:
                    zi = int(np.clip(np.searchsorted(edges, s["z_au"], side="right") - 1, 0, len(edges) - 2))
                    fabs[(c["endpoint"], c["role"], c["band"], zi)] = s["f_abs"]

    def dm_for(e, role, band, zi):
        if band == "J6":
            v = [fabs.get((e, role, b, zi)) for b in P.bands]
            v = [x for x in v if x is not None]
            f = float(np.mean(v)) if v else None
        else:
            f = fabs.get((e, role, band, zi))
        return (None if f is None else float(-2.5 * np.log10(max(1.0 - f, 1e-6)))), f

    out = {"note": "report-level correction; frozen Constraint records unchanged", "families": {}}
    lines = ["---", 'title: "SPHEREx template-absorption-corrected limits (generated by template_absorb.correct)"', "---", "",
             "Frozen injection-calibrated m90 (final-candidate curve, worst temporal model) and the corrected limit "
             "m90 − Δm for a slow source (Δm from the ladder at the cell's m90, µ = 0; joint cell: mean f_abs over "
             "the six detectors). Records unchanged.", ""]
    for fam, rec_path in (("per-detector", REPO / "runs" / "spherex" / "v2" / "records" / "constraint.jsonl"),
                          ("joint6", REPO / "runs" / "spherex" / "joint6" / "records" / "constraint.jsonl")):
        if not rec_path.exists():
            continue
        cons = [c for c in read_records(rec_path) if c["completeness_kind"] == "final_candidate" and c["kind"] == "recovery_curve"]
        rows = []
        for c in cons:
            zi = int(np.argmin([abs(c["z_interval_au"][0] - e) for e in edges[:-1]]))
            dm, f = dm_for(c["endpoint_id"], c["role"], c["band"], zi)
            fl = c["flux_limit"]
            rows.append({"cell": f"{c['endpoint_id']}/{c['role']}/{c['band']}", "set": c["extra"]["set"], "band": c["band"], "z_interval": zi,
                         "z_interval_au": c["z_interval_au"], "m90": fl["value"], "worst_model": fl.get("worst_temporal_model"),
                         "m90_persistent": (fl.get("per_model_m90") or {}).get("persistent"),
                         "f_abs": f, "delta_mag": dm, "m90_corrected": (fl["value"] - dm) if dm is not None else None,
                         "m90_persistent_corrected": ((fl.get("per_model_m90") or {}).get("persistent") - dm)
                         if (dm is not None and (fl.get("per_model_m90") or {}).get("persistent") is not None) else None})
        bands = list(P.bands) if fam == "per-detector" else ["J6"]
        summ = {}
        for s in ("confirmatory", "dev"):
            for b in bands:
                rb = [r for r in rows if r["set"] == s and r["band"] == b and r["delta_mag"] is not None]
                if not rb:
                    continue
                summ[f"{s}/{b}"] = {"n": len(rb), "median_m90": float(np.median([r["m90"] for r in rb])),
                                    "median_m90_corrected": float(np.median([r["m90_corrected"] for r in rb])),
                                    "median_delta_mag": float(np.median([r["delta_mag"] for r in rb])),
                                    "p90_delta_mag": float(np.quantile([r["delta_mag"] for r in rb], 0.9)),
                                    "median_m90_persistent": float(np.median([r["m90_persistent"] for r in rb if r["m90_persistent"] is not None])) if any(r["m90_persistent"] is not None for r in rb) else None,
                                    "median_m90_persistent_corrected": float(np.median([r["m90_persistent_corrected"] for r in rb if r["m90_persistent_corrected"] is not None])) if any(r["m90_persistent_corrected"] is not None for r in rb) else None,
                                    "n_uncorrected": sum(1 for r in rows if r["set"] == s and r["band"] == b and r["delta_mag"] is None)}
        out["families"][fam] = {"summary": summ, "rows": rows}
        lines += [f"## {fam}", "", "| set | band | n | median m90 (worst model) | corrected | median m90 (persistent) | corrected | median Δm | p90 Δm |",
                  "| --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
        f = lambda x, nd=2: "—" if x is None else f"{x:.{nd}f}"  # noqa: E731
        for k, v in summ.items():
            s_, b = k.split("/")
            lines.append(f"| {s_} | {b} | {v['n']} | {f(v['median_m90'])} | {f(v['median_m90_corrected'])} | {f(v['median_m90_persistent'])} | "
                         f"{f(v['median_m90_persistent_corrected'])} | {f(v['median_delta_mag'])} | {f(v['p90_delta_mag'])} |")
        lines.append("")
    res_dir = P.survey_dir / "results"
    (res_dir / "template_absorption_corrected.json").write_text(json.dumps(out, indent=1, default=float))
    (res_dir / "template_absorption_corrected.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["run", "summarise", "correct"])
    ap.add_argument("--set", choices=["dev", "confirmatory", "all"], default="dev")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--no-only-missing", action="store_true")
    a = ap.parse_args()
    if a.mode == "run":
        absorb(a.set, a.workers, not a.no_only_missing)
    elif a.mode == "summarise":
        summarise()
    else:
        correct()
