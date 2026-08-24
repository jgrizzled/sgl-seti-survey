"""v2 engine: image-level injections (profile-parameterised port of
surveys/wise-v2/scripts/inject_pipeline.py). Differences: the PSF
renderer and focal-plane element come from the profile (empirical PRF
grid or Moffat at the frame seeing), flux maps from profile.build_map,
weights use the tensor's per-frame cap, the catalogue from
profile.catalog, magnitudes in the profile's system."""

from __future__ import annotations

import hashlib
import json
import time as _time
import warnings

import numpy as np

warnings.filterwarnings("ignore", message="All-NaN slice")

from sglsurvey import inject, nulls  # noqa: E402
from sglsurvey.photometry import _gaussian_kernel  # noqa: E402
from sglsurvey.vetting import flux_consistency, radial_response_table  # noqa: E402

MASKS = ("primary", "strict", "loose")

Z_WINDOW = 8
DEFAULT_M90 = {"W1": 14.0, "W2": 13.0, "W3": 10.5, "W4": 7.0}
SCHEMA = "v2-injections-v2"
_G: dict = {}


def _v1_m90(P, endpoint, role, band):
    m = P.v1_m90(endpoint, role, band) if P.v1_m90 else None
    return m if m is not None else P.default_m90.get(band, 20.0)


def _envelope_xt(P, endpoint, role):
    env = _G.get("env")
    if env is None:
        p = P.run_dir / "geometry" / "covariance_envelopes.json"
        env = {}
        if p.exists():
            for e in json.loads(p.read_text())["envelopes"]:
                env[(e["endpoint"], e["role"])] = e.get("max_sigma_xt_99") or 0.0
        _G["env"] = env
    return float(env.get((endpoint, role), 0.0))


def oversampled_template(psf, n=641, oversample=8):
    """8x-oversampled unit-sum template of a PSF renderer: the PRF grid's
    mean template, or a synthesised Moffat."""
    if hasattr(psf, "mean_template"):
        return psf.mean_template()
    c = (n - 1) // 2
    yy, xx = np.mgrid[-c:c + 1, -c:c + 1] / oversample
    r2 = xx ** 2 + yy ** 2
    t = (1.0 + r2 / psf.alpha ** 2) ** (-psf.beta)
    return t / t.sum()


def _response_table(P, band, psf, fwhm_arcsec, pix_arcsec):
    key = ("resp", band, round(fwhm_arcsec, 2), round(pix_arcsec, 3))
    if key not in _G:
        fwhm_pix = fwhm_arcsec / pix_arcsec
        k = _gaussian_kernel(fwhm_pix, int(np.ceil(2.5 * fwhm_pix)))
        r, v = radial_response_table(oversampled_template(psf), k)
        _G[key] = (r * pix_arcsec, v)
    return _G[key]


def temporal_on(model, mjd, rng, visit_gap_days=5.0):
    E = len(mjd)
    if model == "persistent":
        return np.ones(E, bool)
    if model == "flicker":
        return rng.random(E) < 0.5
    if model == "visit":
        order = np.argsort(mjd)
        gaps = np.diff(mjd[order]) > visit_gap_days
        vid = np.concatenate([[0], np.cumsum(gaps)])
        on_v = rng.random(vid.max() + 1) < 0.5
        on = np.empty(E, bool)
        on[order] = on_v[vid]
        return on
    if model == "block":
        span = mjd.max() - mjd.min()
        start = mjd.min() + rng.random() * 0.5 * span
        return (mjd >= start) & (mjd <= start + 0.5 * span)
    raise ValueError(model)


def inject_cell(P, d, band, maps, n_inj, rng, catalog):
    """All injections of one cell. ``maps`` is the band's flux-map list
    aligned with the tensor's epochs of that band (None where the
    cutout could not be rebuilt)."""
    endpoint, role = str(d["endpoint"]), str(d["role"])
    b = P.band_idx[band]
    eb = np.where(d["band_idx"] == b)[0]
    if len(eb) < P.min_epochs:
        return None
    fcap = d["frame_cap"][eb]
    f0 = d["f"][0, eb]; v0 = d["v"][0, eb]; g0 = d["g"][0, eb]
    mjd = d["mjd"][eb]; phase = d["phase"][eb]; masks = d["masks"][eb]
    xy0 = d["xy0"][eb]; d0 = d["d0"][eb]; J = d["J"][eb]; p0 = d["p0"][eb]
    origin = d["origin"][eb]; magzp = d["magzp"][eb]
    centre = tuple(d["centre"])
    nz, nm = len(P.z_grid), len(P.mu_grid)
    q_grid = 1.0 / P.z_grid
    E = len(eb)
    psf_cache = {}
    pix_arcsec = float(np.nanmedian(d["pix_arcsec"][eb]))
    fwhm_med = float(np.nanmedian(d["fwhm_arcsec"][eb]))
    # Real S per mask (full grid) and the capped weights, which do not
    # depend on the injection (v, g unchanged): precomputed once for the
    # full, early and late epoch selections so that each injection is a
    # weighted sum over its window.
    real = {}
    f0n = np.nan_to_num(f0.astype(np.float32))
    for mi, m in enumerate(MASKS):
        ok = masks[:, mi]
        if ok.sum() < P.min_epochs:
            real[m] = None
            continue
        S, A, B, n, w = nulls.stack_S(f0, v0, g0, ok, return_parts=True, frame_cap=fcap, clip_sigma=P.clip_sigma)
        split = P.holdout_split_mjd if P.holdout_split_mjd is not None else np.inf
        early = ok & (mjd <= split)
        late = ok & (mjd > split)
        w_e = nulls.stack_weights(v0, g0, early, f=f0, frame_cap=fcap, clip_sigma=P.clip_sigma)[0] if (early.sum() >= P.min_epochs and np.isfinite(split)) else None
        w_l = nulls.stack_weights(v0, g0, late, f=f0, frame_cap=fcap, clip_sigma=P.clip_sigma)[0] if late.sum() >= P.min_epochs else None
        S_e_real = None
        if w_e is not None:
            with np.errstate(all="ignore"):
                Be = w_e.sum(0); S_e_real = (f0n * w_e).sum(0) / np.sqrt(np.where(Be > 0, Be, np.inf))
            S_e_real = np.where((w_e > 0).sum(0) >= P.min_epochs, S_e_real, np.nan)
        real[m] = {"S": S, "ok": ok, "w": w, "early": early, "late": late,
                   "w_e": w_e, "w_l": w_l, "S_e_real": S_e_real}
    T = {m: float(np.nanmax(d["summary"][mi, 1 + np.asarray(d["designated"]), b - 1, 0]))
         for mi, m in enumerate(MASKS)}
    m90 = _v1_m90(P, endpoint, role, band)
    sig_xt = _envelope_xt(P, endpoint, role)
    models = list(P.temporal_models)
    out = []
    dt_yr = (mjd - P.t0_mjd) / 365.25
    cosd = np.cos(np.deg2rad(centre[1]))
    amp_ref = None
    lo_m, hi_m = (P.inj_window(endpoint, role, band) if getattr(P, "inj_window", None)
                  else (m90 - P.inj_mag_halfwidth, m90 + P.inj_mag_halfwidth))
    seed0 = int(hashlib.sha256(f"inj/{endpoint}/{role}/{P.split_seed}".encode()).hexdigest()[:8], 16)
    win_store = {"A": [], "B": [], "n": [], "lo": []}
    for j in range(n_inj):
        model = models[j % len(models)]
        # per-injection seeded draws: the same (z, mu, xt, mag) for the
        # j-th injection of this cell in every archive (joint stage)
        rj = np.random.default_rng([seed0, j])
        z = float(np.exp(rj.uniform(np.log(550.0), np.log(10000.0))))
        mu = rj.uniform(-1.0, 1.0, size=2)
        xt = float(rj.uniform(-1.0, 1.0)) * sig_xt
        mag = float(rj.uniform(lo_m, hi_m))
        on = temporal_on(model, mjd, np.random.default_rng([seed0, j, 1]), P.visit_gap_days)
        amp_ref = 10.0 ** (0.4 * (P.zp_ref - mag))          # tensor units (ZP_REF)
        q = 1.0 / z
        iz = int(np.argmin(np.abs(q_grid - q)))
        lo, hi = max(0, iz - Z_WINDOW), min(nz, iz + Z_WINDOW + 1)
        win = np.arange(lo, hi)
        # per-epoch track offsets (arcsec) at the injection: interpolate d0 in q
        # (d0 is (E, nz, 2), q_grid decreasing in z -> increasing in q order)
        qo = np.argsort(q_grid)
        dz = np.stack([np.array([np.interp(q, q_grid[qo], d0[e, qo, k]) for e in range(E)])
                       for k in (0, 1)], axis=1)                       # (E, 2)
        # local tangent of the locus at iz (for the cross-track direction)
        i1, i2 = max(iz - 1, 0), min(iz + 1, nz - 1)
        tvec = d0[:, i2, :] - d0[:, i1, :]
        nrm = np.linalg.norm(tvec, axis=1, keepdims=True)
        tvec = np.where(nrm > 0, tvec / np.where(nrm > 0, nrm, 1), np.array([[1.0, 0.0]]))
        nvec = np.stack([-tvec[:, 1], tvec[:, 0]], axis=1)
        track = dz + np.stack([mu[0] * dt_yr, mu[1] * dt_yr], axis=1) + xt * nvec   # (E, 2)
        off = np.empty((E, len(win), nm, nm, 2))
        off[..., 0] = d0[:, win, 0][:, :, None, None] + (P.mu_grid * dt_yr[:, None])[:, None, :, None]
        off[..., 1] = d0[:, win, 1][:, :, None, None] + (P.mu_grid * dt_yr[:, None])[:, None, None, :]
        if P.linear_wcs:
            pix_inj = p0 + np.einsum("eij,ej->ei", J, track)                            # (E, 2)
            pix_nodes = p0[:, None, None, None, :] + np.einsum("eij,e...j->e...i", J, off)
        else:
            pix_inj = np.full((E, 2), np.nan); pix_nodes = np.full(off.shape, np.nan)
            for e in range(E):
                fm_ = maps[e]
                if fm_ is None:
                    continue
                ra_q = centre[0] + np.concatenate([[track[e, 0]], off[e, ..., 0].ravel()]) / 3600.0 / cosd
                dec_q = centre[1] + np.concatenate([[track[e, 1]], off[e, ..., 1].ravel()]) / 3600.0
                pp = fm_.world2pix(ra_q, dec_q)
                pix_inj[e] = pp[0]; pix_nodes[e] = pp[1:].reshape(off.shape[1:])
        delta = np.zeros((E, len(win), nm, nm), np.float32)
        resp0 = np.full(E, np.nan)
        for e in range(E):
            fm = maps[e]
            if fm is None or not on[e]:
                continue
            x, y = float(pix_inj[e, 0]), float(pix_inj[e, 1])
            ny, nx = fm.flux.shape
            if not (0 <= x < nx and 0 <= y < ny):
                continue
            psf = psf_cache.get(e)
            if psf is None:
                psf = psf_cache[e] = P.make_psf(band, fm)
            up = int(getattr(fm, "upsample", 1))
            el = P.psf_element(fm, x / up, y / up) if P.psf_element else None
            stamp, ox, oy = psf.render(x / up, y / up, element=el, half=P.stamp_half)
            rw = (P.stamp_response or inject.stamp_response)(fm, stamp, ox, oy)
            resp0[e] = rw.sample(x, y)
            delta[e] = (amp_ref * rw.sample(pix_nodes[e, ..., 0], pix_nodes[e, ..., 1])).astype(np.float32)
        f_inj = f0n[:, win] + delta
        rec = {"j": j, "model": model, "z_au": z, "mu": mu.tolist(), "xt_arcsec": xt, "mag": mag,
               "fnu_jy": P.fnu_from_mag(band, mag),
               "iz": iz, "n_on": int(on.sum()),
               "resp_median": float(np.nanmedian(resp0)) if np.isfinite(resp0).any() else None,
               "masks": {}}
        # distance of every window node from the injection, in nodes
        wz = np.abs(win - iz)[:, None, None]
        wmu = np.maximum(np.abs(np.arange(nm)[:, None] - np.argmin(np.abs(P.mu_grid - mu[0]))),
                         np.abs(np.arange(nm)[None, :] - np.argmin(np.abs(P.mu_grid - mu[1]))))[None]
        near = (wz <= P.recovery_window[0]) & (wmu <= P.recovery_window[1])
        for mi, m in enumerate(MASKS):
            if real[m] is None:
                continue
            R_ = real[m]
            S_real, ok = R_["S"], R_["ok"]
            w_w = R_["w"][:, win]
            A_w = (f_inj * w_w).sum(0); B_w = w_w.sum(0); n_w = (w_w > 0).sum(0)
            if m == "primary":
                pad = 2 * Z_WINDOW + 1 - len(win)
                win_store["A"].append(np.pad(A_w, ((0, pad), (0, 0), (0, 0))).astype(np.float32))
                win_store["B"].append(np.pad(B_w, ((0, pad), (0, 0), (0, 0))).astype(np.float32))
                win_store["n"].append(np.pad(n_w, ((0, pad), (0, 0), (0, 0))).astype(np.int16))
                win_store["lo"].append(lo)
            with np.errstate(all="ignore"):
                S_w = A_w / np.sqrt(np.where(B_w > 0, B_w, np.inf))
            S_w = np.where(n_w >= P.min_epochs, S_w, np.nan)
            S_full = S_real.copy()
            S_full[win] = S_w
            s_max, node = nulls.grid_max(S_full)
            S_near = np.where(near, S_w, np.nan)
            s_peak, pk = nulls.grid_max(S_near)
            r = {"S_peak": s_peak, "S_max": s_max, "T": T[m],
                 "R_peak": (s_peak / T[m]) if T[m] > 0 else None,
                 "argmax_in_window": bool(node is not None and lo <= node[0] < hi and near[node[0] - lo, node[1], node[2]])}
            if pk is None:
                rec["masks"][m] = r
                continue
            pz, pm1, pm2 = pk
            fe = f_inj[:, pz, pm1, pm2].astype(np.float64)
            we = w_w[:, pz, pm1, pm2]
            contrib = fe * we
            Ap = contrib.sum(); Bp = we.sum()
            r["top_share"] = float(contrib.max() / Ap) if Ap > 0 else None
            S_ph = {}
            for p in (0, 1):
                sel = ok & (phase == p)
                Bq = we[sel].sum()
                S_ph[p] = float(contrib[sel].sum() / np.sqrt(Bq)) if Bq > 0 else np.nan
            r["S_phase"] = [S_ph[0], S_ph[1]]
            # held-out-epoch refit on the early epochs (window + real outside)
            hold = None
            if R_["w_e"] is not None:
                w_e = R_["w_e"][:, win]
                A_e = (f_inj * w_e).sum(0); B_e = w_e.sum(0); n_e = (w_e > 0).sum(0)
                with np.errstate(all="ignore"):
                    S_e = A_e / np.sqrt(np.where(B_e > 0, B_e, np.inf))
                S_e = np.where(n_e >= P.min_epochs, S_e, np.nan)
                S_e_full = R_["S_e_real"].copy(); S_e_full[win] = S_e
                s_em, nd = nulls.grid_max(S_e_full)
                if nd is not None and R_["w_l"] is not None:
                    inwin = lo <= nd[0] < hi
                    fl = f_inj[:, nd[0] - lo, nd[1], nd[2]] if inwin else f0n[:, nd[0], nd[1], nd[2]]
                    wl = R_["w_l"][:, nd[0], nd[1], nd[2]]
                    A_l = float((fl * wl).sum()); B_l = float(wl.sum()); n_l = int((wl > 0).sum())
                    fe_ = (A_e[nd[0] - lo, nd[1], nd[2]] / B_e[nd[0] - lo, nd[1], nd[2]]
                           if inwin and B_e[nd[0] - lo, nd[1], nd[2]] > 0 else np.nan)
                    hold = {"S_early": float(s_em), "f_early": float(fe_),
                            "S_late": float(A_l / np.sqrt(B_l)) if B_l > 0 else None,
                            "f_late": float(A_l / B_l) if B_l > 0 else None,
                            "n_late": n_l, "node_in_window": bool(inwin)}
            r["holdout"] = hold
            # flux-consistent static-source test at the peak node
            cat_b = catalog.get(band) if isinstance(catalog, dict) else catalog
            if cat_b is not None and len(cat_b.ra):
                tr = np.empty((E, 2))
                tr[:, 0] = d0[:, win[pz], 0] + P.mu_grid[pm1] * dt_yr
                tr[:, 1] = d0[:, win[pz], 1] + P.mu_grid[pm2] * dt_yr
                ra_t = centre[0] + tr[:, 0] / 3600.0 / cosd
                dec_t = centre[1] + tr[:, 1] / 3600.0
                psf0 = next((p for p in psf_cache.values()), None) or P.make_psf(band, next(m for m in maps if m is not None))
                resp_table = _response_table(P, band, psf0, fwhm_med, pix_arcsec)
                esc = P.static_epoch_scale(d, eb, (win[pz], pm1, pm2)) if P.static_epoch_scale else None
                fc = flux_consistency(cat_b, ra_t, dec_t, we * ok, phase, P.zp_ref,
                                      resp_table[0], resp_table[1], s_peak,
                                      S_by_phase=S_ph, search_arcsec=P.static_search_arcsec,
                                      factor=P.static_factor, min_ndet=P.catalog_min_ndet,
                                      epoch_scale=esc)
                r["static"] = {"S_pred": fc["S_pred"], "consistent": fc["consistent"],
                               "n_sources": fc["n_sources_considered"],
                               "nearest": (fc["sources"][0] if fc["sources"] else None)}
            rec["masks"][m] = r
        out.append(rec)
    cell = {"endpoint": endpoint, "role": role, "band": band, "n_epochs": int(E),
            "m90_v1": m90, "mag_window": [lo_m, hi_m], "sigma_xt_99": sig_xt, "T": T,
            "amp_note": f"tensor units ZP {P.zp_ref}", "injections": out}
    cell["_windows"] = {k: np.array(v) for k, v in win_store.items()}
    return cell


def run_corridor(args):
    corridor, n_inj, only_missing = args
    P = _G["P"]
    t0 = _time.monotonic()
    obs_by_id, manifest = _G["obs_by_id"], _G["manifest"]
    pairs = []
    for e in P.members(corridor):
        for r in ("rx", "tx"):
            for tp in [P.tensor_dir / f"{e}__{r}.npz"] + sorted(P.tensor_dir.glob(f"{e}__{r}__xt*.npz")):
                if tp.exists():
                    pairs.append((e, r, tp))
    P.inj_dir.mkdir(parents=True, exist_ok=True)
    done = []
    for endpoint, role, tpath in pairs:
        out_path = P.inj_dir / (tpath.stem + ".json")
        if only_missing and out_path.exists():
            continue
        d = np.load(tpath)
        centre = tuple(d["centre"])
        catalog = P.catalog(corridor, centre) if P.catalog else None
        result = {"schema": SCHEMA, "survey": P.name, "endpoint": endpoint, "role": role, "corridor": corridor,
                  "xt_arcsec": float(d["xt_arcsec"]) if "xt_arcsec" in d.files else 0.0,
                  "tensor_input_hash": str(d["input_hash"]), "freeze_hash": P.freeze_hash(),
                  "n_per_cell": n_inj,
                  "catalog": (None if catalog is None else (catalog.label if not isinstance(catalog, dict)
                                                          else {k: (v.label if v else None) for k, v in catalog.items()})),
                  "cells": {}}
        seed = int(hashlib.sha256(f"inj/{endpoint}/{role}/{P.split_seed}".encode()).hexdigest()[:8], 16)
        rng = np.random.default_rng(seed)
        result["seed"] = seed
        for band in P.bands:
            b = P.band_idx[band]
            eb = np.where(d["band_idx"] == b)[0]
            if len(eb) < P.min_epochs:
                continue
            maps = []
            for e in eb:
                oid = str(d["oid"][e])
                maps.append(P.build_map(oid, obs_by_id[oid], manifest.get(oid), True, None))
            cell = inject_cell(P, d, band, maps, n_inj, rng, catalog)
            del maps
            if cell:
                wins = cell.pop("_windows")
                np.savez_compressed(P.inj_dir / f"{tpath.stem}__{band}_windows.npz",
                                    A=wins["A"], B=wins["B"], n=wins["n"], lo=wins["lo"],
                                    z_grid=P.z_grid, mu_grid=P.mu_grid, z_window=Z_WINDOW)
                result["cells"][band] = cell
        d.close()
        out_path.write_text(json.dumps(result, default=_jsonable))
        done.append(tpath.stem)
    return f"{corridor}: {len(done)} pairs in {(_time.monotonic() - t0) / 60:.1f} min"


def _jsonable(x):
    if isinstance(x, (np.floating, np.integer)):
        v = x.item()
        return None if isinstance(v, float) and not np.isfinite(v) else v
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, float) and not np.isfinite(x):
        return None
    raise TypeError(type(x))


def _init(P, obs_by_id, manifest):
    _G["P"] = P
    _G["obs_by_id"] = obs_by_id
    _G["manifest"] = manifest


def run(P, set_name, workers=7, corridors=None, only_missing=True, n_per_cell=None):
    n_inj = n_per_cell or P.n_inj_per_cell
    freeze = P.load_freeze()
    if corridors:
        pass
    elif set_name == "all":
        corridors = P.corridors()
    else:
        corridors = freeze["split"]["development" if set_name == "dev" else "confirmatory"]["corridors"]
    obs_by_id, usable, manifest = P.load_inputs()
    # warm the catalogue cache serially (network / snapshot reads) before forking
    for c in corridors:
        for e in P.members(c):
            p = P.tensor_dir / f"{e}__rx.npz"
            if p.exists() and P.catalog:
                with np.load(p) as d:
                    P.catalog(c, tuple(d["centre"]))
                break
    jobs = [(c, n_inj, only_missing) for c in corridors]
    log = P.run_dir / "inject_pipeline.log"
    if workers <= 1:
        _init(P, obs_by_id, manifest)
        for j in jobs:
            msg = run_corridor(j); print(msg, flush=True); open(log, "a").write(msg + "\n")
        return
    import multiprocessing as mp
    with mp.get_context("fork").Pool(workers, initializer=_init, initargs=(P, obs_by_id, manifest)) as pool:
        for msg in pool.imap_unordered(run_corridor, jobs):
            print(msg, flush=True); open(log, "a").write(msg + "\n")
