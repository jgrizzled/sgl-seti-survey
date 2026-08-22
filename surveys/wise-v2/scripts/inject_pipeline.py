"""Step E / G: image-level injections (hypotheses v2.0 §4; v2 plan §4).

Per cell (endpoint x role x band) N_INJ_PER_CELL sources are drawn —
z log-uniform, (mu_a, mu_d) uniform on the L-inf box, cross-track
offset uniform within the propagated 99 % envelope, Vega magnitude
uniform on [m90_v1 - 2, m90_v1 + 2], temporal model from the four
families — and each is added to every cutout of the band as the
frame's empirical PRF template (focal-plane element from the frame
position, sub-pixel phase exact) in DN through the frame MAGZP, before
the matched filter: exactly, via the stamp-response linearity of the
estimator (sglsurvey.inject.stamp_response; invariant-tested against
full re-convolution). The injected track's response is added to the
real-trajectory tensor on a window of z-nodes around the injection,
the window is stacked under each quality mask, and everything the
decision rule and the calibrated vetoes need is recorded per injection:
the peak S near the injection, the full-grid maximum, the argmax node,
the phase split, the top-epoch share, the held-out-epoch refit /
prediction numbers, and the flux-consistent static-source test
against the corridor's CatWISE2020 catalogue.

Recovery decisions (threshold and final-candidate completeness) are
made afterwards by completeness.py with the set's frozen R~_FWER.

Usage: uv run python surveys/wise-v2/scripts/inject_pipeline.py
           --set dev|confirmatory|all [--workers 7] [--only-missing]
           [--n-per-cell 400] [--corridors ...]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time as _time
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v2common as C  # noqa: E402

warnings.filterwarnings("ignore", message="All-NaN slice")

from sglsurvey import inject, nulls  # noqa: E402
from sglsurvey.adapters.irsa_wise import WiseExactFootprint  # noqa: E402
from sglsurvey.corridors import CORRIDOR_OF, MEMBERS  # noqa: E402
from sglsurvey.photometry import _gaussian_kernel, build_flux_map  # noqa: E402
from sglsurvey.records import read_records  # noqa: E402
from sglsurvey.vetting import (flux_consistency, load_catwise_vizier,  # noqa: E402
                               radial_response_table)

MASKS = ("primary", "strict", "loose")
Z_WINDOW = 8
DEFAULT_M90 = {"W1": 14.0, "W2": 13.0, "W3": 10.5, "W4": 7.0}
SCHEMA = "wise-v2-injections-v1"
_G: dict = {}


def _v1_m90(endpoint, role, band):
    d = _G.setdefault("m90", np.load(C.V1_CAL / "m90_curves.npz"))
    key = f"{endpoint}__{role}__{band}__0.5"
    if key in d.files:
        arr = d[key]
        if np.isfinite(arr).any():
            return float(np.nanmedian(arr))
    return DEFAULT_M90[band]


def _envelope_xt(endpoint, role):
    env = _G.get("env")
    if env is None:
        p = C.RUN_DIR / "geometry" / "covariance_envelopes.json"
        env = {}
        if p.exists():
            for e in json.loads(p.read_text())["envelopes"]:
                env[(e["endpoint"], e["role"])] = e.get("max_sigma_xt_99") or 0.0
        _G["env"] = env
    return float(env.get((endpoint, role), 0.0))


def _catalog(corridor, centre):
    key = ("cat", corridor)
    if key not in _G:
        try:
            _G[key] = load_catwise_vizier(centre[0], centre[1], 0.1, C.RUN_DIR / "catwise")
        except Exception as exc:
            print(f"  [{corridor}] CatWISE unavailable: {exc}", flush=True)
            _G[key] = None
    return _G[key]


def _prf(band):
    key = ("prf", band)
    if key not in _G:
        _G[key] = inject.PRFGrid.load(band, C.PRF_DIR / band.lower())
    return _G[key]


def _response_table(band, pix_arcsec):
    key = ("resp", band, round(pix_arcsec, 3))
    if key not in _G:
        prf = _prf(band)
        fwhm_pix = C.PSF_FWHM[band] / pix_arcsec
        k = _gaussian_kernel(fwhm_pix, int(np.ceil(2.5 * fwhm_pix)))
        r, v = radial_response_table(prf.mean_template(), k)
        _G[key] = (r * pix_arcsec, v)
    return _G[key]


def temporal_on(model, mjd, rng):
    E = len(mjd)
    if model == "persistent":
        return np.ones(E, bool)
    if model == "flicker":
        return rng.random(E) < 0.5
    if model == "visit":
        order = np.argsort(mjd)
        gaps = np.diff(mjd[order]) > C.VISIT_GAP_DAYS
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


def inject_cell(d, band, maps, n_inj, rng, catalog, resp_table):
    """All injections of one cell. ``maps`` is the band's flux-map list
    aligned with the tensor's epochs of that band (None where the
    cutout could not be rebuilt)."""
    endpoint, role = str(d["endpoint"]), str(d["role"])
    b = C.BAND_IDX[band]
    eb = np.where(d["band_idx"] == b)[0]
    if len(eb) < C.MIN_EPOCHS:
        return None
    f0 = d["f"][0, eb]; v0 = d["v"][0, eb]; g0 = d["g"][0, eb]
    mjd = d["mjd"][eb]; phase = d["phase"][eb]; masks = d["masks"][eb]
    xy0 = d["xy0"][eb]; d0 = d["d0"][eb]; J = d["J"][eb]; p0 = d["p0"][eb]
    origin = d["origin"][eb]; magzp = d["magzp"][eb]
    centre = tuple(d["centre"])
    nz, nm = len(C.Z_GRID), len(C.MU_GRID)
    q_grid = 1.0 / C.Z_GRID
    E = len(eb)
    prf = _prf(band)
    # Real S per mask (full grid) and the capped weights, which do not
    # depend on the injection (v, g unchanged): precomputed once for the
    # full, early and late epoch selections so that each injection is a
    # weighted sum over its window.
    real = {}
    f0n = np.nan_to_num(f0.astype(np.float32))
    for mi, m in enumerate(MASKS):
        ok = masks[:, mi]
        if ok.sum() < C.MIN_EPOCHS:
            real[m] = None
            continue
        S, A, B, n, w = nulls.stack_S(f0, v0, g0, ok, return_parts=True)
        early = ok & (mjd <= C.HOLDOUT_SPLIT_MJD)
        late = ok & (mjd > C.HOLDOUT_SPLIT_MJD)
        w_e = nulls.stack_weights(v0, g0, early, f=f0)[0] if early.sum() >= C.MIN_EPOCHS else None
        w_l = nulls.stack_weights(v0, g0, late, f=f0)[0] if late.sum() >= C.MIN_EPOCHS else None
        S_e_real = None
        if w_e is not None:
            with np.errstate(all="ignore"):
                Be = w_e.sum(0); S_e_real = (f0n * w_e).sum(0) / np.sqrt(np.where(Be > 0, Be, np.inf))
            S_e_real = np.where((w_e > 0).sum(0) >= C.MIN_EPOCHS, S_e_real, np.nan)
        real[m] = {"S": S, "ok": ok, "w": w, "early": early, "late": late,
                   "w_e": w_e, "w_l": w_l, "S_e_real": S_e_real}
    T = {m: float(np.nanmax(d["summary"][mi, 1 + np.asarray(d["designated"]), b - 1, 0]))
         for mi, m in enumerate(MASKS)}
    m90 = _v1_m90(endpoint, role, band)
    sig_xt = _envelope_xt(endpoint, role)
    models = list(C.TEMPORAL_MODELS)
    out = []
    dt_yr = (mjd - C.T0_MJD) / 365.25
    cosd = np.cos(np.deg2rad(centre[1]))
    amp_ref = None
    for j in range(n_inj):
        model = models[j % len(models)]
        z = float(np.exp(rng.uniform(np.log(550.0), np.log(10000.0))))
        mu = rng.uniform(-1.0, 1.0, size=2)
        xt = float(rng.uniform(-sig_xt, sig_xt)) if sig_xt > 0 else 0.0
        mag = float(rng.uniform(m90 - C.INJ_MAG_HALFWIDTH, m90 + C.INJ_MAG_HALFWIDTH))
        on = temporal_on(model, mjd, rng)
        amp_ref = 10.0 ** (0.4 * (C.ZP_REF - mag))          # tensor units (ZP 20)
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
        pix_inj = p0 + np.einsum("eij,ej->ei", J, track)                            # (E, 2)
        # window node pixel positions (E, nw, nm, nm, 2)
        off = np.empty((E, len(win), nm, nm, 2))
        off[..., 0] = d0[:, win, 0][:, :, None, None] + (C.MU_GRID * dt_yr[:, None])[:, None, :, None]
        off[..., 1] = d0[:, win, 1][:, :, None, None] + (C.MU_GRID * dt_yr[:, None])[:, None, None, :]
        pix_nodes = p0[:, None, None, None, :] + np.einsum("eij,e...j->e...i", J, off)
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
            el = prf.element_of(x + origin[e, 0], y + origin[e, 1])
            stamp, ox, oy = prf.render(x, y, element=el, half=C.STAMP_HALF)
            rw = inject.stamp_response(fm, stamp, ox, oy)
            resp0[e] = rw.sample(x, y)
            delta[e] = (amp_ref * rw.sample(pix_nodes[e, ..., 0], pix_nodes[e, ..., 1])).astype(np.float32)
        f_inj = f0n[:, win] + delta
        rec = {"j": j, "model": model, "z_au": z, "mu": mu.tolist(), "xt_arcsec": xt, "mag": mag,
               "fnu_jy": inject.fnu_from_vega_mag(band, mag, C.SPECTRUM[band]),
               "iz": iz, "n_on": int(on.sum()),
               "resp_median": float(np.nanmedian(resp0)) if np.isfinite(resp0).any() else None,
               "masks": {}}
        # distance of every window node from the injection, in nodes
        wz = np.abs(win - iz)[:, None, None]
        wmu = np.maximum(np.abs(np.arange(nm)[:, None] - np.argmin(np.abs(C.MU_GRID - mu[0]))),
                         np.abs(np.arange(nm)[None, :] - np.argmin(np.abs(C.MU_GRID - mu[1]))))[None]
        near = (wz <= C.RECOVERY_WINDOW[0]) & (wmu <= C.RECOVERY_WINDOW[1])
        for mi, m in enumerate(MASKS):
            if real[m] is None:
                continue
            R_ = real[m]
            S_real, ok = R_["S"], R_["ok"]
            w_w = R_["w"][:, win]
            A_w = (f_inj * w_w).sum(0); B_w = w_w.sum(0); n_w = (w_w > 0).sum(0)
            with np.errstate(all="ignore"):
                S_w = A_w / np.sqrt(np.where(B_w > 0, B_w, np.inf))
            S_w = np.where(n_w >= C.MIN_EPOCHS, S_w, np.nan)
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
                S_e = np.where(n_e >= C.MIN_EPOCHS, S_e, np.nan)
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
            if catalog is not None and len(catalog.ra):
                tr = np.empty((E, 2))
                tr[:, 0] = d0[:, win[pz], 0] + C.MU_GRID[pm1] * dt_yr
                tr[:, 1] = d0[:, win[pz], 1] + C.MU_GRID[pm2] * dt_yr
                ra_t = centre[0] + tr[:, 0] / 3600.0 / cosd
                dec_t = centre[1] + tr[:, 1] / 3600.0
                fc = flux_consistency(catalog, ra_t, dec_t, we * ok, phase, C.ZP_REF,
                                      resp_table[0], resp_table[1], s_peak,
                                      S_by_phase=S_ph, search_arcsec=C.STATIC_SEARCH_ARCSEC,
                                      factor=C.STATIC_FACTOR, min_ndet=C.CATWISE_MIN_NDET)
                r["static"] = {"S_pred": fc["S_pred"], "consistent": fc["consistent"],
                               "n_sources": fc["n_sources_considered"],
                               "nearest": (fc["sources"][0] if fc["sources"] else None)}
            rec["masks"][m] = r
        out.append(rec)
    return {"endpoint": endpoint, "role": role, "band": band, "n_epochs": int(E),
            "m90_v1": m90, "sigma_xt_99": sig_xt, "T": T, "amp_note": "tensor units ZP 20",
            "injections": out}


def run_corridor(args):
    corridor, n_inj, only_missing = args
    t0 = _time.monotonic()
    obs_by_id, manifest = _G["obs_by_id"], _G["manifest"]
    pairs = []
    for e in MEMBERS[corridor]:
        for r in ("rx", "tx"):
            for tp in [C.TENSOR_DIR / f"{e}__{r}.npz"] + sorted(C.TENSOR_DIR.glob(f"{e}__{r}__xt*.npz")):
                if tp.exists():
                    pairs.append((e, r, tp))
    C.INJ_DIR.mkdir(parents=True, exist_ok=True)
    done = []
    for endpoint, role, tpath in pairs:
        out_path = C.INJ_DIR / (tpath.stem + ".json")
        if only_missing and out_path.exists():
            continue
        d = np.load(tpath)
        centre = tuple(d["centre"])
        catalog = _catalog(corridor, centre)
        result = {"schema": SCHEMA, "endpoint": endpoint, "role": role, "corridor": corridor,
                  "xt_arcsec": float(d["xt_arcsec"]) if "xt_arcsec" in d.files else 0.0,
                  "tensor_input_hash": str(d["input_hash"]), "freeze_hash": C.freeze_hash(),
                  "n_per_cell": n_inj, "catalog": None if catalog is None else catalog.label,
                  "cells": {}}
        seed = int(hashlib.sha256(f"inj/{endpoint}/{role}/{C.SPLIT_SEED}".encode()).hexdigest()[:8], 16)
        rng = np.random.default_rng(seed)
        result["seed"] = seed
        for band in ("W1", "W2", "W3", "W4"):
            b = C.BAND_IDX[band]
            eb = np.where(d["band_idx"] == b)[0]
            if len(eb) < C.MIN_EPOCHS:
                continue
            # rebuild the band's flux maps (same code path as the tensor build)
            maps = []
            for e in eb:
                oid = str(d["oid"][e])
                row = manifest.get(oid); obs = obs_by_id[oid]
                fm = None
                if row and "int" in row.get("files", {}):
                    try:
                        fm = build_flux_map(C.CUT_DIR / row["files"]["int"], C.CUT_DIR / row["files"]["unc"],
                                            C.MSK_DIR / obs["products"]["msk"]["url"].rsplit("/", 1)[-1],
                                            band, obs["t_mid_mjd_utc"], WiseExactFootprint.FATAL_MASK,
                                            magzp=obs["quality_flags"].get("magzp"),
                                            keep_inputs=True, pix_scale_from_header=True)
                    except Exception:
                        fm = None
                maps.append(fm)
            pix_arcsec = float(np.median(d["pix_arcsec"][eb]))
            cell = inject_cell(d, band, maps, n_inj, rng, catalog, _response_table(band, pix_arcsec))
            del maps
            if cell:
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


def _init(obs_by_id, manifest):
    _G["obs_by_id"] = obs_by_id
    _G["manifest"] = manifest


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", choices=["dev", "confirmatory", "all"], required=True)
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--only-missing", action="store_true")
    ap.add_argument("--n-per-cell", type=int, default=C.N_INJ_PER_CELL)
    ap.add_argument("--corridors", nargs="*", default=None)
    a = ap.parse_args()
    freeze = C.load_freeze()
    if a.corridors:
        corridors = a.corridors
    elif a.set == "all":
        corridors = sorted(MEMBERS)
    else:
        corridors = freeze["split"]["development" if a.set == "dev" else "confirmatory"]["corridors"]
    obs_by_id = {r["observation_id"]: r for r in read_records(C.COARSE_DIR / "records" / "observation.jsonl")}
    manifest = {}
    for line in open(C.CUT_DIR / "manifest.jsonl"):
        row = json.loads(line); manifest[row["observation_id"]] = row
    # warm the catalogue cache serially (network) before forking
    for c in corridors:
        for e in MEMBERS[c]:
            p = C.TENSOR_DIR / f"{e}__rx.npz"
            if p.exists():
                d = np.load(p); _catalog(c, tuple(d["centre"])); d.close()
                break
    jobs = [(c, a.n_per_cell, a.only_missing) for c in corridors]
    log = C.RUN_DIR / "inject_pipeline.log"
    if a.workers <= 1:
        _init(obs_by_id, manifest)
        for j in jobs:
            msg = run_corridor(j); print(msg, flush=True)
            open(log, "a").write(msg + "\n")
        return
    import multiprocessing as mp
    with mp.get_context("fork").Pool(a.workers, initializer=_init, initargs=(obs_by_id, manifest)) as pool:
        for msg in pool.imap_unordered(run_corridor, jobs):
            print(msg, flush=True)
            open(log, "a").write(msg + "\n")


if __name__ == "__main__":
    main()
