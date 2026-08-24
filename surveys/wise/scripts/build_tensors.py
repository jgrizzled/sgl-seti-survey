"""Step D (part 1): v2 sample tensors and per-trajectory null summaries.

Per corridor, band by band, every usable v1 cutout (loose mask; the
primary / strict masks are epoch subsets applied at stack time) is
turned into a matched-filter flux map once (pixel scale from the
header, spacecraft observer) and sampled on:

  * the real (z, mu) grid (trajectory 0),
  * the 48-offset spatial ring (trajectories 1..48; the 8 designated
    controls are among them, see v2common.DESIGNATED),
  * 50 trajectory-randomised tracks (49..98): another endpoint's
    (z, mu) sky motion re-centred on this corridor at T0.

Per pair writes runs/wise/v2/tensors/<endpoint>__<role>.npz:
  f, v (9, E, 64, 5, 5) float32 / g float16 — real + designated controls
     at common ZP 20 (as v1), for injections, vetting and scrambles;
  summaries for all 99 trajectories x 3 masks x bands: grid maximum,
     argmax node, phase-split S at the argmax, epoch count, top-epoch
     share, held-out-epoch refit/prediction (hypotheses v2.0 §3.7b);
  phase-coherence scramble maxima (200) and per-epoch scramble maxima
     (50, diagnostic) of the real trajectory per mask x band;
  geometry per epoch: node pixel positions at mu = 0, the local
     WCS Jacobian, tangent-plane node offsets, frame origin, MAGZP,
     quality flags, observation ids — what the injection stage needs
     to place a source on any continuous (z, mu, cross-track) track.
  input_hash: hash of the cutout digests used (stale detection).

Usage: uv run python surveys/wise/scripts/build_tensors.py
           [--corridors c1 c2 ...] [--workers 6] [--only-missing]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time as _time
from collections import defaultdict
from pathlib import Path

import warnings

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v2common as C  # noqa: E402

warnings.filterwarnings("ignore", message="All-NaN slice")

from astropy.time import Time  # noqa: E402
from sglseti import Role, adaptive_locus, load_target_registry  # noqa: E402
from sglsurvey import nulls  # noqa: E402
from sglsurvey.adapters.irsa_wise import WiseExactFootprint  # noqa: E402
from sglsurvey.corridors import CORRIDOR_OF, MEMBERS  # noqa: E402
from sglsurvey.geometry import (register_wise_spacecraft_observer,  # noqa: E402
                                wise_v2_context)
from sglsurvey.manifest import HashCache, combined_hash  # noqa: E402
from sglsurvey.photometry import build_flux_map  # noqa: E402
from sglsurvey.records import read_records  # noqa: E402

DONOR_BIN_DAYS = 1.0
SCHEMA = "wise-v2-tensor-v1"
MASKS = ("primary", "strict", "loose")
N_TRAJ = 1 + len(C.RING_OFFSETS) + C.N_DONORS

_G: dict = {}   # per-process globals (registry, context, caches)


# ---------------------------------------------------------------------------
def _init_worker():
    registry = load_target_registry(C.REGISTRY_PATH)
    tab = np.load(C.RUN_DIR / "observer" / "wise_sc_ephemeris.npz")
    ident = json.loads((C.RUN_DIR / "observer" / "summary.json").read_text())["table_sha256"]
    obs = register_wise_spacecraft_observer(tab["mjd_utc"], tab["xyz_au"], ident)
    _G["registry"] = registry
    _G["ctx"] = wise_v2_context(obs)
    _G["zcache"] = {}
    _G["observer_identity"] = ident


def _locus_at_zgrid(endpoint: str, role: str, mjd: float, bin_days: float):
    """(nz, 2) ICRS deg of the endpoint's locus at the Z_GRID nodes, at
    the epoch bin containing ``mjd`` (cached)."""
    key = (endpoint, role, round(mjd / bin_days))
    zc = _G["zcache"]
    if key not in zc:
        ctx = _G["ctx"]
        al = adaptive_locus(
            target=_G["registry"][endpoint], role=Role(role),
            observation_time=Time(key[2] * bin_days, format="mjd"),
            observer=ctx.observer, relay_range=ctx.relay_range,
            tolerance_arcsec=2.0, ephemeris=ctx.ephemeris, model=ctx.model)
        zs = np.array([p.z_au for p in al.points])
        ra = np.array([p.icrs_ra_deg for p in al.points])
        dec = np.array([p.icrs_dec_deg for p in al.points])
        q = 1.0 / zs
        o = np.argsort(q)
        qg = 1.0 / C.Z_GRID
        # unwrap RA across 0/360 before interpolating
        ra_u = np.unwrap(np.deg2rad(ra[o]))
        zc[key] = np.stack([np.rad2deg(np.interp(qg, q[o], ra_u)) % 360.0,
                            np.interp(qg, q[o], dec[o])], axis=1)
    return zc[key]


def _tangent_offsets(pts: np.ndarray, centre: tuple[float, float]) -> np.ndarray:
    """(N, 2) arcsec (dRA cos dec, dDec) of ICRS points from ``centre``."""
    ra0, dec0 = centre
    dra = (pts[:, 0] - ra0 + 180.0) % 360.0 - 180.0
    return np.stack([dra * np.cos(np.deg2rad(dec0)) * 3600.0,
                     (pts[:, 1] - dec0) * 3600.0], axis=1)


def _donors_for(corridor: str, role: str) -> list[str]:
    """Seeded draw of N_DONORS endpoints from other corridors."""
    pool = sorted(e for e, c in CORRIDOR_OF.items() if c != corridor)
    seed = int(hashlib.sha256(f"{corridor}/{role}/{C.SPLIT_SEED}".encode()).hexdigest()[:8], 16)
    rng = np.random.default_rng(seed)
    return sorted(rng.choice(pool, size=C.N_DONORS, replace=False).tolist())


def _holdout(f, v, g, mjd, epoch_ok):
    """Refit on early epochs, forced photometry on late epochs at the
    early argmax node (hypotheses v2.0 §3.7b). Returns a small array."""
    early = epoch_ok & (mjd <= C.HOLDOUT_SPLIT_MJD)
    late = epoch_ok & (mjd > C.HOLDOUT_SPLIT_MJD)
    out = np.full(7, np.nan)   # S_early, f_early, S_late, f_late, n_late, iz, imu
    if early.sum() < C.MIN_EPOCHS:
        return out
    S_e, A_e, B_e, n_e, _ = nulls.stack_S(f, v, g, early, return_parts=True)
    s_max, node = nulls.grid_max(S_e)
    if node is None:
        return out
    out[0] = s_max
    out[1] = A_e[node] / B_e[node] if B_e[node] > 0 else np.nan
    out[5] = node[0]
    out[6] = node[1] * len(C.MU_GRID) + node[2]
    n_l = int(late.sum())
    out[4] = n_l
    if n_l >= C.MIN_EPOCHS:
        S_l, A_l, B_l, n_ln, _ = nulls.stack_S(f, v, g, late, return_parts=True,
                                               min_epochs=1)
        out[2] = S_l[node]
        out[3] = A_l[node] / B_l[node] if B_l[node] > 0 else np.nan
        out[4] = n_ln[node]
    return out


def _summaries(F, V, G, mjd, phase, epoch_ok):
    """Per-trajectory summaries (T, 12) for one mask x band:
    S_max, iz, imu, S_phase0, S_phase1, n_epochs_at_max, top_share,
    then the 7 hold-out numbers... (packed as two arrays)."""
    T = F.shape[0]
    summ = np.full((T, 7), np.nan)
    hold = np.full((T, 7), np.nan)
    for t in range(T):
        S, A, B, n, w = nulls.stack_S(F[t], V[t], G[t], epoch_ok, return_parts=True)
        s_max, node = nulls.grid_max(S)
        if node is None:
            continue
        summ[t, 0] = s_max
        summ[t, 1] = node[0]
        summ[t, 2] = node[1] * len(C.MU_GRID) + node[2]
        summ[t, 5] = n[node]
        fe = np.nan_to_num(F[t][(slice(None),) + node].astype(np.float64))
        we = w[(slice(None),) + node]
        contrib = fe * we
        summ[t, 6] = float(np.max(contrib) / A[node]) if A[node] > 0 else np.nan
        for p in (0, 1):
            sel = epoch_ok & (phase == p)
            if sel.sum() >= 1:
                Bp = we[sel].sum()
                summ[t, 3 + p] = contrib[sel].sum() / np.sqrt(Bp) if Bp > 0 else np.nan
        hold[t] = _holdout(F[t], V[t], G[t], mjd, epoch_ok)
    return summ, hold


# ---------------------------------------------------------------------------
def build_corridor(corridor: str, only_missing: bool) -> str:
    t_start = _time.monotonic()
    registry = _G["registry"]
    obs_by_id = _G["obs_by_id"]
    usable = _G["usable"]
    manifest = _G["manifest"]
    pairs = sorted(p for p in usable if CORRIDOR_OF[p[0]] == corridor)
    # cross-track variants (hypotheses v2.0 §1.4): extra tensors per offset
    xt_cfg = {}
    xt_path = C.CONFIG_DIR / "cross_track_cells.json"
    if xt_path.exists():
        for c in json.loads(xt_path.read_text())["cells"]:
            xt_cfg[(c["endpoint"], c["role"])] = [x for x in c["offsets_arcsec"] if x != 0.0]
    pairs = [(e, r, x) for (e, r) in pairs for x in [0.0] + xt_cfg.get((e, r), [])]
    if only_missing:
        # stale-product detection (v2 plan §8.2): an existing tensor is
        # kept only if its recorded input hash matches the cutouts now
        # on disk; otherwise it is rebuilt (and the event logged).
        fresh = []
        cache = HashCache(C.RUN_DIR / "manifest" / f"cutout_hashes_{corridor}.json")
        for p in pairs:
            out = C.TENSOR_DIR / tensor_name(*p)
            if not out.exists():
                fresh.append(p)
                continue
            try:
                with np.load(out) as d:
                    recorded = str(d["input_hash"])
                    oids = [str(o) for o in d["oid"]]
                digs = []
                for oid in oids:
                    row = manifest[oid]; obs = obs_by_id[oid]
                    digs += [cache.sha256(C.CUT_DIR / row["files"]["int"]),
                             cache.sha256(C.CUT_DIR / row["files"]["unc"]),
                             cache.sha256(C.MSK_DIR / obs["products"]["msk"]["url"].rsplit("/", 1)[-1])]
                cache.flush()
                if combined_hash(digs) != recorded:
                    print(f"  [{corridor}] STALE tensor {out.name}: inputs changed -> rebuilding", flush=True)
                    fresh.append(p)
            except Exception as exc:
                print(f"  [{corridor}] unreadable tensor {out.name} ({exc}) -> rebuilding", flush=True)
                fresh.append(p)
        pairs = fresh
    if not pairs:
        return f"{corridor}: nothing to do"
    frames = sorted({o for p in pairs for o in usable[p[:2]]},
                    key=lambda o: obs_by_id[o]["t_mid_mjd_utc"])
    centre = tuple(manifest[next(o for o in frames if o in manifest)]["center_radec"])
    # phase reference as v1: median day-of-year of the first 50 frames
    doys = np.array([obs_by_id[o]["t_mid_mjd_utc"] % 365.25 for o in frames[:50]])
    phase_ref = float(np.median(doys))

    def phase_of(mjd):
        d = (mjd % 365.25) - phase_ref
        d = (d + 182.625) % 365.25 - 182.625
        return 0 if abs(d) < 91.3 else 1

    donors = {role: _donors_for(corridor, role) for role in ("rx", "tx")}
    cache = HashCache(C.RUN_DIR / "manifest" / f"cutout_hashes_{corridor}.json")
    nz, nm = len(C.Z_GRID), len(C.MU_GRID)
    ring = np.array(C.RING_OFFSETS)            # (48, 2) arcsec
    keep = [0] + [1 + i for i in C.DESIGNATED]
    usable_set = _G["usable_set"]
    store = {p: [] for p in pairs}             # pair -> list of per-band dicts
    digests = {p: [] for p in pairs}
    rngs = {p: np.random.default_rng(int(hashlib.sha256(
        f"{p[0]}/{p[1]}/{C.SPLIT_SEED}".encode()).hexdigest()[:8], 16)) for p in pairs}
    n_maps = 0
    cosd = np.cos(np.deg2rad(centre[1]))
    eps_arcsec = 10.0
    probe = np.array([[centre[0], centre[1]],
                      [centre[0] + eps_arcsec / 3600.0 / cosd, centre[1]],
                      [centre[0], centre[1] + eps_arcsec / 3600.0]])
    for band in ("W1", "W2", "W3", "W4"):
        b = C.BAND_IDX[band]
        # ---- flux maps of this band, cached in memory ----------------------
        maps = []
        for oid in (o for o in frames if obs_by_id[o]["band"] == band):
            row = manifest.get(oid)
            obs = obs_by_id[oid]
            if row is None or "int" not in row.get("files", {}):
                continue
            int_p = C.CUT_DIR / row["files"]["int"]
            unc_p = C.CUT_DIR / row["files"]["unc"]
            msk_p = C.MSK_DIR / obs["products"]["msk"]["url"].rsplit("/", 1)[-1]
            if not (int_p.exists() and unc_p.exists() and msk_p.exists()):
                continue
            try:
                fm = build_flux_map(int_p, unc_p, msk_p, band, obs["t_mid_mjd_utc"],
                                    WiseExactFootprint.FATAL_MASK,
                                    magzp=obs["quality_flags"].get("magzp"),
                                    keep_inputs=True, pix_scale_from_header=True)
            except Exception as exc:
                print(f"  [{corridor}] fluxmap FAIL {oid}: {exc}", flush=True)
                continue
            n_maps += 1
            pp = fm.world2pix(probe[:, 0], probe[:, 1])
            J = np.stack([(pp[1] - pp[0]) / eps_arcsec, (pp[2] - pp[0]) / eps_arcsec], axis=1)
            zp = fm.magzp if fm.magzp else C.ZP_REF
            maps.append({
                "oid": oid, "mjd": obs["t_mid_mjd_utc"], "flux": fm.flux, "var": fm.var,
                "gf": fm.good_frac, "p0": pp[0], "J": J, "magzp": float(zp),
                "scale": 10.0 ** ((C.ZP_REF - zp) / 2.5),
                "origin": np.array(fm.frame_origin, dtype=np.int32),
                "pix_arcsec": float(getattr(fm, "pix_arcsec", 2.75)),
                "quality": [obs["quality_flags"].get(k) for k in
                            ("qual_frame", "qual_scan", "saa_sep", "moon_sep")],
                "mask": [C.quality_ok(obs["quality_flags"], m) for m in MASKS],
                "dig": (cache.sha256(int_p), cache.sha256(unc_p), cache.sha256(msk_p)),
            })
            del fm
        if not maps:
            continue
        # ---- per pair: sample all trajectories, summarise, keep 9 -----------
        for pair in pairs:
            endpoint, role, xt = pair
            sel = [m for m in maps if m["oid"] in usable_set[pair[:2]]]
            if not sel:
                continue
            base_t0 = _locus_at_zgrid(endpoint, role, C.T0_MJD, C.LOCUS_BIN_DAYS)
            d0_t0 = _tangent_offsets(base_t0, centre)
            E = len(sel)
            F_all = np.empty((N_TRAJ, E, nz, nm, nm), dtype=np.float32)
            V_all = np.empty_like(F_all)
            G_all = np.empty((N_TRAJ, E, nz, nm, nm), dtype=np.float16)
            xy0 = np.empty((E, nz, 2), np.float32)
            d0s = np.empty((E, nz, 2), np.float32)
            for ei, m in enumerate(sel):
                mjd = m["mjd"]
                dt_yr = (mjd - C.T0_MJD) / 365.25
                base = _locus_at_zgrid(endpoint, role, mjd, C.LOCUS_BIN_DAYS)
                d0 = _tangent_offsets(base, centre)
                dmu = C.MU_GRID * dt_yr
                if xt != 0.0:
                    # cross-track offset along the local normal of the locus
                    tv = np.gradient(d0, axis=0)
                    nn = np.linalg.norm(tv, axis=1, keepdims=True)
                    tv = np.where(nn > 0, tv / np.where(nn > 0, nn, 1), [[1.0, 0.0]])
                    d0 = d0 + xt * np.stack([-tv[:, 1], tv[:, 0]], axis=1)
                all_off = np.empty((N_TRAJ, nz, nm, nm, 2))
                all_off[0, ..., 0] = d0[:, 0][:, None, None] + dmu[None, :, None]
                all_off[0, ..., 1] = d0[:, 1][:, None, None] + dmu[None, None, :]
                all_off[1:1 + len(ring)] = all_off[0][None] + ring[:, None, None, None, :]
                for k, dn in enumerate(donors[role]):
                    dl_t = _locus_at_zgrid(dn, role, mjd, DONOR_BIN_DAYS)
                    dl_0 = _locus_at_zgrid(dn, role, C.T0_MJD, DONOR_BIN_DAYS)
                    dc = (dl_0[nz // 2, 0], dl_0[nz // 2, 1])
                    delta = _tangent_offsets(dl_t, dc) - _tangent_offsets(dl_0, dc)
                    kk = 1 + len(ring) + k
                    all_off[kk, ..., 0] = (d0_t0[:, 0] + delta[:, 0])[:, None, None] + dmu[None, :, None]
                    all_off[kk, ..., 1] = (d0_t0[:, 1] + delta[:, 1])[:, None, None] + dmu[None, None, :]
                pix = m["p0"][None, None, None, None, :] + np.einsum("ij,...j->...i", m["J"], all_off)
                f, vv, g = _sample_pix(m, pix[..., 0].ravel(), pix[..., 1].ravel())
                F_all[:, ei] = (f.reshape(N_TRAJ, nz, nm, nm) * m["scale"]).astype(np.float32)
                V_all[:, ei] = np.minimum(vv.reshape(N_TRAJ, nz, nm, nm) * m["scale"] ** 2,
                                          1e30).astype(np.float32)
                G_all[:, ei] = g.reshape(N_TRAJ, nz, nm, nm).astype(np.float16)
                xy0[ei] = pix[0, :, nm // 2, nm // 2, :]
                d0s[ei] = d0
                digests[pair].extend(m["dig"])
            mjd_b = np.array([m["mjd"] for m in sel])
            phase_b = np.array([phase_of(m["mjd"]) for m in sel], dtype=np.uint8)
            masks_b = np.array([m["mask"] for m in sel], dtype=bool)
            summ = np.full((len(MASKS), N_TRAJ, 7), np.nan, dtype=np.float32)
            hold = np.full((len(MASKS), N_TRAJ, 7), np.nan, dtype=np.float32)
            scr = np.full((len(MASKS), C.N_SCRAMBLE), np.nan, dtype=np.float32)
            escr = np.full((len(MASKS), C.N_EPOCH_SCRAMBLE), np.nan, dtype=np.float32)
            n_ok = np.zeros(len(MASKS), dtype=np.int32)
            for mi in range(len(MASKS)):
                ok = masks_b[:, mi]
                n_ok[mi] = int(ok.sum())
                if ok.sum() < C.MIN_EPOCHS:
                    continue
                summ[mi], hold[mi] = _summaries(F_all, V_all, G_all, mjd_b, phase_b, ok)
                parts = nulls.phase_parts(F_all[0], V_all[0], G_all[0], phase_b, ok)
                scr[mi] = nulls.phase_scrambled_maxima(parts, C.N_SCRAMBLE, rngs[pair])
                escr[mi] = nulls.epoch_scrambled_maxima(
                    F_all[0], V_all[0], G_all[0], C.N_EPOCH_SCRAMBLE, rngs[pair], ok)
            store[pair].append({
                "band": b, "mjd": mjd_b, "phase": phase_b, "masks": masks_b,
                "f": F_all[keep].copy(), "v": V_all[keep].copy(), "g": G_all[keep].copy(),
                "xy0": xy0, "d0": d0s,
                "J": np.stack([m["J"] for m in sel]).astype(np.float32),
                "p0": np.stack([m["p0"] for m in sel]).astype(np.float32),
                "origin": np.stack([m["origin"] for m in sel]),
                "magzp": np.array([m["magzp"] for m in sel], np.float32),
                "pix_arcsec": np.array([m["pix_arcsec"] for m in sel], np.float32),
                "quality": np.array([[np.nan if q is None else q for q in m["quality"]]
                                     for m in sel], np.float32),
                "oid": np.array([m["oid"] for m in sel]),
                "summary": summ, "holdout": hold, "scramble": scr, "escramble": escr,
                "n_ok": n_ok,
            })
            del F_all, V_all, G_all
        del maps
    cache.flush()

    # ---- assemble and write per pair ----------------------------------------
    written = []
    for pair in pairs:
        endpoint, role, xt = pair
        parts = store[pair]
        if not parts:
            continue
        n_bands = 4
        E = sum(len(q["mjd"]) for q in parts)
        mjd = np.concatenate([q["mjd"] for q in parts])
        order = np.argsort(mjd)
        cat = lambda key: np.concatenate([q[key] for q in parts], axis=0)[order]  # noqa: E731
        catf = lambda key: np.concatenate([q[key] for q in parts], axis=1)[:, order]  # noqa: E731
        band_idx = np.concatenate([np.full(len(q["mjd"]), q["band"], np.uint8) for q in parts])[order]
        summ = np.full((len(MASKS), N_TRAJ, n_bands, 7), np.nan, np.float32)
        hold = np.full((len(MASKS), N_TRAJ, n_bands, 7), np.nan, np.float32)
        scr = np.full((len(MASKS), n_bands, C.N_SCRAMBLE), np.nan, np.float32)
        escr = np.full((len(MASKS), n_bands, C.N_EPOCH_SCRAMBLE), np.nan, np.float32)
        n_ok = np.zeros((len(MASKS), n_bands), np.int32)
        for q in parts:
            bi = q["band"] - 1
            summ[:, :, bi] = q["summary"]
            hold[:, :, bi] = q["holdout"]
            scr[:, bi] = q["scramble"]
            escr[:, bi] = q["escramble"]
            n_ok[:, bi] = q["n_ok"]
        input_hash = combined_hash(digests[pair])
        out = C.TENSOR_DIR / tensor_name(endpoint, role, xt)
        C.TENSOR_DIR.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            out,
            schema=SCHEMA, input_hash=input_hash, endpoint=endpoint, role=role, xt_arcsec=xt,
            corridor=corridor, centre=np.array(centre), phase_ref=phase_ref,
            observer_identity=_G["observer_identity"],
            z_grid=C.Z_GRID, mu_grid=C.MU_GRID, t0_mjd=C.T0_MJD, zp_ref=C.ZP_REF,
            ring_offsets=ring, designated=np.array(C.DESIGNATED),
            kept_trajectories=np.array(keep), donors=np.array(donors[role]),
            mjd=mjd[order], band_idx=band_idx, phase=cat("phase"), masks=cat("masks"),
            mask_names=np.array(MASKS), quality=cat("quality"), magzp=cat("magzp"),
            pix_arcsec=cat("pix_arcsec"), oid=cat("oid"),
            xy0=cat("xy0"), d0=cat("d0"), J=cat("J"), p0=cat("p0"), origin=cat("origin"),
            f=catf("f"), v=catf("v"), g=catf("g"),
            summary=summ, holdout=hold, scramble_max=scr, epoch_scramble_max=escr,
            n_epochs_ok=n_ok,
            summary_fields=np.array(["S_max", "iz", "imu", "S_phase0", "S_phase1",
                                     "n_epochs", "top_share"]),
            holdout_fields=np.array(["S_early", "f_early", "S_late", "f_late",
                                     "n_late", "iz_early", "imu_early"]),
        )
        written.append(f"{out.stem} ({E} epochs)")
    dt = _time.monotonic() - t_start
    return f"{corridor}: {n_maps} maps, wrote {', '.join(written)} in {dt / 60:.1f} min"


def tensor_name(endpoint, role, xt=0.0):
    return f"{endpoint}__{role}.npz" if xt == 0.0 else f"{endpoint}__{role}__xt{xt:+.2f}.npz"


def _sample_pix(m, x, y):
    """Bilinear (flux, var, good_frac) at pixel coordinates (NaN outside)."""
    ny, nx = m["flux"].shape
    out_f = np.full(len(x), np.nan)
    out_v = np.full(len(x), np.nan)
    out_g = np.full(len(x), np.nan)
    ok = (x >= 0) & (x <= nx - 1.001) & (y >= 0) & (y <= ny - 1.001) & np.isfinite(x) & np.isfinite(y)
    if ok.any():
        x0 = np.floor(x[ok]).astype(int)
        y0 = np.floor(y[ok]).astype(int)
        fx, fy = x[ok] - x0, y[ok] - y0

        def bil(a):
            return (a[y0, x0] * (1 - fx) * (1 - fy) + a[y0, x0 + 1] * fx * (1 - fy)
                    + a[y0 + 1, x0] * (1 - fx) * fy + a[y0 + 1, x0 + 1] * fx * fy)

        out_f[ok] = bil(m["flux"])
        out_v[ok] = bil(m["var"])
        out_g[ok] = bil(m["gf"])
    return out_f, out_v, out_g


# ---------------------------------------------------------------------------
def load_inputs():
    obs_by_id = {r["observation_id"]: r for r in read_records(
        C.COARSE_DIR / "records" / "observation.jsonl")}
    usable = defaultdict(list)
    for r in read_records(C.PRECISE_DIR / "records" / "intersection_evaluation.jsonl"):
        if r["stage"] == "precise" and r["usable"] in ("usable", "partial"):
            usable[(r["endpoint_id"], r["role"])].append(r["observation_id"])
    manifest = {}
    for line in open(C.CUT_DIR / "manifest.jsonl"):
        row = json.loads(line)
        manifest[row["observation_id"]] = row
    return obs_by_id, dict(usable), manifest


def _worker_init(obs_by_id, usable, manifest):
    _init_worker()
    _G["obs_by_id"] = obs_by_id
    _G["usable"] = usable
    _G["usable_set"] = {p: set(v) for p, v in usable.items()}
    _G["manifest"] = manifest


def _run(args):
    corridor, only_missing = args
    try:
        return build_corridor(corridor, only_missing)
    except Exception as exc:
        import traceback
        return f"{corridor}: FAILED {exc}\n{traceback.format_exc()}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--corridors", nargs="*", default=None)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--only-missing", action="store_true")
    a = ap.parse_args()
    corridors = a.corridors or sorted(MEMBERS)
    obs_by_id, usable, manifest = load_inputs()
    print(f"{len(corridors)} corridors, {len(usable)} usable pairs, "
          f"{len(obs_by_id)} observations", flush=True)
    (C.RUN_DIR / "manifest").mkdir(parents=True, exist_ok=True)
    log = C.RUN_DIR / "build_tensors.log"
    if a.workers <= 1:
        _worker_init(obs_by_id, usable, manifest)
        for c in corridors:
            msg = _run((c, a.only_missing))
            print(msg, flush=True)
            with open(log, "a") as fh:
                fh.write(msg + "\n")
        return
    import multiprocessing as mp
    ctx = mp.get_context("fork")
    with ctx.Pool(a.workers, initializer=_worker_init,
                  initargs=(obs_by_id, usable, manifest)) as pool:
        for msg in pool.imap_unordered(_run, [(c, a.only_missing) for c in corridors]):
            print(msg, flush=True)
            with open(log, "a") as fh:
                fh.write(msg + "\n")


if __name__ == "__main__":
    main()
