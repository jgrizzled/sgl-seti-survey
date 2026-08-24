"""v2 engine: tensor build with streaming null accumulation.

Per corridor and band, every usable v1 cutout becomes a FluxMap
(profile.build_map) twice: pass 0 measures the frame weight (the
inverse median matched-filter variance over usable pixels) so that the
band's per-frame cap WEIGHT_CAP x median(frame weight) is known; pass 1
samples the real (z, mu) grid, the 48-offset ring and the
trajectory-randomised tracks, applies the cap and accumulates the
stack sums per quality mask — full, per parallax phase (real only),
and early / late epochs for the held-out test — for all trajectories,
so that only the real + 8 designated per-epoch tensors are stored.

Per pair writes <tensor_dir>/<endpoint>__<role>[__xt±X].npz with the
same fields as the WISE v2 tensors plus ``frame_cap`` (E,) and
``cap_kind = "per-frame"``.
"""

from __future__ import annotations

import hashlib
import json
import time as _time
import warnings
from pathlib import Path

import numpy as np
from astropy.time import Time

from sglseti import Role, adaptive_locus, load_target_registry

from sglsurvey import nulls
from sglsurvey.manifest import HashCache, combined_hash
from sglsurvey.v2.profile import ArchiveProfile

warnings.filterwarnings("ignore", message="All-NaN slice")
MASKS = ("primary", "strict", "loose")
SCHEMA = "v2-tensor-v2 (per-frame cap)"
_G: dict = {}


def tensor_name(endpoint, role, xt=0.0):
    return f"{endpoint}__{role}.npz" if xt == 0.0 else f"{endpoint}__{role}__xt{xt:+.2f}.npz"


def init_worker(P: ArchiveProfile, obs_by_id, usable, manifest):
    _G["P"] = P
    _G["registry"] = load_target_registry(P.registry_path)
    _G["ctx"] = P.geometry_context()
    _G["zcache"] = {}
    _G["obs_by_id"] = obs_by_id
    _G["usable"] = usable
    _G["usable_set"] = {p: set(v) for p, v in usable.items()}
    _G["manifest"] = manifest


def locus_at_zgrid(endpoint, role, mjd, bin_days):
    P = _G["P"]
    key = (endpoint, role, round(mjd / bin_days))
    zc = _G["zcache"]
    if key not in zc:
        if len(zc) > 200000:
            zc.clear()
        ctx = _G["ctx"]
        al = adaptive_locus(target=_G["registry"][endpoint], role=Role(role),
                            observation_time=Time(key[2] * bin_days, format="mjd"),
                            observer=ctx.observer, relay_range=ctx.relay_range,
                            tolerance_arcsec=P.locus_tol_arcsec, ephemeris=ctx.ephemeris,
                            model=ctx.model)
        zs = np.array([p.z_au for p in al.points]); ra = np.array([p.icrs_ra_deg for p in al.points])
        dec = np.array([p.icrs_dec_deg for p in al.points])
        q = 1.0 / zs; o = np.argsort(q); qg = 1.0 / P.z_grid
        ra_u = np.unwrap(np.deg2rad(ra[o]))
        zc[key] = np.stack([np.rad2deg(np.interp(qg, q[o], ra_u)) % 360.0,
                            np.interp(qg, q[o], dec[o])], axis=1)
    return zc[key]


def tangent_offsets(pts, centre):
    ra0, dec0 = centre
    dra = (pts[:, 0] - ra0 + 180.0) % 360.0 - 180.0
    return np.stack([dra * np.cos(np.deg2rad(dec0)) * 3600.0, (pts[:, 1] - dec0) * 3600.0], axis=1)


def donors_for(P, corridor, role):
    pool = sorted(e for e in P.endpoints if P.corridor_of(e) != corridor)
    seed = int(hashlib.sha256(f"{corridor}/{role}/{P.split_seed}".encode()).hexdigest()[:8], 16)
    rng = np.random.default_rng(seed)
    n = min(P.n_donors, len(pool))
    return sorted(rng.choice(pool, size=n, replace=False).tolist())


def sample_pix(fm, x, y, aux=False):
    ny, nx = fm.flux.shape
    out = [np.full(len(x), np.nan) for _ in range(4 if aux else 3)]
    ok = (x >= 0) & (x <= nx - 1.001) & (y >= 0) & (y <= ny - 1.001) & np.isfinite(x) & np.isfinite(y)
    if ok.any():
        x0 = np.floor(x[ok]).astype(int); y0 = np.floor(y[ok]).astype(int)
        fx, fy = x[ok] - x0, y[ok] - y0

        def bil(a):
            return (a[y0, x0] * (1 - fx) * (1 - fy) + a[y0, x0 + 1] * fx * (1 - fy)
                    + a[y0 + 1, x0] * (1 - fx) * fy + a[y0 + 1, x0 + 1] * fx * fy)
        out[0][ok] = bil(fm.flux); out[1][ok] = bil(fm.var); out[2][ok] = bil(fm.good_frac)
        if aux:
            out[3][ok] = bil(fm.aux) if fm.aux is not None else np.nan
    return out


class Acc:
    """Streaming stack accumulators for one (pair, band): all trajectories x masks."""

    def __init__(self, T, grid, n_masks):
        shp = (n_masks, T) + grid
        self.A = np.zeros(shp); self.B = np.zeros(shp); self.n = np.zeros(shp, np.int32)
        self.cmax = np.full(shp, -np.inf)
        self.Ae = np.zeros(shp); self.Be = np.zeros(shp); self.ne = np.zeros(shp, np.int32)
        self.Al = np.zeros(shp); self.Bl = np.zeros(shp); self.nl = np.zeros(shp, np.int32)
        # real-only per-phase parts
        self.Ap = np.zeros((n_masks, 2) + grid); self.Bp = np.zeros((n_masks, 2) + grid)
        self.np_ = np.zeros((n_masks, 2) + grid, np.int32)

    def add(self, mi, f, w, phase, early, late):
        c = f * w
        self.A[mi] += c; self.B[mi] += w; self.n[mi] += (w > 0)
        np.maximum(self.cmax[mi], c, out=self.cmax[mi])
        if early:
            self.Ae[mi] += c; self.Be[mi] += w; self.ne[mi] += (w > 0)
        elif late:
            self.Al[mi] += c; self.Bl[mi] += w; self.nl[mi] += (w > 0)
        self.Ap[mi, phase] += c[0]; self.Bp[mi, phase] += w[0]; self.np_[mi, phase] += (w[0] > 0)


def _S(A, B, n, min_epochs):
    with np.errstate(all="ignore"):
        S = A / np.sqrt(np.where(B > 0, B, np.inf))
    return np.where(n >= min_epochs, S, np.nan)


def build_corridor(corridor: str, only_missing: bool) -> str:
    P: ArchiveProfile = _G["P"]
    t_start = _time.monotonic()
    obs_by_id, usable, manifest, usable_set = _G["obs_by_id"], _G["usable"], _G["manifest"], _G["usable_set"]
    pairs = sorted(p for p in usable if P.corridor_of(p[0]) == corridor and p[0] in P.endpoints)
    xt_cfg = {}
    xt_path = P.config_dir / "cross_track_cells.json"
    if xt_path.exists():
        for c in json.loads(xt_path.read_text())["cells"]:
            xt_cfg[(c["endpoint"], c["role"])] = [x for x in c["offsets_arcsec"] if x != 0.0]
    pairs = [(e, r, x) for (e, r) in pairs for x in [0.0] + xt_cfg.get((e, r), [])]
    cache = HashCache(P.run_dir / "manifest" / f"cutout_hashes_{corridor}.json")
    if only_missing:
        fresh = []
        for p in pairs:
            out = P.tensor_dir / tensor_name(*p)
            if not out.exists():
                fresh.append(p); continue
            try:
                with np.load(out) as d:
                    recorded = str(d["input_hash"]); oids = [str(o) for o in d["oid"]]
                digs = []
                for oid in oids:
                    for fp in P.cutout_files(oid, obs_by_id[oid], manifest.get(oid)):
                        digs.append(cache.sha256(fp))
                cache.flush()
                if combined_hash(digs) != recorded:
                    print(f"  [{corridor}] STALE {out.name} -> rebuilding", flush=True); fresh.append(p)
            except Exception as exc:
                print(f"  [{corridor}] unreadable {out.name} ({exc}) -> rebuilding", flush=True); fresh.append(p)
        pairs = fresh
    if not pairs:
        return f"{corridor}: nothing to do"
    frames = sorted({o for p in pairs for o in usable[p[:2]]}, key=lambda o: obs_by_id[o]["t_mid_mjd_utc"])
    centre = P.corridor_centre(corridor, manifest, obs_by_id)
    # phase reference: circular median day-of-year of the corridor's frames
    doys = np.array([obs_by_id[o]["t_mid_mjd_utc"] % 365.25 for o in frames])
    ang = doys / 365.25 * 2 * np.pi
    phase_ref = float((np.arctan2(np.sin(ang).mean(), np.cos(ang).mean()) % (2 * np.pi)) / (2 * np.pi) * 365.25)

    def phase_of(mjd):
        d = (mjd % 365.25) - phase_ref
        d = (d + 182.625) % 365.25 - 182.625
        return 0 if abs(d) < 91.3 else 1

    donors = {role: donors_for(P, corridor, role) for role in ("rx", "tx")}
    nz, nm = len(P.z_grid), len(P.mu_grid)
    grid = (nz, nm, nm)
    ring = np.array(P.ring)
    T = 1 + len(ring) + len(donors["rx"])
    keep = [0] + [1 + i for i in P.designated]
    store = {p: [] for p in pairs}
    digests = {p: [] for p in pairs}
    rngs = {p: np.random.default_rng(int(hashlib.sha256(f"{p[0]}/{p[1]}/{P.split_seed}".encode()).hexdigest()[:8], 16)) for p in pairs}
    cosd = np.cos(np.deg2rad(centre[1]))
    eps = 10.0
    probe = np.array([[centre[0], centre[1]], [centre[0] + eps / 3600.0 / cosd, centre[1]], [centre[0], centre[1] + eps / 3600.0]])
    n_maps = 0
    split = P.holdout_split_mjd
    for band in P.bands:
        b = P.band_idx[band]
        bframes = [o for o in frames if obs_by_id[o]["band"] == band]
        if len(bframes) < P.min_epochs:
            continue
        # ---- pass 0: frame weights -> per-frame cap ------------------------------
        wf = {}
        for oid in bframes:
            fm = P.build_map(oid, obs_by_id[oid], manifest.get(oid), False, None)
            if fm is None or fm.magzp is None:
                continue
            scale = 10.0 ** ((P.zp_ref - fm.magzp) / 2.5)
            sel = np.isfinite(fm.var) & (fm.good_frac >= P.min_good_frac) & (fm.var > 0)
            if sel.sum() < 100:
                continue
            wf[oid] = float(1.0 / (np.median(fm.var[sel]) * scale * scale))
            del fm
        if len(wf) < P.min_epochs:
            continue
        cap = P.weight_cap * float(np.median(list(wf.values())))
        # ---- pass 1: sample and accumulate ------------------------------------------
        accs = {p: Acc(T, grid, len(MASKS)) for p in pairs}
        epochs = {p: [] for p in pairs}
        for oid in bframes:
            if oid not in wf:
                continue
            obs = obs_by_id[oid]; row = manifest.get(oid)
            fm = P.build_map(oid, obs, row, False, None)
            if fm is None or fm.magzp is None:
                continue
            n_maps += 1
            mjd = obs["t_mid_mjd_utc"]; scale = 10.0 ** ((P.zp_ref - fm.magzp) / 2.5)
            dt_yr = (mjd - P.t0_mjd) / 365.25
            ph = phase_of(mjd)
            masks_e = [P.quality_ok(obs["quality_flags"], m) for m in MASKS]
            early = split is not None and mjd <= split
            late = split is not None and mjd > split
            pp = fm.world2pix(probe[:, 0], probe[:, 1])
            J = np.stack([(pp[1] - pp[0]) / eps, (pp[2] - pp[0]) / eps], axis=1); p0 = pp[0]
            digs = [cache.sha256(fp) for fp in P.cutout_files(oid, obs, row)]
            for pair in pairs:
                endpoint, role, xt = pair
                if oid not in usable_set[pair[:2]]:
                    continue
                base = locus_at_zgrid(endpoint, role, mjd, P.locus_bin_days)
                d0 = tangent_offsets(base, centre)
                if xt != 0.0:
                    tv = np.gradient(d0, axis=0); nn = np.linalg.norm(tv, axis=1, keepdims=True)
                    tv = np.where(nn > 0, tv / np.where(nn > 0, nn, 1), [[1.0, 0.0]])
                    d0 = d0 + xt * np.stack([-tv[:, 1], tv[:, 0]], axis=1)
                dmu = P.mu_grid * dt_yr
                all_off = np.empty((T, nz, nm, nm, 2))
                all_off[0, ..., 0] = d0[:, 0][:, None, None] + dmu[None, :, None]
                all_off[0, ..., 1] = d0[:, 1][:, None, None] + dmu[None, None, :]
                all_off[1:1 + len(ring)] = all_off[0][None] + ring[:, None, None, None, :]
                base_t0 = locus_at_zgrid(endpoint, role, P.t0_mjd, P.locus_bin_days)
                d0_t0 = tangent_offsets(base_t0, centre)
                for k, dn in enumerate(donors[role]):
                    dl_t = locus_at_zgrid(dn, role, mjd, P.donor_bin_days)
                    dl_0 = locus_at_zgrid(dn, role, P.t0_mjd, P.donor_bin_days)
                    dc = (dl_0[nz // 2, 0], dl_0[nz // 2, 1])
                    delta = tangent_offsets(dl_t, dc) - tangent_offsets(dl_0, dc)
                    kk = 1 + len(ring) + k
                    all_off[kk, ..., 0] = (d0_t0[:, 0] + delta[:, 0])[:, None, None] + dmu[None, :, None]
                    all_off[kk, ..., 1] = (d0_t0[:, 1] + delta[:, 1])[:, None, None] + dmu[None, None, :]
                if P.linear_wcs:
                    pix = p0[None, None, None, None, :] + np.einsum("ij,...j->...i", J, all_off)
                else:
                    flat = all_off.reshape(-1, 2)
                    ra_q = centre[0] + flat[:, 0] / 3600.0 / cosd
                    dec_q = centre[1] + flat[:, 1] / 3600.0
                    pix = fm.world2pix(ra_q, dec_q).reshape(all_off.shape)
                f, v, g = sample_pix(fm, pix[..., 0].ravel(), pix[..., 1].ravel())
                if P.sample_adjust is not None:
                    f, v = P.sample_adjust(fm, obs, all_off.reshape(-1, 2), centre, f, v)
                F = (f.reshape(T, *grid) * scale).astype(np.float32)
                V = np.minimum(v.reshape(T, *grid) * scale * scale, 1e30).astype(np.float32)
                G = g.reshape(T, *grid).astype(np.float32)
                valid = np.isfinite(F) & np.isfinite(V) & (V > 0) & (G >= P.min_good_frac)
                if P.clip_sigma is not None:
                    with np.errstate(all="ignore"):
                        valid &= np.abs(F) / np.sqrt(np.where(V > 0, V, np.inf)) <= P.clip_sigma
                w = np.where(valid, np.minimum(1.0 / np.where(valid, V, 1.0), cap), 0.0).astype(np.float64)
                Fc = np.where(valid, np.nan_to_num(F), 0.0).astype(np.float64)
                acc = accs[pair]
                for mi in range(len(MASKS)):
                    if masks_e[mi]:
                        acc.add(mi, Fc, w, ph, early, late)
                aux = None
                if fm.aux is not None:
                    aux = sample_pix(fm, pix[0, ..., 0].ravel(), pix[0, ..., 1].ravel(), aux=True)[3]
                    aux = aux.reshape(grid).astype(np.float16)
                extra = P.epoch_extra(fm, obs) if P.epoch_extra else {}
                epochs[pair].append({
                    "mjd": mjd, "phase": ph, "masks": masks_e, "f": F[keep], "v": V[keep],
                    "g": G[keep].astype(np.float16), "aux": aux,
                    "xy0": pix[0, :, nm // 2, nm // 2, :].astype(np.float32), "d0": d0.astype(np.float32),
                    "J": J.astype(np.float32), "p0": p0.astype(np.float32),
                    "origin": np.array(getattr(fm, "frame_origin", (0, 0)), np.int32),
                    "magzp": float(fm.magzp), "pix_arcsec": float(getattr(fm, "pix_arcsec", np.nan)),
                    "fwhm_arcsec": float(getattr(fm, "fwhm_arcsec", np.nan)),
                    "frame_w": wf[oid], "oid": oid, "extra": extra,
                    "quality": obs["quality_flags"],
                })
                digests[pair].extend(digs)
            del fm
        cache.flush()
        # ---- summaries per pair ----------------------------------------------------
        for pair in pairs:
            eps_ = epochs[pair]
            if len(eps_) < P.min_epochs:
                continue
            acc = accs[pair]
            summ = np.full((len(MASKS), T, 7), np.nan, np.float32)
            hold = np.full((len(MASKS), T, 7), np.nan, np.float32)
            scr = np.full((len(MASKS), P.n_scramble), np.nan, np.float32)
            escr = np.full((len(MASKS), P.n_epoch_scramble), np.nan, np.float32)
            n_ok = np.zeros(len(MASKS), np.int32)
            masks_b = np.array([e["masks"] for e in eps_], bool)
            for mi in range(len(MASKS)):
                n_ok[mi] = int(masks_b[:, mi].sum())
                if n_ok[mi] < P.min_epochs:
                    continue
                S = _S(acc.A[mi], acc.B[mi], acc.n[mi], P.min_epochs)
                for t in range(T):
                    s_max, node = nulls.grid_max(S[t])
                    if node is None:
                        continue
                    summ[mi, t, 0] = s_max; summ[mi, t, 1] = node[0]; summ[mi, t, 2] = node[1] * nm + node[2]
                    summ[mi, t, 5] = acc.n[mi][t][node]
                    summ[mi, t, 6] = acc.cmax[mi][t][node] / acc.A[mi][t][node] if acc.A[mi][t][node] > 0 else np.nan
                    if t == 0:
                        for p in (0, 1):
                            Bp = acc.Bp[mi, p][node]
                            summ[mi, 0, 3 + p] = acc.Ap[mi, p][node] / np.sqrt(Bp) if Bp > 0 else np.nan
                    if split is not None and acc.ne[mi][t].max() >= P.min_epochs:
                        Se = _S(acc.Ae[mi][t], acc.Be[mi][t], acc.ne[mi][t], P.min_epochs)
                        s_em, nd = nulls.grid_max(Se)
                        if nd is not None:
                            hold[mi, t, 0] = s_em
                            hold[mi, t, 1] = acc.Ae[mi][t][nd] / acc.Be[mi][t][nd]
                            hold[mi, t, 5] = nd[0]; hold[mi, t, 6] = nd[1] * nm + nd[2]
                            nl = acc.nl[mi][t][nd]; hold[mi, t, 4] = nl
                            if nl >= P.min_epochs:
                                Bl = acc.Bl[mi][t][nd]
                                hold[mi, t, 2] = acc.Al[mi][t][nd] / np.sqrt(Bl) if Bl > 0 else np.nan
                                hold[mi, t, 3] = acc.Al[mi][t][nd] / Bl if Bl > 0 else np.nan
                parts = {p: (acc.Ap[mi, p], acc.Bp[mi, p], acc.np_[mi, p]) for p in (0, 1) if acc.np_[mi, p].max() > 0}
                if parts:
                    scr[mi] = nulls.phase_scrambled_maxima(parts, P.n_scramble, rngs[pair])
                f0 = np.stack([e["f"][0] for e in eps_]); v0 = np.stack([e["v"][0] for e in eps_])
                g0 = np.stack([e["g"][0] for e in eps_]).astype(np.float32)
                fc = np.full(len(eps_), cap, np.float32)
                escr[mi] = nulls.epoch_scrambled_maxima(f0, v0, g0, P.n_epoch_scramble, rngs[pair],
                                                        masks_b[:, mi], frame_cap=fc, clip_sigma=P.clip_sigma)
            store[pair].append({"band": b, "eps": eps_, "summary": summ, "holdout": hold, "scramble": scr,
                                "escramble": escr, "n_ok": n_ok, "cap": cap,
                                "accA": acc.A.astype(np.float32), "accB": acc.B.astype(np.float32),
                                "accn": acc.n.astype(np.int32), "accAp": acc.Ap.astype(np.float32),
                                "accBp": acc.Bp.astype(np.float32), "accnp": acc.np_.astype(np.int32)})
            del accs[pair]
        del accs, epochs
    # ---- write per pair ----------------------------------------------------------
    written = []
    for pair in pairs:
        parts = store[pair]
        if not parts:
            continue
        endpoint, role, xt = pair
        eps_all = [(q["band"], e, q["cap"]) for q in parts for e in q["eps"]]
        order = np.argsort([e["mjd"] for _, e, _ in eps_all])
        eps_all = [eps_all[i] for i in order]
        E = len(eps_all)
        nb = len(P.bands)
        summ = np.full((len(MASKS), T, nb, 7), np.nan, np.float32)
        hold = np.full((len(MASKS), T, nb, 7), np.nan, np.float32)
        scr = np.full((len(MASKS), nb, P.n_scramble), np.nan, np.float32)
        escr = np.full((len(MASKS), nb, P.n_epoch_scramble), np.nan, np.float32)
        n_ok = np.zeros((len(MASKS), nb), np.int32)
        for q in parts:
            bi = q["band"] - 1
            summ[:, :, bi] = q["summary"]; hold[:, :, bi] = q["holdout"]
            scr[:, bi] = q["scramble"]; escr[:, bi] = q["escramble"]; n_ok[:, bi] = q["n_ok"]
        has_aux = any(e["aux"] is not None for _, e, _ in eps_all)
        arr = lambda key, dt=None: np.array([e[key] for _, e, _ in eps_all], dtype=dt)  # noqa: E731
        out = P.tensor_dir / tensor_name(endpoint, role, xt)
        P.tensor_dir.mkdir(parents=True, exist_ok=True)
        extra_keys = sorted({k for _, e, _ in eps_all for k in e["extra"]})
        # accumulated stack sums per (mask, trajectory, band, grid) for the joint stage
        acc_fields = {}
        for key in ("accA", "accB", "accn", "accAp", "accBp", "accnp"):
            shp = parts[0][key].shape
            arrs = np.zeros((len(MASKS), shp[1], nb) + tuple(shp[2:]), dtype=parts[0][key].dtype)
            for q in parts:
                arrs[:, :, q["band"] - 1] = q[key]
            acc_fields[key] = arrs
        np.savez_compressed(
            out, schema=SCHEMA, cap_kind="per-frame", input_hash=combined_hash(digests[pair]),
            endpoint=endpoint, role=role, xt_arcsec=xt, corridor=corridor, centre=np.array(centre),
            phase_ref=phase_ref, observer_identity=P.observer_identity,
            z_grid=P.z_grid, mu_grid=P.mu_grid, t0_mjd=P.t0_mjd, zp_ref=P.zp_ref,
            ring_offsets=ring, designated=np.array(P.designated), kept_trajectories=np.array(keep),
            donors=np.array(donors[role]), bands=np.array(P.bands),
            mjd=arr("mjd"), band_idx=np.array([b for b, _, _ in eps_all], np.uint8),
            phase=arr("phase", np.uint8), masks=arr("masks", bool), mask_names=np.array(MASKS),
            frame_cap=np.array([c for _, _, c in eps_all], np.float32), frame_w=arr("frame_w", np.float32),
            magzp=arr("magzp", np.float32), pix_arcsec=arr("pix_arcsec", np.float32),
            fwhm_arcsec=arr("fwhm_arcsec", np.float32), oid=arr("oid"),
            quality=np.array([json.dumps(e["quality"], default=float) for _, e, _ in eps_all]),
            xy0=np.stack([e["xy0"] for _, e, _ in eps_all]), d0=np.stack([e["d0"] for _, e, _ in eps_all]),
            J=np.stack([e["J"] for _, e, _ in eps_all]), p0=np.stack([e["p0"] for _, e, _ in eps_all]),
            origin=np.stack([e["origin"] for _, e, _ in eps_all]),
            f=np.stack([e["f"] for _, e, _ in eps_all], axis=1), v=np.stack([e["v"] for _, e, _ in eps_all], axis=1),
            g=np.stack([e["g"] for _, e, _ in eps_all], axis=1),
            aux=(np.stack([e["aux"] if e["aux"] is not None else np.full(grid, np.nan, np.float16)
                           for _, e, _ in eps_all]) if has_aux else np.zeros(0, np.float16)),
            extra_keys=np.array(extra_keys),
            extra=np.array([[float(e["extra"].get(k, np.nan)) for k in extra_keys] for _, e, _ in eps_all], np.float32),
            summary=summ, holdout=hold, scramble_max=scr, epoch_scramble_max=escr, n_epochs_ok=n_ok,
            **acc_fields,
            summary_fields=np.array(["S_max", "iz", "imu", "S_phase0", "S_phase1", "n_epochs", "top_share"]),
            holdout_fields=np.array(["S_early", "f_early", "S_late", "f_late", "n_late", "iz_early", "imu_early"]),
        )
        written.append(f"{out.stem} ({E} epochs)")
    return f"{corridor}: {n_maps} maps, wrote {', '.join(written)} in {(_time.monotonic() - t_start) / 60:.1f} min"


def run(P: ArchiveProfile, corridors=None, workers=7, only_missing=True):
    obs_by_id, usable, manifest = P.load_inputs()
    corridors = corridors or P.corridors()
    print(f"[{P.name}] {len(corridors)} corridors, {len(usable)} usable pairs", flush=True)
    (P.run_dir / "manifest").mkdir(parents=True, exist_ok=True)
    log = P.run_dir / "build_tensors.log"

    def _run(c):
        try:
            return build_corridor(c, only_missing)
        except Exception as exc:
            import traceback
            return f"{c}: FAILED {exc}\n{traceback.format_exc()}"
    P.only_missing_flag = only_missing
    if workers <= 1:
        init_worker(P, obs_by_id, usable, manifest)
        for c in corridors:
            msg = _run(c); print(msg, flush=True); open(log, "a").write(msg + "\n")
        return
    import multiprocessing as mp
    with mp.get_context("fork").Pool(workers, initializer=init_worker,
                                     initargs=(P, obs_by_id, usable, manifest)) as pool:
        for msg in pool.imap_unordered(_run_global, corridors):
            print(msg, flush=True); open(log, "a").write(msg + "\n")


def _run_global(c):
    try:
        return build_corridor(c, getattr(_G["P"], "only_missing_flag", True))
    except Exception as exc:
        import traceback
        return f"{c}: FAILED {exc}\n{traceback.format_exc()}"
