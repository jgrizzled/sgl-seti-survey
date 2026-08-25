"""Injection completeness (threshold freeze v1.1 §completeness).

Exact stamp-response injections through the identical confirmatory
chain (same fitted Gaussian kernel, same per-cadence matched-filter
maps, same track sampling, same baseline/k/WEIGHT_CAP statistic
recipe), with the SPOC per-camera/CCD PRF as the injection PSF
(archive.stsci.edu prf_fitsfiles start_s0004, nearest 5x5 grid point
to each cutout's CCD position; 13 x 13 px at 9x oversampling,
sha-recorded). The matched filter is linear, so the response of every
statistic to an added F x PRF stamp is computed locally and exactly
(sglsurvey.inject.stamp_response on the per-cadence denom/good/kernel).

Per B unit: magnitude grid T 10-18 (0.5 steps), 200 draws per temporal
model (freeze requires >= 100). A draw picks the injected z from the
recoverable z grid, a shared per-cadence background realization from
the unit's common off-window sample pool (cadences valid and unclipped
at every z node, so cross-z correlation of the nested max is
preserved), and for the pulsed model a period from the unit's frozen
6-point log grid with a uniform phase; the boxcar duty is d = 0.1 with
sub-cadence overlap integrated (x d scaling below one cadence, per
hypotheses v1.0). Both frozen statistics are evaluated for both
temporal models (2 x 2 recovery matrix); recovery is against the
unit's actual frozen threshold max(T, 0) from the confirmatory run.
The quoted m90 pairs each temporal model with its natural statistic
(chord -> S_c, pulse -> S_p). m90 at 18.0 is a grid-censored lower
limit; a unit whose recovery never reaches 90% even at T = 10 is
reported insensitive (m90 None).

Channel-A reference rows (teegarden s71, van-maanen s43): the same
injection machinery on the on-star per-cadence series, but
THRESHOLD-FREE by freeze — depths are quoted against a labeled
S = 5 reference level (no controls exist for zero-trial rows), with
the pulse period grid built by the same 6-point rule (2 x cadence ->
window). Flux scale: B = field-star ZP re-measured here through the
UNIT kernel (finding C1: dev_validation's zp_check hardcoded a
1.5 px kernel, not the per-cube fitted kernel the freeze requires;
the wolf-359 3.0 px maps therefore need their own measurement —
flux-scale-only correction, no search statistic touched); A =
target-star self-calibration against the TIC Tmag (van-maanen
carries the declared +/-0.3 mag scale caveat).

Because the star ZP is measured through the chain, it already
absorbs the kernel<->PRF throughput; every stamp response is
therefore normalized by the same-cadence response at the calibration
reference (the cutout centre — the response is position-independent
to well under the ZP scatter), so a magnitude-m injection produces
exactly the chain-seen flux of a real magnitude-m star.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from astropy.io import fits

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import sglsurvey.photometry as phot
from sglsurvey.inject import OversampledPSF, stamp_response
from sglsurvey.photometry import (_gaussian_kernel, _matched_filter_full,
                                  flux_map_from_arrays)
from tesscut_lib import PRODUCTS, load_cube

import confirmatory_search as cs
import dev_validation as dv

phot.PSF_FWHM_ARCSEC.setdefault("T", 1.5 * 21.0)

D = REPO / "surveys" / "tess-crossings"
OUT = D / "results"
PRF_DIR = REPO / "runs" / "tess-crossings" / "products" / "prf"
PRF_FILES = {
    (1, 1): "tess2019107181900-prf-1-1-row1536-col2092.fits",
    (2, 1): "tess2019107181901-prf-2-1-row0513-col2092.fits",
    (1, 2): "tess2019107181901-prf-1-2-row1536-col0045.fits",
    (1, 4): "tess2019107181901-prf-1-4-row1536-col0045.fits",
}
MAGS = np.arange(10.0, 18.01, 0.5)
NDRAW = 200
DUTY = 0.1
NSUB = 32                      # sub-cadence samples for boxcar overlap
S_REF_A = 5.0                  # labeled reference level, A rows
ZP_NOMINAL = 20.44
RNG = np.random.default_rng(20260824)

FREEZE = cs.FREEZE
Z_GRID = cs.Z_GRID
TMAG = {t["target_id"]: t["tmag"] for t in json.loads(
    (OUT / "tic_cut_v1.json").read_text())["targets"]}
CONF = json.loads((OUT / "confirmatory_v1.json").read_text())
DEVVAL = json.loads((OUT / "dev_validation_v1.json").read_text())


def load_prf(cam, ccd):
    path = PRF_DIR / PRF_FILES[(cam, ccd)]
    plane = np.asarray(fits.getdata(path), float)
    import hashlib
    sha = "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()
    return OversampledPSF(plane=plane, oversample=9, centre=58,
                          band="T"), path.name, sha


def build_maps_full(cube, fwhm_pix):
    """cs.build_maps plus the per-cadence denom and good-pixel mask
    (what stamp_response needs)."""
    kern = _gaussian_kernel(fwhm_pix,
                            max(3, int(np.ceil(2 * fwhm_pix))))
    idx = np.where((cube.quality == 0)
                   & np.isfinite(cube.mjd_utc))[0]
    n = len(idx)
    ny, nx = cube.flux.shape[1:]
    F = np.empty((n, ny, nx), np.float32)
    V = np.empty((n, ny, nx), np.float32)
    G = np.empty((n, ny, nx), np.float32)
    DEN = np.empty((n, ny, nx), np.float32)
    GM = np.empty((n, ny, nx), bool)
    for j, i in enumerate(idx):
        img = np.asarray(cube.flux[i], float)
        err = np.asarray(cube.flux_err[i], float)
        good = np.isfinite(img) & np.isfinite(err) & (err > 0)
        var = np.where(good, err ** 2, 1e30)
        f, v, g, den, _, _ = _matched_filter_full(img, var, good,
                                                  fwhm_pix, True, kern)
        F[j], V[j], G[j], DEN[j], GM[j] = f, v, g, den, good
    return idx, F, V, G, DEN, GM, kern


def clip_keep(fo):
    keep = np.ones(len(fo), bool)
    for _ in range(3):
        med = np.median(fo[keep])
        sd = 1.4826 * np.median(np.abs(fo[keep] - med)) or 1.0
        keep &= np.abs(fo - med) <= 3 * sd
    return keep


def bilinear2d(A, x, y):
    ny, nx = A.shape
    if not (0 <= x <= nx - 1.001 and 0 <= y <= ny - 1.001):
        return np.nan
    xi, yi = int(np.floor(x)), int(np.floor(y))
    fx, fy = x - xi, y - yi
    return float(A[yi, xi] * (1 - fx) * (1 - fy)
                 + A[yi, xi + 1] * fx * (1 - fy)
                 + A[yi + 1, xi] * (1 - fx) * fy
                 + A[yi + 1, xi + 1] * fx * fy)


def zp_unit(cube, kern, fwhm_pix):
    """Field-star ZP through the UNIT kernel (finding C1), sub-pixel
    bilinear sampling of the median-stack matched-filter map."""
    ok = (cube.quality == 0) & np.isfinite(cube.mjd_utc)
    stack = np.nanmedian(cube.flux[ok][::10], axis=0)
    err = np.nanmedian(cube.flux_err[ok][::10], axis=0)
    good = np.isfinite(stack) & np.isfinite(err) & (err > 0)
    var = np.where(good, err ** 2, 1e30)
    f, v, g, *_ = _matched_filter_full(stack, var, good, fwhm_pix,
                                       True, kern)
    ny, nx = stack.shape
    ctr = cube.wcs.wcs_pix2world([[nx / 2, ny / 2]], 0)[0]
    rows = dv.tic_cone(float(ctr[0]), float(ctr[1]),
                       max(nx, ny) * cs.PIX_ARCSEC / 3600.0 / 1.4)
    zps = []
    for r in rows:
        t = r.get("Tmag")
        if t is None or not (9.0 <= float(t) <= 15.0):
            continue
        x, y = cube.wcs.wcs_world2pix([[float(r["ra"]),
                                        float(r["dec"])]], 0)[0]
        if not (2 <= x <= nx - 3 and 2 <= y <= ny - 3):
            continue
        fl = bilinear2d(f, x, y)
        vv = bilinear2d(v, x, y)
        if np.isfinite(fl) and np.isfinite(vv) and fl > 5 * np.sqrt(vv):
            zps.append(float(t) + 2.5 * np.log10(fl))
    zps = np.asarray(zps)
    if len(zps) < 3:
        return {"n_stars": int(len(zps)), "status": "sparse"}
    mad = 1.4826 * np.median(np.abs(zps - np.median(zps)))
    return {"n_stars": int(len(zps)),
            "zp_median": round(float(np.median(zps)), 3),
            "zp_scatter_mad": round(float(mad), 3),
            "offset_vs_nominal": round(float(np.median(zps)
                                             - ZP_NOMINAL), 3),
            "status": "pass" if mad <= 0.2 else "FAIL"}


def m90_from(recovery):
    r = np.asarray(recovery)
    if r[0] < 0.9:
        return None, False
    below = np.where(r < 0.9)[0]
    if not below.size:
        return float(MAGS[-1]), True
    i = below[0]
    frac = (r[i - 1] - 0.9) / max(r[i - 1] - r[i], 1e-9)
    return float(MAGS[i - 1] + frac * (MAGS[i] - MAGS[i - 1])), False


def boxcar_frac(mjd, cadence_s, period_s, phase_s):
    """Overlap fraction of the d=0.1 boxcar train with each cadence."""
    half = cadence_s / 2.0 / 86400.0
    sub = mjd[:, None] + np.linspace(-half, half, NSUB)[None, :]
    ph = ((sub * 86400.0 - phase_s) % period_s) / period_s
    return (ph < DUTY).mean(axis=1)


def period_grid(cadence_s, window_days):
    lo, hi = 2.0 * cadence_s, window_days * 86400.0
    return np.exp(np.linspace(np.log(lo), np.log(hi), 6))


def complete_b_cube(tid, prefix):
    cube = load_cube(sorted(sorted(PRODUCTS.glob(f"{prefix}*"))[0]
                            .glob("*.fits"))[0], tid)
    idx0 = np.where((cube.quality == 0) & np.isfinite(cube.mjd_utc))[0]
    stack = np.nanmedian(cube.flux[idx0][::10], axis=0)
    fwhm, fwhm_src = cs.fit_fwhm_field(stack)
    print(f"[{tid}] fwhm={fwhm:.2f}px ({fwhm_src}); maps + responses",
          flush=True)
    idx, F, V, G, DEN, GM, kern = build_maps_full(cube, fwhm)
    chk = zp_unit(cube, kern, fwhm)
    print(f"[{tid}] unit-kernel zp_check: {chk}", flush=True)
    if chk.get("status") == "sparse":
        zp_info = {"zp": ZP_NOMINAL, "zp_check": chk,
                   "source": f"nominal 20.44 fallback (unit-kernel "
                             f"zp_check sparse, n={chk['n_stars']})"}
    else:
        zp_info = {
            "zp": chk["zp_median"], "zp_check": chk,
            "source": (f"field-star ZP through the unit kernel "
                       f"(finding C1; n={chk['n_stars']}, scatter "
                       f"{chk['zp_scatter_mad']}"
                       + ("; FAILS the 0.2 mag gate - depths carry a "
                          f"declared +/-{chk['zp_scatter_mad']} mag "
                          "scale caveat"
                          if chk["status"] == "FAIL" else "; pass")
                       + ")")}
    mjd = cube.mjd_utc[idx]
    units = [u for u in FREEZE["search_units"]["units"]
             if u["target_id"] == tid]
    ev = cs.event_row(units[0]["event_id"])
    tracks = cs.relay_tracks(ev, mjd)                # (nz, n, 2)
    cosd = np.cos(np.radians(tracks[..., 1]))
    nz = len(Z_GRID)

    # real-track samples per z
    pix = np.empty((nz, len(mjd), 2))
    fS = np.empty((nz, len(mjd)))
    vS = np.empty((nz, len(mjd)))
    gS = np.empty((nz, len(mjd)))
    for zi in range(nz):
        p = cube.wcs.wcs_world2pix(
            np.column_stack([tracks[zi, :, 0], tracks[zi, :, 1]]), 0)
        pix[zi] = p
        fS[zi] = cs.bilinear(F, p[:, 0], p[:, 1])
        vS[zi] = cs.bilinear(V, p[:, 0], p[:, 1])
        gS[zi] = cs.bilinear(G, p[:, 0], p[:, 1])
    valid = (np.isfinite(fS) & np.isfinite(vS) & (vS > 0)
             & np.isfinite(gS) & (gS >= cs.GOOD_MIN))

    prf, prf_name, prf_sha = load_prf(cube.camera, cube.ccd)
    zp = zp_info["zp"]

    rows = []
    for u in units:
        hd = u["window_days"] / 2.0
        lo, hi = u["t_ca_mjd"] - hd, u["t_ca_mjd"] + hd
        inw = (mjd >= lo) & (mjd <= hi)

        # per-z baseline, k, clip mask on off-window samples
        base = np.full(nz, np.nan)
        kz = np.ones(nz)
        offkeep = np.zeros((nz, len(mjd)), bool)
        for zi in range(nz):
            om = valid[zi] & ~inw
            oi = np.where(om)[0]
            keep = clip_keep(fS[zi, oi])
            base[zi] = float(np.median(fS[zi, oi[keep]]))
            if keep.sum() >= 6:
                kz[zi] = max(1.0, float(np.median(
                    (fS[zi, oi[keep]] - base[zi]) ** 2
                    / vS[zi, oi[keep]]) / 0.4549))
            offkeep[zi, oi[keep]] = True

        # common off-window pool (valid + unclipped at every z)
        pool_idx = np.where(offkeep.all(axis=0))[0]
        pool_source = "common_all_z"
        if len(pool_idx) < 50:
            pool_idx = np.where(offkeep[0])[0]
            pool_source = "z550_only_fallback"

        # in-window cadence sets per z + injection responses
        in_z = [np.where(valid[zi] & inw)[0] for zi in range(nz)]
        union = sorted(set().union(*[set(i) for i in in_z]))
        resp = {}                       # cad index -> (nz_inj, nz_meas)
        ny_c, nx_c = F.shape[1:]
        for i in union:
            shim = SimpleNamespace(denom=DEN[i], kernel=kern,
                                   good=GM[i])
            # calibration-reference response: PRF stamp at the cutout
            # centre through the same cadence's chain (the ZP already
            # absorbs the kernel<->PRF throughput; see docstring)
            stamp_c, oxc, oyc = prf.render(nx_c / 2.0, ny_c / 2.0,
                                           half=6)
            r_ref = float(stamp_response(shim, stamp_c, oxc, oyc)
                          .sample([nx_c / 2.0], [ny_c / 2.0])[0])
            if not np.isfinite(r_ref) or r_ref <= 1e-6:
                r_ref = 1.0
            R = np.zeros((nz, nz))
            for zi in range(nz):
                x, y = pix[zi, i]
                if not (0 <= x <= nx_c - 1 and 0 <= y <= ny_c - 1):
                    continue
                stamp, ox, oy = prf.render(x, y, half=6)
                rw = stamp_response(shim, stamp, ox, oy)
                R[zi] = rw.sample(pix[:, i, 0], pix[:, i, 1]) / r_ref
            resp[i] = R
        alive = [zi for zi in range(nz)
                 if any(resp[i][zi].max() > 1e-6 for i in union)]
        if not alive:
            rows.append({"target_id": tid, "radius_au": u["radius_au"],
                         "status": "all_z_masked"})
            continue

        # per-z precomputations for the chord statistic
        chord = []                # per z: dict(w, sw, mjd, kv, R, loc)
        for zi in range(nz):
            ii = in_z[zi]
            w = 1.0 / (vS[zi, ii] * kz[zi])
            w = np.minimum(w, cs.WEIGHT_CAP * np.median(w))
            chord.append({
                "ii": ii, "w": w, "sw": float(np.sqrt(w.sum())),
                "kv": np.sqrt(vS[zi, ii] * kz[zi]),
                "R": np.array([[resp[i][zj, zi] for i in ii]
                               for zj in range(nz)]),  # (nz_inj, n_in)
            })

        conf = next(r for r in CONF["units"]
                    if r["target_id"] == tid
                    and r["radius_au"] == u["radius_au"])
        th_c = max(conf["T_c"], 0.0)
        th_p = max(conf["T_p"], 0.0)
        pgrid = np.asarray(u["pulse_periods_s"])

        rec = {m: {"S_c": np.zeros(len(MAGS)),
                   "S_p": np.zeros(len(MAGS))}
               for m in ("chord", "pulse")}
        FLUX = 10 ** (0.4 * (zp - MAGS))
        for _ in range(NDRAW):
            zi_inj = alive[int(RNG.integers(len(alive)))]
            draw = RNG.choice(pool_idx, size=len(union))
            j_of = {c: draw[a] for a, c in enumerate(union)}
            P = float(pgrid[int(RNG.integers(len(pgrid)))])
            phase = float(RNG.uniform(0.0, P))
            fr_pulse = {c: f for c, f in zip(union, boxcar_frac(
                mjd[np.asarray(union)], u["cadence_s"], P, phase))}
            for model in ("chord", "pulse"):
                a_all, b_all = [], []
                Sc = np.full((nz, len(MAGS)), -np.inf)
                for zm in range(nz):
                    ch = chord[zm]
                    ii = ch["ii"]
                    if not len(ii):
                        continue
                    bg = fS[zm, [j_of[c] for c in ii]] - base[zm]
                    fr = (np.ones(len(ii)) if model == "chord" else
                          np.asarray([fr_pulse[c] for c in ii]))
                    sig = ch["R"][zi_inj] * fr
                    Sc[zm] = (float(np.sum(ch["w"] * bg)) / ch["sw"]
                              + FLUX * float(np.sum(ch["w"] * sig))
                              / ch["sw"])
                    a_all.append(bg / ch["kv"])
                    b_all.append(sig / ch["kv"])
                Scm = Sc.max(axis=0)
                a = np.concatenate(a_all)
                b = np.concatenate(b_all)
                Spm = (a[None, :] + FLUX[:, None] * b[None, :]).max(
                    axis=1)
                rec[model]["S_c"] += Scm > th_c
                rec[model]["S_p"] += Spm > th_p
        for m in rec:
            for s in rec[m]:
                rec[m][s] /= NDRAW
        m90_c, cen_c = m90_from(rec["chord"]["S_c"])
        m90_p, cen_p = m90_from(rec["pulse"]["S_p"])
        row = {
            "target_id": tid, "event_id": u["event_id"],
            "radius_au": u["radius_au"], "b_rsun": u["b_rsun"],
            "sector": u["sector"], "window_days": u["window_days"],
            "status": "ok",
            "threshold_c": round(th_c, 3), "threshold_p": round(th_p, 3),
            "zp_used": round(zp, 3), "zp_source": zp_info["source"],
            "m90_chord": m90_c, "m90_chord_censored": cen_c,
            "m90_pulse": m90_p, "m90_pulse_censored": cen_p,
            "recovery": {m: {s: [round(float(x), 3)
                                 for x in rec[m][s]]
                             for s in rec[m]} for m in rec},
            "z_alive": [float(Z_GRID[z]) for z in alive],
            "z_dead": [float(Z_GRID[z]) for z in range(nz)
                       if z not in alive],
            "n_in_union": len(union), "n_off_pool": len(pool_idx),
            "pool_source": pool_source,
            "median_unit_response": round(float(np.median(
                [resp[i][zi, zi] for i in union
                 for zi in range(nz)])), 3),
            "prf": prf_name, "prf_sha256": prf_sha,
        }
        rows.append(row)
        print(json.dumps({k: row[k] for k in
                          ("target_id", "radius_au", "m90_chord",
                           "m90_pulse", "z_dead")}), flush=True)
    return rows


def complete_a_row(row):
    tid, sec = row["target_id"], row["sector"]
    cube = dv.cube_for(f"A-{tid}-s{sec:04d}")
    ok = np.where((cube.quality == 0) & np.isfinite(cube.mjd_utc))[0]
    stack = np.nanmedian(cube.flux[ok][::10], axis=0)
    ny, nx = stack.shape
    cx, cy = nx // 2, ny // 2
    fwhm = dv.fit_fwhm(stack, cx, cy)
    kern = _gaussian_kernel(fwhm, max(3, int(np.ceil(2 * fwhm))))
    prf, prf_name, prf_sha = load_prf(cube.camera, cube.ccd)
    stamp, ox, oy = prf.render(float(cx), float(cy), half=6)
    f_out, v_out, r_out = [], [], []
    for i in ok:
        img = np.asarray(cube.flux[i], float)
        err = np.asarray(cube.flux_err[i], float)
        good = np.isfinite(img) & np.isfinite(err) & (err > 0)
        var = np.where(good, err ** 2, 1e30)
        fm = flux_map_from_arrays(img, var, good, cube.wcs, "T",
                                  float(cube.mjd_utc[i]), kernel=kern,
                                  pix_arcsec=cs.PIX_ARCSEC,
                                  keep_inputs=True)
        f_out.append(float(fm.flux[cy, cx]))
        v_out.append(float(fm.var[cy, cx]))
        r_out.append(float(stamp_response(fm, stamp, ox, oy)
                           .sample([cx], [cy])[0]))
    mjd = cube.mjd_utc[ok]
    f = np.asarray(f_out)
    v = np.asarray(v_out)
    # injection position == calibration position (the star), so the
    # normalized response is exactly 1; keep the raw response as a
    # throughput diagnostic only
    R_raw = np.asarray(r_out)
    R = np.ones_like(R_raw)
    hd = row["window_days"] / 2.0
    lo, hi = row["t_ca_mjd"] - hd, row["t_ca_mjd"] + hd
    inw = (mjd >= lo) & (mjd <= hi)
    oi = np.where(~inw)[0]
    keep = clip_keep(f[oi])
    base = float(np.median(f[oi[keep]]))
    k = max(1.0, float(np.median((f[oi[keep]] - base) ** 2
                                 / v[oi[keep]]) / 0.4549)) \
        if keep.sum() >= 6 else 1.0
    pool = oi[keep]
    ii = np.where(inw)[0]
    w = 1.0 / (v[ii] * k)
    w = np.minimum(w, cs.WEIGHT_CAP * np.median(w))
    sw = float(np.sqrt(w.sum()))
    kv = np.sqrt(v[ii] * k)
    zp = TMAG[tid] + 2.5 * np.log10(max(base, 1e-3))
    pgrid = period_grid(row["cadence_s"], row["window_days"])
    FLUX = 10 ** (0.4 * (zp - MAGS))
    rec = {m: {"S_c": np.zeros(len(MAGS)), "S_p": np.zeros(len(MAGS))}
           for m in ("chord", "pulse")}
    for _ in range(NDRAW):
        draw = RNG.choice(pool, size=len(ii))
        bg = f[draw] - base
        P = float(pgrid[int(RNG.integers(len(pgrid)))])
        phase = float(RNG.uniform(0.0, P))
        frp = boxcar_frac(mjd[ii], row["cadence_s"], P, phase)
        for model, fr in (("chord", np.ones(len(ii))), ("pulse", frp)):
            sig = R[ii] * fr
            Sc = (float(np.sum(w * bg)) / sw
                  + FLUX * float(np.sum(w * sig)) / sw)
            Sp = ((bg / kv)[None, :]
                  + FLUX[:, None] * (sig / kv)[None, :]).max(axis=1)
            rec[model]["S_c"] += Sc > S_REF_A
            rec[model]["S_p"] += Sp > S_REF_A
    for m in rec:
        for s in rec[m]:
            rec[m][s] /= NDRAW
    m90_c, cen_c = m90_from(rec["chord"]["S_c"])
    m90_p, cen_p = m90_from(rec["pulse"]["S_p"])
    out = {"target_id": tid, "sector": sec,
           "event_id": row["event_id"], "b_rsun": row["b_rsun"],
           "window_days": row["window_days"],
           "depth_kind": f"reference_S{S_REF_A:.0f}_threshold_free",
           "zp_selfcal": round(float(zp), 3),
           "tic_tmag": TMAG[tid],
           "star_flux_e_s": round(base, 1), "k": round(k, 2),
           "m90_chord": m90_c, "m90_chord_censored": cen_c,
           "m90_pulse": m90_p, "m90_pulse_censored": cen_p,
           "recovery": {m: {s: [round(float(x), 3) for x in rec[m][s]]
                            for s in rec[m]} for m in rec},
           "pulse_period_grid_s": [round(float(p), 1) for p in pgrid],
           "median_raw_response": round(float(np.median(R_raw[ii])),
                                        3),
           "n_in": int(len(ii)), "n_off_pool": int(len(pool)),
           "prf": prf_name, "prf_sha256": prf_sha}
    if tid == "van-maanen":
        out["scale_caveat"] = ("self-cal offset -0.31 mag at FWHM-fit "
                               "ceiling (dev validation): depths "
                               "carry +/-0.3 mag")
    print(json.dumps({k: out[k] for k in
                      ("target_id", "m90_chord", "m90_pulse",
                       "zp_selfcal")}), flush=True)
    return out


def main():
    out = {"freeze": FREEZE["freeze_content_hash"],
           "mag_grid": [float(m) for m in MAGS], "n_draws": NDRAW,
           "duty": DUTY, "seed": 20260824,
           "z_grid_au": list(Z_GRID),
           "reference_level_a": S_REF_A,
           "note": ("2x2 recovery per unit: temporal model (chord d=1"
                    " / pulse d=0.1 boxcar) x statistic (S_c / S_p); "
                    "quoted m90 pairs chord->S_c, pulse->S_p; "
                    "backgrounds are iid draws from the common "
                    "off-window pool, so the sector-scale drift enters"
                    " through the frozen thresholds (real-control "
                    "maxima), not the synthetic background - the "
                    "programme's established construction"),
           "finding_c1": ("dev_validation's zp_check used a hardcoded "
                          "1.5 px kernel; the freeze requires the "
                          "identical (per-cube fitted) kernel, so the "
                          "flux scale is re-measured here through the "
                          "unit kernel with responses normalized at "
                          "the calibration reference. Flux-scale-only "
                          "correction - no search statistic touched. "
                          "Dev-stage value for the record: "
                          "B-wolf-359-s42 zp 20.251 (1.5 px kernel)"),
           "B": [], "A": []}
    for tid, prefix in cs.CUBE_PREFIX.items():
        out["B"].extend(complete_b_cube(tid, prefix))
    for row in FREEZE["a_reference_rows"]["rows"]:
        out["A"].append(complete_a_row(row))

    (OUT / "completeness_v1.json").write_text(
        json.dumps(out, indent=1, default=float) + "\n")
    print("\nsummary:")
    for r in out["B"]:
        print(" B", r["target_id"], r["radius_au"],
              "m90_c:", r.get("m90_chord"),
              "m90_p:", r.get("m90_pulse"))
    for r in out["A"]:
        print(" A", r["target_id"], "m90_c:", r.get("m90_chord"),
              "m90_p:", r.get("m90_pulse"), f"({r['depth_kind']})")


if __name__ == "__main__":
    main()
