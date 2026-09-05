"""Spectral-archive family v1 search (thresholds §2-§4).

usage: search.py dev|confirmatory [--checks]

Per star x instrument combo used by the family: build the null Ensemble
(seeded draw from units_v1.json), null statistics (LOO), thresholds
T_line (per cell) and T_coadd (per unit multiplicity, 60 seeded
complement-template draws); then the in-window spectra of every unit:
S_line, S_coadd, exceedances, automatic adjudication ladder.

Outputs: results/<family>_search_v1.json, results/<family>_v1.md,
runs/spectral-archives/v1/ensembles/<combo>.npz (template, sigma,
null stats) for completeness.py.
"""
from __future__ import annotations

import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import yaml
from astropy.time import Time
from scipy.optimize import curve_fit

sys.path.insert(0, str(Path(__file__).resolve().parent))
from spec_lib import CFG, C_KMS, STELLAR_REGIONS_NM, Ensemble, Grid, despike, gauss_kernel, load_spectrum, normalise, stellar_line_mask, to_grid, xcorr_shift  # noqa: E402

HERE = Path(__file__).resolve().parents[1]
REPO = HERE.parents[1]
RUN = REPO / "runs" / "spectral-archives" / "v1"
DATA = RUN / "data"
ENS = RUN / "ensembles"
units = json.load(open(HERE / "results" / "units_v1.json"))
REG = yaml.safe_load(open(REPO / "registries" / "pilot_wise_2026.yaml"))["targets"]
STAR_RV = {t: float(REG[t]["state"]["astrometry"]["radial_velocity_km_s"]) for t in ("wolf-359", "ross-128", "teegarden", "ross-154", "gj-908")}
SKY_LINES_NM_VAC = {  # compact observer-frame list (vacuum): airglow + Calar Alto / La Silla / Maunakea lamps
    "OI_5577": 557.89, "OI_6300": 630.20, "OI_6364": 636.56, "NaI_5890": 589.16, "NaI_5896": 589.76,
    "Hg_4047": 404.77, "Hg_4358": 435.96, "Hg_5461": 546.23, "Hg_5770": 577.12, "Hg_5791": 579.23,
    "HeI_10830_sky": 1083.33, "OH_P1_10?": None,
}
SKY_LINES_NM_VAC = {k: v for k, v in SKY_LINES_NM_VAC.items() if v}
STELLAR = CFG["stellar_lines_nm_vacuum"]


def load_norm(inst, grid, rec):
    path = DATA / inst / rec["file"]
    if not path.exists():
        return None
    try:
        sp = load_spectrum(inst, path)
    except Exception as exc:
        return {"error": repr(exc)}
    f, e, m = to_grid(sp, grid)
    nf, ne, nm = normalise(f, e, m)
    with np.errstate(invalid="ignore", divide="ignore"):
        snr = np.where(m, np.nan, f / e)
    cell_snr = {cl: float(np.nanmedian(snr[grid.cell_slice(cl)])) if np.isfinite(snr[grid.cell_slice(cl)]).any() else 0.0
                for cl in CFG["cells_nm"]}
    return {"F": nf, "E": ne, "M": nm, "mjd_mid": sp["mjd_mid"], "meta": sp["meta"], "raw_flux": f, "cell_snr": cell_snr}


def berv_sign_check(inst, grid, recs):
    """D8 (v1.1 form): stellar-region cross-correlation between the two nulls
    with the largest |dBERV| under three wavelength treatments -- the applied
    barycentric shift must bring the stellar features to ~0 px; 'none' must
    show ~dBERV (or 0 if the product is already barycentric)."""
    loaded = []
    for r in recs[:40]:
        sp_path = DATA / inst / r["file"]
        if sp_path.exists():
            try:
                sp = load_spectrum(inst, sp_path)
            except Exception as exc:
                print(f"   [berv check] skip {sp_path.name}: {exc!r}", flush=True); continue
            if sp["meta"].get("berv") is not None:
                loaded.append(sp)
    if len(loaded) < 2:
        return {"skipped": "too few / no BERV"}
    bervs = np.array([s["meta"]["berv"] for s in loaded])
    i, j = int(np.argmin(bervs)), int(np.argmax(bervs))
    db = float(bervs[j] - bervs[i])
    if abs(db) < 15:
        return {"skipped": f"max |dBERV| {abs(db):.1f} km/s < 15"}
    already_bary = loaded[0]["meta"].get("specsys") == "BARYCENT"
    def variant(sp, mode):
        berv = sp["meta"]["berv"]
        if already_bary:
            fac = {"applied": 1.0, "none": 1 / (1 + berv / C_KMS), "flipped": (1 - berv / C_KMS)}[mode]
        else:
            fac = {"applied": 1.0, "none": 1 / (1 + berv / C_KMS), "flipped": (1 - berv / C_KMS) / (1 + berv / C_KMS)}[mode]
        s2 = dict(sp); s2["orders"] = [(w * fac, f, e, m) for (w, f, e, m) in sp["orders"]]
        return normalise(*to_grid(s2, grid))[0]
    out = {"berv_a": float(bervs[i]), "berv_b": float(bervs[j]), "dberv_px": db / C_KMS * 3 * grid.R, "regions": {}}
    for mode in ("none", "applied", "flipped"):
        fa, fb = variant(loaded[i], mode), variant(loaded[j], mode)
        out["regions"][mode] = {f"{lo}-{hi}": xcorr_shift(fa[int(grid.idx(lo)):int(grid.idx(hi))], fb[int(grid.idx(lo)):int(grid.idx(hi))])
                                for lo, hi in STELLAR_REGIONS_NM.get(inst, [(grid.w[grid.n // 3], grid.w[grid.n // 3] * 1.03)])}
    shifts_applied = np.abs(list(out["regions"]["applied"].values()))
    out["verdict"] = "applied sign correct" if shifts_applied.max() <= 3 else "CHECK: stellar features not aligned"
    return out


def gauss(x, a, mu, sig, c):
    return a * np.exp(-0.5 * ((x - mu) / sig) ** 2) + c


def shape_fit(z, peak, fwhm_px):
    lo, hi = max(0, peak - 8), min(len(z), peak + 9)
    x = np.arange(lo, hi); y = z[lo:hi]
    ok = np.isfinite(y)
    if ok.sum() < 6:
        return None
    try:
        p, _ = curve_fit(gauss, x[ok], y[ok], p0=[np.nanmax(y), peak, fwhm_px / 2.3548, 0.0], maxfev=4000)
        return {"amp_z": float(p[0]), "mu_px": float(p[1]), "fwhm_px": float(abs(p[2]) * 2.3548), "fwhm_ratio": float(abs(p[2]) * 2.3548 / fwhm_px)}
    except Exception:
        return None


def adjudicate(ens, grid, unit, rec, S, peak, T, others_S, berv):
    """Frozen ladder (hypotheses §8). Returns disposition + evidence."""
    lam = float(grid.w[peak])
    ev = {"lambda_bary_vac_nm": lam}
    z = ens.z_single(rec["F"], e=rec["E"])
    fit = shape_fit(z, peak, grid.fwhm_px)   # re-fit for the record (statistic already PSF-consistent)
    ev["shape"] = fit
    lo_r, hi_r = CFG["shape_fwhm_bounds"]
    if fit is None or not (lo_r <= fit["fwhm_ratio"] <= hi_r):
        return "defect/cosmic", ev
    if ens.mask_frac[peak] >= 0.2:
        ev["mask_frac"] = float(ens.mask_frac[peak]); return "defect/cosmic", ev
    rv = STAR_RV[unit["target"]]
    tol = CFG["stellar_line_tolerance_km_s"]
    for name, l0 in STELLAR.items():
        dv = (lam / (l0 * (1 + rv / C_KMS)) - 1) * C_KMS
        if abs(dv) <= tol:
            ev["stellar_line"] = {"line": name, "dv_km_s": float(dv)}; return "stellar_flare", ev
    lam_obs = lam / (1 + (berv or 0.0) / C_KMS)
    res_el = lam_obs / grid.R
    for name, l0 in SKY_LINES_NM_VAC.items():
        if abs(lam_obs - l0) <= res_el:
            ev["sky_line"] = {"line": name, "lambda_obs_nm": lam_obs}; return "telluric", ev
    sl = slice(max(0, peak - 5), peak + 6)
    med_sig = np.nanmedian(ens.sigma[ens.pix_ok])
    if np.nanmax(ens.sigma[sl]) > 3 * med_sig:
        ev["sigma_ratio"] = float(np.nanmax(ens.sigma[sl]) / med_sig); return "telluric", ev
    # recurrence within the unit
    rec_hits = [k for k, So in others_S.items() if np.nanmax(So[max(0, peak - 2):peak + 3]) > T]
    ev["recurrence_in_unit"] = rec_hits
    return "retained-ambiguous", ev


def main():
    fam = sys.argv[1]
    do_checks = "--checks" in sys.argv
    ENS.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(CFG["seed"])
    frozen_sha_file = HERE / "results" / ("threshold_freeze_v1_1_sha.txt" if (HERE / "results" / "threshold_freeze_v1_1_sha.txt").exists() else "threshold_freeze_v1_sha.txt")
    fam_units = [u for u in units["units"] if u["family"] == fam]
    combos = sorted(set(f"{u['target']}|{u['instrument']}" for u in fam_units))
    results = {"family": fam, "frozen_sha": open(frozen_sha_file).read().split()[0], "frozen_config": frozen_sha_file.name,
               "run_utc": Time.now().isot, "combos": {}, "units": [], "machinery_checks": {}}
    t0 = time.time()
    for key in combos:
        tid, inst = key.split("|")
        c = units["combos"][key]
        grid = Grid(inst)
        print(f"\n== {key}: grid {grid.n} px, nulls {len(c['null_draw'])}", flush=True)
        if do_checks:
            chk = berv_sign_check(inst, grid, c["null_draw"])
            results["machinery_checks"][key] = {"berv_sign": chk}
            print("   BERV check:", chk, flush=True)
        # ---- null ensemble
        NF, NE, NM, nmeta = [], [], [], []
        for r in c["null_draw"]:
            L = load_norm(inst, grid, r)
            if L is None or "error" in L:
                nmeta.append({"file": r["file"], "status": "load_failed" if L else "missing", "error": (L or {}).get("error")}); continue
            frac_ok = np.isfinite(L["F"]).mean()
            if frac_ok < 0.5:
                nmeta.append({"file": r["file"], "status": "too_masked", "frac_ok": float(frac_ok)}); continue
            if L["cell_snr"]["generic"] < CFG.get("pixel_snr_min", 5.0):
                nmeta.append({"file": r["file"], "status": "low_snr", "cell_snr": L["cell_snr"]}); continue
            NF.append(L["F"]); NE.append(L["E"]); NM.append(L["M"])
            nmeta.append({"file": r["file"], "status": "ok", "mjd_mid": L["mjd_mid"], "snr": L["meta"].get("snr"), "berv": L["meta"].get("berv")})
        nmeta_ok = [m for m in nmeta if m["status"] == "ok"]
        n_ok = len(NF)
        if n_ok < CFG["min_null_pool"]:
            results["combos"][key] = {"status": "insufficient_nulls", "n_ok": n_ok, "nulls": nmeta}
            print(f"   only {n_ok} usable nulls -> combo not searchable", flush=True)
            continue
        ens = Ensemble(grid, NF, NE, NM, extra_mask=stellar_line_mask(grid, STAR_RV[tid]))
        # v1.1 A1: wavelength-solution gate on the nulls (stellar-region shift vs template <= 3 px)
        regions = STELLAR_REGIONS_NM.get(inst, [])
        keep = []
        for i in range(ens.N):
            sh, _ = ens.wavelength_shift_px(ens.F[i], regions)
            nmeta_ok[i]["wave_shift_px"] = sh
            if abs(sh) <= CFG.get("wave_shift_max_px", 3):
                keep.append(i)
            else:
                nmeta_ok[i]["status"] = "wavelength_solution_failed"
        if len(keep) < ens.N:
            print(f"   dropped {ens.N - len(keep)} null(s) failing the wavelength-solution gate", flush=True)
            ens = Ensemble(grid, [NF[i] for i in keep], [NE[i] for i in keep], [NM[i] for i in keep], extra_mask=stellar_line_mask(grid, STAR_RV[tid]))
            n_ok = len(keep)
            if n_ok < CFG["min_null_pool"]:
                results["combos"][key] = {"status": "insufficient_nulls_after_gate", "n_ok": n_ok, "nulls": nmeta}
                continue
        cells = next(u["cells"] for u in fam_units if u["instrument"] == inst and u["target"] == tid)
        nstats, npeaks = ens.null_stats(cells)
        T_line = {cl: float(np.nanmax(nstats[cl])) for cl in cells}
        null_dist = {cl: {"median": float(np.nanmedian(nstats[cl])), "p90": float(np.nanpercentile(nstats[cl], 90)), "max": T_line[cl],
                          "argmax_file": nmeta_ok[int(np.nanargmax(nstats[cl]))]["file"] if np.isfinite(nstats[cl]).any() else None,
                          "argmax_lambda_nm": float(grid.w[npeaks[cl][int(np.nanargmax(nstats[cl]))]]) if np.isfinite(nstats[cl]).any() else None} for cl in cells}
        cell_mask_frac = {cl: float(np.nanmean(ens.mask_frac[grid.cell_slice(cl)])) for cl in cells}
        phot_ratio = {cl: float(np.nanmedian((ens.sigma / np.nanmedian(ens.E, axis=0))[grid.cell_slice(cl)])) for cl in cells}
        results["combos"][key] = {"status": "ok", "n_nulls": n_ok, "grid_px": grid.n, "T_line": T_line, "null_S_distribution": null_dist,
                                  "stellar_mask_frac": float(ens.extra_mask.mean()),
                                  "null_S_line": {cl: [None if not np.isfinite(v) else float(v) for v in nstats[cl]] for cl in cells},
                                  "cell_mask_frac": cell_mask_frac, "sigma_over_photon": phot_ratio, "nulls": nmeta}
        np.savez_compressed(ENS / f"{tid}_{inst}.npz", template=ens.template, sigma=ens.sigma, pix_ok=ens.pix_ok, mad=ens.mad, extra_mask=ens.extra_mask,
                            F=ens.F, E=ens.E, T_line=json.dumps(T_line), cells=json.dumps(cells), grid_w=grid.w)
        print(f"   nulls ok {n_ok}; T_line {T_line}; null dist {null_dist}; sigma/photon {phot_ratio}; mask {cell_mask_frac}  [{time.time()-t0:.0f}s]", flush=True)

        # ---- units of this combo
        for u in [u for u in fam_units if u["instrument"] == inst and u["target"] == tid]:
            ures = {"unit_id": u["unit_id"], "target": tid, "instrument": inst, "event_id": u["event_id"], "t_ca_utc": u["t_ca_utc"],
                    "b_rsun": u["b_rsun"], "rungs": u["rungs"], "cells": cells, "spectra": [], "S_line": {}, "T_line": T_line,
                    "exceedances": [], "trials": 0}
            loaded = []
            for s in u["spectra"]:
                L = load_norm(inst, grid, s)
                srec = {"file": s["file"], "utc": s["utc"], "dt_days": s["nearest_dt_days"], "in_window": s["in_window"]}
                if L is None or "error" in L:
                    srec["status"] = "load_failed"; ures["spectra"].append(srec); continue
                if np.isfinite(L["F"]).mean() < 0.5:
                    srec["status"] = "too_masked"; ures["spectra"].append(srec); continue
                srec["cell_snr"] = L["cell_snr"]
                if L["cell_snr"]["generic"] < CFG.get("pixel_snr_min", 5.0):
                    srec["status"] = "low_snr"; ures["spectra"].append(srec); continue
                srec["mjd_mid_header"] = L["mjd_mid"]
                sh, _ = ens.wavelength_shift_px(L["F"], STELLAR_REGIONS_NM.get(inst, []))
                srec["wave_shift_px"] = sh
                if abs(sh) > CFG.get("wave_shift_max_px", 3):
                    srec["status"] = "wavelength_solution_failed"; ures["spectra"].append(srec); continue
                z = ens.z_single(L["F"], e=L["E"])
                S = ens.matched(z)
                srec["status"] = "ok"; srec["snr"] = L["meta"].get("snr"); srec["berv"] = L["meta"].get("berv")
                srec["S_cell"] = {}
                for cl in cells:
                    v, p = ens.cell_max(S, cl, z)
                    srec["S_cell"][cl] = {"S": v, "peak_px": p, "lambda_nm": None if p is None else float(grid.w[p])}
                srec["n_spikes_removed"] = int(getattr(ens, "last_spikes", 0))
                loaded.append((s, L, S, srec))
                ures["spectra"].append(srec)
            n_ok_u = len(loaded)
            ures["n_usable"] = n_ok_u
            if n_ok_u == 0:
                ures["status"] = "no_usable_spectra"; results["units"].append(ures); continue
            ures["status"] = "searched"
            # S_line
            for cl in cells:
                vals = [(sr["S_cell"][cl]["S"], i) for i, (_, _, _, sr) in enumerate(loaded) if sr["S_cell"][cl]["S"] is not None and np.isfinite(sr["S_cell"][cl]["S"])]
                if not vals:
                    ures["S_line"][cl] = None; continue
                v, i = max(vals)
                ures["S_line"][cl] = {"S": v, "spectrum": loaded[i][3]["file"], "peak_px": loaded[i][3]["S_cell"][cl]["peak_px"],
                                      "lambda_nm": loaded[i][3]["S_cell"][cl]["lambda_nm"], "T": T_line[cl], "exceeds": bool(v > T_line[cl])}
                ures["trials"] += 1
                if v > T_line[cl] and ures["S_line"][cl]["peak_px"] is not None:
                    others = {loaded[j][3]["file"]: loaded[j][2] for j in range(n_ok_u) if j != i}
                    disp, ev = adjudicate(ens, grid, u, loaded[i][1], loaded[i][2], ures["S_line"][cl]["peak_px"], T_line[cl], others, loaded[i][3].get("berv"))
                    ures["exceedances"].append({"statistic": "S_line", "cell": cl, "S": v, "T": T_line[cl], "spectrum": loaded[i][3]["file"],
                                                "utc": loaded[i][3]["utc"], "disposition": disp, "evidence": ev})
            # S_coadd
            if n_ok_u >= 2:
                zc = ens.z_mean([L["F"] for _, L, _, _ in loaded], [L["E"] for _, L, _, _ in loaded])
                Sc = ens.matched(zc)
                # null draws with complement templates
                draws = {cl: [] for cl in cells}
                n_draw = CFG["n_coadd_draws"]
                for d in range(n_draw):
                    if n_ok_u <= ens.N // 2:
                        idx = rng.choice(ens.N, size=n_ok_u, replace=False)
                    else:
                        idx = rng.choice(ens.N, size=n_ok_u, replace=True)
                    comp = np.setdiff1d(np.arange(ens.N), idx)
                    tmpl = np.nanmedian(ens.F[comp], axis=0) if len(comp) >= 5 else ens.template
                    Fm = ens.F[idx]
                    n_eff = np.isfinite(Fm).sum(axis=0)
                    sig_d = np.fmax(ens.mad, np.sqrt(np.nanmean(ens.E[idx] ** 2, axis=0)))
                    zd = (np.nanmean(Fm, axis=0) - tmpl) / (sig_d / np.sqrt(np.fmax(n_eff, 1)))
                    zd[(n_eff < max(1, n_ok_u // 2)) | ~ens.pix_ok | ens.extra_mask] = np.nan
                    zd, _ = despike(zd)
                    Sd = ens.matched(zd)
                    for cl in cells:
                        draws[cl].append(ens.cell_max(Sd, cl, zd)[0])
                ures["S_coadd"] = {}
                for cl in cells:
                    v, p = ens.cell_max(Sc, cl, zc)
                    Tc = float(np.nanmax(draws[cl]))
                    ures["S_coadd"][cl] = {"S": v, "T": Tc, "peak_px": p, "lambda_nm": None if p is None else float(grid.w[p]),
                                           "exceeds": bool(np.isfinite(v) and v > Tc), "n": n_ok_u, "null_draws": [None if not np.isfinite(x) else float(x) for x in draws[cl]]}
                    ures["trials"] += 1
                    if np.isfinite(v) and v > Tc and p is not None:
                        # adjudicate on the spectrum with the largest single S at the coadd peak
                        best = max(range(n_ok_u), key=lambda j: np.nanmax(loaded[j][2][max(0, p - 2):p + 3]))
                        others = {loaded[j][3]["file"]: loaded[j][2] for j in range(n_ok_u) if j != best}
                        disp, ev = adjudicate(ens, grid, u, loaded[best][1], loaded[best][2], p, Tc, others, loaded[best][3].get("berv"))
                        ures["exceedances"].append({"statistic": "S_coadd", "cell": cl, "S": v, "T": Tc, "spectrum": loaded[best][3]["file"],
                                                    "disposition": disp, "evidence": ev})
            # v1.1 A7: S_line compares the max over n in-window spectra with the max over N nulls -> P = n/(N+n); coadd 1/(draws+1)
            n_line_trials = sum(1 for cl in cells if ures["S_line"].get(cl))
            n_coadd_trials = len(cells) if n_ok_u >= 2 else 0
            ures["expected_exceedances"] = n_line_trials * n_ok_u / (ens.N + n_ok_u) + n_coadd_trials / (CFG["n_coadd_draws"] + 1)
            ures["n_nulls"] = ens.N
            results["units"].append(ures)
            print(f"   unit {u['unit_id']}: usable {n_ok_u}/{u['n_spectra']}; S_line " +
                  ", ".join(f"{cl}={ures['S_line'][cl]['S']:.2f}/{T_line[cl]:.2f}" for cl in cells if ures["S_line"].get(cl)) +
                  (" ; S_coadd " + ", ".join(f"{cl}={ures['S_coadd'][cl]['S']:.2f}/{ures['S_coadd'][cl]['T']:.2f}" for cl in cells) if n_ok_u >= 2 else "") +
                  f" ; exceedances {len(ures['exceedances'])}", flush=True)

    # ---- summary
    n_tr = sum(u.get("trials", 0) for u in results["units"])
    n_ex = sum(len(u.get("exceedances", [])) for u in results["units"])
    exp = sum(u.get("expected_exceedances", 0) for u in results["units"])
    results["summary"] = {"units": len(results["units"]), "searched": sum(1 for u in results["units"] if u.get("status") == "searched"),
                          "trials": n_tr, "exceedances": n_ex, "expected_exceedances": exp,
                          "dispositions": dict(defaultdict(int, {d: sum(1 for u in results["units"] for e in u.get("exceedances", []) if e["disposition"] == d)
                                                                 for d in set(e["disposition"] for u in results["units"] for e in u.get("exceedances", []))}))}
    def clean(o):
        if isinstance(o, dict): return {k: clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)): return [clean(v) for v in o]
        if isinstance(o, (np.floating,)): return None if not np.isfinite(o) else float(o)
        if isinstance(o, (np.integer,)): return int(o)
        if isinstance(o, np.bool_): return bool(o)
        if isinstance(o, float) and not np.isfinite(o): return None
        return o
    out = HERE / "results" / f"{fam}_search_v1.json"
    out.write_text(json.dumps(clean(results), indent=1))
    print("\nSUMMARY", results["summary"], f"[{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
