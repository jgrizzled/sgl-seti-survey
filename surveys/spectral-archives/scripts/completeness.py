"""Injection completeness and power conversion (hypotheses D7, thresholds §6).

usage: completeness.py dev|confirmatory

Per searched unit x cell: inject Gaussian lines (instrumental FWHM) at
50 seeded random unmasked positions per amplitude step (12 log steps,
0.02..2.0 x local pseudo-continuum) into each usable in-window
spectrum; recovered iff the local matched-filter S > T_line. A90 per
spectrum by linear interpolation of the recovery curve; unit A90 = the
best spectrum (continuous-emission hypothesis: recovery in any spectrum
is recovery). Line flux F = 1.0645 A90 Flambda(lambda) FWHM_lambda with
Flambda from Gaia DR3 + 2MASS photometry (log-log interpolation);
P = F pi b_rung^2 for each rung of the unit.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import convolve1d

sys.path.insert(0, str(Path(__file__).resolve().parent))
from spec_lib import CFG, Grid, gauss_fit, gauss_kernel, load_spectrum, normalise, to_grid  # noqa: E402

HERE = Path(__file__).resolve().parents[1]
REPO = HERE.parents[1]
RUN = REPO / "runs" / "spectral-archives" / "v1"
DATA = RUN / "data"; ENS = RUN / "ensembles"
units = json.load(open(HERE / "results" / "units_v1.json"))
phot = json.load(open(HERE / "results" / "photometry_v1.json"))
AMPS = np.logspace(np.log10(0.02), np.log10(2.0), 12)
N_POS = 50
# Vega zero points F_lambda(m=0) [erg/s/cm2/A] and effective wavelengths [nm]
ZP = {"BP": (511.0, 4.084e-9), "G": (621.7, 2.495e-9), "RP": (777.0, 1.269e-9),
      "J": (1235.0, 3.129e-10), "H": (1662.0, 1.133e-10), "K": (2159.0, 4.283e-11)}


def nirps_sed(target):
    """Median absolute-flux stellar SED (erg/s/cm2/A vs nm) from the target's
    NIRPS phase-3 spectra in hand (v1.1: NIRPS calibration checked against
    2MASS J/H to 0.12-0.17 mag). Telluric-corrected flux, 10-nm bins each
    needing >= 200 pixels with transmission > 0.95; opaque bands are bridged
    by log interpolation. None if no NIRPS files for the target."""
    import glob
    from astropy.io import fits
    files = sorted(glob.glob(str(DATA / "NIRPS" / "*.fits")))
    if not files:
        return None
    edges = np.arange(965, 1930, 10.0); centres = 0.5 * (edges[1:] + edges[:-1])
    acc = []
    for f in files[:40]:
        h = fits.open(f); obj = str(h[0].header.get("OBJECT", "")).upper()
        if target == "wolf-359" and not any(k in obj for k in ("406", "WOLF")):
            continue
        if target == "teegarden" and not any(k in obj for k in ("TEEG", "GAT", "SO0253", "J02530")):
            continue
        if target not in ("wolf-359", "teegarden"):
            continue
        d = h[1].data; w = d["WAVE"][0] / 10
        fl = d["FLUX_TELL_CAL"][0] if "FLUX_TELL_CAL" in d.columns.names else d["FLUX_CAL"][0]
        tr = d["ATM_TRANSM"][0] if "ATM_TRANSM" in d.columns.names else np.ones_like(fl)
        ok = np.isfinite(fl) & (fl > 0) & np.isfinite(tr) & (tr > 0.95)
        row = np.full(len(centres), np.nan)
        idx = np.digitize(w[ok], edges) - 1
        for i in range(len(centres)):
            v = fl[ok][idx == i]
            if len(v) >= 200:
                row[i] = np.median(v)
        acc.append(row)
    if not acc:
        return None
    med = np.nanmedian(np.asarray(acc), axis=0)
    good = np.isfinite(med) & (med > 0.2 * np.nanmedian(med))     # drop opaque-band bins the transmission model lets through
    if good.sum() < 10:
        return None
    filled = np.exp(np.interp(centres, centres[good], np.log(med[good])))
    return centres, filled


# SED-shape correction for NIR cells where no measured SED exists: ratio of the
# measured wolf-359 NIRPS SED to the log-log photometric interpolation
# (2.7 at 1060 nm, 1.5 at 1550 nm); applied to the other stars as a conservative
# (limit-weakening) factor. Optical cells: no correction, +-30 % systematic declared.
SED_CORR = {"1064": 2.7, "1550": 1.5, "generic_nir": 2.0}


def flam_interp(target):
    p = phot[target]
    mags = {"G": float(p["gaia"]["phot_g_mean_mag"]), "BP": float(p["gaia"]["phot_bp_mean_mag"]), "RP": float(p["gaia"]["phot_rp_mean_mag"]),
            "J": float(p["twomass"]["j_m"]), "H": float(p["twomass"]["h_m"]), "K": float(p["twomass"]["k_m"])}
    pts = sorted((ZP[b][0], ZP[b][1] * 10 ** (-0.4 * m)) for b, m in mags.items())
    lw = np.log([x[0] for x in pts]); lf = np.log([x[1] for x in pts])
    sed = nirps_sed(target)
    def f(lam_nm, cell=None):
        if sed is not None and 966 <= lam_nm <= 1923:
            return float(np.interp(lam_nm, sed[0], sed[1]))
        base = float(np.exp(np.interp(np.log(lam_nm), lw, lf)))
        if lam_nm > 950:
            return base * SED_CORR.get(cell if cell in ("1064", "1550") else "generic_nir", 2.0)
        return base
    f.source = "NIRPS FLUX_CAL median SED (measured)" if sed is not None else "Gaia+2MASS log-log interpolation (+SED_CORR in the NIR)"
    return f, mags


def local_S(z, sigma, k, pos, A, fwhm_px=2.5):
    """Matched-filter S at pos after injecting amplitude A (normalised units);
    A3: the injected peak must pass the same Gaussian-shape test as the
    statistic, else it is not recovered (returns -inf)."""
    half = len(k) // 2
    lo, hi = pos - 2 * half - 2, pos + 2 * half + 3
    if lo < 0 or hi > len(z):
        return np.nan
    zz = z[lo:hi].copy()
    inj = np.zeros_like(zz); c = pos - lo
    inj[c - half:c + half + 1] = A * k / k.max()
    zz = zz + inj / sigma[lo:hi]
    fin = np.isfinite(zz)
    if fin.mean() < 0.6:
        return np.nan
    zf = np.where(fin, zz, 0.0); wgt = fin.astype(float)
    num = convolve1d(zf, k, mode="constant"); den = np.sqrt(convolve1d(wgt, k ** 2, mode="constant"))
    with np.errstate(invalid="ignore", divide="ignore"):
        S = num / den
    loc = S[c - 2:c + 3]
    if not np.isfinite(loc).any():
        return -np.inf
    pk = c - 2 + int(np.nanargmax(loc))
    fit = gauss_fit(zz, pk, fwhm_px)
    lo_r, hi_r = CFG["shape_fwhm_bounds"]
    if fit is None or not (lo_r <= fit["fwhm_ratio"] <= hi_r) or abs(fit["mu_px"] - pk) > 2.0:
        return -np.inf
    return float(loc.max() if np.isfinite(loc).all() else np.nanmax(loc))


def main():
    fam = sys.argv[1]
    res = json.load(open(HERE / "results" / f"{fam}_search_v1.json"))
    rng = np.random.default_rng(CFG["seed"] + 1)
    out = {"family": fam, "amplitudes": AMPS.tolist(), "n_pos": N_POS, "units": []}
    for u in res["units"]:
        if u.get("status") != "searched":
            continue
        inst = u["instrument"]; tid = u["target"]
        grid = Grid(inst); k = gauss_kernel(grid.fwhm_px)
        E = np.load(ENS / f"{tid}_{inst}.npz")
        template, sigma, pix_ok = E["template"], E["sigma"], E["pix_ok"]
        mad = E["mad"] if "mad" in E.files else sigma
        stellar_mask = E["extra_mask"] if "extra_mask" in E.files else np.zeros_like(pix_ok)
        T_line = json.loads(str(E["T_line"]))
        flam, mags = flam_interp(tid)
        ures = {"unit_id": u["unit_id"], "rungs": u["rungs"], "cells": {}, "photometry": mags}
        # load usable in-window spectra
        specs = []
        for s in u["spectra"]:
            if s.get("status") != "ok":
                continue
            sp = load_spectrum(inst, DATA / inst / s["file"])
            F, Er, M = to_grid(sp, grid); nf, ne, nm = normalise(F, Er, M)
            sig_i = np.fmax(mad, np.where(np.isfinite(ne), ne, np.nan))
            z = (nf - template) / sig_i; z[~pix_ok | stellar_mask] = np.nan
            specs.append((s["file"], z, nf, sig_i))
        for cl in u["cells"]:
            sl = grid.cell_slice(cl)
            cand = np.arange(sl.start, sl.stop)
            T = T_line[cl]
            per_spec = []
            for fn, z, nf, sig_i in specs:
                okpos = cand[np.isfinite(z[cand])]
                if len(okpos) < 100:
                    per_spec.append({"file": fn, "A90": None, "reason": "too few unmasked px"}); continue
                rec = []
                for A in AMPS:
                    pos = rng.choice(okpos, size=N_POS, replace=len(okpos) < N_POS)
                    hits = [local_S(z, sig_i, k, int(p), A, grid.fwhm_px) > T for p in pos]
                    rec.append(float(np.mean(hits)))
                rec = np.array(rec)
                # first crossing of 0.9 (monotone-ish); linear interpolation in log A
                A90 = None
                if rec[-1] >= 0.9:
                    idx = int(np.argmax(rec >= 0.9))
                    if idx == 0:
                        A90 = float(AMPS[0])
                    else:
                        a0, a1 = np.log(AMPS[idx - 1]), np.log(AMPS[idx]); r0, r1 = rec[idx - 1], rec[idx]
                        A90 = float(np.exp(a0 + (0.9 - r0) * (a1 - a0) / max(r1 - r0, 1e-9)))
                per_spec.append({"file": fn, "recovery": rec.tolist(), "A90": A90})
            a90s = [p["A90"] for p in per_spec if p.get("A90")]
            cres = {"T_line": T, "per_spectrum": per_spec, "A90_best": min(a90s) if a90s else None,
                    "A90_median": float(np.median(a90s)) if a90s else None, "n_spectra": len(per_spec)}
            # power conversion at the cell's reference wavelength (geometric centre) and at the cell edges
            iv = CFG["cells_nm"][cl]
            lam_ref = float(np.sqrt(iv[0] * iv[1])) if iv else float(np.sqrt(grid.w[sl.start] * grid.w[sl.stop - 1]))
            lams = [lam_ref] + ([iv[0], iv[1]] if iv else [float(grid.w[sl.start]), float(grid.w[sl.stop - 1])])
            cres["lambda_ref_nm"] = lam_ref
            cres["Flambda_ref_erg_s_cm2_A"] = flam(lam_ref, cl)
            cres["Flambda_source"] = flam.source
            if cres["A90_best"] is not None:
                cres["power_W"] = {}
                for rung in u["rungs"]:
                    b = CFG["rung_radius_m"][rung]
                    P = {}
                    for lam in lams:
                        fwhm_A = 10.0 * lam / grid.R          # Angstrom
                        Fline = 1.0645 * cres["A90_best"] * flam(lam, cl) * fwhm_A     # erg/s/cm2
                        P[f"{lam:.0f}"] = Fline * 1e-3 * np.pi * b ** 2               # W
                    cres["power_W"][rung] = {"at_ref": P[f"{lam_ref:.0f}"], "range_over_cell": [min(P.values()), max(P.values())]}
            ures["cells"][cl] = cres
        out["units"].append(ures)
        print(f"{u['unit_id']:36s} " + "  ".join(f"{cl}: A90 {ures['cells'][cl]['A90_best'] if ures['cells'][cl]['A90_best'] is None else round(ures['cells'][cl]['A90_best'],3)}"
                                              + (" P=" + ", ".join(f"{r} {v['at_ref']:.3g} W" for r, v in ures['cells'][cl].get('power_W', {}).items()) if ures['cells'][cl].get('power_W') else "")
                                              for cl in u["cells"]), flush=True)
    (HERE / "results" / f"completeness_{fam}_v1.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
