"""Static-sky template for the SPHEREx pilot (the substitute for ZTF's
reference-subtracted difference images).

SPHEREx has no difference products, and static stars below the single-
epoch clip stack coherently along every trajectory (real and control)
over hundreds of epochs, so the stack would be confusion-limited at
~16-17 AB. A relay at 550-10,000 AU moves 20-375" per year and spends
only a few percent of the baseline at any fixed sky node, whereas a
static source is always there, with a flux that depends on the
exposure wavelength (linear-variable filter) but not on time. We
therefore fit, per corridor x detector and per sky node on a 3" tangent-
plane grid, f(lambda) = a + b (lambda - lambda0) across epochs by
weighted least squares (two passes, 3-sigma clipped), and the sampling
stage subtracts the template from each epoch's matched-filter flux.
Nodes with fewer than MIN_EPOCHS epochs get no template (flagged).

Writes runs/spherex/calib_v1/templates/<corridor>__<band>.npz.

Usage: uv run python surveys/spherex/scripts/static_template.py
"""

from __future__ import annotations

import json
import os
import sys
import time as _time
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.io import fits

from sglsurvey.adapters.irsa_spherex import (SpherexExactFootprint,
                                             wavelength_at)
from sglsurvey.photometry import build_flux_map_spherex
from sglsurvey.records import read_records

sys.path.insert(0, str(Path(__file__).resolve().parent))
from spherex_corridors import CORRIDOR_OF  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
COARSE_DIR = REPO / "runs" / "spherex" / "coarse_v1"
PRECISE_DIR = REPO / "runs" / "spherex" / "precise_v1"
OUT_DIR = REPO / "runs" / "spherex" / os.environ.get("SPHEREX_CALIB_RUN", "calib_v1") / "templates"

SPACING_ARCSEC = 3.0
MIN_EPOCHS = 3       # median-equivalent floor for any subtraction
MIN_EPOCHS_SLOPE = 8  # below this, fit the constant term only
CLIP_SIGMA = 3.0
MIN_GOOD_FRAC = 0.7


def tangent_grid(ra_c, dec_c, half_extent_arcsec, spacing):
    n = int(np.ceil(2 * half_extent_arcsec / spacing)) + 1
    xi = (np.arange(n) - (n - 1) / 2) * spacing
    XI, ETA = np.meshgrid(xi, xi)
    return xi, XI, ETA


def inverse_gnomonic(xi_arcsec, eta_arcsec, ra_c, dec_c):
    x = np.deg2rad(np.asarray(xi_arcsec) / 3600.0)
    y = np.deg2rad(np.asarray(eta_arcsec) / 3600.0)
    ra0, dec0 = np.deg2rad(ra_c), np.deg2rad(dec_c)
    rho = np.hypot(x, y)
    c = np.arctan(rho)
    with np.errstate(invalid="ignore", divide="ignore"):
        dec = np.arcsin(np.cos(c) * np.sin(dec0) + np.where(rho > 0, y * np.sin(c) * np.cos(dec0) / rho, 0))
        ra = ra0 + np.arctan2(x * np.sin(c), rho * np.cos(dec0) * np.cos(c) - y * np.sin(dec0) * np.sin(c))
    return np.rad2deg(ra) % 360.0, np.rad2deg(dec)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    obs_by_id = {r["observation_id"]: r for r in read_records(
        COARSE_DIR / "records" / "observation.jsonl")}
    usable_obs = set()
    for r in read_records(PRECISE_DIR / "records" / "intersection_evaluation.jsonl"):
        if r["stage"] == "precise" and r["usable"] in ("usable", "partial"):
            usable_obs.add(r["observation_id"])
    groups = defaultdict(list)
    with (PRECISE_DIR / "records" / "cutout_index.jsonl").open() as fh:
        for line in fh:
            rec = json.loads(line)
            if rec["available"] and rec["observation_id"] in usable_obs:
                groups[(CORRIDOR_OF[rec["endpoints"][0]], rec["band"])].append(rec)

    for (corridor, band), recs in sorted(groups.items()):
        t0 = _time.monotonic()
        ra_c = float(np.median([r["center_ra_deg"] for r in recs]))
        dec_c = float(np.median([r["center_dec_deg"] for r in recs]))
        cosd = np.cos(np.deg2rad(dec_c))
        dra = np.array([(r["center_ra_deg"] - ra_c) * cosd * 3600 for r in recs])
        ddec = np.array([(r["center_dec_deg"] - dec_c) * 3600 for r in recs])
        half = float(np.max(np.hypot(dra, ddec))) + 0.5 * 100 * 6.15 + 10
        xi, XI, ETA = tangent_grid(ra_c, dec_c, half, SPACING_ARCSEC)
        RA, DEC = inverse_gnomonic(XI.ravel(), ETA.ravel(), ra_c, dec_c)
        nn = RA.size
        # --- pass 1 ----------------------------------------------------
        epochs = []
        for r in recs:
            obs = obs_by_id[r["observation_id"]]
            path = REPO / r["path"]
            try:
                fm = build_flux_map_spherex(path, SpherexExactFootprint.FATAL_MASK,
                                            obs["band"], obs["t_mid_mjd_utc"])
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    with fits.open(path) as hdul:
                        hdr = hdul["IMAGE"].header
                        ny, nx = hdul["IMAGE"].data.shape
                        lam, _ = wavelength_at(hdul["WCS-WAVE"].data, hdr, nx / 2, ny / 2)
            except Exception as exc:
                print(f"  fluxmap FAIL {r['observation_id']}: {exc}", flush=True)
                continue
            f, v, g = fm.sample(RA, DEC)
            ok = np.isfinite(f) & np.isfinite(v) & (v > 0) & (g >= MIN_GOOD_FRAC)
            w = np.where(ok, 1.0 / np.where(v > 0, v, 1.0), 0.0)
            epochs.append((lam, np.where(ok, f, 0.0).astype(np.float32), w.astype(np.float32)))
        if not epochs:
            continue
        lam0 = float(np.median([e[0] for e in epochs]))

        def fit(mask_fn=None):
            S0 = np.zeros(nn); S1 = np.zeros(nn); S2 = np.zeros(nn)
            F0 = np.zeros(nn); F1 = np.zeros(nn); N = np.zeros(nn, dtype=np.int32)
            for k, (lam, f, w) in enumerate(epochs):
                wk = w if mask_fn is None else w * mask_fn(k, lam, f)
                dl = lam - lam0
                S0 += wk; S1 += wk * dl; S2 += wk * dl * dl
                F0 += wk * f; F1 += wk * f * dl; N += (wk > 0)
            with np.errstate(invalid="ignore", divide="ignore"):
                det = S0 * S2 - S1 * S1
                b = np.where(det > 0, (S0 * F1 - S1 * F0) / det, 0.0)
                a_slope = np.where(det > 0, (S2 * F0 - S1 * F1) / det, np.nan)
                a_const = np.where(S0 > 0, F0 / S0, np.nan)
            a = np.where(N >= MIN_EPOCHS_SLOPE, a_slope, a_const)
            b = np.where(N >= MIN_EPOCHS_SLOPE, b, 0.0)
            a[N < MIN_EPOCHS] = np.nan
            b[N < MIN_EPOCHS] = 0.0
            return a, b, N

        a, b, N = fit()

        def clip_mask(k, lam, f):
            lam_, f_, w_ = epochs[k]
            with np.errstate(invalid="ignore"):
                resid = (f_ - (a + b * (lam_ - lam0))) * np.sqrt(w_)
            return (np.abs(np.nan_to_num(resid)) <= CLIP_SIGMA).astype(np.float32)

        a, b, N = fit(clip_mask)
        # residual chi2 per node for diagnostics
        chi2 = np.zeros(nn); Nc = np.zeros(nn)
        for lam, f, w in epochs:
            with np.errstate(invalid="ignore"):
                r_ = (f - (a + b * (lam - lam0))) ** 2 * w
            ok = np.isfinite(r_) & (w > 0)
            chi2 += np.where(ok, r_, 0.0); Nc += ok
        with np.errstate(invalid="ignore", divide="ignore"):
            rchi2 = np.where(Nc > 2, chi2 / np.maximum(Nc - 2, 1), np.nan)
        side = XI.shape[0]
        np.savez_compressed(
            OUT_DIR / f"{corridor}__{band}.npz",
            ra_c=ra_c, dec_c=dec_c, spacing_arcsec=SPACING_ARCSEC, xi=xi,
            lam0=lam0, a=a.reshape(side, side).astype(np.float32),
            b=b.reshape(side, side).astype(np.float32),
            n=N.reshape(side, side), rchi2=rchi2.reshape(side, side).astype(np.float32),
            n_epochs=len(epochs))
        fin = np.isfinite(a)
        print(f"[{corridor}/{band}] {len(epochs)} epochs, grid {side}x{side}, "
              f"template defined on {fin.mean():.2f} of nodes "
              f"(median N={np.median(N[fin]) if fin.any() else 0:.0f}), "
              f"median reduced chi2 {np.nanmedian(rchi2):.2f} "
              f"({_time.monotonic() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
