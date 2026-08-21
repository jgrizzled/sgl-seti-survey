"""Positive control for the SPHEREx pilot: flux-scale recovery of
catalogued stars through the same matched-filter estimator used on the
SGL tracks (plan §6 step 6; the lesson of the ZTF asteroid control).

For each corridor: snapshot 2MASS PSC, CatWISE2020 and Gaia DR3 sources
within the cutout field; select isolated, well-measured stars
(2MASS ph_qual AAA, cc_flg 000, 11 <= J <= 14.5, no 2MASS neighbour
within 15"); for every usable cutout, measure the matched-filter flux
at each control star's position and compare with the catalogue SED
interpolated (log F_nu vs log lambda) to the exposure's wavelength at
that pixel. Gaia RP (0.78 um), J, H, Ks, W1, W2 bracket all six
detectors. Writes runs/spherex/control_v1/{stars.jsonl, measurements
.jsonl, summary.json}; the summary's throughput_offset_mag is consumed
by injection_calibrate.py in the conservative direction.

Usage: uv run python surveys/spherex/scripts/star_control.py [--max-cutouts N]
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import sys
import warnings
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import requests
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning

from sglsurvey.adapters.irsa_spherex import (SpherexExactFootprint,
                                             wavelength_at)
from sglsurvey.photometry import build_flux_map_spherex
from sglsurvey.records import read_records
from sglsurvey.snapshots import SnapshotStore

sys.path.insert(0, str(Path(__file__).resolve().parent))
from spherex_corridors import CORRIDOR_OF  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
COARSE_DIR = REPO / "runs" / "spherex" / "coarse_v1"
PRECISE_DIR = REPO / "runs" / "spherex" / "precise_v1"
RUN_DIR = REPO / "runs" / "spherex" / "control_v1"
TAP_SYNC = "https://irsa.ipac.caltech.edu/TAP/sync"

SEARCH_RADIUS_DEG = 0.12
J_RANGE = (11.0, 14.5)
ISOLATION_ARCSEC = 15.0
MIN_GOOD_FRAC = 0.95
MAX_PER_CORRIDOR_DEFAULT = 400
#: Vega zero points (Jy) and effective wavelengths (um).
ZP = {"RP": (2111.0, 0.783), "J": (1594.0, 1.235), "H": (1024.0, 1.662),
      "K": (666.7, 2.159), "W1": (309.54, 3.368), "W2": (171.79, 4.618)}


def tap(session, store, query):
    request_utc = datetime.now(timezone.utc).isoformat()
    r = session.get(TAP_SYNC, params={"QUERY": query, "FORMAT": "CSV"},
                    timeout=900)
    r.raise_for_status()
    if r.text.lstrip().startswith("<"):
        raise RuntimeError(r.text[:300])
    rows = list(csv.DictReader(io.StringIO(r.text)))
    store.store(service_url=TAP_SYNC, query=query, request_utc=request_utc,
                response_bytes=r.content, row_count=len(rows),
                http_status=r.status_code)
    return rows


def sep_arcsec(ra1, dec1, ra2, dec2):
    c = np.cos(np.deg2rad(dec1))
    return 3600.0 * np.hypot((np.asarray(ra2) - ra1) * c, np.asarray(dec2) - dec1)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-cutouts", type=int, default=MAX_PER_CORRIDOR_DEFAULT)
    args = ap.parse_args()
    session = requests.Session()
    store = SnapshotStore(RUN_DIR)
    RUN_DIR.mkdir(parents=True, exist_ok=True)

    obs_by_id = {r["observation_id"]: r for r in read_records(
        COARSE_DIR / "records" / "observation.jsonl")}
    cut_by_corridor = defaultdict(list)
    with (PRECISE_DIR / "records" / "cutout_index.jsonl").open() as fh:
        for line in fh:
            rec = json.loads(line)
            if rec["available"]:
                cut_by_corridor[CORRIDOR_OF[rec["endpoints"][0]]].append(rec)

    stars_out = (RUN_DIR / "stars.jsonl").open("w")
    meas_out = (RUN_DIR / "measurements.jsonl").open("w")
    all_dm = defaultdict(list)
    for corridor in sorted(cut_by_corridor):
        cuts = cut_by_corridor.get(corridor, [])
        if not cuts:
            continue
        cra = float(np.median([c["center_ra_deg"] for c in cuts]))
        cdec = float(np.median([c["center_dec_deg"] for c in cuts]))
        circle = f"CIRCLE('ICRS',{cra:.6f},{cdec:.6f},{SEARCH_RADIUS_DEG})"
        tm = tap(session, store,
                 "SELECT designation, ra, dec, j_m, h_m, k_m, ph_qual, cc_flg "
                 f"FROM fp_psc WHERE CONTAINS(POINT('ICRS',ra,dec),{circle})=1")
        cw = tap(session, store,
                 "SELECT source_name, ra, dec, w1mpro, w2mpro, w1sigmpro, "
                 "w2sigmpro, cc_flags, ab_flags FROM catwise_2020 "
                 f"WHERE CONTAINS(POINT('ICRS',ra,dec),{circle})=1")
        ga = tap(session, store,
                 "SELECT source_id, ra, dec, phot_rp_mean_mag, phot_g_mean_mag, "
                 "pmra, pmdec FROM gaia_dr3_source "
                 f"WHERE CONTAINS(POINT('ICRS',ra,dec),{circle})=1 "
                 "AND phot_rp_mean_mag < 16")
        print(f"[{corridor}] 2MASS {len(tm)}, CatWISE {len(cw)}, Gaia {len(ga)} "
              f"around ({cra:.4f}, {cdec:.4f})", flush=True)
        tra = np.array([float(r["ra"]) for r in tm])
        tdec = np.array([float(r["dec"]) for r in tm])
        cwra = np.array([float(r["ra"]) for r in cw]) if cw else np.zeros(0)
        cwdec = np.array([float(r["dec"]) for r in cw]) if cw else np.zeros(0)
        gra = np.array([float(r["ra"]) for r in ga]) if ga else np.zeros(0)
        gdec = np.array([float(r["dec"]) for r in ga]) if ga else np.zeros(0)
        stars = []
        for i, r in enumerate(tm):
            try:
                j, h, k = float(r["j_m"]), float(r["h_m"]), float(r["k_m"])
            except ValueError:
                continue
            if r["ph_qual"] != "AAA" or r["cc_flg"] != "000":
                continue
            if not (J_RANGE[0] <= j <= J_RANGE[1]):
                continue
            d = sep_arcsec(tra[i], tdec[i], tra, tdec)
            d[i] = np.inf
            if d.min() < ISOLATION_ARCSEC:
                continue
            sed = {"J": j, "H": h, "K": k}
            if len(cw):
                dc = sep_arcsec(tra[i], tdec[i], cwra, cwdec)
                jc = int(np.argmin(dc))
                if dc[jc] < 3.0:
                    try:
                        w1, w2 = float(cw[jc]["w1mpro"]), float(cw[jc]["w2mpro"])
                        if cw[jc]["cc_flags"].strip("0") == "":
                            sed["W1"], sed["W2"] = w1, w2
                    except ValueError:
                        pass
            if len(ga):
                dg = sep_arcsec(tra[i], tdec[i], gra, gdec)
                jg = int(np.argmin(dg))
                if dg[jg] < 2.0:
                    try:
                        sed["RP"] = float(ga[jg]["phot_rp_mean_mag"])
                    except ValueError:
                        pass
            lam = np.array([ZP[b][1] for b in sed])
            fnu = np.array([ZP[b][0] * 10 ** (-0.4 * sed[b]) for b in sed])
            o = np.argsort(lam)
            star = {"corridor": corridor, "designation": r["designation"],
                    "ra": tra[i], "dec": tdec[i], "sed_mag": sed,
                    "lam_um": lam[o].tolist(), "fnu_jy": fnu[o].tolist()}
            stars.append(star)
            stars_out.write(json.dumps(star) + "\n")
        print(f"[{corridor}] {len(stars)} control stars", flush=True)
        if not stars:
            continue
        sra = np.array([s["ra"] for s in stars])
        sdec = np.array([s["dec"] for s in stars])
        rng = np.random.default_rng(1)
        sel = cuts if len(cuts) <= args.max_cutouts else list(
            rng.choice(cuts, size=args.max_cutouts, replace=False))
        n_meas = 0
        for ci in sel:
            path = REPO / ci["path"]
            obs = obs_by_id[ci["observation_id"]]
            try:
                fm = build_flux_map_spherex(path, SpherexExactFootprint.FATAL_MASK,
                                            obs["band"], obs["t_mid_mjd_utc"])
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", AstropyWarning)
                    with fits.open(path) as hdul:
                        hdr = hdul["IMAGE"].header
                        wave_tab = hdul["WCS-WAVE"].data
            except Exception as exc:
                print(f"  fluxmap FAIL {ci['observation_id']}: {exc}")
                continue
            pix = fm.world2pix(sra, sdec)
            f, v, g = fm.sample(sra, sdec)
            for k in range(len(stars)):
                if not (np.isfinite(f[k]) and g[k] >= MIN_GOOD_FRAC and v[k] > 0):
                    continue
                x, y = pix[k]
                ny, nx = fm.flux.shape
                if not (6 <= x <= nx - 7 and 6 <= y <= ny - 7):
                    continue
                lam, bw = wavelength_at(wave_tab, hdr, x, y)
                s = stars[k]
                lamg, fg = np.array(s["lam_um"]), np.array(s["fnu_jy"])
                if lam < lamg[0] * 0.9 or lam > lamg[-1] * 1.1:
                    continue
                pred_ujy = 1e6 * 10 ** np.interp(np.log10(lam), np.log10(lamg),
                                                 np.log10(fg))
                snr = f[k] / np.sqrt(v[k])
                if snr < 10:
                    continue
                dm = -2.5 * np.log10(max(f[k], 1e-9) / pred_ujy)
                rec = {"corridor": corridor, "observation_id": ci["observation_id"],
                       "band": obs["band"], "designation": s["designation"],
                       "wave_um": round(lam, 4), "pred_ujy": round(pred_ujy, 1),
                       "meas_ujy": round(float(f[k]), 1),
                       "snr": round(float(snr), 1), "dmag": round(float(dm), 3),
                       "var_scale": round(fm.var_scale, 3),
                       "pred_mag_ab": round(23.9 - 2.5 * np.log10(pred_ujy), 2)}
                meas_out.write(json.dumps(rec) + "\n")
                all_dm[obs["band"]].append(dm)
                n_meas += 1
        print(f"[{corridor}] {n_meas} star measurements from {len(sel)} cutouts",
              flush=True)
    stars_out.close()
    meas_out.close()

    summary = {"bands": {}, "selection": {
        "j_range": J_RANGE, "isolation_arcsec": ISOLATION_ARCSEC,
        "min_good_frac": MIN_GOOD_FRAC, "min_snr": 10,
        "sed_model": "log-log interpolation of Gaia RP, 2MASS JHKs, CatWISE W1W2 (Vega ZPs)"}}
    meds = []
    for band in sorted(all_dm):
        arr = np.array(all_dm[band])
        med = float(np.median(arr))
        mad = float(1.4826 * np.median(np.abs(arr - med)))
        summary["bands"][band] = {"n": int(len(arr)), "median_dmag": round(med, 3),
                                  "robust_scatter_mag": round(mad, 3),
                                  "mean_dmag_clipped": round(float(np.mean(
                                      arr[np.abs(arr - med) < 3 * max(mad, 0.05)])), 3)}
        meds.append(med)
    overall = float(np.median(meds)) if meds else 0.0
    summary["throughput_offset_mag"] = round(max(0.0, overall), 3)
    summary["median_dmag_all_bands"] = round(overall, 3)
    summary["note"] = ("dmag = -2.5 log10(measured/predicted); positive = "
                       "pipeline recovers less flux than the catalogue SED "
                       "predicts; applied to depths only when positive. "
                       "Catalogue SED interpolation is itself uncertain at "
                       "the ~0.1 mag level between bands.")
    (RUN_DIR / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
