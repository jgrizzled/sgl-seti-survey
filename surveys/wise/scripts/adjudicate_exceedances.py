"""Full-exceedance census and adjudication (stage 7b).

For every endpoint x role x band whose real-trajectory stack maximum
exceeds the 8-control threshold (runs/wise/calib_v1/threshold_report.json),
recompute the diagnostics that earlier batches applied by hand and
emit Candidate records with a verdict:

  phase split      S restricted to parallax phase 0 / phase 1 epochs at
                   the peak cell. A persistent relay contributes at both
                   alternating phases; a static star at one (screen_v1
                   lesson). Single-phase if the weaker phase carries
                   < PHASE_MIN of the stronger.
  epoch dominance  largest single-epoch share of the stack numerator at
                   the peak (v0.2.0 weight cap should hold this < ~0.1;
                   > EPOCH_MAX is a transient/artifact veto).
  colour           S in the companion band (W1<->W2) at the same cell;
                   S_W1/S_W2 >~ 1.5 with both positive is star-like.
  cryo visits      W3/W4 cells with < CRYO_MIN_EPOCHS epochs are
                   single-visit cryo cells: a coincident W1 static
                   source at the peak position vetoes them.
  static source    nearest CatWISE2020 source (IRSA TAP) to the peak
                   trajectory position at each phase's median epoch;
                   a source within STATIC_R of one phase position is
                   the static-star mechanism.
  bright halo      a CatWISE source brighter than BRIGHT_W1 within
                   BRIGHT_R of a phase position: the trajectory rides
                   the star's PSF wings (luyten-star/tx/W1 lesson —
                   S rises and decays over years as the fitted
                   residual motion slides across the halo; a
                   station-kept relay is persistent by hypothesis).
  marginal         excess < MARGINAL over a contamination-limited null
                   (control max S >> photon-noise S) with balanced
                   phases: batch-4 "dense-corridor" class; the catalog
                   screening null is the backstop (reported, not
                   asserted here).

Cells that fail no veto are written with status "retained" and MUST be
reviewed by hand; the script never invents a veto.

Usage:
  uv run python surveys/wise/scripts/adjudicate_exceedances.py \
      [--endpoints e1 e2 ...] [--run-id run-xxxx] [--no-tap] [--write]
Without --write, prints the table only.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import time
from pathlib import Path

import numpy as np
import requests
from astropy.time import Time
from sglseti import Role, adaptive_locus, load_target_registry, stable_id
from sglsurvey.geometry import GeometryContext
from sglsurvey.records import Candidate, append_records, read_records

REPO = Path(__file__).resolve().parents[3]
CAL_DIR = REPO / "runs" / "wise" / "calib_v1"
TENSOR_DIR = CAL_DIR / "tensors"
REGISTRY_PATH = REPO / "registries" / "pilot_wise_2026.yaml"
TAP = "https://irsa.ipac.caltech.edu/TAP/sync"

MIN_GOOD_FRAC = 0.7
WEIGHT_CAP_FACTOR = 20.0
BAND_NAME = {1: "W1", 2: "W2", 3: "W3", 4: "W4"}
PHASE_MIN = 0.25
EPOCH_MAX = 0.5
CRYO_MIN_EPOCHS = 60
STATIC_R = 6.0        # arcsec
BRIGHT_R = 15.0       # arcsec: PSF-wing / halo reach of a bright star
BRIGHT_W1 = 10.0      # CatWISE W1 brighter than this = halo hazard
MARGINAL = 0.05
LOCUS_BIN_DAYS = 0.5


def stack_S(f, v, g, sel):
    """S grid (64,5,5) over epochs `sel`, v0.2.0 weighting; also the
    per-epoch numerator contributions."""
    fb, vb, gb = f[0][sel], v[0][sel], g[0][sel]
    valid = np.isfinite(fb) & np.isfinite(vb) & (vb > 0) & (gb >= MIN_GOOD_FRAC)
    w = np.where(valid, 1.0 / np.where(vb > 0, vb, 1.0), 0.0)
    with np.errstate(all="ignore"):
        wmed = np.nanmedian(np.where(w > 0, w, np.nan), axis=0)
    cap = WEIGHT_CAP_FACTOR * np.nan_to_num(wmed, nan=np.inf, posinf=np.inf)
    w = np.minimum(w, cap[None])
    contrib = np.where(valid, fb, 0.0) * w
    A, B = contrib.sum(axis=0), w.sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        S = A / np.sqrt(np.where(B > 0, B, np.inf))
    return S, contrib, A


def peak_position(registry, ctx, endpoint, role, z, mu, t0, mjd):
    al = adaptive_locus(target=registry[endpoint], role=Role(role),
                        observation_time=Time(mjd, format="mjd"),
                        observer=ctx.observer, relay_range=ctx.relay_range,
                        tolerance_arcsec=2.0, ephemeris=ctx.ephemeris,
                        model=ctx.model)
    zs = np.array([p.z_au for p in al.points])
    ra = np.array([p.icrs_ra_deg for p in al.points])
    dec = np.array([p.icrs_dec_deg for p in al.points])
    o = np.argsort(1 / zs)
    ra0 = np.interp(1 / z, (1 / zs)[o], ra[o])
    dec0 = np.interp(1 / z, (1 / zs)[o], dec[o])
    dt = (mjd - t0) / 365.25
    return (ra0 + mu[0] * dt / 3600.0 / np.cos(np.deg2rad(dec0)),
            dec0 + mu[1] * dt / 3600.0)


def catwise_nearest(session, ra, dec, r_arcsec=15.0):
    """(sep, w1, w2) of the nearest source, plus the brightest source
    within the cone as a 4th element (sep, w1)."""
    q = ("SELECT ra, dec, w1mpro, w2mpro FROM catwise_2020 WHERE "
         f"CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.6f},{dec:.6f},"
         f"{r_arcsec / 3600:.6f}))=1")
    for k in range(4):
        try:
            r = session.get(TAP, params={"QUERY": q, "FORMAT": "CSV"}, timeout=120)
            if r.status_code == 200:
                rows = list(csv.DictReader(io.StringIO(r.text)))
                break
        except requests.RequestException:
            pass
        time.sleep(10 * (k + 1))
    else:
        return None
    best, bright = None, None
    for row in rows:
        d = np.hypot((float(row["ra"]) - ra) * np.cos(np.deg2rad(dec)),
                     float(row["dec"]) - dec) * 3600
        w1 = float(row["w1mpro"] or "nan")
        if best is None or d < best[0]:
            best = (d, w1, float(row["w2mpro"] or "nan"))
        if np.isfinite(w1) and (bright is None or w1 < bright[1]):
            bright = (d, w1)
    return best + (bright,) if best else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--endpoints", nargs="*", default=None)
    ap.add_argument("--run-id", default=None)
    ap.add_argument("--no-tap", action="store_true")
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()

    rep = json.load(open(CAL_DIR / "threshold_report.json"))
    runs = read_records(CAL_DIR / "records" / "analysis_run.jsonl")
    run_id = a.run_id or runs[-1]["analysis_run_id"]
    registry = load_target_registry(REGISTRY_PATH)
    ctx = GeometryContext.wise_coarse_default()
    session = requests.Session()

    keys = [k for k, r in rep.items() if r["exceeds"]
            and (a.endpoints is None or k.split("/")[0] in a.endpoints)]
    print(f"{len(keys)} exceeding cells; run {run_id}\n")
    out, cands = [], []
    for key in sorted(keys, key=lambda k: -rep[k]["real_max_S"] / rep[k]["T"]):
        endpoint, role, band = key.split("/")
        b = {v: k for k, v in BAND_NAME.items()}[band]
        d = np.load(TENSOR_DIR / f"{endpoint}__{role}.npz")
        z_grid, mu_grid, t0 = d["z_grid"], d["mu_grid"], float(d["t0_mjd"])
        mjd, band_idx, phase = d["mjd"], d["band_idx"], d["phase"]
        f, v, g = d["f"], np.asarray(d["v"]), np.asarray(d["g"], dtype=np.float32)
        eb = band_idx == b
        S, contrib, A = stack_S(f, v, g, eb)
        pk = np.unravel_index(np.nanargmax(S), S.shape)
        S_real = float(S[pk])
        ph = phase[eb]
        S0 = float(stack_S(f, v, g, np.where(eb)[0][ph == 0])[0][pk]) if (ph == 0).any() else 0.0
        S1 = float(stack_S(f, v, g, np.where(eb)[0][ph == 1])[0][pk]) if (ph == 1).any() else 0.0
        c = contrib[(slice(None),) + pk]
        top_frac = float(np.max(c) / A[pk]) if A[pk] > 0 else 0.0
        n_ep = int(eb.sum())
        # companion band
        other = {1: 2, 2: 1, 3: 1, 4: 1}[b]
        eo = band_idx == other
        S_other = float(stack_S(f, v, g, eo)[0][pk]) if eo.any() else float("nan")
        T = rep[key]["T"]
        z, mu = float(z_grid[pk[0]]), (float(mu_grid[pk[1]]), float(mu_grid[pk[2]]))
        # static-source check at each phase's median epoch
        static = []
        if not a.no_tap:
            for p in (0, 1):
                sel = np.where(eb)[0][ph == p]
                if len(sel) == 0:
                    static.append(None)
                    continue
                m = float(np.median(mjd[sel]))
                ra, dec = peak_position(registry, ctx, endpoint, role, z, mu, t0, m)
                static.append(catwise_nearest(session, ra, dec))
        d.close()

        # ---- verdict ----------------------------------------------
        reasons = []
        lo, hi = sorted([abs(S0), abs(S1)])
        single_phase = (hi > 0 and lo / hi < PHASE_MIN) or (S0 <= 0) != (S1 <= 0)
        if single_phase:
            reasons.append(f"single-phase (S0={S0:.1f}, S1={S1:.1f})")
        if top_frac > EPOCH_MAX:
            reasons.append(f"single-epoch dominance ({top_frac:.0%} of stack)")
        if b in (1, 2) and np.isfinite(S_other) and S_other > 0:
            ratio = (S_real / S_other) if b == 1 else (S_other / S_real)
            if ratio > 1.5:
                reasons.append(f"W1:W2 significance ratio {ratio:.1f} star-like")
        near = [s for s in static if s and s[0] < STATIC_R]
        halo = [s[3] for s in static if s and s[3] and s[3][1] < BRIGHT_W1
                and s[3][0] < BRIGHT_R]
        if halo:
            reasons.append("bright-star halo: W1=%.1f star %.1f\" from a phase "
                           "position (trajectory rides its PSF wings)"
                           % (halo[0][1], halo[0][0]))
        if b in (3, 4) and n_ep < CRYO_MIN_EPOCHS:
            reasons.append(f"single-visit cryo cell ({n_ep} epochs)"
                           + (f"; CatWISE source {near[0][0]:.1f}\" W1={near[0][1]:.1f}" if near else ""))
        elif near:
            desc = "; ".join(f"phase {i}: {s[0]:.1f}\" W1={s[1]:.1f}"
                             for i, s in enumerate(static) if s and s[0] < STATIC_R)
            reasons.append(f"static CatWISE source within {STATIC_R:.0f}\" ({desc})"
                           + (" — distinct sources at both phases (dense-corridor class)"
                              if len(near) == 2 else ""))
        excess = S_real / T - 1
        if not reasons and excess < MARGINAL:
            reasons.append(f"+{excess:.0%} over contamination-limited null "
                           f"(T={T:.0f}); balanced phases from dense-field "
                           "static sources; screening null is the backstop")
        status = "vetoed" if reasons else "retained"
        row = dict(key=key, S=S_real, T=T, z=z, mu=mu, S0=S0, S1=S1,
                   top_frac=top_frac, n_ep=n_ep, S_other=S_other,
                   static=static, status=status, reason="; ".join(reasons))
        out.append(row)
        print(f"{key:22s} S={S_real:7.1f} T={T:7.1f} z={z:6.0f} "
              f"ph=({S0:6.1f},{S1:6.1f}) top={top_frac:4.0%} n={n_ep:4d} "
              f"Sb={S_other:6.1f} {status:8s} {row['reason']}")
        cands.append(Candidate(
            candidate_id=stable_id("cnd", {"endpoint": endpoint, "role": role,
                                           "band": band, "run": run_id,
                                           "source": "exceedance-census"}),
            analysis_run_id=run_id, endpoint_id=endpoint, role=role,
            observation_ids=(), fitted_z_au=z,
            fitted_residual_motion={"mu_ra": mu[0], "mu_dec": mu[1]},
            model_comparison={"real_max_S": round(S_real, 1),
                              "threshold_8_controls": round(T, 1),
                              "phase_split_S": [round(S0, 1), round(S1, 1)],
                              "top_epoch_fraction": round(top_frac, 3),
                              "companion_band_S": (round(S_other, 1)
                                                   if np.isfinite(S_other) else None),
                              "catwise_nearest": static},
            status=status, veto_reason=row["reason"] or None,
            extra={"band": band, "origin": "full-calibration exceedance census",
                   "n_epochs": n_ep}))
    n_ret = sum(1 for r in out if r["status"] == "retained")
    print(f"\n{len(out)} exceedances: {len(out) - n_ret} vetoed, {n_ret} RETAINED for review")
    if a.write:
        append_records(CAL_DIR / "records" / "candidate.jsonl", cands)
        (CAL_DIR / f"exceedance_census_{run_id}.json").write_text(
            json.dumps(out, indent=1, default=str))
        print("candidate records appended")


if __name__ == "__main__":
    main()
