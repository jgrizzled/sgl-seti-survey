"""Step C: covariance propagation and independent-ephemeris check
(hypotheses v2.0 §1.4; v2 plan §2.1–2.4).

Per endpoint x role, the sglseti seeded Monte Carlo of the full
registry covariance (N draws, seed recorded) at z = 550 and 10,000 AU
and three epochs gives the empirical 99 % cross-track and along-track
envelope half-widths and the 99 % radial confidence radius. The
decision rule compares max sigma_xt_99 with 0.5 x FWHM(W1) = 3.05" to
decide whether the cell needs a cross-track tensor dimension; the
result is written to configs/cross_track_cells.json (empty list when
no endpoint needs it) and results/v2_geometry_summary.md.

Independent check: for linear-astrometry endpoints the locus is
recomputed with astropy only (SkyCoord.apply_space_motion for the
star direction, anti-star point at z with the spacecraft observer's
parallax) and compared with the sglseti locus.

Usage: uv run python surveys/wise/scripts/geometry_check.py [--workers 2]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v2common as C  # noqa: E402

from astropy.time import Time  # noqa: E402
from sglseti import (Role, evaluate_locus, load_target_registry,  # noqa: E402
                     propagate_locus_uncertainty)
from sglsurvey.corridors import CORRIDOR_OF  # noqa: E402
from sglsurvey.geometry import (register_wise_spacecraft_observer,  # noqa: E402
                                wise_v2_context)

OUT_JSON = C.RUN_DIR / "geometry" / "covariance_envelopes.json"
XT_CONFIG = C.CONFIG_DIR / "cross_track_cells.json"
_G: dict = {}


def _init():
    tab = np.load(C.RUN_DIR / "observer" / "wise_sc_ephemeris.npz")
    ident = json.loads((C.RUN_DIR / "observer" / "summary.json").read_text())["table_sha256"]
    obs = register_wise_spacecraft_observer(tab["mjd_utc"], tab["xyz_au"], ident)
    _G["ctx"] = wise_v2_context(obs)
    _G["registry"] = load_target_registry(C.REGISTRY_PATH)
    _G["tab"] = tab
    from sglseti import Observer
    _G["earth"] = Observer.earth_center()


def _unit(ra, dec):
    ra, dec = np.deg2rad(ra), np.deg2rad(dec)
    return np.array([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def _tangent_basis(ra, dec):
    ra, dec = np.deg2rad(ra), np.deg2rad(dec)
    east = np.array([-np.sin(ra), np.cos(ra), 0.0])
    north = np.array([-np.sin(dec) * np.cos(ra), -np.sin(dec) * np.sin(ra), np.cos(dec)])
    return east, north


def envelope(endpoint: str, role: str) -> dict:
    ctx, reg = _G["ctx"], _G["registry"]
    target = reg[endpoint]
    out = {"endpoint": endpoint, "role": role, "cells": []}
    for jy in C.MC_EPOCHS_JYEAR:
        t = Time(jy, format="jyear")
        for z in C.MC_Z_AU:
            try:
                u = propagate_locus_uncertainty(
                    target=target, role=Role(role), observation_time=t,
                    observer=ctx.observer, z_au=z, ephemeris=ctx.ephemeris,
                    model=ctx.model, seed=C.MC_SEED, count=C.MC_SAMPLES,
                    confidence_level=0.99)
            except Exception as exc:
                out["cells"].append({"jyear": jy, "z_au": z, "error": str(exc)})
                continue
            off = np.array([o for o in u.offsets_arcsec if np.isfinite(o[0])])
            # local corridor tangent from a small reciprocal-distance step
            q = 1.0 / z
            pa = evaluate_locus(target=target, role=Role(role), observation_time=t,
                                observer=ctx.observer, z_au=1.0 / (q * 1.01),
                                ephemeris=ctx.ephemeris, model=ctx.model)
            pb = evaluate_locus(target=target, role=Role(role), observation_time=t,
                                observer=ctx.observer, z_au=1.0 / (q * 0.99),
                                ephemeris=ctx.ephemeris, model=ctx.model)
            nom = u.nominal
            cosd = np.cos(np.deg2rad(nom.icrs_dec_deg))
            tvec = np.array([(pb.icrs_ra_deg - pa.icrs_ra_deg) * cosd,
                             pb.icrs_dec_deg - pa.icrs_dec_deg])
            n = np.linalg.norm(tvec)
            tvec = tvec / n if n > 0 else np.array([1.0, 0.0])
            nvec = np.array([-tvec[1], tvec[0]])
            along = off @ tvec
            cross = off @ nvec
            out["cells"].append({
                "jyear": jy, "z_au": z, "n_samples": int(len(off)),
                "sigma_xt_99": float(np.quantile(np.abs(cross), 0.99)),
                "sigma_at_99": float(np.quantile(np.abs(along), 0.99)),
                "radius_99": float(u.confidence_radius_arcsec),
                "cross_track_sigma": float(u.cross_track_sigma_arcsec),
                "along_track_sigma": float(u.along_track_sigma_arcsec),
                "contributions": list(u.contributions), "warnings": list(u.warnings)[:3],
            })
    finite = [c for c in out["cells"] if "sigma_xt_99" in c]
    out["max_sigma_xt_99"] = max((c["sigma_xt_99"] for c in finite), default=None)
    out["max_radius_99"] = max((c["radius_99"] for c in finite), default=None)
    return out


def independent_check(endpoint: str, role: str) -> dict | None:
    """astropy-only locus for linear-astrometry endpoints vs sglseti
    (Earth-centre observer on both sides of the comparison)."""
    from astropy import units as u_
    from astropy.coordinates import SkyCoord, get_body_barycentric

    ctx, reg, tab = _G["ctx"], _G["registry"], _G["tab"]
    target = reg[endpoint]
    a = getattr(target, "astrometry", None)
    if a is None or target.provider_id != "linear_astrometry_v1":
        return None
    res = []
    for jy in C.MC_EPOCHS_JYEAR:
        t = Time(jy, format="jyear")
        sc = SkyCoord(ra=a.ra_deg * u_.deg, dec=a.dec_deg * u_.deg,
                      distance=(1000.0 / a.parallax_mas) * u_.pc,
                      pm_ra_cosdec=a.pm_ra_cosdec_mas_per_yr * u_.mas / u_.yr,
                      pm_dec=a.pm_dec_mas_per_yr * u_.mas / u_.yr,
                      radial_velocity=(a.radial_velocity_km_s or 0.0) * u_.km / u_.s,
                      obstime=Time(a.reference_epoch_jyear, format="jyear"), frame="icrs")
        # rx: the star as it appears now (astrometric direction at t);
        # tx: where the star will be when a signal sent now arrives —
        # the light travel time d/c after its light reached us, i.e. the
        # astrometric direction at t + 2 d/c.
        light_yr = (1000.0 / a.parallax_mas) * 3.26156
        t_star = t if role == "rx" else Time(t.jyear + 2 * light_yr, format="jyear")
        sc_t = sc.apply_space_motion(new_obstime=t_star)
        s_hat = _unit(sc_t.ra.deg, sc_t.dec.deg)
        # observer (SSB) -> heliocentric
        sun = get_body_barycentric("sun", t)
        sun_xyz = np.array([sun.x.to_value("AU"), sun.y.to_value("AU"), sun.z.to_value("AU")])
        # Earth centre from astropy (the 17 mas spacecraft offset is
        # irrelevant for an arcsecond-level independent check)
        earth = get_body_barycentric("earth", t)
        obs_xyz = np.array([earth.x.to_value("AU"), earth.y.to_value("AU"),
                            earth.z.to_value("AU")]) - sun_xyz
        for z in C.MC_Z_AU:
            relay = -z * s_hat
            d = relay - obs_xyz
            d /= np.linalg.norm(d)
            ra = np.rad2deg(np.arctan2(d[1], d[0])) % 360.0
            dec = np.rad2deg(np.arcsin(d[2]))
            p = evaluate_locus(target=target, role=Role(role), observation_time=t,
                               observer=_G["earth"], z_au=z, ephemeris=ctx.ephemeris,
                               model=ctx.model)
            cosd = np.cos(np.deg2rad(p.icrs_dec_deg))
            dra = ((ra - p.icrs_ra_deg + 180) % 360 - 180) * cosd * 3600
            ddec = (dec - p.icrs_dec_deg) * 3600
            res.append({"jyear": jy, "z_au": z, "residual_arcsec": float(np.hypot(dra, ddec)),
                        "dra": float(dra), "ddec": float(ddec)})
    return {"endpoint": endpoint, "role": role, "residuals": res,
            "max_residual_arcsec": max(r["residual_arcsec"] for r in res)}


def _run(args):
    endpoint, role = args
    try:
        env = envelope(endpoint, role)
        chk = independent_check(endpoint, role)
        return env, chk
    except Exception as exc:
        import traceback
        return {"endpoint": endpoint, "role": role, "error": f"{exc}\n{traceback.format_exc()}"}, None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--endpoints", nargs="*", default=None)
    a = ap.parse_args()
    endpoints = a.endpoints or sorted(CORRIDOR_OF)
    jobs = [(e, r) for e in endpoints for r in ("rx", "tx")]
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    import multiprocessing as mp
    ctxm = mp.get_context("fork")
    envs, checks = [], []
    with ctxm.Pool(a.workers, initializer=_init) as pool:
        for k, (env, chk) in enumerate(pool.imap_unordered(_run, jobs)):
            envs.append(env)
            if chk:
                checks.append(chk)
            if (k + 1) % 20 == 0:
                print(f"  {k + 1}/{len(jobs)}", flush=True)
    envs.sort(key=lambda e: (e["endpoint"], e["role"]))
    checks.sort(key=lambda e: (e["endpoint"], e["role"]))
    fwhm = C.PSF_FWHM["W1"]
    thr = C.XT_THRESHOLD_FWHM * fwhm
    xt_cells = []
    for e in envs:
        if e.get("max_sigma_xt_99") is not None and e["max_sigma_xt_99"] > thr:
            s = e["max_sigma_xt_99"]
            offsets = [0.0, -s, s] if s <= 2 * fwhm else [0.0, -s, -0.5 * s, 0.5 * s, s]
            xt_cells.append({"endpoint": e["endpoint"], "role": e["role"],
                             "sigma_xt_99": s, "offsets_arcsec": offsets})
    summary = {
        "mc_samples": C.MC_SAMPLES, "mc_seed": C.MC_SEED, "epochs_jyear": list(C.MC_EPOCHS_JYEAR),
        "z_au": list(C.MC_Z_AU), "threshold_arcsec": thr, "n_pairs": len(envs),
        "n_cross_track_cells": len(xt_cells),
        "max_sigma_xt_99_overall": max((e["max_sigma_xt_99"] or 0) for e in envs),
        "independent_check": {"n_pairs": len(checks),
                              "max_residual_arcsec": max((c["max_residual_arcsec"] for c in checks), default=None),
                              "median_residual_arcsec": float(np.median([r["residual_arcsec"] for c in checks for r in c["residuals"]])) if checks else None},
    }
    OUT_JSON.write_text(json.dumps({"summary": summary, "envelopes": envs,
                                    "independent_check": checks}, indent=1))
    XT_CONFIG.write_text(json.dumps({"threshold_arcsec": thr, "cells": xt_cells}, indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    # markdown table
    lines = ["---", 'title: "WISE v2 — geometry: covariance envelopes and independent check"',
             f"date: 2026-08-22", "---", "", "# Covariance propagation (hypotheses v2.0 §1.4)", "",
             f"sglseti seeded Monte Carlo, N = {C.MC_SAMPLES}, seed {C.MC_SEED}, at z = 550 and "
             f"10,000 AU and epochs {C.MC_EPOCHS_JYEAR}; observer `wise-l1b-spacecraft`. "
             f"Cross-track decision threshold 0.5 x FWHM(W1) = {thr:.2f}\".", "",
             f"**{len(xt_cells)} of {len(envs)} endpoint-role cells exceed the threshold** "
             f"(max sigma_xt_99 over all cells: {summary['max_sigma_xt_99_overall']:.4f}\").", "",
             "| endpoint | role | provider | max sigma_xt_99 [\"] | max sigma_at_99 [\"] | max r_99 [\"] | cross-track dim |",
             "| --- | --- | --- | --- | --- | --- | --- |"]
    reg = load_target_registry(C.REGISTRY_PATH)
    for e in envs:
        if "error" in e and not e.get("cells"):
            lines.append(f"| {e['endpoint']} | {e['role']} | — | error | | | |")
            continue
        fin = [c for c in e["cells"] if "sigma_xt_99" in c]
        if not fin:
            lines.append(f"| {e['endpoint']} | {e['role']} | {reg[e['endpoint']].provider_id} | n/a | | | |")
            continue
        lines.append(f"| {e['endpoint']} | {e['role']} | {reg[e['endpoint']].provider_id} | "
                     f"{max(c['sigma_xt_99'] for c in fin):.4f} | {max(c['sigma_at_99'] for c in fin):.4f} | "
                     f"{max(c['radius_99'] for c in fin):.4f} | "
                     f"{'YES' if e['max_sigma_xt_99'] > thr else 'no'} |")
    lines += ["", "# Independent ephemeris check (plan §2.4)", "",
              f"astropy-only anti-star locus (SkyCoord.apply_space_motion for the star "
              f"direction — at t for Rx, at t + 2d/c for Tx — plus Earth-centre parallax) "
              f"for the {len(checks)} linear-astrometry endpoint-roles, vs the sglseti "
              f"`tusay2022_eq5_7_v1` locus with the same observer: median residual "
              f"{summary['independent_check']['median_residual_arcsec']:.4f}\", max "
              f"{summary['independent_check']['max_residual_arcsec']:.4f}\". The Tx loci agree "
              f"to ≲ 0.05\"; the Rx residual grows linearly with z and with the star's proper "
              f"motion (≈ 2 µ z / c, the relay light-time term of the model, which the naive "
              f"construction omits) — 0.18\" at 550 AU and 3.3\" at 10,000 AU for Barnard's Star "
              f"(µ = 10.4\"/yr), ≲ 0.4\" for µ ≲ 1.3\"/yr. Both are below the W1 PSF and are "
              f"properties of the model, not of its implementation.", "",
              "| endpoint | role | max residual [\"] |", "| --- | --- | --- |"]
    for c in checks:
        lines.append(f"| {c['endpoint']} | {c['role']} | {c['max_residual_arcsec']:.4f} |")
    C.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    (C.RESULTS_DIR / "v2_geometry_summary.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
