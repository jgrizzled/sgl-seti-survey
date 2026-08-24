"""v2 engine: covariance envelopes and the cross-track decision
(hypotheses §1.4 of each survey; WISE v2_plan §2). Same Monte Carlo as
surveys/wise-v2/scripts/geometry_check.py with the profile's observer,
epochs and per-band nominal PSF (threshold = xt_threshold_fwhm x the
smallest band FWHM)."""

from __future__ import annotations

import json

import numpy as np
from astropy.time import Time

from sglseti import Role, evaluate_locus, load_target_registry, propagate_locus_uncertainty

_G: dict = {}


def _init(P):
    from sglseti import Observer
    _G["P"] = P
    ctx = P.geometry_context()
    # The envelope is a property of the target astrometry; the observer
    # changes the locus by <= 0.02" for terrestrial sites, so the Monte
    # Carlo uses the Earth centre (site transforms cost ~10x per draw).
    _G["ctx"] = type(ctx)(model=ctx.model, ephemeris=ctx.ephemeris, observer=Observer.earth_center(),
                          relay_range=ctx.relay_range, tolerance_arcsec=ctx.tolerance_arcsec,
                          padding_arcsec=ctx.padding_arcsec, confidence_level=ctx.confidence_level)
    _G["registry"] = load_target_registry(P.registry_path)


def envelope(endpoint, role):
    P, ctx, reg = _G["P"], _G["ctx"], _G["registry"]
    target = reg[endpoint]
    out = {"endpoint": endpoint, "role": role, "cells": []}
    for jy in P.mc_epochs_jyear:
        t = Time(jy, format="jyear")
        for z in P.mc_z_au:
            try:
                u = propagate_locus_uncertainty(target=target, role=Role(role), observation_time=t,
                                                observer=ctx.observer, z_au=z, ephemeris=ctx.ephemeris,
                                                model=ctx.model, seed=P.mc_seed, count=P.mc_samples,
                                                confidence_level=0.99)
            except Exception as exc:
                out["cells"].append({"jyear": jy, "z_au": z, "error": str(exc)}); continue
            off = np.array([o for o in u.offsets_arcsec if np.isfinite(o[0])])
            q = 1.0 / z
            pa = evaluate_locus(target=target, role=Role(role), observation_time=t, observer=ctx.observer,
                                z_au=1.0 / (q * 1.01), ephemeris=ctx.ephemeris, model=ctx.model)
            pb = evaluate_locus(target=target, role=Role(role), observation_time=t, observer=ctx.observer,
                                z_au=1.0 / (q * 0.99), ephemeris=ctx.ephemeris, model=ctx.model)
            cosd = np.cos(np.deg2rad(u.nominal.icrs_dec_deg))
            tvec = np.array([(pb.icrs_ra_deg - pa.icrs_ra_deg) * cosd, pb.icrs_dec_deg - pa.icrs_dec_deg])
            n = np.linalg.norm(tvec); tvec = tvec / n if n > 0 else np.array([1.0, 0.0])
            nvec = np.array([-tvec[1], tvec[0]])
            out["cells"].append({"jyear": jy, "z_au": z, "n_samples": int(len(off)),
                                 "sigma_xt_99": float(np.quantile(np.abs(off @ nvec), 0.99)),
                                 "sigma_at_99": float(np.quantile(np.abs(off @ tvec), 0.99)),
                                 "radius_99": float(u.confidence_radius_arcsec)})
    fin = [c for c in out["cells"] if "sigma_xt_99" in c]
    out["max_sigma_xt_99"] = max((c["sigma_xt_99"] for c in fin), default=None)
    out["max_radius_99"] = max((c["radius_99"] for c in fin), default=None)
    return out


def _run(job):
    try:
        return envelope(*job)
    except Exception as exc:
        return {"endpoint": job[0], "role": job[1], "error": str(exc), "cells": []}


def run(P, workers=2, endpoints=None):
    endpoints = endpoints or P.endpoints
    jobs = [(e, r) for e in endpoints for r in ("rx", "tx")]
    import multiprocessing as mp
    envs = []
    with mp.get_context("fork").Pool(workers, initializer=_init, initargs=(P,)) as pool:
        for k, env in enumerate(pool.imap_unordered(_run, jobs)):
            envs.append(env)
            if (k + 1) % 20 == 0:
                print(f"  {k + 1}/{len(jobs)}", flush=True)
    envs.sort(key=lambda e: (e["endpoint"], e["role"]))
    fwhm = min(P.psf_fwhm_nominal.values())
    thr = P.xt_threshold_fwhm * fwhm
    xt = []
    for e in envs:
        s = e.get("max_sigma_xt_99")
        if s is not None and s > thr:
            offsets = [0.0, -s, s] if s <= 2 * fwhm else [0.0, -s, -0.5 * s, 0.5 * s, s]
            xt.append({"endpoint": e["endpoint"], "role": e["role"], "sigma_xt_99": s, "offsets_arcsec": offsets})
    summary = {"mc_samples": P.mc_samples, "mc_seed": P.mc_seed, "epochs_jyear": list(P.mc_epochs_jyear),
               "z_au": list(P.mc_z_au), "threshold_arcsec": thr, "fwhm_arcsec": fwhm, "n_pairs": len(envs),
               "n_cross_track_cells": len(xt),
               "max_sigma_xt_99_overall": max((e.get("max_sigma_xt_99") or 0) for e in envs)}
    (P.run_dir / "geometry").mkdir(parents=True, exist_ok=True)
    (P.run_dir / "geometry" / "covariance_envelopes.json").write_text(
        json.dumps({"summary": summary, "envelopes": envs}, indent=1))
    P.config_dir.mkdir(parents=True, exist_ok=True)
    (P.config_dir / "cross_track_cells.json").write_text(json.dumps({"threshold_arcsec": thr, "cells": xt}, indent=1) + "\n")
    P.results_dir.mkdir(parents=True, exist_ok=True)
    lines = ["---", f'title: "{P.name} v2 — geometry: covariance envelopes"', "date: 2026-08-22", "---", "",
             f"sglseti seeded Monte Carlo, N = {P.mc_samples}, seed {P.mc_seed}, z = {P.mc_z_au} AU, epochs "
             f"{P.mc_epochs_jyear}; observer {P.observer_identity}. Cross-track threshold {P.xt_threshold_fwhm} x "
             f"FWHM ({fwhm}\") = {thr:.2f}\". **{len(xt)} of {len(envs)} endpoint-roles exceed it** "
             f"(max {summary['max_sigma_xt_99_overall']:.3f}\").", "",
             "| endpoint | role | max σ_xt,99 [\"] | max r_99 [\"] | cross-track |", "| --- | --- | --- | --- | --- |"]
    for e in envs:
        s = e.get("max_sigma_xt_99")
        lines.append(f"| {e['endpoint']} | {e['role']} | {'—' if s is None else f'{s:.3f}'} | "
                     f"{'—' if e.get('max_radius_99') is None else f'{e['max_radius_99']:.3f}'} | "
                     f"{'YES' if s is not None and s > thr else 'no'} |")
    (P.results_dir / "v2_geometry_summary.md").write_text("\n".join(lines) + "\n")
    print(json.dumps(summary, indent=1))
