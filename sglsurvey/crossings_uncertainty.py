"""Propagated impact-parameter uncertainty for a crossing-list product
(project plan §5.27 standing item; sglseti roadmap §3.3).

Runs :func:`sglseti.crossing_uncertainty` over every valid event of one
``crossings/<product>/`` list: each event's closest approach is
re-minimized for ``--count`` Monte Carlo draws of the target state
(positions, proper motions, parallax, radial velocity and — for the
orbit providers — the orbital elements, all from the schema-v2 registry
provenance uncertainties), inside a ``--window-days`` span centred on
the nominal ``t_ca``. The observer table, the ephemeris and the
geometry-model floor are declared ``not_propagated`` on every row, as
the library labels them; this product answers one question only — how
much the *target-state* knowledge moves each crossing's ``b_min``,
``t_ca``, ``v_perp`` and side-of-axis.

Outputs, next to the list's ``events.ecsv``:

* ``uncertainty.ecsv`` — one row per propagated event (quantile bounds
  at the declared confidence level, sigmas, side consistency, window-
  edge count, validity, warnings), joinable on ``event_id``;
* ``uncertainty_samples.npz`` — the per-event ``b_min`` sample vectors
  (``event_id``, ``b_min_samples_au[n_event, count]``), so the
  empirical distributions are reproducible without a re-run;
* ``uncertainty_summary.json`` — run parameters, input hashes
  (events.ecsv, registry, observer table), per-target maxima and the
  programme-level census (how many events could change rung
  membership at the 95 % bound).

Per-event seeds derive from the run seed and the ``event_id`` so the
product is independent of worker scheduling. Invalid status rows (no
nominal closest approach) are skipped and counted.

Usage::

    python -m sglsurvey.crossings_uncertainty --product universal_v1
    python -m sglsurvey.crossings_uncertainty --product universal_v1 \
        --targets sirius-a ez-aqr --count 32 --workers 4
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import sys
import time
from datetime import datetime, timezone
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from astropy.table import Table

from sglseti import Observer, load_target_registry
from sglseti.ephemeris import AstropyEphemeris
from sglseti.geometry import Tusay2022Eq57V1
from sglseti.models import (BeamSide, CrossingEvent, LinkDirection, Role,
                            UncertaintyMethod, Validity)
from sglseti.uncertainty import crossing_uncertainty

from sglsurvey.crossings import OBSERVER_DEFAULTS, REPO, resolve_observer
from sglsurvey.manifest import file_sha256

R_SUN_AU = 0.00465047
DEFAULT_SEED = 20260906
#: rung radii (AU) whose membership the census tests against the 95 %
#: bound: the grazing family (~1 R_sun beam at the 550 AU focus) and the
#: 0.1 / 1 AU sunward rungs used by the crossings surveys
CENSUS_RUNGS_AU = {"1.2_rsun": 1.2 * R_SUN_AU, "2.5_rsun": 2.5 * R_SUN_AU,
                   "0.1_au": 0.1, "1.0_au": 1.0}

_ENUM_FIELDS = {
    "link_direction": LinkDirection, "role": Role, "validity": Validity,
    "uncertainty_method": UncertaintyMethod,
}


def event_from_row(row) -> CrossingEvent:
    """Rebuild a :class:`CrossingEvent` from one ``events.ecsv`` row
    (windows are not carried; ``crossing_uncertainty`` does not read
    them)."""
    kw = {}
    for f in dataclasses.fields(CrossingEvent):
        if f.name in ("warnings", "windows"):
            continue
        v = row[f.name]
        if f.name in _ENUM_FIELDS:
            v = _ENUM_FIELDS[f.name](str(v))
        elif f.name == "side":
            v = BeamSide(str(v)) if str(v) else None
        elif f.name == "minimum_index":
            v = int(v)
        elif isinstance(v, (np.floating,)):
            v = float(v)
        elif isinstance(v, (np.str_,)):
            v = str(v)
        kw[f.name] = v
    w = str(row["warnings"])
    kw["warnings"] = tuple(w.split(";")) if w else ()
    return CrossingEvent(**kw)


def event_seed(run_seed: int, event_id: str) -> int:
    h = hashlib.sha256(f"{run_seed}:{event_id}".encode()).digest()
    return int.from_bytes(h[:4], "big")


# ---------------------------------------------------------------------------
# worker
# ---------------------------------------------------------------------------

_W: dict = {}


def _init_worker(registry_path: str, observer_name: str) -> None:
    _W["registry"] = load_target_registry(registry_path)
    _W["observer"], _ = resolve_observer(observer_name)
    _W["ephemeris"] = AstropyEphemeris()
    _W["model"] = Tusay2022Eq57V1()


def _run_one(args):
    row_dict, run_seed, count, confidence, window_days, tol = args
    ev = event_from_row(row_dict)
    t0 = time.time()
    try:
        u = crossing_uncertainty(
            event=ev, target=_W["registry"][ev.target_id],
            observer=_W["observer"], ephemeris=_W["ephemeris"],
            model=_W["model"], seed=event_seed(run_seed, ev.event_id),
            count=count, confidence_level=confidence,
            window_days=window_days, refine_tolerance_s=tol)
    except Exception as exc:  # surfaced as a failed row, never dropped
        return {"event_id": ev.event_id, "target_id": ev.target_id,
                "error": f"{type(exc).__name__}: {exc}",
                "elapsed_s": time.time() - t0}
    d = dataclasses.asdict(u)
    for k in ("link_direction", "method", "validity"):
        d[k] = d[k].value
    d["side_nominal"] = u.side_nominal.value if u.side_nominal else ""
    d["warnings"] = ";".join(u.warnings)
    d["contributions"] = ";".join(u.contributions)
    d["nominal_validity"] = ev.validity.value
    d["nominal_warnings"] = ";".join(ev.warnings)
    d["elapsed_s"] = time.time() - t0
    return d


# ---------------------------------------------------------------------------
# product assembly
# ---------------------------------------------------------------------------

ROW_COLUMNS = (
    "event_id", "target_id", "link_direction", "observer_id", "z_au",
    "method", "distribution", "sample_count", "seed", "confidence_level",
    "b_min_nominal_au", "b_min_lower_au", "b_min_median_au",
    "b_min_upper_au", "b_min_sigma_au",
    "b_min_sigma_solar_radii", "b_min_halfwidth95_solar_radii",
    "b_min_median_shift_au",
    "t_ca_nominal_tdb_jd", "t_ca_lower_tdb_jd", "t_ca_median_tdb_jd",
    "t_ca_upper_tdb_jd", "t_ca_sigma_days", "v_perp_sigma_km_s",
    "side_nominal", "side_consistency_fraction", "window_days",
    "refine_tolerance_s", "window_edge_count", "axis_model_id",
    "axis_model_version", "model_id", "model_version", "ephemeris_id",
    "contributions", "validity", "warnings",
    "nominal_validity", "nominal_warnings", "elapsed_s",
)


def build_table(results: list[dict]) -> Table:
    cols = {c: [] for c in ROW_COLUMNS}
    for d in results:
        d = dict(d)
        d["b_min_sigma_solar_radii"] = d["b_min_sigma_au"] / R_SUN_AU
        d["b_min_halfwidth95_solar_radii"] = (
            0.5 * (d["b_min_upper_au"] - d["b_min_lower_au"]) / R_SUN_AU)
        d["b_min_median_shift_au"] = d["b_min_median_au"] - d["b_min_nominal_au"]
        for c in ROW_COLUMNS:
            cols[c].append(d[c])
    tab = Table(cols)
    units = {"z_au": "AU", "b_min_nominal_au": "AU", "b_min_lower_au": "AU",
             "b_min_median_au": "AU", "b_min_upper_au": "AU",
             "b_min_sigma_au": "AU", "b_min_median_shift_au": "AU",
             "b_min_sigma_solar_radii": "solRad",
             "b_min_halfwidth95_solar_radii": "solRad",
             "t_ca_nominal_tdb_jd": "d", "t_ca_lower_tdb_jd": "d",
             "t_ca_median_tdb_jd": "d", "t_ca_upper_tdb_jd": "d",
             "t_ca_sigma_days": "d", "v_perp_sigma_km_s": "km / s",
             "window_days": "d", "refine_tolerance_s": "s", "elapsed_s": "s"}
    for c, u in units.items():
        tab[c].unit = u
    return tab


def census(tab: Table, events: Table) -> dict:
    """Programme-level summary: per-target maxima and rung-membership
    stability at the 95 % bound."""
    out: dict = {}
    nom_ok = tab["nominal_validity"] == "valid"
    per_target = {}
    for tid in sorted(set(tab["target_id"])):
        m = tab["target_id"] == tid
        mv = m & nom_ok
        per_target[tid] = {
            "n_events": int(m.sum()),
            "n_nominal_degraded": int((m & ~nom_ok).sum()),
            "max_b_sigma_solar_radii": float(tab["b_min_sigma_solar_radii"][m].max()),
            "max_b_halfwidth95_solar_radii": float(
                tab["b_min_halfwidth95_solar_radii"][m].max()),
            # shift census on nominally valid events only: for the list's
            # interval-boundary rows the refinement walks to the true
            # minimum outside the interval, a correction not an uncertainty
            "max_abs_median_shift_solar_radii": float(
                np.abs(tab["b_min_median_shift_au"][mv]).max() / R_SUN_AU)
            if mv.any() else float("nan"),
            "max_t_ca_sigma_s": float(tab["t_ca_sigma_days"][m].max() * 86400.0),
            "max_v_perp_sigma_km_s": float(tab["v_perp_sigma_km_s"][m].max()),
            "min_side_consistency": float(tab["side_consistency_fraction"][m].min()),
            "n_degraded": int((tab["validity"][m] == "degraded").sum()),
        }
    out["per_target"] = per_target
    rung = {}
    for name, r_au in CENSUS_RUNGS_AU.items():
        nominal_in = tab["b_min_nominal_au"] <= r_au
        lower_in = tab["b_min_lower_au"] <= r_au
        upper_in = tab["b_min_upper_au"] <= r_au
        rung[name] = {
            "radius_au": r_au,
            "nominal_inside": int(nominal_in.sum()),
            # membership could change within the 95 % interval
            "ambiguous_at_95": int((lower_in & ~upper_in).sum()),
            "nominal_inside_but_upper_outside": int((nominal_in & ~upper_in).sum()),
            "nominal_outside_but_lower_inside": int((~nominal_in & lower_in).sum()),
        }
    out["rung_membership"] = rung
    bnd = tab[~nom_ok]
    out["nominal_boundary_rows"] = {
        "n": len(bnd),
        "note": "events the list flagged degraded (interval-boundary "
                "minima); the propagated median is the refined minimum "
                "inside +/- window_days/2 and may lie outside the list's "
                "interval — reported separately from the uncertainty census",
        "max_abs_median_shift_solar_radii": float(
            np.abs(bnd["b_min_median_shift_au"]).max() / R_SUN_AU) if len(bnd) else 0.0,
        "max_abs_t_ca_shift_days": float(
            np.abs(bnd["t_ca_median_tdb_jd"] - bnd["t_ca_nominal_tdb_jd"]).max())
        if len(bnd) else 0.0,
        "n_still_at_window_edge": int((bnd["window_edge_count"] > 0).sum()),
    }
    valid = tab[nom_ok]
    out["global"] = {
        "n_rows": len(tab),
        "n_nominal_valid": len(valid),
        "max_abs_median_shift_solar_radii_nominal_valid": float(
            np.abs(valid["b_min_median_shift_au"]).max() / R_SUN_AU),
        "max_abs_t_ca_shift_s_nominal_valid": float(
            np.abs(valid["t_ca_median_tdb_jd"] - valid["t_ca_nominal_tdb_jd"]).max() * 86400.0),
        "n_degraded": int((tab["validity"] == "degraded").sum()),
        "n_side_flip_any": int((tab["side_consistency_fraction"] < 1.0).sum()),
        "max_b_sigma_solar_radii": float(tab["b_min_sigma_solar_radii"].max()),
        "median_b_sigma_solar_radii": float(np.median(tab["b_min_sigma_solar_radii"])),
        "max_b_halfwidth95_solar_radii": float(tab["b_min_halfwidth95_solar_radii"].max()),
        "max_t_ca_sigma_s": float(tab["t_ca_sigma_days"].max() * 86400.0),
        "median_t_ca_sigma_s": float(np.median(tab["t_ca_sigma_days"]) * 86400.0),
        "max_v_perp_sigma_km_s": float(tab["v_perp_sigma_km_s"].max()),
    }
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--product", default="universal_v1",
                    help="crossings/<product> directory (default universal_v1)")
    ap.add_argument("--observer", default=None,
                    help="observer name as in sglsurvey.crossings "
                         "(default: inferred from the product name)")
    ap.add_argument("--registry", type=Path, default=None,
                    help="default: the registry recorded in the product "
                         "manifest (its sha256 must still match)")
    ap.add_argument("--count", type=int, default=128)
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    ap.add_argument("--confidence", type=float, default=0.95)
    ap.add_argument("--window-days", type=float, default=90.0)
    ap.add_argument("--refine-tolerance-s", type=float, default=1.0)
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--targets", nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=None,
                    help="first N eligible events only (smoke tests)")
    ap.add_argument("--out", type=Path, default=None,
                    help="output directory (default: the product dir)")
    a = ap.parse_args(argv)

    product = REPO / "crossings" / a.product
    manifest = json.loads((product / "manifest.json").read_text())
    recorded = manifest["science_inputs"]["input_file_hashes"]
    reg_rel = next(k for k in recorded if k.startswith("registries/"))
    registry_path = a.registry or (REPO / reg_rel)
    reg_sha = "sha256:" + file_sha256(registry_path)
    if a.registry is None and reg_sha != recorded[reg_rel]:
        raise SystemExit(f"{reg_rel} no longer matches the product manifest "
                         f"({reg_sha[:23]}… vs {recorded[reg_rel][:23]}…); "
                         "rebuild the list or pass --registry explicitly")
    observer_name = a.observer
    if observer_name is None:
        observer_name = next(k for k, v in OBSERVER_DEFAULTS.items()
                             if v[2] == a.product or
                             a.product.startswith(v[2].split("_")[0] + "_"))
    observer, obs_notes = resolve_observer(observer_name)

    events = Table.read(product / "events.ecsv")
    ok = np.isfinite(events["t_ca_tdb_jd"]) & (events["validity"] != "invalid")
    n_invalid = int((~ok).sum())
    sel = events[ok]
    if a.targets:
        sel = sel[np.isin(sel["target_id"], a.targets)]
    if a.limit:
        sel = sel[: a.limit]
    print(f"{a.product}: {len(events)} events, {n_invalid} invalid skipped, "
          f"{len(sel)} to propagate; count={a.count} seed={a.seed} "
          f"window={a.window_days} d tol={a.refine_tolerance_s} s "
          f"observer={observer.observer_id}", file=sys.stderr)

    jobs = [({c: sel[c][i] for c in sel.colnames}, a.seed, a.count,
             a.confidence, a.window_days, a.refine_tolerance_s)
            for i in range(len(sel))]
    t0 = time.time()
    results, failures = [], []
    with Pool(a.workers, initializer=_init_worker,
              initargs=(str(registry_path), observer_name)) as pool:
        for k, d in enumerate(pool.imap_unordered(_run_one, jobs, chunksize=4)):
            (failures if "error" in d else results).append(d)
            if (k + 1) % 500 == 0 or k + 1 == len(jobs):
                el = time.time() - t0
                print(f"  {k + 1}/{len(jobs)} {el / 60:.1f} min "
                      f"(eta {el / (k + 1) * (len(jobs) - k - 1) / 60:.1f} min)",
                      file=sys.stderr)
    elapsed = time.time() - t0
    order = {eid: i for i, eid in enumerate(sel["event_id"])}
    results.sort(key=lambda d: order[d["event_id"]])

    out = a.out or product
    out.mkdir(parents=True, exist_ok=True)
    tab = build_table(results)
    generated = datetime.now(timezone.utc).isoformat(timespec="seconds")
    tab.meta = {
        "crossings_id": manifest["science_inputs"]["crossings_id"],
        "product": a.product, "generated_utc": generated,
        "uncertainty_schema_version": 1,
        "contributions_note": "target_state propagated; observer_state, "
                              "ephemeris and geometry_model_floor not "
                              "propagated (sglseti uncertainty v1.1)",
    }
    tab.write(out / "uncertainty.ecsv", format="ascii.ecsv", overwrite=True)
    np.savez_compressed(
        out / "uncertainty_samples.npz",
        event_id=np.array([d["event_id"] for d in results]),
        b_min_samples_au=np.array([d["b_min_samples_au"] for d in results]),
        seed=np.array([d["seed"] for d in results]))

    input_hashes = {
        str((product / "events.ecsv").relative_to(REPO)):
            "sha256:" + file_sha256(product / "events.ecsv"),
        reg_rel: reg_sha,
    }
    if "observer_table" in obs_notes:
        input_hashes[obs_notes["observer_table"]] = obs_notes["observer_table_sha256"]
    summary = {
        "product": a.product,
        "crossings_id": manifest["science_inputs"]["crossings_id"],
        "generated_utc": generated,
        "observer_id": observer.observer_id,
        "parameters": {"count": a.count, "seed": a.seed,
                       "confidence_level": a.confidence,
                       "window_days": a.window_days,
                       "refine_tolerance_s": a.refine_tolerance_s,
                       "targets": a.targets, "limit": a.limit},
        "input_file_hashes": input_hashes,
        "output_files": {
            "uncertainty_ecsv": "sha256:" + file_sha256(out / "uncertainty.ecsv"),
            "uncertainty_samples_npz":
                "sha256:" + file_sha256(out / "uncertainty_samples.npz")},
        "n_events_in_list": len(events), "n_invalid_skipped": n_invalid,
        "n_propagated": len(results), "n_failed": len(failures),
        "failures": failures[:50],
        "elapsed_s": elapsed, "workers": a.workers,
        "versions": manifest["run"]["versions"],
        "census": census(tab, events),
        **obs_notes,
    }
    (out / "uncertainty_summary.json").write_text(
        json.dumps(summary, indent=1) + "\n")
    g = summary["census"]["global"]
    print(json.dumps({k: summary[k] for k in ("n_propagated", "n_failed",
                                              "elapsed_s")}, indent=1))
    print(json.dumps(g, indent=1))
    print(json.dumps(summary["census"]["rung_membership"], indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
