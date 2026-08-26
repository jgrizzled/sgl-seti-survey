"""Universal historical beam-crossing list (Pipeline B, project plan §3.5).

Computes, once per observer, every local minimum of the impact parameter
b(t) between the observer and each registry endpoint's Sun-star beam axis
(``sglseti.crossings``, ``sun_star_axis_v1``) over a generous historical
window, for both link directions. The product is survey-independent:
each ``surveys/<name>-crossings/`` sub-project intersects it with that
archive's actual exposure coverage and applies beam radius / wavelength /
duty-cycle hypotheses afterwards. Accordingly no beam radii and no
``report_max_b_au`` cut are applied here — every minimum is kept, and
ephemeris-coverage failures are kept as ``invalid`` rows rather than
dropped.

Earth-center is the canonical observer (ground surveys, LEO). Spacecraft
observers are derivative runs keyed by ``observer_id``, produced from the
same request with a per-epoch observer (plan §11.2 step 4):

* ``--observer wise`` — the WISE v2 spacecraft table
  (``runs/wise/v2/observer/wise_sc_ephemeris.npz``, built from L1b frame
  headers; exact at frame epochs, Earth-center fallback in gaps > 2 d
  incl. the 2011–2013 hibernation; LEO bound ≤ ~7,100 km ≈ 0.010 R☉).
* ``--observer tess`` — a JPL Horizons SSB vector table (spacecraft
  -95) fetched by ``--fetch-observer tess`` into
  ``crossings/observers/`` (raw response + sha pinned). TESS's HEO
  (apogee ~0.54 R☉) is the case where the Earth-center list is invalid
  at grazing impact parameters.
* ``--observer spherex`` — Earth center over the SPHEREx era with the
  LEO offset declared as a budget (a coarse-sampled table of a 95-min
  orbit would be no more accurate than Earth center; the per-frame
  route via L2 headers is a survey-stage task).
* ``--observer soho`` — a JPL Horizons SSB vector table (spacecraft
  -21) fetched by ``--fetch-observer soho``. SOHO's L1 halo orbit has
  ~0.9 R_sun transverse amplitude (measured directly in the LASCO
  recon star check), so the Earth-center list is invalid at grazing
  impact parameters, as for TESS. Era stop 2026-10-01 is limited by
  the Horizons SPK end (2026-10-05); extend at the yearly refresh.

Usage::

    python -m sglsurvey.crossings --observer earth \
        --start 1980-01-01 --stop 2028-01-01 --out crossings/universal_v1
    python -m sglsurvey.crossings --fetch-observer tess
    python -m sglsurvey.crossings --observer tess --out crossings/tess_v1
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

from astropy.time import Time

from sglseti import (CrossingsRequest, LinkDirection, Observer, TimeInterval,
                     find_crossings, load_target_registry, write_crossings_products)

REPO = Path(__file__).resolve().parents[1]
DEFAULT_REGISTRY = REPO / "registries" / "pilot_wise_2026.yaml"
MODEL_ID = "tusay2022_eq5_7_v1"
RELAY_DISTANCE_AU = 550.0   # representative; the axis is independent of it
OBSERVER_TABLE_DIR = REPO / "crossings" / "observers"
WISE_TABLE = REPO / "runs" / "wise" / "v2" / "observer" / "wise_sc_ephemeris.npz"

#: per-observer default eras and product dirs
OBSERVER_DEFAULTS = {
    "earth": ("1980-01-01", "2028-01-01", "universal_v1"),
    # L1b frames span MJD 55203-60523 (2010-01-07 .. 2024-08-01)
    "wise": ("2010-01-01", "2024-08-01", "wise_v1"),
    # science operations after commissioning; Horizons TESS_merged SPK
    "tess": ("2018-07-04", "2026-08-24", "tess_v1"),
    # launch 2025-03-11; era extends to the universal window end so the
    # list survives quick-release growth
    "spherex": ("2025-03-15", "2028-01-01", "spherex_v1"),
    # LASCO science era; stop bounded by the Horizons -21 SPK end
    # (2026-10-05), not the universal 2028 window — extend at refresh
    "soho": ("1996-01-01", "2026-10-01", "soho_v1"),
}
HORIZONS_IDS = {"tess": "-95", "soho": "-21"}


def _sha256(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def fetch_horizons_table(observer: str, start: str, stop: str,
                         step: str = "6h") -> Path:
    """Fetch an SSB ICRF position table for a spacecraft from JPL
    Horizons into ``crossings/observers/`` (raw response retained, npz
    keyed by MJD UTC). Returns the npz path."""
    import numpy as np
    import requests

    OBSERVER_TABLE_DIR.mkdir(parents=True, exist_ok=True)
    resp = requests.get(
        "https://ssd.jpl.nasa.gov/api/horizons.api",
        params={"format": "text", "COMMAND": f"'{HORIZONS_IDS[observer]}'",
                "EPHEM_TYPE": "VECTORS", "CENTER": "'500@0'",
                "START_TIME": f"'{start}'", "STOP_TIME": f"'{stop}'",
                "STEP_SIZE": f"'{step}'", "REF_PLANE": "'FRAME'",
                "OUT_UNITS": "'AU-D'", "VEC_TABLE": "'1'",
                "CSV_FORMAT": "'YES'"}, timeout=600)
    resp.raise_for_status()
    txt = resp.text
    if "$$SOE" not in txt:
        raise RuntimeError("Horizons returned no vector block:\n"
                           + txt[:1000])
    raw = OBSERVER_TABLE_DIR / f"{observer}_horizons_response.txt"
    raw.write_text(txt)
    rows = []
    for line in txt.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines():
        parts = [p.strip() for p in line.split(",")]
        if len(parts) < 5:
            continue
        rows.append((float(parts[0]), float(parts[2]), float(parts[3]),
                     float(parts[4])))
    arr = np.array(rows)
    t = Time(arr[:, 0], format="jd", scale="tdb")
    mjd_utc = t.utc.mjd
    npz = OBSERVER_TABLE_DIR / f"{observer}_sc_ephemeris.npz"
    np.savez_compressed(npz, mjd_utc=mjd_utc, xyz_au=arr[:, 1:4])
    meta = {"observer": observer, "horizons_command": HORIZONS_IDS[observer],
            "center": "500@0 (SSB)", "ref_plane": "FRAME (ICRF)",
            "start": start, "stop": stop, "step": step,
            "n_rows": int(len(rows)),
            "mjd_utc_range": [float(mjd_utc.min()), float(mjd_utc.max())],
            "raw_response": raw.name, "raw_sha256": _sha256(raw),
            "table_sha256": _sha256(npz)}
    (OBSERVER_TABLE_DIR / f"{observer}_table_meta.json").write_text(
        json.dumps(meta, indent=1) + "\n")
    print(json.dumps(meta, indent=1))
    return npz


def resolve_observer(name: str):
    """Returns (Observer, note_dict) for a --observer name."""
    import numpy as np

    from sglsurvey.geometry import register_spacecraft_table_observer

    if name == "earth":
        return Observer.earth_center(), {}
    if name == "spherex":
        return Observer.earth_center(), {
            "observer_applicability": "SPHEREx (LEO ~650 km): Earth-center "
            "geometry with the spacecraft offset carried as a budget — "
            "|Δb| <= ~7,100 km = 0.010 R_sun = 4.7e-5 AU, |Δt_ca| <= ~4 min "
            "at v_perp = 30 km/s. Matches the SPHEREx pipeline-A "
            "observer convention."}
    if name == "wise":
        tab = np.load(WISE_TABLE)
        ident = _sha256(WISE_TABLE)
        obs = register_spacecraft_table_observer(
            "wise-l1b-spacecraft", tab["mjd_utc"], tab["xyz_au"], ident)
        return obs, {
            "observer_table": str(WISE_TABLE.relative_to(REPO)),
            "observer_table_sha256": ident,
            "observer_accuracy": "exact at L1b frame epochs; linear "
            "interpolation between frames <= 2 d apart; Earth-center "
            "fallback in larger gaps (incl. the 2011-02..2013-12 "
            "hibernation) — bounded everywhere by the geocentric radius "
            "~7,000 km = 0.010 R_sun"}
    if name == "tess":
        npz = OBSERVER_TABLE_DIR / "tess_sc_ephemeris.npz"
        if not npz.exists():
            raise SystemExit("no TESS table; run --fetch-observer tess first")
        tab = np.load(npz)
        ident = _sha256(npz)
        obs = register_spacecraft_table_observer(
            "tess-spacecraft", tab["mjd_utc"], tab["xyz_au"], ident,
            max_gap_days=1.0)
        return obs, {
            "observer_table": str(npz.relative_to(REPO)),
            "observer_table_sha256": ident,
            "observer_accuracy": "JPL Horizons -95 SSB vectors at 6 h "
            "sampling; linear-interpolation error <= ~1,400 km near "
            "perigee (0.002 R_sun), smaller elsewhere. TESS HEO apogee "
            "~376,000 km = 0.54 R_sun: the Earth-center list is invalid "
            "for TESS at grazing impact parameters, which this run "
            "corrects"}
    if name == "soho":
        npz = OBSERVER_TABLE_DIR / "soho_sc_ephemeris.npz"
        if not npz.exists():
            raise SystemExit("no SOHO table; run --fetch-observer soho first")
        tab = np.load(npz)
        ident = _sha256(npz)
        obs = register_spacecraft_table_observer(
            "soho-spacecraft", tab["mjd_utc"], tab["xyz_au"], ident,
            max_gap_days=1.0)
        return obs, {
            "observer_table": str(npz.relative_to(REPO)),
            "observer_table_sha256": ident,
            "observer_accuracy": "JPL Horizons -21 SSB vectors at 6 h "
            "sampling; the ~178 d L1 halo period makes linear-"
            "interpolation error negligible (<< 0.001 R_sun). SOHO's "
            "halo cross-track amplitude ~0.9 R_sun (measured in the "
            "LASCO recon star check) invalidates the Earth-center list "
            "at grazing impact parameters, which this run corrects"}
    raise ValueError(name)


def build_request(target_ids, start: str, stop: str, observer: Observer,
                  coarse_step_days: float = 10.0) -> CrossingsRequest:
    return CrossingsRequest(
        target_ids=tuple(target_ids),
        link_directions=(LinkDirection.INBOUND, LinkDirection.OUTBOUND),
        intervals=(TimeInterval(interval_id="historical",
                                start=Time(start, scale="utc"),
                                stop=Time(stop, scale="utc")),),
        observer=observer,
        relay_distance_au=RELAY_DISTANCE_AU,
        model_id=MODEL_ID,
        beam_radii_au=(),
        report_max_b_au=None,
        coarse_step_days=coarse_step_days,
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    ap.add_argument("--observer", default="earth",
                    choices=sorted(OBSERVER_DEFAULTS))
    ap.add_argument("--fetch-observer", choices=sorted(HORIZONS_IDS),
                    help="fetch the Horizons table for this observer "
                         "and exit")
    ap.add_argument("--start", default=None)
    ap.add_argument("--stop", default=None)
    ap.add_argument("--coarse-step-days", type=float, default=10.0)
    ap.add_argument("--targets", nargs="*", help="subset of registry IDs (default: all)")
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)

    d_start, d_stop, d_dir = OBSERVER_DEFAULTS[a.observer]
    start, stop = a.start or d_start, a.stop or d_stop
    if a.fetch_observer:
        fs, fe, _ = OBSERVER_DEFAULTS[a.fetch_observer]
        fetch_horizons_table(a.fetch_observer, a.start or fs, a.stop or fe)
        return 0
    out = a.out or (REPO / "crossings" / d_dir)

    registry = load_target_registry(a.registry)
    ids = a.targets or sorted(registry.ids)
    observer, obs_notes = resolve_observer(a.observer)
    request = build_request(ids, start, stop, observer, a.coarse_step_days)
    print(f"{len(ids)} targets x 2 link directions, {start} -> {stop}, "
          f"observer={observer.observer_id}", file=sys.stderr)
    result = find_crossings(request, registry)

    generated = datetime.now(timezone.utc).isoformat(timespec="seconds")
    input_hashes = {str(a.registry.relative_to(REPO)): _sha256(a.registry)}
    if "observer_table" in obs_notes:
        input_hashes[obs_notes["observer_table"]] = \
            obs_notes["observer_table_sha256"]
    written = write_crossings_products(
        result, out, generated_utc=generated,
        input_file_hashes=input_hashes)
    n_inv = sum(1 for e in result.events if e.validity.value == "invalid")
    summary = {
        "crossings_id": result.crossings_id, "generated_utc": generated,
        "observer_id": observer.observer_id, "n_targets": len(ids),
        "n_events": len(result.events), "n_invalid": n_inv,
        "warnings": list(result.warnings),
        "files": {k: str(v) for k, v in written.items()},
        **obs_notes,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
