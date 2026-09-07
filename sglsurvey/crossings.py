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
* ``--observer kepler`` — a JPL Horizons SSB vector table (spacecraft
  -227) fetched by ``--fetch-observer kepler``; Earth-trailing
  heliocentric orbit, up to ~1 AU from Earth (plan §5.8 item 8).
* ``--observer stereoa`` — a JPL Horizons SSB vector table (spacecraft
  -234) fetched by ``--fetch-observer stereoa``; heliocentric ~0.96 AU
  orbit leading Earth, so Earth-center is invalid at every rung (plan
  §5.15 O2, STEREO-A HI-1 sunward survey).
* ``--observer psp`` — a JPL Horizons SSB vector table (spacecraft
  -96) fetched by ``--fetch-observer psp`` at **10-min** sampling in
  yearly chunks (the 0.046–0.7 AU orbit sweeps ~34° of heliocentric
  longitude per 6 h at perihelion, so the 6 h default is useless
  there). The observer's star-side and anti-star-side minima can be
  only ~2 d apart around perihelion, so the list is built with
  ``--coarse-step-days 0.5`` (plan §5.15 O5, PSP/WISPR geometry pass).
* ``--observer solo`` — a JPL Horizons SSB vector table (spacecraft
  -144, Solar Orbiter) fetched by ``--fetch-observer solo`` at 10-min
  sampling like PSP; heliocentric 0.28-1.0 AU (period ~150-180 d,
  Venus-resonant), so the same fast-observer settings apply
  (``--coarse-step-days 0.5``; plan §5.23, SoloHI sunward survey).
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
    # Kepler science era: Q0 start 2009-05-02 .. K2 C19 end 2018-09-26
    # (spacecraft retired 2018-10-30); Earth-trailing heliocentric orbit
    # drifts ~0.1 AU/yr from Earth, so Earth-center is invalid at every
    # rung (plan §5.8 item 8 footprint intersect)
    "kepler": ("2009-05-01", "2018-11-01", "kepler_v1"),
    # STEREO-A SECCHI/HI science era (HI-1A first light 2006-12, routine
    # synoptic imaging from 2007-01); stop bounded by the Horizons -234
    # predicted-trajectory end (2026-12-18) — extend at refresh. The
    # spacecraft drifts ~22 deg/yr ahead of Earth on a 346-d heliocentric
    # orbit (0.96 AU), so Earth-center is invalid at every rung
    # (plan §5.15 O2)
    "stereoa": ("2007-01-01", "2026-12-01", "stereoa_v1"),
    # Parker Solar Probe: launch 2018-08-12, WISPR first light 2018-09,
    # first perihelion (E1) 2018-11-06; stop bounded by the Horizons -96
    # SPK end — extend at refresh. Heliocentric 0.046-0.73 AU orbit
    # (period 88-150 d): Earth-center is meaningless at every rung
    # (plan §5.15 O5)
    "psp": ("2018-08-15", "2026-12-01", "psp_v1"),
    # Solar Orbiter: launch 2020-02-10, SoloHI first light 2020-05
    # (cruise-phase remote-sensing checkouts), nominal mission from the
    # 2021-11-27 Earth flyby; the Horizons -144 SPK runs to 2030-11-20
    # so the era extends to the universal window end. Heliocentric
    # 0.28-1.0 AU, period ~150-180 d (Venus-resonant, inclination rising
    # to ~24 deg by 2029): Earth-center is meaningless at every rung
    # (plan §5.23)
    "solo": ("2020-05-01", "2028-01-01", "solo_v1"),
}
HORIZONS_IDS = {"tess": "-95", "soho": "-21", "kepler": "-227",
                "stereoa": "-234", "psp": "-96", "solo": "-144"}
#: per-observer Horizons sampling step (default 6 h); the inner-
#: heliosphere observers (PSP, Solar Orbiter) need 10 min
OBSERVER_FETCH_STEP = {"psp": "10m", "solo": "10m"}
#: longest single Horizons request (rows are capped server-side at ~90k)
FETCH_CHUNK_DAYS = 365
#: raw responses larger than this are stored gzipped (sha of the .gz)
RAW_GZIP_BYTES = 20_000_000


def _sha256(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


def _horizons_vectors(command: str, start: str, stop: str, step: str) -> str:
    import requests

    resp = requests.get(
        "https://ssd.jpl.nasa.gov/api/horizons.api",
        params={"format": "text", "COMMAND": f"'{command}'",
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
    return txt


def fetch_horizons_table(observer: str, start: str, stop: str,
                         step: str | None = None) -> Path:
    """Fetch an SSB ICRF position table for a spacecraft from JPL
    Horizons into ``crossings/observers/`` (raw response retained, npz
    keyed by MJD UTC). Requests longer than ``FETCH_CHUNK_DAYS`` are
    split into consecutive chunks (Horizons caps a response at ~90k
    rows) and concatenated; the raw file keeps every chunk's response.
    Returns the npz path."""
    import numpy as np

    step = step or OBSERVER_FETCH_STEP.get(observer, "6h")
    OBSERVER_TABLE_DIR.mkdir(parents=True, exist_ok=True)
    t0, t1 = Time(start, scale="utc"), Time(stop, scale="utc")
    edges = np.arange(t0.mjd, t1.mjd, FETCH_CHUNK_DAYS)
    edges = np.append(edges, t1.mjd)
    chunks, rows = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        a = Time(lo, format="mjd", scale="utc").isot[:19]
        b = Time(hi, format="mjd", scale="utc").isot[:19]
        txt = _horizons_vectors(HORIZONS_IDS[observer], a, b, step)
        chunks.append(txt)
        n0 = len(rows)
        for line in txt.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines():
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 5:
                continue
            jd = float(parts[0])
            if rows and jd <= rows[-1][0]:
                continue      # chunk boundary epoch repeated
            rows.append((jd, float(parts[2]), float(parts[3]),
                         float(parts[4])))
        print(f"  {a} -> {b}: {len(rows) - n0} rows", file=sys.stderr)
    raw = OBSERVER_TABLE_DIR / f"{observer}_horizons_response.txt"
    text = "\n".join(chunks)
    if len(text) > RAW_GZIP_BYTES:       # 10-min tables run to ~50 MB
        import gzip
        raw = raw.with_suffix(".txt.gz")
        with gzip.open(raw, "wt", compresslevel=9) as fh:
            fh.write(text)
    else:
        raw.write_text(text)
    arr = np.array(rows)
    t = Time(arr[:, 0], format="jd", scale="tdb")
    mjd_utc = t.utc.mjd
    npz = OBSERVER_TABLE_DIR / f"{observer}_sc_ephemeris.npz"
    np.savez_compressed(npz, mjd_utc=mjd_utc, xyz_au=arr[:, 1:4])
    meta = {"observer": observer, "horizons_command": HORIZONS_IDS[observer],
            "center": "500@0 (SSB)", "ref_plane": "FRAME (ICRF)",
            "start": start, "stop": stop, "step": step,
            "n_chunks": len(chunks), "n_rows": int(len(rows)),
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
    if name == "kepler":
        npz = OBSERVER_TABLE_DIR / "kepler_sc_ephemeris.npz"
        if not npz.exists():
            raise SystemExit("no Kepler table; run --fetch-observer kepler first")
        tab = np.load(npz)
        ident = _sha256(npz)
        obs = register_spacecraft_table_observer(
            "kepler-spacecraft", tab["mjd_utc"], tab["xyz_au"], ident,
            max_gap_days=1.0)
        return obs, {
            "observer_table": str(npz.relative_to(REPO)),
            "observer_table_sha256": ident,
            "observer_accuracy": "JPL Horizons -227 SSB vectors at 6 h "
            "sampling; heliocentric Earth-trailing orbit (period 372.5 d) "
            "makes linear-interpolation error negligible. The spacecraft "
            "was 0.04 AU (2009) to 1.14 AU (2018) from Earth, so the "
            "Earth-center list is invalid for Kepler at every rung; the "
            "crossing epochs shift by the trailing angle (weeks to months)"}
    if name == "stereoa":
        npz = OBSERVER_TABLE_DIR / "stereoa_sc_ephemeris.npz"
        if not npz.exists():
            raise SystemExit("no STEREO-A table; run --fetch-observer "
                             "stereoa first")
        tab = np.load(npz)
        ident = _sha256(npz)
        obs = register_spacecraft_table_observer(
            "stereoa-spacecraft", tab["mjd_utc"], tab["xyz_au"], ident,
            max_gap_days=1.0)
        return obs, {
            "observer_table": str(npz.relative_to(REPO)),
            "observer_table_sha256": ident,
            "observer_accuracy": "JPL Horizons -234 SSB vectors at 6 h "
            "sampling; heliocentric orbit at ~0.96 AU (period ~346 d) "
            "leading Earth by ~22 deg/yr (full circuit 2006-2023), so "
            "linear-interpolation error is negligible and the "
            "Earth-center list is invalid for STEREO-A at every rung; "
            "crossing epochs shift by the leading angle (weeks to months) "
            "and the Sun-star axis is sampled at a different heliocentric "
            "radius"}
    if name == "psp":
        npz = OBSERVER_TABLE_DIR / "psp_sc_ephemeris.npz"
        if not npz.exists():
            raise SystemExit("no PSP table; run --fetch-observer psp first")
        tab = np.load(npz)
        ident = _sha256(npz)
        obs = register_spacecraft_table_observer(
            "psp-spacecraft", tab["mjd_utc"], tab["xyz_au"], ident,
            max_gap_days=0.5)
        return obs, {
            "observer_table": str(npz.relative_to(REPO)),
            "observer_table_sha256": ident,
            "observer_accuracy": "JPL Horizons -96 SSB vectors at 10 min "
            "sampling (yearly chunks); at the 0.046 AU perihelion the "
            "chord sagitta over one step is ~2e-6 AU = 0.0003 R_sun, "
            "negligible. Heliocentric 0.046-0.73 AU orbit (period "
            "88-150 d, Venus-flyby shrinking): the Earth-center list is "
            "meaningless for PSP at every rung; the observer crosses "
            "each Sun-star axis twice per orbit at whatever heliocentric "
            "radius the orbit has at the star's (anti-)longitude"}
    if name == "solo":
        npz = OBSERVER_TABLE_DIR / "solo_sc_ephemeris.npz"
        if not npz.exists():
            raise SystemExit("no Solar Orbiter table; run --fetch-observer "
                             "solo first")
        tab = np.load(npz)
        ident = _sha256(npz)
        obs = register_spacecraft_table_observer(
            "solo-spacecraft", tab["mjd_utc"], tab["xyz_au"], ident,
            max_gap_days=0.5)
        return obs, {
            "observer_table": str(npz.relative_to(REPO)),
            "observer_table_sha256": ident,
            "observer_accuracy": "JPL Horizons -144 SSB vectors at 10 min "
            "sampling (yearly chunks); at the 0.28 AU perihelion the "
            "chord sagitta over one step is ~5e-8 AU, negligible. "
            "Heliocentric 0.28-1.0 AU orbit (period ~150-180 d, "
            "Venus-resonant): the Earth-center list is meaningless for "
            "Solar Orbiter at every rung; the observer crosses each "
            "Sun-star axis twice per orbit at whatever heliocentric "
            "radius the orbit has at the star's (anti-)longitude"}
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
    ap.add_argument("--coarse-step-days", type=float, default=10.0,
                    help="coarse-scan step of the minimum bracketing; "
                         "10 d for 1-AU-class observers, 0.5 d for PSP and Solar Orbiter")
    ap.add_argument("--step", default=None,
                    help="Horizons sampling step for --fetch-observer "
                         "(default 6h; psp and solo 10m)")
    ap.add_argument("--targets", nargs="*", help="subset of registry IDs (default: all)")
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)

    d_start, d_stop, d_dir = OBSERVER_DEFAULTS[a.observer]
    start, stop = a.start or d_start, a.stop or d_stop
    if a.fetch_observer:
        fs, fe, _ = OBSERVER_DEFAULTS[a.fetch_observer]
        fetch_horizons_table(a.fetch_observer, a.start or fs, a.stop or fe,
                             a.step)
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
