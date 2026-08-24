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
observers (WISE, SPHEREx, TESS, ...) are derivative runs keyed by
``observer_id`` and should be produced from the same request with the
per-epoch observer in :mod:`sglsurvey.geometry` when a survey needs them.

Usage::

    python -m sglsurvey.crossings --observer earth \
        --start 1980-01-01 --stop 2028-01-01 --out crossings/universal_v1
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


def _sha256(path: Path) -> str:
    return "sha256:" + hashlib.sha256(path.read_bytes()).hexdigest()


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
    ap.add_argument("--observer", default="earth", choices=["earth"])
    ap.add_argument("--start", default="1980-01-01")
    ap.add_argument("--stop", default="2028-01-01")
    ap.add_argument("--coarse-step-days", type=float, default=10.0)
    ap.add_argument("--targets", nargs="*", help="subset of registry IDs (default: all)")
    ap.add_argument("--out", type=Path, default=REPO / "crossings" / "universal_v1")
    a = ap.parse_args(argv)

    registry = load_target_registry(a.registry)
    ids = a.targets or sorted(registry.ids)
    observer = Observer.earth_center()
    request = build_request(ids, a.start, a.stop, observer, a.coarse_step_days)
    print(f"{len(ids)} targets x 2 link directions, {a.start} -> {a.stop}, "
          f"observer={observer.observer_id}", file=sys.stderr)
    result = find_crossings(request, registry)

    generated = datetime.now(timezone.utc).isoformat(timespec="seconds")
    written = write_crossings_products(
        result, a.out, generated_utc=generated,
        input_file_hashes={str(a.registry.relative_to(REPO)): _sha256(a.registry)})
    n_inv = sum(1 for e in result.events if e.validity.value == "invalid")
    summary = {
        "crossings_id": result.crossings_id, "generated_utc": generated,
        "observer_id": observer.observer_id, "n_targets": len(ids),
        "n_events": len(result.events), "n_invalid": n_inv,
        "warnings": list(result.warnings),
        "files": {k: str(v) for k, v in written.items()},
    }
    (a.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
