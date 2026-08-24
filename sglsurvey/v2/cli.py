"""Command-line driver for the v2 engine: ``python run.py <stage> [opts]``
from a survey's ``-v2`` directory, where ``run.py`` imports its profile.
Stages: freeze, build, geometry, nulls, inject, completeness,
adjudicate, report, control."""

from __future__ import annotations

import argparse
import os


def main(P, argv=None):
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    ap = argparse.ArgumentParser(prog=f"{P.name}-v2")
    ap.add_argument("stage", choices=["freeze", "build", "geometry", "nulls", "inject",
                                      "completeness", "adjudicate", "report", "control"])
    ap.add_argument("--set", choices=["dev", "confirmatory", "all"], default=None)
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--corridors", nargs="*", default=None)
    ap.add_argument("--endpoints", nargs="*", default=None)
    ap.add_argument("--no-only-missing", action="store_true")
    ap.add_argument("--n-per-cell", type=int, default=None)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args(argv)
    if a.stage == "freeze":
        from sglsurvey.v2 import freeze
        freeze.run(P)
    elif a.stage == "build":
        from sglsurvey.v2 import build
        build.run(P, corridors=a.corridors, workers=a.workers, only_missing=not a.no_only_missing)
    elif a.stage == "geometry":
        from sglsurvey.v2 import geometry_stage
        geometry_stage.run(P, workers=a.workers, endpoints=a.endpoints)
    elif a.stage == "nulls":
        from sglsurvey.v2 import nulls_stage
        nulls_stage.run(P, a.set or "dev", force=a.force)
    elif a.stage == "inject":
        from sglsurvey.v2 import inject_stage
        inject_stage.run(P, a.set or "dev", workers=a.workers, corridors=a.corridors,
                         only_missing=not a.no_only_missing, n_per_cell=a.n_per_cell)
    elif a.stage == "completeness":
        from sglsurvey.v2 import completeness_stage
        completeness_stage.run(P, a.set or "dev")
    elif a.stage == "adjudicate":
        from sglsurvey.v2 import adjudicate_stage
        adjudicate_stage.run(P, a.set or "dev")
    elif a.stage == "report":
        from sglsurvey.v2 import report_stage
        raise SystemExit(report_stage.run(P))
    elif a.stage == "control":
        raise SystemExit("the positive control is a survey script (scripts/asteroid_control_v2.py)")
