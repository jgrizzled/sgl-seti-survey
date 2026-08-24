#!/usr/bin/env bash
# WISE v2 programme runner (v2_plan.md §10 steps A–I). Every stage is
# idempotent / resumable; the invariant tests run before any batch.
# Usage: surveys/wise/scripts/run_v2.sh <stage>
#   A  freeze        hypothesis v2.0 freeze + stratified corridor hold-out
#   B  tests         pytest invariants (also run before every other stage)
#   C  geometry      observer table, covariance envelopes, asteroid control
#   D  build         v2 tensors + null summaries (all corridors)
#   D2 nulls-dev     null ensemble, FWER, mask sensitivity on the dev set
#   E  inject-dev    image-level injections + completeness on the dev set
#   F  adjudicate-dev
#   G  confirmatory  nulls -> injections -> completeness -> adjudication, once
#   I  report        tables from the ledger
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
REPO=$(cd "$(dirname "$0")/../../.." && pwd)
S=$REPO/surveys/wise/scripts
cd "$REPO"
stage=${1:-}
run_tests() { uv run pytest tests -q; }
case "$stage" in
  A|freeze)      uv run python $S/freeze_v2.py ;;
  B|tests)       run_tests ;;
  C|geometry)    run_tests
                 uv run python $S/build_observer_table.py
                 uv run python $S/geometry_check.py --workers 2
                 uv run python $S/asteroid_control.py --asteroid 14000 ;;
  D|build)       run_tests
                 uv run python $S/build_tensors.py --workers 7 --only-missing ;;
  D2|nulls-dev)  uv run python $S/null_ensemble.py --set dev ;;
  E|inject-dev)  run_tests
                 uv run python $S/inject_pipeline.py --set dev --workers 7 --only-missing
                 uv run python $S/completeness.py --set dev ;;
  F|adjudicate-dev) uv run python $S/adjudicate_v2.py --set dev ;;
  G|confirmatory) run_tests
                 uv run python $S/null_ensemble.py --set confirmatory
                 uv run python $S/inject_pipeline.py --set confirmatory --workers 7 --only-missing
                 uv run python $S/completeness.py --set confirmatory
                 uv run python $S/adjudicate_v2.py --set confirmatory ;;
  I|report)      uv run python $S/report_tables.py ;;
  *) echo "usage: $0 {A|B|C|D|D2|E|F|G|I}"; exit 2 ;;
esac
