#!/usr/bin/env bash
# PS1 v2: per batch of corridors, re-fetch the purged warp cutouts with
# the v1 fetcher, build the v2 tensors, run the image-level injections
# (all corridors; recovery decisions are made later per hold-out set),
# and purge the cutouts again. Requires the freeze and the geometry
# check to exist. Usage: nohup bash surveys/panstarrs/scripts/run_batches.sh > runs/panstarrs/v2/batches.log 2>&1 &
set -u -o pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$(dirname "$0")/../../.."
S1=surveys/panstarrs/scripts
R=surveys/panstarrs/run.py
BATCH=${BATCH:-4}
WORKERS=${WORKERS:-6}
ALL=$(uv run python -c "import sys; sys.path.insert(0,'surveys/panstarrs'); from profile import PROFILE as P; print(' '.join(P.corridors()))")
set -- $ALL
mkdir -p runs/panstarrs/v2
while [ $# -gt 0 ]; do
  b=("${@:1:$BATCH}"); shift $(( $# < BATCH ? $# : BATCH ))
  echo "=== batch ${b[*]} $(date)"
  uv run python $S1/fetch_cutouts.py --workers 6 --corridors "${b[@]}" || { echo "fetch failed"; exit 1; }
  uv run python $R build --workers $WORKERS --corridors "${b[@]}" 2>&1 | grep -v Warn || { echo "build failed"; exit 1; }
  uv run python $R inject --set all --workers $WORKERS --corridors "${b[@]}" 2>&1 | grep -v Warn || { echo "inject failed"; exit 1; }
  PURGE_TENSOR_DIR=runs/panstarrs/v2/tensors uv run python $S1/purge_products.py "${b[@]}"
  df -h . | tail -1
done
echo "=== DONE $(date)"
