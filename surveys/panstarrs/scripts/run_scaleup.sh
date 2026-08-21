#!/bin/bash
# PS1 scale-up driver: global coarse/precise/screen stages over every
# overlay-visible corridor, then per-batch cutouts -> tensors -> purge
# (disk budget), then calibration + figures.
# Usage: nohup bash surveys/panstarrs/scripts/run_scaleup.sh [from] > runs/panstarrs/scaleup.log 2>&1 &
set -u
cd "$(dirname "$0")/../../.."
S=surveys/panstarrs/scripts
FROM=${1:-coarse}   # resume point: coarse|precise|screen|batches|calibrate
BATCH=${BATCH:-4}
stage() { case "$FROM" in
  coarse) return 0;; precise) [ "$1" != coarse ];; screen) [[ "$1" != coarse && "$1" != precise ]];;
  batches) [[ "$1" == batches || "$1" == calibrate ]];; calibrate) [ "$1" == calibrate ];; esac; }
stage coarse  && { echo "=== coarse $(date)";   uv run python $S/coarse_discovery.py || exit 1; }
stage precise && { echo "=== precise $(date)";  uv run python $S/precise_pass.py --workers 8 || exit 1; }
stage screen  && { echo "=== screen $(date)";   uv run python $S/catalog_screen.py || exit 1;
  echo "=== recurrence $(date)"; uv run python $S/screen_recurrence.py > runs/panstarrs/screen_v1/recurrence_stdout.txt 2>&1 || echo "recurrence failed"; }
if stage batches; then
  QUEUE=$(uv run python -c "import sys; sys.path.insert(0,'$S'); import ps1_corridors as p; print(' '.join(c for c in p.QUEUE))")
  set -- $QUEUE
  while [ $# -gt 0 ]; do
    b=("${@:1:$BATCH}"); shift $(( $# < BATCH ? $# : BATCH ))
    echo "=== batch ${b[*]} $(date)"
    uv run python $S/fetch_cutouts.py --workers 6 --corridors "${b[@]}" || exit 1
    uv run python $S/sample_tensor.py --only-missing --corridors "${b[@]}" || exit 1
    uv run python $S/purge_products.py "${b[@]}"
    df -h . | tail -1
  done
fi
echo "=== calibrate $(date)"; rm -rf runs/panstarrs/calib_v1/records; uv run python $S/injection_calibrate.py || exit 1
echo "=== figures $(date)";  uv run python $S/make_figures.py
echo "=== DONE $(date)"
