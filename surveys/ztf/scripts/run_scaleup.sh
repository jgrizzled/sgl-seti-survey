#!/bin/bash
# ZTF scale-up driver: runs the full chain over every overlay-visible
# corridor (pilot corridors are re-evaluated cheaply from cache).
# Usage: nohup bash surveys/ztf/scripts/run_scaleup.sh > runs/ztf/scaleup.log 2>&1 &
set -u
cd "$(dirname "$0")/../../.."
S=surveys/ztf/scripts
FROM=${1:-coarse}   # resume point: coarse|precise|screen|cutouts|tensors|calibrate
stage() { case "$FROM" in
  coarse) return 0;; precise) [ "$1" != coarse ];; screen) [[ "$1" != coarse && "$1" != precise ]];;
  cutouts) [[ "$1" == cutouts || "$1" == tensors || "$1" == calibrate ]];;
  tensors) [[ "$1" == tensors || "$1" == calibrate ]];; calibrate) [ "$1" == calibrate ];; esac; }
stage coarse  && { echo "=== coarse $(date)";   uv run python $S/coarse_discovery.py || exit 1; }
stage precise && { echo "=== precise $(date)";  uv run python $S/precise_pass.py --workers 6 || exit 1; }
stage screen  && { echo "=== screen $(date)";   uv run python $S/catalog_screen.py || exit 1;
  echo "=== recurrence $(date)"; uv run python $S/screen_recurrence.py > runs/ztf/screen_v1/recurrence_stdout.txt 2>&1 || echo "recurrence failed"; }
stage cutouts && { echo "=== cutouts $(date)";  uv run python $S/fetch_cutouts.py --workers 6 || exit 1; }
stage tensors && { echo "=== tensors $(date)";  uv run python $S/sample_tensor.py --only-missing || exit 1; }
echo "=== calibrate $(date)"; rm -rf runs/ztf/calib_v1/records; uv run python $S/injection_calibrate.py || exit 1
echo "=== figures $(date)";  uv run python $S/make_figures.py
echo "=== DONE $(date)"
