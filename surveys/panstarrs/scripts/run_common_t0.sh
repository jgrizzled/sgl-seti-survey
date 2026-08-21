#!/bin/bash
# Rebuild the PS1 sample tensors on the common mu reference epoch
# T0 = 59800 (ZTF's), so the joint PS1+ZTF stage can calibrate the full
# (z, mu) family. Re-fetches the purged image products per batch, writes
# tensors to runs/panstarrs/calib_t0_59800/tensors, purges again.
# Usage: nohup bash surveys/panstarrs/scripts/run_common_t0.sh > runs/panstarrs/common_t0.log 2>&1 &
set -u
cd "$(dirname "$0")/../../.."
S=surveys/panstarrs/scripts
T0=59800
OUT=runs/panstarrs/calib_t0_$T0/tensors
BATCH=${BATCH:-4}
ALL=$(uv run python -c "import sys; sys.path.insert(0,'$S'); import ps1_corridors as p; print(' '.join(p.ALL_CORRIDORS))")
set -- $ALL
while [ $# -gt 0 ]; do
  b=("${@:1:$BATCH}"); shift $(( $# < BATCH ? $# : BATCH ))
  echo "=== batch ${b[*]} $(date)"
  uv run python $S/fetch_cutouts.py --workers 6 --corridors "${b[@]}" || exit 1
  uv run python $S/sample_tensor.py --only-missing --t0 $T0 --out-dir $OUT --corridors "${b[@]}" || exit 1
  PURGE_TENSOR_DIR=$OUT uv run python $S/purge_products.py "${b[@]}"
  df -h . | tail -1
done
echo "=== DONE $(date)"
