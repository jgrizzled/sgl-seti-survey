#!/usr/bin/env bash
# Generic v2 chain for a survey directory with run.py: waits for its
# geometry check, then build -> dev (nulls, inject, completeness,
# adjudicate) -> confirmatory (once) -> report.
# Usage: bash sglsurvey/v2/chain.sh surveys/ztf [workers] [skip-build]
set -u
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$(dirname "$0")/../.."
D=$1; W=${2:-7}; SKIP=${3:-}
R=$D/run.py
until [ -f "$D/configs/cross_track_cells.json" ]; do sleep 60; done
echo "geometry ready $(date)"
if [ -z "$SKIP" ]; then
  uv run python $R build --workers $W 2>&1 | grep -v Warn | tail -80
fi
uv run python $R nulls --set dev 2>&1 | grep -v Warn | tail -4
uv run python $R inject --set dev --workers $W 2>&1 | grep -v Warn | tail -40
uv run python $R completeness --set dev 2>&1 | grep -v Warn | tail -30
uv run python $R adjudicate --set dev 2>&1 | grep -v Warn | tail -20
uv run python $R nulls --set confirmatory 2>&1 | grep -v Warn | tail -4
uv run python $R inject --set confirmatory --workers $W 2>&1 | grep -v Warn | tail -60
uv run python $R completeness --set confirmatory 2>&1 | grep -v Warn | tail -30
uv run python $R adjudicate --set confirmatory 2>&1 | grep -v Warn | tail -20
uv run python $R report 2>&1 | tail -5
echo "chain done $(date)"
