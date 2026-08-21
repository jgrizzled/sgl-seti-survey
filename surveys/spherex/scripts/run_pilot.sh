#!/usr/bin/env bash
# Run the SPHEREx pilot chain end to end (plan §6). Idempotent: every
# stage skips work already recorded under runs/spherex/.
set -euo pipefail
cd "$(dirname "$0")/../../.."
S=surveys/spherex/scripts
uv run python $S/coarse_discovery.py "$@"
uv run python $S/precise_pass.py --workers "${WORKERS:-10}"
uv run python $S/star_control.py
uv run python $S/static_template.py
uv run python $S/track_screen.py
rm -rf runs/spherex/calib_v1/tensors runs/spherex/calib_v1/records
uv run python $S/sample_tensor.py
uv run python $S/injection_calibrate.py
