#!/bin/bash
# PTF corridor survey: post-freeze chain (hypotheses v1.0). Runs the
# remaining Pipeline A stages over every searchable corridor, then the
# v2 engine on the development set, and — only if the dev stage needs
# no amendment (checked by the operator between the two halves) — the
# blind confirmatory set once.
#   stage 1: bash run_chain.sh data     (fetch + zeropoints, all corridors)
#   stage 2: bash run_chain.sh dev      (freeze -> build dev -> nulls/inject/completeness/adjudicate)
#   stage 3: bash run_chain.sh conf     (build conf -> nulls/inject/completeness/adjudicate -> report)
set -eo pipefail
cd "$(dirname "$0")/../../.."
PY=.venv/bin/python
W=${WORKERS:-7}
filt() { grep -vE "IERS|polar|WARNING|FITSFixedWarning|non-standard|invalid value|RuntimeWarning|AstropyWarning|warnings.warn" || true; }
DEV=$($PY -c "import json;print(' '.join(json.load(open('surveys/ptf/configs/v2_freeze.json'))['split']['development']['corridors']))" 2>/dev/null || true)
CONF=$($PY -c "import json;print(' '.join(json.load(open('surveys/ptf/configs/v2_freeze.json'))['split']['confirmatory']['corridors']))" 2>/dev/null || true)

case "$1" in
  data)
    echo "=== fetch cutouts (all searchable corridors) ==="
    $PY surveys/ptf/scripts/fetch_cutouts.py --workers 6 2>&1 | filt | tail -3
    echo "=== zeropoints ==="
    $PY surveys/ptf/scripts/calibrate_zeropoints.py --workers $W 2>&1 | filt | tail -2
    $PY surveys/ptf/scripts/flux_scale_check.py 2>&1 | filt
    ;;
  dev)
    echo "=== freeze ==="
    $PY surveys/ptf/run.py freeze 2>&1 | filt
    DEV=$($PY -c "import json;print(' '.join(json.load(open('surveys/ptf/configs/v2_freeze.json'))['split']['development']['corridors']))")
    echo "=== build dev: $DEV ==="
    $PY surveys/ptf/run.py build --corridors $DEV --workers $W 2>&1 | filt | tail -25
    for st in nulls inject completeness adjudicate; do
      echo "=== $st dev ==="
      $PY surveys/ptf/run.py $st --set dev --workers $W 2>&1 | filt | tail -30
    done
    echo DEV-CHAIN-DONE
    ;;
  conf)
    echo "=== build confirmatory: $CONF ==="
    $PY surveys/ptf/run.py build --corridors $CONF --workers $W 2>&1 | filt | tail -50
    for st in nulls inject completeness adjudicate; do
      echo "=== $st confirmatory ==="
      $PY surveys/ptf/run.py $st --set confirmatory --workers $W 2>&1 | filt | tail -40
    done
    echo "=== report tables ==="
    $PY surveys/ptf/run.py report 2>&1 | filt | tail -20
    echo CONF-CHAIN-DONE
    ;;
  *) echo "usage: run_chain.sh {data|dev|conf}"; exit 2;;
esac
