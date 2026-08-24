#!/bin/bash
# Pipeline A + zeropoints over the confirmatory corridors (freeze
# split, hypotheses v1.0). One corridor at a time: coarse -> precise
# -> screen -> fetch -> zeropoints. Frozen exposure selection is
# applied in precise_pass; the deepest-per-night cap at tensor build.
set -e
cd "$(dirname "$0")/../../.."

endpoints_of() {
  case "$1" in
    ross248)  echo "ross-248";;
    cyg61)    echo "61-cyg-a 61-cyg-b";;
    grb34)    echo "groombridge-34-a groombridge-34-b";;
    gj1221)   echo "gj-1221";;
    gj338)    echo "gj-338-a gj-338-b";;
    gj625)    echo "gj-625";;
    gj687)    echo "gj-687";;
    gj251)    echo "gj-251";;
    wolf1069) echo "wolf-1069";;
    gj3512)   echo "gj-3512";;
    gj13157)  echo "gj-13157";;
  esac
}

ORDER="ross248 cyg61 grb34 gj338 gj625 gj251 wolf1069 gj3512 gj13157 gj687 gj1221"
for c in $ORDER; do
  eps=$(endpoints_of $c)
  echo "=== [$c] coarse ==="
  uv run python surveys/decam/scripts/coarse_discovery.py $eps 2>&1 | grep -vE "IERS|polar|WARNING" | tail -4
  echo "=== [$c] precise ==="
  uv run python surveys/decam/scripts/precise_pass.py $eps 2>&1 | grep -vE "IERS|polar|WARNING|invalid|_oi\[|non-standard" | tail -4
  echo "=== [$c] screen ==="
  uv run python surveys/decam/scripts/catalog_screen.py $c 2>&1 | grep -vE "IERS|polar|WARNING" | tail -4
  echo "=== [$c] fetch ==="
  uv run python surveys/decam/scripts/fetch_cutouts.py --corridors $c --workers 4 2>&1 | tail -3
  echo "=== [$c] zeropoints ==="
  uv run python surveys/decam/scripts/calibrate_zeropoints.py --corridors $c 2>&1 | grep -vE "IERS|polar|WARNING|invalid|_oi\[" | tail -2
done
echo CONFIRMATORY-PIPELINE-DONE
