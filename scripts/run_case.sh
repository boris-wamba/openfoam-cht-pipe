#!/bin/sh
# Generate, run and post-process one case. An existing case is never overwritten
# unless FORCE=1 is set (a regenerated case erases its results).
# Usage: scripts/run_case.sh <mesh> <Re> <wall> [np]
#   e.g. scripts/run_case.sh medium 25000 steel 4
set -e
MESH=${1:?mesh: coarse|medium|fine}
RE=${2:?Reynolds number}
WALL=${3:-steel}
NP=${4:-}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT"
CASE="runs/${MESH}_Re${RE}_${WALL}"
if [ -d "$CASE" ] && [ -z "$FORCE" ]
then
    echo "SKIP: $CASE already exists (use FORCE=1 to regenerate it and erase its results)"
    exit 0
fi
python3 scripts/make_case.py --mesh "$MESH" --re "$RE" --wall "$WALL" --overwrite ${NP:+--np "$NP"}
"$CASE/Allrun"
scripts/post_case.sh "$CASE"
