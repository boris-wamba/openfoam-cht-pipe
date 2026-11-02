#!/bin/sh
# OpenFOAM-side post-processing of one finished case, then the Python analysis.
# Usage: scripts/post_case.sh runs/<case>
# Every step prints [ok] or [FAILED] (with the end of its log); nothing is silent.
CASE=${1:?usage: post_case.sh runs/<case>}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$CASE" || exit 1

step() {            # step <name> <command...>
    name=$1; shift
    if "$@" > "log.post.$name" 2>&1; then
        echo "[ok]     $name"
    else
        echo "[FAILED] $name   (full log: $CASE/log.post.$name)"
        tail -n 15 "log.post.$name"
        return 1
    fi
}

echo "case: $CASE"
LAST=$(ls -d [0-9]* 2>/dev/null | sort -g | tail -1)
echo "last time directory: ${LAST:-none}"
echo "fluid fields there : $(ls "$LAST/fluid" 2>/dev/null | tr '\n' ' ')"
[ -f "$LAST/fluid/phi" ] || echo "WARNING: no phi in $LAST/fluid (needed by outletTb)"

# bulk temperature at the outlet, for every written time except 0 (steady-state check)
step outletTb       postProcess -noZero -region fluid -func outletTb
# fields and profiles at the last time
step sampleWall     postProcess -region fluid -latestTime -func sampleWall
step sampleProfiles postProcess -region fluid -latestTime -func sampleProfiles
step sampleU        postProcess -region fluid -latestTime -func sampleU
step sampleLine     postProcess -region fluid -latestTime -func sampleLine
# (the yPlus, wallHeatFlux and wallShearStress functions are not used: with this multi-region setup they write
#  zeros, no rows, or fail; y+ is computed from the pressure-gradient friction velocity in the Python analysis)

echo "--- data rows in outletTb output:"
grep -vc '^#' postProcessing/fluid/outletTb/*/*.dat 2>/dev/null || echo "no file"

cd "$ROOT"
python3 scripts/analyse_case.py "$CASE" --out "results/$(basename "$CASE")"
