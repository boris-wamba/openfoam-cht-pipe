#!/bin/sh
# Run all six cases, cheapest first. Existing cases are skipped (see run_case.sh).
# Usage: NP=4 scripts/run_all.sh
cd "$(dirname "$0")/.." || exit
NP=${NP:-4}
for c in "coarse 25000 steel 1" "medium 25000 steel" "medium 10000 steel" "medium 50000 steel" "medium 25000 copper" "fine 25000 steel"
do
    scripts/run_case.sh $c ${NP} || echo "FAILED: $c"
done
python3 scripts/summarize.py results
