#!/bin/sh
# Run the remaining medium-mesh cases one after the other (existing cases are skipped).
# Usage: NP=4 scripts/run_queue.sh            (NP = number of MPI processes, default 4)
#        MPI_EXTRA_OPTS="--oversubscribe" NP=8 scripts/run_queue.sh
cd "$(dirname "$0")/.." || exit
NP=${NP:-4}
for c in "medium 25000 steel" "medium 10000 steel" "medium 50000 steel" "medium 25000 copper"
do
    echo "=== $c  ($(date))"
    scripts/run_case.sh $c "$NP" || echo "FAILED: $c"
done
python3 scripts/summarize.py results
