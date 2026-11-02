#!/usr/bin/env python3
"""Create one case (or the whole case matrix) under runs/.

Examples
--------
    python scripts/make_case.py --mesh medium --re 25000 --wall steel
    python scripts/make_case.py --all
"""
import argparse
import sys

from chtpipe.casegen import make_case
from chtpipe.config import CASE_MATRIX, MESHES, WALLS, CaseSpec


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mesh", choices=sorted(MESHES))
    ap.add_argument("--re", type=float)
    ap.add_argument("--wall", choices=sorted(WALLS), default="steel")
    ap.add_argument("--np", type=int, help="number of MPI processes (default depends on the mesh)")
    ap.add_argument("--all", action="store_true", help="generate the six cases of the study")
    ap.add_argument("--out", default="runs")
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args(argv)

    if a.all:
        specs = [CaseSpec(m, r, w, a.np) for m, r, w in CASE_MATRIX]
    elif a.mesh and a.re:
        specs = [CaseSpec(a.mesh, a.re, a.wall, a.np)]
    else:
        ap.error("give --mesh and --re, or --all")
    for s in specs:
        print(make_case(s, a.out, a.overwrite))


if __name__ == "__main__":
    sys.exit(main())
