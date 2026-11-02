#!/usr/bin/env python3
"""Analyse one finished case: python scripts/analyse_case.py runs/medium_Re25000_steel"""
import argparse
import json

from chtpipe.analysis import analyse_case

ap = argparse.ArgumentParser()
ap.add_argument("case")
ap.add_argument("--out", help="output directory (default: results/<case name>)")
ap.add_argument("--no-plots", action="store_true")
a = ap.parse_args()
print(json.dumps(analyse_case(a.case, a.out, not a.no_plots), indent=2))
