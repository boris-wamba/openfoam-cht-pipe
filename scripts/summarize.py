#!/usr/bin/env python3
"""Cross-case summary (GCI, Reynolds sweep, wall material): python scripts/summarize.py [results]"""
import sys
from pathlib import Path

from chtpipe.summary import load_summaries, markdown_report, plots

results = Path(sys.argv[1] if len(sys.argv) > 1 else "results")
s = load_summaries(results)
if not s:
    sys.exit(f"no summary.json found under {results}/*/")
md = markdown_report(s)
(results / "summary.md").write_text(md)
plots(s, results)
print(md)
