"""Cross-case summary: grid convergence (GCI), Reynolds sweep, wall-material effect."""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from .analysis import gci
from .config import PR_F
from .correlations import dittus_boelter, gnielinski


def load_summaries(results: str | Path) -> dict:
    out = {}
    for f in sorted(Path(results).glob("*/summary.json")):
        s = json.loads(f.read_text())
        out[s["name"]] = s
    return out


def mesh_study(summaries: dict, re_: float = 25000, wall: str = "steel") -> dict | None:
    """GCI of Nu and f on coarse / medium / fine at fixed Re and wall."""
    sel = {s["mesh"]: s for s in summaries.values() if s["Re"] == re_ and s["wall"] == wall}
    if not all(k in sel for k in ("coarse", "medium", "fine")):
        return None
    n1, n2, n3 = (sel[k]["n_cells"] for k in ("fine", "medium", "coarse"))
    r21 = math.sqrt(n1 / n2)
    r32 = math.sqrt(n2 / n3)
    res = {"r21": r21, "r32": r32, "cells": {k: sel[k]["n_cells"] for k in sel}}
    res["Nu_fd"] = gci(sel["fine"]["Nu_fd"], sel["medium"]["Nu_fd"], sel["coarse"]["Nu_fd"], r21, r32)
    if all("f_cfd" in sel[k] for k in sel):
        res["f"] = gci(sel["fine"]["f_cfd"], sel["medium"]["f_cfd"], sel["coarse"]["f_cfd"], r21, r32)
    return res


def reynolds_sweep(summaries: dict, mesh: str = "medium", wall: str = "steel") -> list:
    rows = [s for s in summaries.values() if s["mesh"] == mesh and s["wall"] == wall]
    return sorted(rows, key=lambda s: s["Re"])


def markdown_report(summaries: dict) -> str:
    L = ["# Summary of the study", ""]
    L += ["## All cases", "",
          "| case | cells | y+ (est.) | Nu CFD | Nu Gnielinski | err % | Nu D-B | err % | f CFD | f Petukhov | err % | energy bal. % |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for s in sorted(summaries.values(), key=lambda s: s["name"]):
        L.append("| {name} | {n_cells} | {yplus_est:.2f} | {Nu_fd:.1f} | {Nu_gnielinski:.1f} | {err_vs_gnielinski_pct:+.1f} | "
                 "{Nu_dittus_boelter:.1f} | {err_vs_dittus_boelter_pct:+.1f} | {f} | {fp} | {ef} | {eb} |".format(
                     f=f"{s['f_cfd']:.4f}" if "f_cfd" in s else "-",
                     fp=f"{s['f_petukhov']:.4f}" if "f_petukhov" in s else "-",
                     ef=f"{s['err_f_pct']:+.1f}" if "err_f_pct" in s else "-",
                     eb=f"{s['energy_balance_err_pct']:+.2f}" if "energy_balance_err_pct" in s else "-",
                     **s))
    ms = mesh_study(summaries)
    L += ["", "## Grid convergence (Re = 25000, steel)", ""]
    if ms is None:
        L.append("coarse, medium and fine cases are not all available.")
    else:
        for key, label in (("Nu_fd", "Nu"), ("f", "f")):
            if key in ms:
                g = ms[key]
                if g.get("status") == "ok":
                    L.append(f"- {label}: p = {g['p']:.2f}, extrapolated value = {g['phi_ext21']:.5g}, "
                             f"GCI(fine) = {100 * g['gci_fine']:.2f} %"
                             + (" (oscillatory convergence)" if g.get("oscillatory") else ""))
                else:
                    L.append(f"- {label}: {g.get('status')}")
    return "\n".join(L) + "\n"


def plots(summaries: dict, results: str | Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    results = Path(results)
    sweep = reynolds_sweep(summaries)
    if len(sweep) >= 2:
        re_a = np.array([s["Re"] for s in sweep])
        fig, ax = plt.subplots(figsize=(6, 4))
        rr = np.logspace(math.log10(re_a.min() * 0.7), math.log10(re_a.max() * 1.4), 50)
        ax.loglog(re_a, [s["Nu_fd"] for s in sweep], "o", label="OpenFOAM (medium mesh)")
        ax.loglog(rr, [gnielinski(r, PR_F) for r in rr], "k--", label="Gnielinski")
        ax.loglog(rr, [dittus_boelter(r, PR_F) for r in rr], "k:", label="Dittus-Boelter")
        ax.set_xlabel("Re")
        ax.set_ylabel("Nu")
        ax.legend()
        fig.tight_layout()
        fig.savefig(results / "Nu_vs_Re.png", dpi=150)
        plt.close(fig)

    mats = {s["wall"]: s for s in summaries.values() if s["mesh"] == "medium" and s["Re"] == 25000}
    if len(mats) == 2:
        fig, ax = plt.subplots(figsize=(7, 4))
        for wall, s in mats.items():
            d = np.genfromtxt(results / s["name"] / "Nu_x.csv", delimiter=",", names=True)
            ax.plot(d["xh_over_D"], d["Nu"], label=wall)
        ax.set_xlim(0, 15)
        ax.set_xlabel("(x - L_entry) / D")
        ax.set_ylabel("Nu")
        ax.legend()
        fig.tight_layout()
        fig.savefig(results / "Nu_wall_material.png", dpi=150)
        plt.close(fig)

