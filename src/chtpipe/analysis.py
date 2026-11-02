"""Post-processing of one case: Nusselt number, friction factor, bulk temperature, GCI.

Inputs are the files written by ``scripts/post_case.sh`` (OpenFOAM ``postProcess``
function objects) and ``case_meta.json`` written by the case generator.
"""
from __future__ import annotations

import json
import math
import re
from pathlib import Path

import numpy as np

from .correlations import dittus_boelter, gnielinski, petukhov_friction

# ---------------------------------------------------------------- file helpers


def latest_time_dir(parent: Path) -> Path:
    """Directory named by the largest number under ``parent`` (searched at any depth).

    Sampling output is ``postProcessing/<func>/<time>/`` for surfaces and
    ``postProcessing/<func>/<region>/<time>/`` for sets.
    """
    dirs = []
    for d in Path(parent).rglob("*"):
        if d.is_dir():
            try:
                dirs.append((float(d.name), len(d.parts), d))
            except ValueError:
                continue
    if not dirs:
        raise FileNotFoundError(f"no time directory in {parent}")
    best = max(t for t, _, _ in dirs)
    return min((n, d) for t, n, d in dirs if t == best)[1]


def find_output_dir(case: Path, func: str) -> Path:
    """``postProcessing/[region/]func`` (the region level is present for multi-region runs)."""
    hits = sorted(Path(case, "postProcessing").glob(f"**/{func}"))
    hits = [h for h in hits if h.is_dir()]
    if not hits:
        raise FileNotFoundError(f"postProcessing/**/{func} not found in {case}")
    return hits[0]


def read_table(path: Path) -> tuple[np.ndarray, list[str] | None]:
    """Numeric table, ignoring comment lines; also returns the last comment as column names."""
    header = None
    rows = []
    for line in Path(path).read_text().splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith("#"):
            header = s.lstrip("#").split()
            continue
        s = s.replace("(", " ").replace(")", " ")
        try:
            rows.append([float(v) for v in s.split()])
        except ValueError:
            continue
    if not rows:
        raise ValueError(f"no numeric data in {path}")
    n = min(len(r) for r in rows)
    return np.array([r[:n] for r in rows]), header


def sampled_file(case: Path, func: str, key: str) -> Path:
    """File of the latest time written by sampling function ``func`` whose name contains ``key``."""
    tdir = latest_time_dir(find_output_dir(case, func))
    files = sorted(f for f in tdir.iterdir() if f.is_file() and key in f.name)
    if not files:
        raise FileNotFoundError(f"no file containing {key!r} in {tdir}")
    return files[0]


def load_meta(case: Path) -> dict:
    return json.loads(Path(case, "case_meta.json").read_text())


# ---------------------------------------------------------------- physics


def bulk_temperature(x: np.ndarray, m: dict) -> np.ndarray:
    """Bulk temperature from the global energy balance (uniform inner-wall heat flux)."""
    area_rate = 2.0 * math.pi * m["R_IN"]            # heated inner perimeter [m]
    xh = np.clip(np.asarray(x) - m["L_ENTRY"], 0.0, m["L_HEATED"])
    return m["T_IN"] + m["Q_INNER"] * area_rate * xh / (m["MDOT"] * m["CP_F"])


def nusselt_profile(x: np.ndarray, t_wall: np.ndarray, m: dict) -> dict:
    """Local Nu(x) = q_i D / (kappa (Tw - Tb)) for x in the heated section."""
    x = np.asarray(x, dtype=float)
    mask = x > m["L_ENTRY"]
    tb = bulk_temperature(x, m)
    dt = np.asarray(t_wall, dtype=float) - tb
    nu = np.full_like(x, np.nan)
    nu[mask] = m["Q_INNER"] * m["D"] / (m["K_F"] * dt[mask])
    return {"x": x, "xh_over_D": (x - m["L_ENTRY"]) / m["D"], "Tw": np.asarray(t_wall), "Tb": tb, "Nu": nu}


def window_mean(x: np.ndarray, y: np.ndarray, x0: float, x1: float) -> float:
    """Length-weighted mean of y over [x0, x1]."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    sel = (x >= x0) & (x <= x1) & np.isfinite(y)
    if sel.sum() < 2:
        raise ValueError("fewer than two points in the averaging window")
    xs, ys = x[sel], y[sel]
    trap = getattr(np, "trapezoid", None) or np.trapz
    return float(trap(ys, xs) / (xs[-1] - xs[0]))


def friction_factor(x: np.ndarray, p: np.ndarray, m: dict, x0: float, x1: float) -> dict:
    """Darcy friction factor from the linear fit of p(x) on [x0, x1]."""
    x = np.asarray(x, dtype=float)
    p = np.asarray(p, dtype=float)
    sel = (x >= x0) & (x <= x1)
    slope = float(np.polyfit(x[sel], p[sel], 1)[0])    # dp/dx [Pa/m]
    f = -slope * m["D"] / (0.5 * m["RHO_F"] * m["UB"] ** 2)
    return {"dpdx": slope, "f": f}


def gci(phi1: float, phi2: float, phi3: float, r21: float, r32: float, fs: float = 1.25) -> dict:
    """Grid convergence index (Celik et al., 2008). phi1: finest, phi3: coarsest."""
    e21, e32 = phi2 - phi1, phi3 - phi2
    out = {"phi1": phi1, "phi2": phi2, "phi3": phi3, "r21": r21, "r32": r32}
    if e21 == 0.0 or e32 == 0.0:
        out.update(p=float("nan"), status="no_difference")
        return out
    s = math.copysign(1.0, e32 / e21)
    if s < 0:
        out["oscillatory"] = True
    p = abs(math.log(abs(e32 / e21))) / math.log(r21)
    for _ in range(200):
        q = math.log((r21 ** p - s) / (r32 ** p - s))
        p_new = abs(math.log(abs(e32 / e21)) + q) / math.log(r21)
        if abs(p_new - p) < 1e-10:
            p = p_new
            break
        p = p_new
    phi_ext = (r21 ** p * phi1 - phi2) / (r21 ** p - 1.0)
    out.update(
        p=p,
        phi_ext21=phi_ext,
        e_a21=abs((phi1 - phi2) / phi1),
        e_ext21=abs((phi_ext - phi1) / phi_ext),
        gci_fine=fs * abs((phi1 - phi2) / phi1) / (r21 ** p - 1.0),
        status="ok",
    )
    return out


# ---------------------------------------------------------------- one case


def read_outlet_tb(case: Path) -> tuple[np.ndarray, np.ndarray]:
    """Time and outlet bulk temperature; files without data rows (old/failed runs) are skipped."""
    d = find_output_dir(case, "outletTb")
    rows = []
    for f in sorted(d.glob("*/surfaceFieldValue*.dat")):
        try:
            data, _ = read_table(f)
        except ValueError:
            continue
        rows.append(data[:, :2])
    if not rows:
        raise FileNotFoundError(f"no data rows in {d}/*/surfaceFieldValue*.dat (rerun postProcess -noZero)")
    arr = np.vstack(rows)
    arr = arr[np.argsort(arr[:, 0])]
    _, idx = np.unique(arr[:, 0], return_index=True)
    arr = arr[idx]
    return arr[:, 0], arr[:, 1]


def analyse_case(case: str | Path, outdir: str | Path | None = None, make_plots: bool = True) -> dict:
    case = Path(case)
    m = load_meta(case)
    out = Path(outdir) if outdir else Path("results") / m["name"]
    out.mkdir(parents=True, exist_ok=True)
    x0 = m["L_ENTRY"] + m["FD_WINDOW_D"][0] * m["D"]
    x1 = m["L_ENTRY"] + m["FD_WINDOW_D"][1] * m["D"]
    re_ = m["re"]

    # wall temperature -> Nusselt
    data, _ = read_table(sampled_file(case, "sampleWall", ""))
    order = np.argsort(data[:, 0])
    x, tw = data[order, 0], data[order, -1]
    prof = nusselt_profile(x, tw, m)
    nu_fd = window_mean(prof["x"], prof["Nu"], x0, x1)
    nu_gn = gnielinski(re_, m["PR_F"])
    nu_db = dittus_boelter(re_, m["PR_F"])

    summary = {
        "name": m["name"], "mesh": m["mesh"], "Re": re_, "wall": m["wall"],
        "n_cells": (m["NX_E"] + m["NX_H"]) * (m["NR_F"] + m["NR_S"]),
        "n_cells_fluid": (m["NX_E"] + m["NX_H"]) * m["NR_F"],
        "yplus_est": m["YPLUS_EST"],
        "window_xh_over_D": m["FD_WINDOW_D"],
        "Nu_fd": nu_fd, "Nu_gnielinski": nu_gn, "Nu_dittus_boelter": nu_db,
        "err_vs_gnielinski_pct": 100.0 * (nu_fd - nu_gn) / nu_gn,
        "err_vs_dittus_boelter_pct": 100.0 * (nu_fd - nu_db) / nu_db,
        "Tw_mean_fd": window_mean(prof["x"], prof["Tw"], x0, x1),
    }

    # pressure gradient -> friction factor
    try:
        pdat, _ = read_table(sampled_file(case, "sampleLine", "mid"))
        fr = friction_factor(pdat[:, 0], pdat[:, -1], m, m["L_ENTRY"] + 5 * m["D"], m["L_TOTAL"] - 5 * m["D"])
        f_ref = petukhov_friction(re_)
        summary.update(f_cfd=fr["f"], dpdx=fr["dpdx"], f_petukhov=f_ref,
                       err_f_pct=100.0 * (fr["f"] - f_ref) / f_ref)
        u_tau = m["UB"] * math.sqrt(fr["f"] / 8.0)
        nu_f = m["MU_F"] / m["RHO_F"]
        summary.update(u_tau=u_tau, yplus_first_cell=0.5 * m["FIRST_CELL"] * u_tau / nu_f)
    except FileNotFoundError:
        u_tau = None

    # outlet bulk temperature: steady state and energy balance
    try:
        t, tb_out = read_outlet_tb(case)
        tb_th = float(bulk_temperature(np.array([m["L_TOTAL"]]), m)[0])
        rise = tb_th - m["T_IN"]
        summary.update(
            Tb_out_cfd=float(tb_out[-1]), Tb_out_theory=tb_th,
            energy_balance_err_pct=100.0 * (tb_out[-1] - tb_th) / rise,
            drift_last_interval_pct=(100.0 * (tb_out[-1] - tb_out[-2]) / rise) if len(tb_out) > 1 else None,
            n_written_times=int(len(t)),
        )
    except FileNotFoundError:
        t = tb_out = None

    pd_csv = out / "Nu_x.csv"
    with open(pd_csv, "w") as fh:
        fh.write("x_m,xh_over_D,Tw_K,Tb_K,Nu\n")
        for row in zip(prof["x"], prof["xh_over_D"], prof["Tw"], prof["Tb"], prof["Nu"]):
            fh.write(",".join(f"{v:.8g}" for v in row) + "\n")
    (out / "summary.json").write_text(json.dumps(summary, indent=2))

    if make_plots:
        _plots(case, out, m, prof, summary, t, tb_out, u_tau)
    return summary


def _plots(case, out, m, prof, s, t, tb_out, u_tau):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Nu(x)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(prof["xh_over_D"], prof["Nu"], label="OpenFOAM (SST, CHT)")
    ax.axhline(s["Nu_gnielinski"], ls="--", c="k", label="Gnielinski")
    ax.axhline(s["Nu_dittus_boelter"], ls=":", c="gray", label="Dittus-Boelter")
    ax.axvspan(*m["FD_WINDOW_D"], color="0.9", label="averaging window")
    ax.set_xlim(0, m["L_HEATED"] / m["D"])
    ax.set_ylim(0, 2.0 * s["Nu_gnielinski"])
    ax.set_xlabel("(x - L_entry) / D")
    ax.set_ylabel("Nu")
    ax.set_title(f"{m['name']}: local Nusselt number")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "Nu_x.png", dpi=150)
    plt.close(fig)

    # Tw, Tb
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(prof["x"], prof["Tw"], label="wall (fluid side)")
    ax.plot(prof["x"], prof["Tb"], label="bulk (energy balance)")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("T [K]")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "T_wall_bulk.png", dpi=150)
    plt.close(fig)

    # convergence of the outlet bulk temperature
    if t is not None:
        fig, ax = plt.subplots(figsize=(6, 3.5))
        ax.plot(t, tb_out, "o-")
        ax.axhline(s["Tb_out_theory"], ls="--", c="k", label="energy balance")
        ax.set_xlabel("time [s]")
        ax.set_ylabel("outlet bulk T [K]")
        ax.legend()
        fig.tight_layout()
        fig.savefig(out / "outlet_Tb_history.png", dpi=150)
        plt.close(fig)

    # radial temperature profiles
    try:
        tdir = latest_time_dir(find_output_dir(case, "sampleProfiles"))
        fig, ax = plt.subplots(figsize=(6, 4))
        files = sorted(tdir.glob("prof_h*")) or sorted(f for f in tdir.iterdir() if f.is_file())
        for f in files:
            d, _ = read_table(f)
            ax.plot(d[:, 1] * 1e3, d[:, -1], label=f.name.split("_T")[0])
        ax.set_xlabel("y [mm]")
        ax.set_ylabel("T [K]")
        ax.legend()
        fig.tight_layout()
        fig.savefig(out / "T_profiles.png", dpi=150)
        plt.close(fig)
    except FileNotFoundError:
        pass

    # law of the wall
    if u_tau:
        try:
            d, _ = read_table(sampled_file(case, "sampleU", "uprof"))
            y_wall = m["R_IN"] * math.cos(math.radians(m["HALF_ANGLE_DEG"]))
            dist = y_wall - d[:, 1]
            yp = dist * u_tau / (m["MU_F"] / m["RHO_F"])
            up = d[:, 3] / u_tau
            sel = yp > 0.5
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.semilogx(yp[sel], up[sel], label="OpenFOAM")
            yy = np.logspace(0, math.log10(max(yp[sel].max(), 10)), 100)
            ax.semilogx(yy[yy < 11], yy[yy < 11], "k:", label="u+ = y+")
            ax.semilogx(yy[yy > 11], np.log(yy[yy > 11]) / 0.41 + 5.2, "k--", label="log law")
            ax.set_xlabel("y+")
            ax.set_ylabel("u+")
            ax.legend()
            fig.tight_layout()
            fig.savefig(out / "u_plus.png", dpi=150)
            plt.close(fig)
        except FileNotFoundError:
            pass
