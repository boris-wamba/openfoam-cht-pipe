"""Post-processing checks on synthetic OpenFOAM-like output with a known answer."""
import json
import math

import numpy as np
import pytest

from chtpipe.analysis import (analyse_case, bulk_temperature, friction_factor, gci,
                              nusselt_profile, read_table, window_mean)
from chtpipe.casegen import make_case
from chtpipe.config import CaseSpec, derive
from chtpipe.correlations import dittus_boelter, gnielinski, petukhov_friction
from chtpipe.summary import markdown_report, mesh_study


def test_correlations_reference_values():
    # Gnielinski, Re = 1e4, Pr = 6.13 ; textbook range 75-80
    assert 70 < gnielinski(1e4, 6.13) < 82
    assert petukhov_friction(1e4) == pytest.approx((0.79 * math.log(1e4) - 1.64) ** -2) == pytest.approx(0.0315, abs=1e-4)
    assert dittus_boelter(1e4, 6.13) == pytest.approx(0.023 * 1e4 ** 0.8 * 6.13 ** 0.4)


def test_bulk_temperature_energy_balance():
    m = derive(CaseSpec("medium", 25000, "steel"))
    x = np.array([0.0, m["L_ENTRY"], m["L_TOTAL"]])
    tb = bulk_temperature(x, m)
    assert tb[0] == tb[1] == m["T_IN"]
    # Q = q * pi D Lh  must equal mdot cp dT
    q_tot = m["Q_INNER"] * math.pi * m["D"] * m["L_HEATED"]
    assert (tb[2] - tb[1]) * m["MDOT"] * m["CP_F"] == pytest.approx(q_tot, rel=1e-12)


def test_outer_flux_is_consistent_with_inner_flux():
    m = derive(CaseSpec("medium", 25000, "steel"))
    assert m["Q_OUTER"] * m["R_OUT"] == pytest.approx(m["Q_INNER"] * m["R_IN"], rel=1e-12)
    assert m["GRAD_OUTER"] == pytest.approx(m["Q_OUTER"] / m["KAPPA_S"])
    assert m["GRAD_OUTER"] > 0     # heating: boundary hotter than the interior (outward normal gradient)


def test_nusselt_roundtrip():
    m = derive(CaseSpec("medium", 25000, "steel"))
    x = np.linspace(0.0, m["L_TOTAL"], 400)
    nu_true = 150.0
    tw = bulk_temperature(x, m) + m["Q_INNER"] * m["D"] / (m["K_F"] * nu_true)
    prof = nusselt_profile(x, tw, m)
    assert np.all(np.isnan(prof["Nu"][x <= m["L_ENTRY"]]))
    assert np.nanmean(prof["Nu"]) == pytest.approx(nu_true)


def test_window_mean_and_friction():
    x = np.linspace(0, 1, 101)
    assert window_mean(x, 2 * x, 0.2, 0.8) == pytest.approx(1.0)
    m = derive(CaseSpec("medium", 25000, "steel"))
    f_true = 0.025
    dpdx = -f_true * 0.5 * m["RHO_F"] * m["UB"] ** 2 / m["D"]
    xx = np.linspace(0, m["L_TOTAL"], 200)
    r = friction_factor(xx, 1e5 + dpdx * xx, m, 0.3, 1.2)
    assert r["f"] == pytest.approx(f_true, rel=1e-9)


def test_gci_recovers_second_order():
    exact, c, p = 100.0, 4.0, 2.0
    h = {1: 1.0, 2: 2.0, 3: 4.0}
    phi = {k: exact + c * v ** p for k, v in h.items()}
    g = gci(phi[1], phi[2], phi[3], 2.0, 2.0)
    assert g["status"] == "ok"
    assert g["p"] == pytest.approx(2.0, rel=1e-6)
    assert g["phi_ext21"] == pytest.approx(exact, rel=1e-9)
    assert g["gci_fine"] == pytest.approx(1.25 * abs((phi[1] - phi[2]) / phi[1]) / 3.0, rel=1e-6)


def test_gci_flags_oscillatory_and_zero_difference():
    assert gci(1.0, 1.0, 2.0, 2.0, 2.0)["status"] == "no_difference"
    assert gci(1.0, 1.2, 1.1, 2.0, 2.0).get("oscillatory") is True


def test_read_table_handles_headers_and_vectors(tmp_path):
    f = tmp_path / "a.xy"
    f.write_text("# x y z T\n1 2 3 4\n5 6 7 8\n")
    d, h = read_table(f)
    assert d.shape == (2, 4) and h == ["x", "y", "z", "T"]
    f.write_text("1 2 3 (4 5 6)\n")
    assert read_table(f)[0].shape == (1, 6)


def _fake_run(tmp_path, nu_true=140.0, f_true=0.026, mesh="medium"):
    case = make_case(CaseSpec(mesh, 25000, "steel"), tmp_path / "runs")
    m = json.loads((case / "case_meta.json").read_text())
    pp = case / "postProcessing" / "fluid"
    # wall temperature
    n = m["NX_E"] + m["NX_H"]
    x = (np.arange(n) + 0.5) * m["L_TOTAL"] / n
    tw = bulk_temperature(x, m) + m["Q_INNER"] * m["D"] / (m["K_F"] * nu_true)
    d = pp / "sampleWall" / "3.5"
    d.mkdir(parents=True)
    with open(d / "wall.raw", "w") as fh:
        fh.write("# x y z T\n")
        for xi, ti in zip(x, tw):
            fh.write(f"{xi} 0.00999 0 {ti}\n")
    # pressure line
    dpdx = -f_true * 0.5 * m["RHO_F"] * m["UB"] ** 2 / m["D"]
    d = pp / "sampleLine" / "3.5"
    d.mkdir(parents=True)
    xs = np.linspace(1e-4, m["L_TOTAL"] - 1e-4, 750)
    np.savetxt(d / "mid_p.xy", np.c_[xs, 0.005 * np.ones_like(xs), np.zeros_like(xs), 1e5 + dpdx * xs])
    # outlet bulk temperature history
    d = pp / "outletTb" / "0"
    d.mkdir(parents=True)
    tb_th = float(bulk_temperature(np.array([m["L_TOTAL"]]), m)[0])
    with open(d / "surfaceFieldValue.dat", "w") as fh:
        fh.write("# Time weightedAverage(T)\n")
        for t, v in ((1.0, tb_th - 2.0), (2.0, tb_th - 0.01), (3.5, tb_th - 0.001)):
            fh.write(f"{t} {v}\n")
    return case, m, tb_th


def test_analyse_case_on_synthetic_output(tmp_path):
    case, m, tb_th = _fake_run(tmp_path, nu_true=140.0, f_true=0.026)
    s = analyse_case(case, tmp_path / "res")
    assert s["Nu_fd"] == pytest.approx(140.0, rel=1e-6)
    assert s["Nu_gnielinski"] == pytest.approx(gnielinski(25000, m["PR_F"]))
    assert s["f_cfd"] == pytest.approx(0.026, rel=1e-6)
    u_tau = m["UB"] * math.sqrt(0.026 / 8.0)
    assert s["yplus_first_cell"] == pytest.approx(0.5 * m["FIRST_CELL"] * u_tau * m["RHO_F"] / m["MU_F"], rel=1e-6)
    assert s["Tb_out_theory"] == pytest.approx(tb_th)
    assert s["energy_balance_err_pct"] == pytest.approx(-100 * 0.001 / (tb_th - m["T_IN"]), rel=1e-6)
    assert s["n_written_times"] == 3
    for name in ("Nu_x.csv", "summary.json", "Nu_x.png", "T_wall_bulk.png", "outlet_Tb_history.png"):
        assert (tmp_path / "res" / name).exists(), name


def test_mesh_study_and_report():
    def mk(mesh, n, nu, f):
        return {"name": f"{mesh}_Re25000_steel", "mesh": mesh, "Re": 25000, "wall": "steel", "n_cells": n,
                "yplus_est": 1.0, "Nu_fd": nu, "Nu_gnielinski": 150.0, "Nu_dittus_boelter": 160.0,
                "err_vs_gnielinski_pct": 1.0, "err_vs_dittus_boelter_pct": -3.0, "f_cfd": f, "f_petukhov": 0.025,
                "err_f_pct": 1.0}
    s = {k["name"]: k for k in (mk("coarse", 1000, 160.0, 0.0270), mk("medium", 4000, 154.0, 0.0258),
                                 mk("fine", 16000, 152.5, 0.02545))}
    ms = mesh_study(s)
    assert ms["r21"] == pytest.approx(2.0)
    assert ms["Nu_fd"]["p"] == pytest.approx(2.0, rel=1e-6)
    assert "GCI(fine)" in markdown_report(s)


def test_stale_header_only_outlet_files_are_ignored(tmp_path):
    case, m, tb_th = _fake_run(tmp_path)
    stale = case / "postProcessing" / "fluid" / "outletTb" / "0" / "surfaceFieldValue_0.dat"
    stale.write_text("# Region type : patch outlet\n# Time weightedAverage(T)\n")
    s = analyse_case(case, tmp_path / "res")
    assert s["n_written_times"] == 3


def test_nested_region_layout_of_set_output(tmp_path):
    """Sets are written to postProcessing/<func>/<region>/<time>/<set>_<field>.xy."""
    case, m, _ = _fake_run(tmp_path)
    pp = case / "postProcessing"
    # move the pressure line to the nested layout used by sets
    src = pp / "fluid" / "sampleLine" / "3.5"
    dst = pp / "sampleLine" / "fluid" / "3.5"
    dst.parent.mkdir(parents=True)
    src.rename(dst)
    (pp / "fluid" / "sampleLine").rmdir()
    s = analyse_case(case, tmp_path / "res")
    assert s["f_cfd"] == pytest.approx(0.026, rel=1e-6)
