import json
import re
from pathlib import Path

import pytest

from chtpipe.casegen import FUNCTION_FILES, make_case
from chtpipe.config import CASE_MATRIX, CaseSpec


@pytest.mark.parametrize("mesh,re_,wall", CASE_MATRIX)
def test_every_case_generates_without_placeholders(tmp_path, mesh, re_, wall):
    case = make_case(CaseSpec(mesh, re_, wall), tmp_path)
    for f in case.rglob("*"):
        if f.is_file():
            assert "@@" not in f.read_text(errors="ignore"), f
    assert (case / "system" / "blockMeshDict").exists()
    for name in FUNCTION_FILES:
        assert (case / "system" / name).read_text().startswith("FoamFile")
        assert (case / "system" / "fluid" / name).exists()
    assert (case / "Allrun").stat().st_mode & 0o111


def test_expected_files_and_values(tmp_path):
    case = make_case(CaseSpec("medium", 25000, "copper"), tmp_path)
    meta = json.loads((case / "case_meta.json").read_text())
    assert meta["UB"] == pytest.approx(25000 * 8.9e-4 / (997 * 0.02))
    t_solid = (case / "0.orig" / "solid" / "T").read_text()
    assert f"gradient        uniform {meta['GRAD_OUTER']:.10g};" in t_solid
    assert "kappa       400;" in (case / "constant" / "solid" / "thermophysicalProperties").read_text()
    assert "numberOfSubdomains  4;" in (case / "system" / "decomposeParDict").read_text()
    for region, files in {"fluid": ["U", "p", "p_rgh", "T", "k", "omega", "nut", "alphat"], "solid": ["T", "p"]}.items():
        for f in files:
            assert (case / "0.orig" / region / f).exists()


def test_sampling_points_are_inside_the_wedge(tmp_path):
    case = make_case(CaseSpec("coarse", 25000, "steel"), tmp_path)
    meta = json.loads((case / "case_meta.json").read_text())
    txt = (case / "system" / "sampleProfiles").read_text()
    ends = [float(m) for m in re.findall(r"end\s+\(\S+ (\S+) 0\)", txt)]
    assert len(ends) == 3 and all(0 < y < meta["R_IN"] for y in ends)


def test_overwrite_protection(tmp_path):
    make_case(CaseSpec("coarse", 10000, "steel"), tmp_path)
    with pytest.raises(FileExistsError):
        make_case(CaseSpec("coarse", 10000, "steel"), tmp_path)
    make_case(CaseSpec("coarse", 10000, "steel"), tmp_path, overwrite=True)


def test_unknown_mesh_or_wall_rejected(tmp_path):
    with pytest.raises(ValueError):
        make_case(CaseSpec("huge", 1e4, "steel"), tmp_path)
    with pytest.raises(ValueError):
        make_case(CaseSpec("coarse", 1e4, "gold"), tmp_path)


def test_gravity_files_are_global_and_regional(tmp_path):
    """chtMultiRegionFoam (ESI v2512) reads constant/g from the case root."""
    case = make_case(CaseSpec("coarse", 25000, "steel"), tmp_path)
    assert (case / "constant" / "g").exists()
    assert (case / "constant" / "hRef").exists()


def test_turbulence_dictionaries_present(tmp_path):
    case = make_case(CaseSpec("coarse", 25000, "steel"), tmp_path)
    for name in ("turbulenceProperties", "momentumTransport"):
        assert "kOmegaSST" in (case / "constant" / "fluid" / name).read_text()


def test_no_catch_all_regex_in_boundary_fields(tmp_path):
    """OpenFOAM tries regex keys last-to-first, so a '.*' key would also capture the wedge patches."""
    case = make_case(CaseSpec("coarse", 25000, "steel"), tmp_path)
    for f in (case / "0.orig").rglob("*"):
        if f.is_file():
            assert '".*"' not in f.read_text(), f


@pytest.mark.parametrize("mesh,re_,wall", CASE_MATRIX)
def test_last_write_time_coincides_with_end_time(mesh, re_, wall):
    from chtpipe.config import derive
    p = derive(CaseSpec(mesh, re_, wall))
    assert p["END_TIME"] == pytest.approx(20 * p["WRITE_INTERVAL"], rel=1e-9)
    assert p["END_TIME"] >= 8 * p["L_TOTAL"] / p["UB"] - 1e-9
