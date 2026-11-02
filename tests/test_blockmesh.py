"""Geometric consistency of the generated blockMeshDict (no OpenFOAM needed)."""
import itertools
import math

import numpy as np
import pytest

from chtpipe import blockmesh
from chtpipe.config import CaseSpec, MESHES, derive, first_cell_height


@pytest.fixture(params=sorted(MESHES))
def geo(request):
    p = derive(CaseSpec(request.param, 25000, "steel"))
    return p, blockmesh.build(p)


GAUSS = [(1 - 1 / math.sqrt(3)) / 2, (1 + 1 / math.sqrt(3)) / 2]


def _map(corner, u, v, w):
    """Trilinear map of the 8 hex corners (OpenFOAM ordering)."""
    c = np.asarray(corner, dtype=float)
    return ((1 - u) * (1 - v) * (1 - w) * c[0] + u * (1 - v) * (1 - w) * c[1]
            + u * v * (1 - w) * c[2] + (1 - u) * v * (1 - w) * c[3]
            + (1 - u) * (1 - v) * w * c[4] + u * (1 - v) * w * c[5]
            + u * v * w * c[6] + (1 - u) * v * w * c[7])


def block_volume(corner):
    """Volume and minimum Jacobian over the 2x2x2 Gauss points (exact for a trilinear map)."""
    vol, jmin, h = 0.0, math.inf, 1e-6
    for u, v, w in itertools.product(GAUSS, repeat=3):
        du = (_map(corner, u + h, v, w) - _map(corner, u - h, v, w)) / (2 * h)
        dv = (_map(corner, u, v + h, w) - _map(corner, u, v - h, w)) / (2 * h)
        dw = (_map(corner, u, v, w + h) - _map(corner, u, v, w - h)) / (2 * h)
        j = float(np.dot(du, np.cross(dv, dw)))
        jmin = min(jmin, j)
        vol += j / 8.0
    return vol, jmin


def test_blocks_are_right_handed_and_volumes_match_wedge(geo):
    p, g = geo
    th = math.radians(p["HALF_ANGLE_DEG"])
    fluid = solid = 0.0
    for b in g["blocks"]:
        corners = [g["vertices"][i] for i in b["hex"]]
        vol, jmin = block_volume(corners)
        assert jmin > 0, f"{b['name']} has a non-positive Jacobian"
        if b["zone"] == "fluid":
            fluid += vol
        else:
            solid += vol
    sc = math.sin(th) * math.cos(th)
    assert fluid == pytest.approx(p["L_TOTAL"] * p["R_IN"] ** 2 * sc, rel=1e-6)
    assert solid == pytest.approx(p["L_TOTAL"] * (p["R_OUT"] ** 2 - p["R_IN"] ** 2) * sc, rel=1e-6)


def _dedup(face):
    out = []
    for v in face:
        if not out or out[-1] != v:
            out.append(v)
    if len(out) > 1 and out[0] == out[-1]:
        out.pop()
    return out


def _normal(face, verts):
    pts = np.array([verts[i] for i in _dedup(face)])
    n = np.zeros(3)
    for k in range(len(pts)):
        n += np.cross(pts[k], pts[(k + 1) % len(pts)])
    return n / 2.0


def test_patch_faces_point_outward(geo):
    p, g = geo
    verts = g["vertices"]
    for name, _, faces in g["patches"]:
        for f in faces:
            n = _normal(f, verts)
            if name == "axis":
                assert np.linalg.norm(n) == 0     # collapsed: zero area
                continue
            assert np.linalg.norm(n) > 0, f"degenerate face in {name}: {f}"
            if name in ("inlet", "solid_inlet"):
                assert n[0] < 0
            elif name in ("outlet", "solid_outlet"):
                assert n[0] > 0
            elif name.startswith("outer"):
                assert n[1] > 0
            elif name == "front":
                assert n[2] > 0
            elif name == "back":
                assert n[2] < 0


def test_collapsed_faces_keep_repeated_labels(geo):
    """Axis-side faces keep four labels with a repeated pair (what blockMesh matches)."""
    p, g = geo
    by_name = {n: f for n, _, f in g["patches"]}
    for name in ("inlet", "outlet"):
        f = by_name[name][0]
        assert len(f) == 4 and len(set(f)) == 3
    for f in by_name["axis"]:
        assert len(f) == 4 and len(set(f)) == 2
    for name, _, faces in g["patches"]:
        if name not in ("inlet", "outlet", "axis"):
            assert all(len(set(f)) == 4 for f in faces), name


def test_interface_faces_match_between_fluid_and_solid(geo):
    _, g = geo
    from chtpipe.blockmesh import HEX_FACES
    def face_set(block, name):
        return frozenset(g["blocks"][block]["hex"][k] for k in HEX_FACES[name])
    for i in range(2):
        assert face_set(i, "ymax") == face_set(2 + i, "ymin")
    # axial junctions
    assert face_set(0, "xmax") == face_set(1, "xmin")
    assert face_set(2, "xmax") == face_set(3, "xmin")


def test_cell_counts_consistent_across_shared_faces(geo):
    _, g = geo
    b = g["blocks"]
    assert b[0]["cells"][0] == b[2]["cells"][0] and b[1]["cells"][0] == b[3]["cells"][0]
    assert b[0]["cells"][1] == b[1]["cells"][1]
    assert b[2]["cells"][1] == b[3]["cells"][1]


def test_every_boundary_face_is_listed_once(geo):
    """Count block boundary faces vs listed patch faces (axis collapsed faces excluded)."""
    _, g = geo
    from chtpipe.blockmesh import HEX_FACES
    count = {}
    for bi, b in enumerate(g["blocks"]):
        for name, loc in HEX_FACES.items():
            labels = frozenset(b["hex"][k] for k in loc)
            count[labels] = count.get(labels, 0) + 1
    boundary = {k for k, v in count.items() if v == 1}
    listed = [frozenset(f) for _, _, faces in g["patches"] for f in faces]
    assert len(listed) == len(set(listed))
    assert set(listed) == boundary


def test_first_cell_height_and_total_length():
    for n in (30, 60, 120):
        h = first_cell_height(0.01, n, 30.0)
        q = 30.0 ** (1.0 / (n - 1))
        total = sum(h * q ** k for k in range(n))
        assert total == pytest.approx(0.01, rel=1e-9)
    assert first_cell_height(0.01, 10, 1.0) == pytest.approx(0.001)


def test_systematic_refinement_halves_first_cell():
    h = [first_cell_height(0.01, MESHES[m]["nr_f"], 30.0) for m in ("coarse", "medium", "fine")]
    assert h[0] / h[1] == pytest.approx(2.0, rel=0.1)
    assert h[1] / h[2] == pytest.approx(2.0, rel=0.1)
