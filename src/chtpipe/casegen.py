"""Generate a ready-to-run OpenFOAM case directory from the template."""
from __future__ import annotations

import json
import math
import re
import shutil
from pathlib import Path

from . import blockmesh
from .config import CaseSpec, derive

TEMPLATE = Path(__file__).resolve().parents[2] / "template"
PLACEHOLDER = re.compile(r"@@([A-Z0-9_]+)@@")
FUNCTION_FILES = ("outletTb", "sampleWall", "sampleProfiles", "sampleU", "sampleLine")
HEADER = "FoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       dictionary;\n    object      %s;\n}\n\n"


def _fmt(v) -> str:
    if isinstance(v, float):
        return f"{v:.10g}"
    return str(v)


def render(text: str, p: dict, where: str = "") -> str:
    """Replace every ``@@KEY@@`` by ``p[KEY]``; fail loudly on unknown keys."""
    def sub(m):
        key = m.group(1)
        if key not in p:
            raise KeyError(f"unknown placeholder @@{key}@@ in {where}")
        return _fmt(p[key])
    return PLACEHOLDER.sub(sub, text)


def _set_block(name: str, start, end, npts: int) -> str:
    return (
        f"    {name}\n    {{\n        type    uniform;\n        axis    xyz;\n"
        f"        start   ({_fmt(start[0])} {_fmt(start[1])} {_fmt(start[2])});\n"
        f"        end     ({_fmt(end[0])} {_fmt(end[1])} {_fmt(end[2])});\n"
        f"        nPoints {npts};\n    }}"
    )


def sampling_sets(p: dict) -> dict:
    """Text blocks for the sets used by the sampling dictionaries."""
    y_wall = p["R_IN"] * math.cos(math.radians(p["HALF_ANGLE_DEG"]))
    y0, y1 = 1.0e-6, y_wall - 1.0e-7
    profs = []
    for d in p["PROFILE_STATIONS_D"]:
        x = p["L_ENTRY"] + d * p["D"]
        profs.append(_set_block(f"prof_h{d}", (x, y0, 0.0), (x, y1, 0.0), 400))
    xu = p["L_ENTRY"] + p["U_STATION_D"] * p["D"]
    uset = _set_block(f"uprof_h{p['U_STATION_D']}", (xu, y0, 0.0), (xu, y1, 0.0), 1200)
    line = _set_block("mid", (1.0e-4, 0.5 * y_wall, 0.0), (p["L_TOTAL"] - 1.0e-4, 0.5 * y_wall, 0.0), 750)
    return {"PROFILE_SETS": "\n".join(profs), "U_SET": uset, "LINE_SET": line, "Y_WALL": y_wall}


def make_case(spec: CaseSpec, out_root: str | Path = "runs", overwrite: bool = False) -> Path:
    p = derive(spec)
    sets = sampling_sets(p)
    p.update({k: v for k, v in sets.items() if k != "Y_WALL"})
    p["Y_WALL"] = sets["Y_WALL"]

    out = Path(out_root) / p["name"]
    if out.exists():
        if not overwrite:
            raise FileExistsError(f"{out} exists (use --overwrite)")
        shutil.rmtree(out)
    shutil.copytree(TEMPLATE, out)

    for f in sorted(out.rglob("*")):
        if f.is_file():
            text = f.read_text()
            if PLACEHOLDER.search(text):
                f.write_text(render(text, p, str(f.relative_to(out))))

    # function-object dictionaries: add a header and make them visible for -region fluid
    for name in FUNCTION_FILES:
        f = out / "system" / name
        body = f.read_text()
        if not body.startswith("FoamFile"):
            body = HEADER % name + body
            f.write_text(body)
        shutil.copy(f, out / "system" / "fluid" / name)

    blockmesh.write(p, out / "system" / "blockMeshDict")

    meta = {k: v for k, v in p.items() if k not in ("PROFILE_SETS", "U_SET", "LINE_SET")}
    (out / "case_meta.json").write_text(json.dumps(meta, indent=2))

    leftovers = [f for f in out.rglob("*") if f.is_file() and PLACEHOLDER.search(f.read_text(errors="ignore"))]
    if leftovers:
        raise RuntimeError(f"unresolved placeholders in {leftovers}")
    return out
