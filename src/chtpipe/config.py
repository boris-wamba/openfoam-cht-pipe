"""Case definition: geometry, materials, mesh levels and derived quantities."""
from __future__ import annotations

import math
from dataclasses import dataclass, asdict

from .correlations import gnielinski

# --- geometry (SI) ---------------------------------------------------------
D = 0.020                 # inner diameter [m]
R_IN = D / 2.0            # inner radius [m]
T_WALL = 0.002            # wall thickness [m]
R_OUT = R_IN + T_WALL     # outer radius [m]
L_ENTRY = 30.0 * D        # unheated entry length [m]
L_HEATED = 45.0 * D       # heated length [m]
L_TOTAL = L_ENTRY + L_HEATED
HALF_ANGLE_DEG = 2.5      # half opening angle of the axisymmetric wedge

# --- fluid: water, constant properties (about 25 degC) -----------------------
RHO_F = 997.0
CP_F = 4180.0
MU_F = 8.9e-4
K_F = 0.607
PR_F = CP_F * MU_F / K_F
PRT = 0.85
T_IN = 300.0

# --- thermal load -------------------------------------------------------------
Q_INNER = 30000.0         # heat flux through the inner wall [W/m2]

# --- turbulence at the inlet -------------------------------------------------
INTENSITY = 0.05
MIX_LENGTH = 0.07 * D
CMU = 0.09

# --- wall materials (constant properties, typical handbook values) -----------
WALLS = {
    "steel": {"kappa": 16.0, "rho": 8000.0, "cp": 500.0},
    "copper": {"kappa": 400.0, "rho": 8960.0, "cp": 385.0},
}

# --- mesh levels -------------------------------------------------------------
# nx_e / nx_h: axial cells in the unheated / heated part; nr_f / nr_s: radial cells
# in the fluid / in the wall. Each level doubles every count (refinement ratio 2).
MESHES = {
    "coarse": {"nx_e": 150, "nx_h": 225, "nr_f": 30, "nr_s": 8, "np": 1},
    "medium": {"nx_e": 300, "nx_h": 450, "nr_f": 60, "nr_s": 16, "np": 4},
    "fine": {"nx_e": 600, "nx_h": 900, "nr_f": 120, "nr_s": 32, "np": 4},
}
EXPANSION_TOTAL = 30.0    # ratio of the largest to the smallest radial fluid cell

# stations (distance from the start of heating, in diameters)
PROFILE_STATIONS_D = (10, 25, 40)
U_STATION_D = 30
FD_WINDOW_D = (25.0, 43.0)   # fully developed window for the averages (from heating start)

# case matrix: (mesh, Re, wall)
CASE_MATRIX = [
    ("coarse", 25000, "steel"),
    ("medium", 25000, "steel"),
    ("fine", 25000, "steel"),
    ("medium", 10000, "steel"),
    ("medium", 50000, "steel"),
    ("medium", 25000, "copper"),
]


@dataclass
class CaseSpec:
    mesh: str
    re: float
    wall: str
    np: int | None = None

    @property
    def name(self) -> str:
        return f"{self.mesh}_Re{int(self.re)}_{self.wall}"


def first_cell_height(length: float, n: int, er: float) -> float:
    """Height of the smallest cell for a geometric distribution of n cells.

    ``er`` is the ratio of the largest to the smallest cell.
    """
    if n < 1:
        raise ValueError("n must be >= 1")
    if n == 1 or abs(er - 1.0) < 1e-12:
        return length / n
    q = er ** (1.0 / (n - 1))
    return length * (q - 1.0) / (q ** n - 1.0)


def derive(spec: CaseSpec) -> dict:
    """All numbers needed to write a case, as a flat dictionary."""
    if spec.mesh not in MESHES:
        raise ValueError(f"unknown mesh level {spec.mesh!r}")
    if spec.wall not in WALLS:
        raise ValueError(f"unknown wall material {spec.wall!r}")
    m = MESHES[spec.mesh]
    w = WALLS[spec.wall]
    re = float(spec.re)

    ub = re * MU_F / (RHO_F * D)
    mdot = RHO_F * ub * math.pi * R_IN ** 2
    q_outer = Q_INNER * R_IN / R_OUT
    grad_outer = q_outer / w["kappa"]      # > 0: wall hotter than the solid interior, heat enters
    k_in = 1.5 * (INTENSITY * ub) ** 2
    omega_in = math.sqrt(k_in) / (CMU ** 0.25 * MIX_LENGTH)

    transit = L_TOTAL / ub
    h_est = gnielinski(re, PR_F) * K_F / D
    tau_wall = w["rho"] * w["cp"] * T_WALL / h_est
    end_req = max(8.0 * transit, 10.0 * tau_wall)
    # write 20 times, the last one exactly at endTime (writeInterval rounded UP to 3 significant digits)
    wi = end_req / 20.0
    step = 10.0 ** (math.floor(math.log10(wi)) - 2)
    write_interval = float(f"{math.ceil(round(wi / step, 6)) * step:.6g}")
    end_time = float(f"{20.0 * write_interval:.6g}")

    dr = first_cell_height(R_IN, m["nr_f"], EXPANSION_TOTAL)
    u_tau_est = ub * math.sqrt(0.0791 * re ** -0.25 / 2.0)   # Blasius (Fanning f = 0.0791 Re^-1/4)
    yplus_est = 0.5 * dr * u_tau_est / (MU_F / RHO_F)
    nproc = spec.np if spec.np is not None else m["np"]

    p = dict(asdict(spec))
    p.update(
        name=spec.name, np=nproc, NP=nproc,
        D=D, R_IN=R_IN, R_OUT=R_OUT, T_WALL=T_WALL, L_ENTRY=L_ENTRY, L_HEATED=L_HEATED,
        L_TOTAL=L_TOTAL, HALF_ANGLE_DEG=HALF_ANGLE_DEG,
        RHO_F=RHO_F, CP_F=CP_F, MU_F=MU_F, K_F=K_F, PR_F=PR_F, PRT=PRT, T_IN=T_IN,
        Q_INNER=Q_INNER, Q_OUTER=q_outer, GRAD_OUTER=grad_outer,
        KAPPA_S=w["kappa"], RHO_S=w["rho"], CP_S=w["cp"], WALL_NAME=spec.wall,
        UB=ub, MDOT=mdot, K_IN=k_in, OMEGA_IN=omega_in, L_MIX=MIX_LENGTH,
        END_TIME=end_time, WRITE_INTERVAL=write_interval,
        MAX_CO=2.0, MAX_DT=0.05,
        EXPANSION_TOTAL=EXPANSION_TOTAL, FIRST_CELL=dr, YPLUS_EST=yplus_est,
        PROFILE_STATIONS_D=list(PROFILE_STATIONS_D), U_STATION_D=U_STATION_D,
        FD_WINDOW_D=list(FD_WINDOW_D),
        **{k.upper(): v for k, v in m.items() if k != "np"},
    )
    return p
