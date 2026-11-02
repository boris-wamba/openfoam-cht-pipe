#!/usr/bin/env python3
"""Figures of the README: 3D view of the geometry and a schematic of the expected fields.

Both are illustrations. The geometry view is drawn from the case dimensions; the field schematic is
reconstructed from the computed wall temperature of the coarse case (Nu_x.csv) and from the shape of the
radial profile of the same case, not rendered from the solution fields.

Usage: python scripts/make_readme_figures.py [results/coarse_Re25000_steel/Nu_x.csv] [out_dir]
"""
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

CSV = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("results/coarse_Re25000_steel/Nu_x.csv")
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("docs/figures")
OUT.mkdir(parents=True, exist_ok=True)

RI, RO, L, LE = 10.0, 12.0, 1.5, 0.6          # mm, mm, m, m
UB, N = 1.116, 8.0
DT_SOLID = 3.42

# ============================================================== 3D geometry view
def geometry_3d(path):
    EX = 8.0                                   # radial exaggeration of the drawing
    ri, ro = RI * 1e-3 * EX, RO * 1e-3 * EX
    nx_e, nx_h, nth = 12, 20, 60
    xs = np.concatenate([np.linspace(0, LE, nx_e + 1), np.linspace(LE, L, nx_h + 1)[1:]])
    th = np.radians(np.linspace(180, 450, nth + 1))      # open quadrant: theta in [90, 180] deg

    polys, cols = [], []

    def P(x, r, t):
        return (x, r * math.cos(t), r * math.sin(t))

    steel_u, steel_h, interface = "#b7c0ca", "#e4673d", "#7db9e8"
    for i in range(len(xs) - 1):
        heated = xs[i] >= LE - 1e-9
        for j in range(nth):
            # outer surface
            polys.append([P(xs[i], ro, th[j]), P(xs[i + 1], ro, th[j]), P(xs[i + 1], ro, th[j + 1]), P(xs[i], ro, th[j + 1])])
            cols.append(steel_h if heated else steel_u)
            # interface (inner surface of the steel = pipe wall seen by the water)
            polys.append([P(xs[i], ri, th[j]), P(xs[i + 1], ri, th[j]), P(xs[i + 1], ri, th[j + 1]), P(xs[i], ri, th[j + 1])])
            cols.append(interface)
    # end rings (steel) and discs (water) at x = 0 and x = L
    for xe, cs in ((0.0, "#8d98a5"), (L, "#8d98a5")):
        for j in range(nth):
            polys.append([P(xe, ri, th[j]), P(xe, ro, th[j]), P(xe, ro, th[j + 1]), P(xe, ri, th[j + 1])])
            cols.append(cs)
            polys.append([(xe, 0, 0), P(xe, ri, th[j]), P(xe, ri, th[j + 1])])
            cols.append("#d4e9f8")
    # cut planes at theta = 180 deg and theta = 90 deg (steel in dark grey, water in light blue)
    for tcut in (math.pi, math.pi / 2):
        for i in range(len(xs) - 1):
            polys.append([P(xs[i], 0.0, tcut), P(xs[i + 1], 0.0, tcut), P(xs[i + 1], ri, tcut), P(xs[i], ri, tcut)])
            cols.append("#d4e9f8")
            polys.append([P(xs[i], ri, tcut), P(xs[i + 1], ri, tcut), P(xs[i + 1], ro, tcut), P(xs[i], ro, tcut)])
            cols.append("#5f6c79")
    # the simulated sector (angle exaggerated) inside the open quadrant
    tc, tw = math.radians(135), math.radians(10)
    for i in range(len(xs) - 1):
        for (ra, rb, c) in ((0.0, ri, "#1f78c8"), (ri, ro, "#b2481e")):
            polys.append([P(xs[i], ra, tc - tw), P(xs[i + 1], ra, tc - tw), P(xs[i + 1], rb, tc - tw), P(xs[i], rb, tc - tw)])
            cols.append(c)
            polys.append([P(xs[i], ra, tc + tw), P(xs[i + 1], ra, tc + tw), P(xs[i + 1], rb, tc + tw), P(xs[i], rb, tc + tw)])
            cols.append(c)
    for xe in (0.0, L):
        polys.append([(xe, 0, 0), P(xe, ri, tc - tw), P(xe, ri, tc + tw)])
        cols.append("#1f78c8")
        polys.append([P(xe, ri, tc - tw), P(xe, ro, tc - tw), P(xe, ro, tc + tw), P(xe, ri, tc + tw)])
        cols.append("#b2481e")

    fig = plt.figure(figsize=(13, 6.2))
    ax = fig.add_axes([-0.05, -0.12, 1.1, 1.25], projection="3d")
    pc = Poly3DCollection(polys, facecolors=cols, edgecolors="none", linewidths=0)
    ax.add_collection3d(pc)
    ax.set_xlim(-0.15, L + 0.15)
    ax.set_ylim(-ro * 1.4, ro * 1.4)
    ax.set_zlim(-ro * 1.4, ro * 1.4)
    ax.set_box_aspect((1.8, 0.28, 0.28), zoom=1.25)
    ax.view_init(elev=26, azim=-62)
    ax.set_axis_off()

    def proj(p):
        x2, y2, _ = proj3d.proj_transform(p[0], p[1], p[2], ax.get_proj())
        return x2, y2

    def note(text, target, offset, color="k", ha="center"):
        x2, y2 = proj(target)
        ax.annotate(text, xy=(x2, y2), xytext=offset, textcoords="offset points", color=color, ha=ha, fontsize=9,
                    arrowprops=dict(arrowstyle="->", color=color, lw=1.0), annotation_clip=False)

    note("Inlet: water, T = 300 K,\nuniform velocity", P(0.0, 0.0, 0), (-70, 60), "#1f4e8c")
    note("Outlet: fixed pressure", P(L, 0.0, 0), (55, 55), "#1f4e8c")
    note("Unheated entry\n30 D = 0.6 m", P(0.3, ro, math.radians(300)), (-20, -105), "#4d5964")
    note("Heated length 45 D = 0.9 m\n(imposed heat flux on the outer wall)", P(1.05, ro, math.radians(300)), (55, -80), "#b2481e")
    note("Steel wall, 2 mm\n(second region, solid)", P(0.9, 0.5 * (ri + ro), math.radians(90)), (0, 85), "#4d5964")
    note("Water, D = 20 mm\n(first region, fluid)", P(0.45, 0.4 * ri, math.radians(180)), (-150, -25), "#1f4e8c")
    note("Simulated sector: 5° wedge, 1 cell thick\n(angle exaggerated in this drawing)", P(1.2, 0.5 * ri, tc + tw), (60, 100), "#0b3d73")
    fig.text(0.5, 0.04, "Radial dimensions are exaggerated ×8 in this drawing (the real pipe is 75 times longer than wide).",
             ha="center", fontsize=8.5, style="italic")
    fig.savefig(path, dpi=170, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    from PIL import Image, ImageChops
    im = Image.open(path).convert("RGB")
    bg = Image.new("RGB", im.size, (255, 255, 255))
    box = ImageChops.difference(im, bg).getbbox()
    if box:
        pad = 20
        im.crop((max(box[0] - pad, 0), max(box[1] - pad, 0), min(box[2] + pad, im.size[0]), min(box[3] + pad, im.size[1]))).save(path)


# ============================================================== field schematic
def field_schematic(csv, path):
    d = np.genfromtxt(csv, delimiter=",", names=True)
    X, TW, TB = d["x_m"], d["Tw_K"], d["Tb_K"]
    eta_p = np.array([0, .4, .6, .8, .9, .95, .98, .99, 1.0])
    f_p = np.array([0, .014, .04, .082, .142, .245, .47, .73, 1.0])

    nx, nr = 700, 480
    x = np.linspace(0, L, nx)
    r = np.linspace(0, RO, nr)
    T = np.zeros((nr, nx))
    tw = np.interp(x, X, TW); tb = np.interp(x, X, TB)
    for i, xi in enumerate(x):
        heated = xi > LE
        g = np.clip(((xi - LE) / 0.4) ** 0.5, 0.15, 1.0) if heated else 1.0
        tc = tb[i] - (0.38 * g if heated else 0.0)
        for j, rj in enumerate(r):
            if rj <= RI:
                eta = rj / RI
                s = (1 - eta) / g
                f = np.interp(1 - s, eta_p, f_p) if s < 1 else 0.0
                T[j, i] = tc + (tw[i] - tc) * f
            else:
                T[j, i] = tw[i] + DT_SOLID * (rj - RI) / (RO - RI) * (1.0 if heated else 0.0)
    xu = np.linspace(0, L, nx); ru = np.linspace(0, RI, nr); eta = ru / RI
    U = np.zeros((nr, nx))
    trapz = getattr(np, "trapezoid", None) or np.trapz
    for i, xi in enumerate(xu):
        dl = min(1.0, 0.02 + xi / 0.55)
        prof = np.where(1 - eta < dl, ((1 - eta) / dl).clip(1e-6) ** (1 / N), 1.0)
        U[:, i] = prof * UB / trapz(prof * 2 * eta, eta)

    fig = plt.figure(figsize=(13, 13.2))
    gs = GridSpec(4, 3, height_ratios=[1.35, 1.35, 1.6, 0.05], hspace=0.62, wspace=0.3,
                  left=0.06, right=0.97, top=0.895, bottom=0.04)
    fig.suptitle("Conjugate heat transfer in a heated pipe: typical fields (schematic, coarse mesh, Re = 25,000, steel wall)",
                 fontsize=13.5, fontweight="bold", y=0.992)

    ax = fig.add_subplot(gs[0, :])
    im = ax.imshow(T, origin="lower", extent=[0, L, 0, RO], aspect="auto", cmap="inferno", vmin=300, vmax=311)
    ax.axhline(RI, color="w", lw=0.8, ls="--"); ax.axvline(LE, color="c", lw=1.0, ls=":")
    ax.set_xlim(0, L); ax.set_ylim(0, RO)
    ax.set_title("A. Temperature T (water + steel wall) at the final time, 15.4 s; radial axis stretched ×15; scale 300 → 311 K",
                 fontsize=10.5, loc="left", pad=44)
    ax.set_xlabel("x [m]"); ax.set_ylabel("r [mm]")
    ax.text(0.30, 3.6, "water: T = 300 K\nunheated entry: nothing to see", color="w", ha="center", va="center", fontsize=9)
    ax.text(1.05, 3.2, "water core: 300.3 → 301 K (nearly uniform)", color="w", ha="center", fontsize=9)
    ax.text(0.30, 11.0, "steel wall (r = 10 to 12 mm)", color="w", ha="center", va="center", fontsize=8.5)
    ax.text(LE + 0.01, 0.6, "start of heating, x = 0.6 m", color="c", fontsize=8.5, va="bottom")
    ax.text(1.0, 8.6, "thin hot layer along the wall\n(307 K at the outlet, < 1 mm)", color="w", ha="center", fontsize=9)
    ax.annotate("symmetry axis (r = 0)", xy=(1.35, 0.0), xytext=(1.30, 1.8), color="w", fontsize=8.5,
                arrowprops=dict(arrowstyle="->", color="w"))
    ax.annotate("", xy=(0.14, 6.5), xytext=(0.0, 6.5), arrowprops=dict(arrowstyle="-|>", color="w", lw=1.6))
    ax.text(0.01, 8.0, "inlet, x = 0\nT = 300 K, U = 1.116 m/s uniform", color="w", ha="left", va="bottom", fontsize=8.5)
    ax.text(L + 0.01, 6.0, "outlet\nTb = 301.16 K\nTw = 307.6 K", ha="left", va="center", fontsize=8.5)
    for xx in np.linspace(0.7, 1.4, 8):
        ax.annotate("", xy=(xx, RO), xytext=(xx, RO + 1.6), annotation_clip=False,
                    arrowprops=dict(arrowstyle="-|>", color="firebrick", lw=1.3))
    ax.text(1.05, RO + 2.0, "imposed flux 25 kW/m² on the outer wall (= 30 kW/m² on the water side), heated length 0.6 → 1.5 m",
            color="firebrick", ha="center", fontsize=9)
    cb = fig.colorbar(im, ax=ax, pad=0.12, fraction=0.025); cb.set_label("T [K]")

    ax = fig.add_subplot(gs[1, :])
    im2 = ax.imshow(U, origin="lower", extent=[0, L, 0, RI], aspect="auto", cmap="viridis", vmin=0, vmax=1.4)
    ax.set_title("B. Axial velocity U_x (water only), scale 0 → 1.4 m/s", fontsize=10.5, loc="left")
    ax.set_xlabel("x [m]"); ax.set_ylabel("r [mm]")
    ax.text(0.04, 5, "flat\nprofile", color="w", ha="center", va="center", fontsize=9)
    ax.text(0.65, 6.0, "the boundary layer grows: the profile develops over ≈ 25 to 30 D (0.5 – 0.6 m)", color="k", ha="center", fontsize=9)
    ax.text(1.1, 2.3, "developed flow: centreline U ≈ 1.33 m/s (1.19 × Ub)", color="k", ha="center", fontsize=9.5)
    ax.text(1.3, 9.2, "velocity → 0 at the wall", color="w", ha="center", fontsize=9)
    cb = fig.colorbar(im2, ax=ax, pad=0.12, fraction=0.025); cb.set_label("U_x [m/s]")

    ax = fig.add_subplot(gs[2, 0])
    xz = (x > 0.55) & (x < 0.8); rz = (r > 8.5)
    ax.imshow(T[np.ix_(rz, xz)], origin="lower", extent=[0.55, 0.8, r[rz][0], RO], aspect="auto",
              cmap="inferno", vmin=300, vmax=308)
    ax.axhline(RI, color="w", lw=0.8, ls="--"); ax.axvline(LE, color="c", lw=1, ls=":")
    ax.set_title("C. Zoom on the start of heating\nx = 0.55–0.8 m, r = 8.5–12 mm\n(colour scale 300–308 K)", fontsize=9.5, loc="left")
    ax.set_xlabel("x [m]"); ax.set_ylabel("r [mm]")
    ax.text(0.72, 11.0, "steel: hotter than\nthe water by ≈ 3.4 K", color="k", ha="center", fontsize=8)
    ax.text(0.68, 9.0, "hot layer: 300 → 306 K\nwithin ≈ 5 D (0.1 m)", color="w", ha="center", fontsize=8)
    ax.text(0.603, 8.65, "x = 0.6 m", color="c", ha="left", fontsize=8)
    ax.text(0.5525, 11.2, "before:\nthe steel warms\nslightly (0.02 K)", color="w", fontsize=7, va="center")

    ax = fig.add_subplot(gs[2, 1])
    q = 30.0 ** (1 / 29)
    h1 = 0.01 * (q - 1) / (q ** 30 - 1) * 1e6
    hf = [h1 * q ** k for k in range(8)]
    cw, ncol = 120.0, 4
    edges = np.concatenate([[0], np.cumsum(hf)])
    for k in range(len(edges)):
        ax.plot([0, ncol * cw], [-edges[k]] * 2, color="k", lw=0.8)
    for a_ in range(ncol + 1):
        ax.plot([a_ * cw] * 2, [-edges[-1], 0], color="k", lw=0.8)
    for k in range(4):
        ax.plot([0, ncol * cw], [250 * k] * 2, color="saddlebrown", lw=0.8)
    for a_ in range(ncol + 1):
        ax.plot([a_ * cw] * 2, [0, 750], color="saddlebrown", lw=0.8)
    ax.axhline(0, color="b", lw=2)
    ax.set_xlim(-190, ncol * cw + 540); ax.set_ylim(-edges[-1] - 120, 800)
    ax.set_xticks([])
    ax.set_ylabel("distance to the wall [µm] (true scale in r)")
    ax.set_title("D. Mesh next to the wall\n(strong zoom)", fontsize=9.5, loc="left")
    ax.text(ncol * cw + 20, 25, "fluid / steel interface\n(patch fluid_to_solid)", color="b", fontsize=8, va="bottom")
    ax.text(ncol * cw + 20, -60, "1st fluid cell:\n38 µm (y+ ≈ 1.3)", fontsize=8, va="top")
    ax.text(ncol * cw + 20, -edges[-1] * 0.72, "fluid cells: +12 %\nper cell towards the axis", fontsize=8, va="center")
    ax.text(ncol * cw + 20, 375, "steel: 250 µm cells\n(8 in total)", fontsize=8, color="saddlebrown", va="center")
    ax.text(-170, -edges[-1] - 95, "width in x not to scale (4 mm cells)", fontsize=7.5, style="italic")

    ax = fig.add_subplot(gs[2, 2])
    e = np.linspace(0, 0.9999, 400)
    tcore, tfirst = 300.27, 306.1
    ax.plot(e * RI, tcore + (tfirst - tcore) * np.interp(e, eta_p, f_p), color="C1", lw=2)
    ax.axhline(300.65, color="gray", ls="--", lw=1)
    ax.text(3.2, 299.95, "bulk T (energy balance) 300.65 K ---", fontsize=8, color="gray")
    ax.set_title("E. Radial profile of T at x = 1.1 m\n(from the axis to the first wall cell)", fontsize=9.5, loc="left")
    ax.set_xlabel("r [mm]"); ax.set_ylabel("T [K]")
    ax.annotate("core ≈ 300.3 K", xy=(1.0, 300.3), xytext=(1.2, 302.5), fontsize=8.5, arrowprops=dict(arrowstyle="->"))
    ax.annotate("steep rise: 5 K in\nless than 0.5 mm", xy=(9.85, 304.2), xytext=(4.3, 304.4), fontsize=8.5, arrowprops=dict(arrowstyle="->"))
    ax.annotate("first cell ≈ 306.1 K\n(wall value: 307.1 K)", xy=(9.99, 306.0), xytext=(5.0, 306.3), fontsize=8.5, arrowprops=dict(arrowstyle="->"))
    ax.set_ylim(299.8, 307.4)

    fig.text(0.5, 0.012, "Schematic reconstructed from the computed wall temperature and radial profile of the coarse case; "
             "it is not a rendering of the solution fields.", ha="center", fontsize=8.5, style="italic")
    fig.savefig(path, dpi=140, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    geometry_3d(OUT / "geometry_3d.png")
    field_schematic(CSV, OUT / "fields_schematic.png")
    print("written to", OUT)
