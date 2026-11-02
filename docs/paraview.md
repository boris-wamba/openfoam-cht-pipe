# ParaView checklist

Open `runs/<case>/case.foam`. In the Properties panel, under *Mesh Regions*, tick
`fluid/internalMesh` and `solid/internalMesh` (and the patches you need), then Apply.
The last time step is the converged solution.

Wedge geometry: the domain is a 5 degree sector, 1.5 m long and 10 mm in radius. For readable
images, apply *Transform* to each region with scale (1, 15, 1) for the whole pipe, or crop
to a region of interest (Clip / Extract Subset).

## Figures to produce

1. **Mesh near the wall at the start of heating** (Surface With Edges, both regions, zoom
   at x = 0.6 m, y between 9 and 12 mm). Shows the first cells and the fluid/solid interface.
2. **Temperature field**, fluid and solid together, range 300 K to the maximum. One image on the
   whole pipe (scaled) and one zoom on x = 0.55-0.8 m.
3. **Axial velocity** `U_X` in the fluid, whole pipe (scaled) to show the development.
4. **Turbulent viscosity ratio**: Calculator `nut/(mu/rho)` with mu/rho = 8.93e-7 m2/s.
5. **Wall shear stress**: after `post_case.sh`, the field `wallShearStress` exists on the interface
   patch (optional step). Use it as an independent check of the friction velocity
   `u_tau = sqrt(|tau_w| / rho)` that the analysis derives from the pressure gradient.
6. **Wall y+**: the first-cell y+ is reported in `summary.json` (`yplus_first_cell`); the OpenFOAM
   `yPlus` function returns zeros with this setup and is not used.
7. **Wall temperature, steel versus copper**, zoom on the first 10 D of heating: open the two
   cases side by side with the same colour range.
8. **Streamlines or vectors** are not informative here (straight flow): skip.

## Quantitative data

The numbers used in the report (Nu, f, GCI, profiles) come from the OpenFOAM
`postProcess` outputs and `scripts/analyse_case.py`, not from screenshots. Use ParaView for
qualitative figures and for checks, for instance a *Plot Over Line* of `T` between
(1.2, 0, 0) and (1.2, 0.01, 0), compared with `results/<case>/T_profiles.png`.

## Export

File > Save Screenshot at 1920 x 1080 or larger, PNG, white background, colour bar with units.
Name files `<case>_<quantity>.png`.
