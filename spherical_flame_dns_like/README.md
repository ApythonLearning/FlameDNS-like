# Low-Mach Axisymmetric Lean H2/O2 Spherical Flame DNS-like Framework

This project is a runnable Python framework for low-Mach, two-dimensional
axisymmetric lean hydrogen/oxygen spherical expanding flame studies.

It is **DNS-like**, not a strict production DNS solver. The first implemented
stage uses Cantera one-dimensional freely propagating premixed flames, maps the
resulting profiles onto an axisymmetric `r-z` mesh, creates an initial spherical
hot kernel, and runs a post-processing workflow for flame radius, burning
velocity, stretch rate, and Markstein length. The two-dimensional reacting-flow
solver interfaces are present as extension points.

## Quick Start

```bash
cd spherical_flame_dns_like
python -m pip install -r requirements.txt
python scripts/run_sweep.py
python scripts/plot_results.py
```

If Cantera is unavailable or a very lean case fails to converge, the framework
falls back to an analytic placeholder flame profile and records this in the case
summary. This keeps the project runnable while preserving the Cantera path. The
fallback laminar-speed model is a smooth `S_L(phi)` curve with a peak near
`phi = 1.1`; the default 6%-14% H2/O2 cases correspond to `phi ~= 0.032-0.081`,
so they lie on the far-lean branch where `S_L` should increase with H2 content.
For the Stage 1 expanding-flame diagnostics, `free_flame.speed_calibration`
stores the experimental-style no-stretch intercept `S_b0` from the `S_b-kappa`
plot and the density ratio `rho_b/rho_u`. The scalar laminar burning velocity
is then computed as `S_L = (rho_b/rho_u) S_b0`. The 1D Cantera or fallback
profiles still provide the flame structure.

## Implemented Stage 1

- Compute 1D freely propagating premixed H2/O2 flames with Cantera.
- Extract laminar burning velocity `S_L`, flame thickness `delta_f`, temperature
  profiles, and species profiles.
- Check that the sweep-level `S_L` trend is non-decreasing with H2 content for
  the configured far-lean cases.
- Map 1D profiles to a 2D axisymmetric `r-z` grid.
- Build a spherical high-temperature ignition kernel.
- Extract flame radius from a temperature isosurface or maximum gradient.
- Extract 2D axisymmetric flame-front contours from the
  `T = (Tu + Tb) / 2` isotherm.
- Compute expanding flame speed `S_b = dR_f/dt`.
- Compute stretch rate `kappa = 2 S_b / R_f`.
- Fit `S_b = S_b0 - L_b kappa` for Markstein length.
- Save all case data as CSV.
- Save plots as PNG and SVG.

The Stage 1 expanding-flame history is controlled until the 2D solver is added.
Its configured Markstein-number trend is real-valued and increases with H2
volume fraction. The default lean H2/O2 cases use negative Markstein numbers
that become less negative as H2 content increases, so the plotted
`S_b-kappa` relation has a positive fitted slope. The code reports both the
direct slope `dS_b/dkappa` and the conventional `L_b = -dS_b/dkappa`. Complex
growth rates or complex eigenvalues belong to the reserved cellular-instability
analysis interface, not to the scalar `S_b-kappa` Markstein fit.

## Reserved Stage 2 Interfaces

- 2D axisymmetric low-Mach reacting-flow solver. `src/solver.py` now contains
  an explicit predictor/RHS framework for continuity, momentum, temperature,
  and species equations; pressure projection and full chemistry coupling remain
  extension points.
- Gravity and buoyancy source terms.
- Local flame-front curvature.
- Upward/downward local flame speeds.
- Flame rise velocity.
- Cellular instability perturbation analysis.

## Model Scope Versus Strict DNS

A strict DNS would directly resolve the full unsteady reacting Navier-Stokes or
low-Mach equations, detailed transport, detailed chemistry, pressure projection,
species diffusion, viscous stresses, heat release, and all relevant physical
time and length scales in the 2D/3D domain. This framework currently uses
resolved 1D Cantera flame structure and DNS-style grid/diagnostic checks, then
generates a controlled expanding-flame history for post-processing. It is meant
for research workflow development, sensitivity studies, initialization, and
diagnostic prototyping before a full 2D solver is added.

## Units

SI units are used throughout:

- Length: m
- Time: s
- Temperature: K
- Pressure: Pa
- Speed: m/s
- Stretch rate: 1/s
- Markstein length: m
- Markstein number: dimensionless, `Ma = L_b / delta_f`

## Main Outputs

For each case under `results/<case_name>/`:

- `free_flame_profile.csv`
- `initial_field.csv`
- `flame_history.csv`
- `flame_front_history.csv` for the representative plotted case
- `case_summary.csv`

Sweep-level outputs under `results/`:

- `sweep_summary.csv`
- `figures/Rf_t.png`, `Rf_t.svg`
- `figures/Sb_kappa.png`, `Sb_kappa.svg`
- `figures/Sb_Rf.png`, `Sb_Rf.svg`
- `figures/Ma_h2.png`, `Ma_h2.svg`
- `figures/flame_front_evolution_<case>.png`,
  `figures/flame_front_evolution_<case>.svg`
