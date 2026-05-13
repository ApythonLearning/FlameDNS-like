# Low-Mach Axisymmetric H2/O2 Spherical Flame Workflow

This project prepares and postprocesses a PeleLMeX low-Mach, two-dimensional
axisymmetric premixed H2/O2 spherical-flame calculation.

It does not fabricate DNS results. The Python-only synthetic expansion scripts
are retained only for reduced-order diagnostics and are not a substitute for
PeleLMeX plotfiles.

## Dependency Check

Run:

```powershell
cd spherical_flame_dns_like
python scripts/check_dependencies.py
```

Current local check from this workspace:

- Cantera: found, version 3.2.0 in the active `flame_dns` conda environment.
- PeleLMeX executable: not found.
- `PELE_HOME`: not set.
- `AMREX_HOME`: not set.
- `yt`: installed after dependency setup.
- `h5py`: installed but currently has a Windows DLL import failure in this
  environment; repair with conda-forge or a clean environment before requiring
  HDF5 output.

Python dependencies:

```powershell
python -m pip install -r requirements.txt
```

If `h5py` fails to import on Windows, prefer a conda-forge install in the active
environment:

```powershell
conda install -c conda-forge h5py hdf5
```

## PeleLMeX / AMReX Install Notes

PeleLMeX and AMReX are not installed in the current machine environment. A
typical Linux or WSL setup is:

```bash
git clone --recursive https://github.com/AMReX-Combustion/PeleLMeX.git
cd PeleLMeX
git submodule update --init --recursive
export PELE_HOME=$PWD
```

Install a C++ compiler, MPI, CMake, Python, and the PeleLMeX third-party
libraries required by your target build. From a case directory:

```bash
make -j4 TPL
make -j4
./PeleLMeX2d.gnu.MPI.ex inputs.axisym_h2o2
```

On native Windows, WSL is the more practical path for PeleLMeX/AMReX builds.

## Step 1: 1D Cantera Premixed Flame

Compute the 1D freely propagating H2/O2 flame:

```powershell
python scripts/compute_1d_h2o2_flame.py --phi 0.6 --output-dir profiles/h2o2_phi_060
```

Outputs:

- `h2o2_1d_summary.csv`
- `h2o2_1d_summary.json`
- `h2o2_1d_profile.csv`
- `h2o2_1d_profile.h5`

The summary contains `S_L_m_s`, `delta_T_m`, `Tb_K`, `rho_u_kg_m3`, and
`rho_b_kg_m3`. The profile file stores temperature, velocity, density, mass
fractions, and mole fractions along the 1D flame coordinate.

Very lean H2/O2, for example `--h2-volume-fraction 0.10` (`phi ~= 0.0556`), may
fail Cantera's freely propagating flame solve. The script intentionally fails in
that case instead of substituting synthetic flame data.

## Step 2: PeleLMeX Case Template

The case template is in:

```text
examples/pelelmex_axisymmetric_h2o2/
```

It includes:

- `GNUmakefile`
- `inputs.axisym_h2o2`
- `pelelmex_prob.H`
- `pelelmex_prob.cpp`

The inputs expose:

- `prob.phi` or `prob.h2_volume_fraction`
- `prob.T0`
- `prob.P_mean`
- `prob.ignition_radius`
- `prob.ignition_temperature`
- `prob.gravity_magnitude`
- `prob.gravity_direction`
- `prob.domain_size`
- `prob.grid_resolution`
- `prob.max_step`
- `prob.cfl`
- `prob.plot_int`

`pelelmex_prob.H` initializes a spherical high-temperature kernel in the RZ
domain. The generated 1D profiles are referenced by `prob.profile_csv` for the
next initialization refinement step.

## Step 3: Plotfile Postprocessing

Postprocessing requires real PeleLMeX/AMReX plotfiles and `yt`:

```powershell
python scripts/postprocess_pelelmex_plotfiles.py plt00020 plt00040 --summary profiles/h2o2_phi_060/h2o2_1d_summary.json --output-dir postprocess_pelelmex
```

The script reads the true temperature field, extracts the
`T = (Tu + Tb) / 2` contour, and writes:

- `flame_radii.csv`
- `flame_radii.png`
- `flame_contours.png`

It reports upper, lower, mean, and maximum flame radius. It does not create
placeholder plotfile data.

## Step 4: DNS Resolution Check

Check the DNS mesh and CFL criteria with:

```powershell
python scripts/check_dns_resolution.py --summary profiles/h2o2_phi_060/h2o2_1d_summary.json --dx 1.875e-5 --cfl 0.3
```

The required criterion is:

```text
dx <= delta_T / 10
CFL < 0.5
```

The recommended criterion is:

```text
dx <= delta_T / 20
```

## Legacy Reduced-Order Scripts

The older `run_case.py`, `run_sweep.py`, and plotting scripts below are
diagnostic prototypes. They may use analytic fallback profiles and synthetic
expansion histories to keep the reduced-order workflow runnable. Do not treat
their outputs as DNS or PeleLMeX results. For DNS postprocessing, use
`postprocess_pelelmex_plotfiles.py` on real plotfiles only.

The legacy fallback laminar-speed model is a smooth `S_L(phi)` curve with a
peak near `phi = 1.1`; the default 6%-16% H2/O2 cases correspond to
`phi ~= 0.032-0.095`, so they lie on the far-lean branch where `S_L` should
increase with H2 content.
For the Stage 1 expanding-flame diagnostics, `free_flame.speed_calibration`
stores the experimental-style no-stretch intercept `S_b0` from the `S_b-kappa`
plot and the density ratio `rho_b/rho_u`. The scalar laminar burning velocity
is then computed as `S_L = (rho_b/rho_u) S_b0`. The 1D Cantera or fallback
profiles still provide the flame structure.

## Implemented Stage 1

- Compute 1D freely propagating premixed H2/O2 flames with Cantera.
- Extract laminar burning velocity `S_L`, preheat-zone flame thickness
  `delta_f`, temperature profiles, and species profiles.
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
- Evaluate a linear DL/TD/RT cellular-instability diagnostic from a Law-model
  dispersion relation.
- Integrate a weakly nonlinear cellular-amplitude model for post-linear onset
  and saturation diagnostics.
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

## Cellular Instability Diagnostic

The current framework evaluates a low-order planar dispersion diagnostic for
the combined Darrieus-Landau (DL), thermal-diffusive/Markstein (TD), and
Rayleigh-Taylor (RT) effects discussed by Zheng et al. (2024). The implemented
Law-model form is

```text
omega = sigma S_L / (1 + sigma)
        [-k + sqrt((sigma^2 + sigma - 1) k^2 / sigma
        - 2 (sigma + 1) (1 - Ma) delta_f k^3
        - ((sigma^2 - 1) / sigma^2) (g / S_L^2) k)]
```

where `sigma = rho_u / rho_b`, `k` is the perturbation wave number, and
`Ma = L_b / delta_f`. The square root is allowed to be complex; CSV output
therefore stores both `omega_real_1_s` and `omega_imag_1_s`. The real part is
used for the reported growth/decay rate.

The sign convention follows the paper's Law-model discussion:
`cellular_instability.rt_mode: unstable` uses `g < 0`, `stable` uses `g > 0`,
and `neutral` sets `g = 0`. By default `cellular_instability.rt_modes` evaluates
all three conditions, RT-Unstable, RT-neutral, and RT-Stable, for the same H2
case list.

The output also records the deficient-reactant weighted effective Lewis number
used for labeling cellular-instability comparisons:

```text
Le_eff = (Le_O + A Le_F) / (1 + A),   phi < 1
A = 1 + beta (phi^-1 - 1)
```

For unity-Lewis runs this value is forced to `1`; otherwise the configurable
defaults are `Le_F = 0.32`, `Le_O = 1.15`, and `beta = 6.0`. These are
diagnostic defaults from the H2/air stability literature and should be replaced
if a more appropriate H2/O2 transport fit is available. This diagnostic is not
a nonlinear cellular DNS; it is a configurable stability-analysis layer used to
identify the most amplified wavelength, cutoff wave number, and relative
DL/TD/RT contributions.

The linear dispersion layer is coupled to an optional weakly nonlinear
Stuart-Landau amplitude model:

```text
dA / dt = omega_eff A [1 - (A / A_sat)^2]
omega_eff = omega_linear - C S_b / R_f
```

This adds a finite-amplitude post-linear phase, onset time, saturation
amplitude, and a simple geometric cellular-speed enhancement factor. It is
still a reduced-order diagnostic, not a resolved nonlinear cellular DNS. It is
intended to compare RT-Unstable, RT-neutral, and RT-Stable trends before the
full two-dimensional reacting-flow solver is coupled to flame-front
perturbations. The plotting script also creates a Fig. 4-style front-location
visualization for RT-Stable, RT-neutral, and RT-Unstable conditions. Its
`amplitude_visual_scale` only magnifies the plotted front shape; it does not
modify the saved physical perturbation amplitude.

## Flame Thickness Definition

For laminar hydrogen flames, the default flame thickness is the preheat-zone
definition

```text
delta_f = (k / cp)_{T0} / (rho_u S_L)
        = alpha(T0) / S_L
```

where `k` is thermal conductivity, `cp` is mass specific heat, `rho_u` is the
unburned-mixture density, and `S_L` is the laminar burning velocity. The
characteristic inner-layer temperature `T0` is configured in
`config/base.yaml` as `free_flame.flame_thickness.inner_layer_temperature_K`.
Replace the default value with the gas-specific value from the cited reference
for final studies. The older thermal-gradient thickness is still written to
`delta_f_gradient_m` in `case_summary.csv` for diagnostics.

## Reserved Stage 2 Interfaces

- 2D axisymmetric low-Mach reacting-flow solver. `src/solver.py` now contains
  an explicit predictor/RHS framework for continuity, momentum, temperature,
  and species equations; pressure projection and full chemistry coupling remain
  extension points.
- Gravity and buoyancy source terms.
- Local flame-front curvature.
- Upward/downward local flame speeds.
- Flame rise velocity.
- Resolved nonlinear cellular instability perturbation analysis.

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
- `cellular_instability.csv`
- `nonlinear_cellular.csv`
- `flame_front_history.csv` for the representative plotted case
- `case_summary.csv`

Sweep-level outputs under `results/`:

- `sweep_summary.csv`
- `figures/Rf_t.png`, `Rf_t.svg`
- `figures/Sb_kappa.png`, `Sb_kappa.svg`
- `figures/Sb_Rf.png`, `Sb_Rf.svg`
- `figures/Ma_h2.png`, `Ma_h2.svg`
- `figures/cellular_dispersion_rt_modes.png`,
  `figures/cellular_dispersion_rt_modes.svg`
- `figures/cellular_dispersion_2x4_rt_modes.png`,
  `figures/cellular_dispersion_2x4_rt_modes.svg`
- `figures/cellular_growth_rt_modes.png`,
  `figures/cellular_growth_rt_modes.svg`
- `figures/cellular_lambda_rt_modes.png`,
  `figures/cellular_lambda_rt_modes.svg`
- `figures/nonlinear_cellular_amplitude_2x4_rt_modes.png`,
  `figures/nonlinear_cellular_amplitude_2x4_rt_modes.svg`
- `figures/nonlinear_cellular_onset_rt_modes.png`,
  `figures/nonlinear_cellular_onset_rt_modes.svg`
- `figures/nonlinear_cellular_speed_factor_rt_modes.png`,
  `figures/nonlinear_cellular_speed_factor_rt_modes.svg`
- `figures/nonlinear_front_locations_<case>.png`,
  `figures/nonlinear_front_locations_<case>.svg`
- `figures/flame_front_evolution_<case>.png`,
  `figures/flame_front_evolution_<case>.svg`
