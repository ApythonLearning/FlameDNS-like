#ifndef PELELMEX_PROB_H
#define PELELMEX_PROB_H

#include <AMReX_Geometry.H>
#include <AMReX_REAL.H>
#include <AMReX_Array4.H>
#include <AMReX_Math.H>

#include "mechanism.H"

struct MyProbParm : public ProbParmDefault
{
  amrex::Real T0 = 298.0;
  amrex::Real phi = 0.05555555555555556;
  amrex::Real h2_volume_fraction = 0.10;
  int use_phi = 0;
  amrex::Real ignition_radius = 7.5e-4;
  amrex::Real initial_flame_radius = 1.5e-3;
  amrex::Real ignition_temperature = 1800.0;
  amrex::Real gravity_magnitude = 0.0;
  amrex::Real gravity_direction[2] = {0.0, -1.0};
  amrex::Real domain_size[2] = {0.006, 0.012};
  int grid_resolution[2] = {320, 640};
  int max_step = 200;
  amrex::Real cfl = 0.3;
  int plot_int = 20;
  char profile_csv[256] = "../../profiles/h2o2_phi_auto/h2o2_1d_profile.csv";
};

struct MyProblemSpecificFunctions : public DefaultProblemSpecificFunctions
{
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  void initdata(
    int i,
    int j,
    int k,
    int /*is_incompressible*/,
    amrex::Array4<amrex::Real> const& state,
    amrex::Array4<amrex::Real> const& /*aux*/,
    amrex::GeometryData const& geomdata,
    MyProbParm const& prob_parm,
    pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/) const noexcept
  {
    const auto plo = geomdata.ProbLo();
    const auto dx = geomdata.CellSize();
    const amrex::Real r = plo[0] + (static_cast<amrex::Real>(i) + 0.5) * dx[0];
    const amrex::Real z = plo[1] + (static_cast<amrex::Real>(j) + 0.5) * dx[1];
    const amrex::Real radius = amrex::Math::sqrt(r * r + z * z);

    const amrex::Real x_h2 = prob_parm.use_phi != 0
                                ? 2.0 * prob_parm.phi / (1.0 + 2.0 * prob_parm.phi)
                                : prob_parm.h2_volume_fraction;
    const amrex::Real x_o2 = 1.0 - x_h2;
    const amrex::Real T = radius <= prob_parm.ignition_radius ? prob_parm.ignition_temperature : prob_parm.T0;

    amrex::Real molefrac[NUM_SPECIES] = {0.0};
    molefrac[H2_ID] = x_h2;
    molefrac[O2_ID] = x_o2;

    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    eos.X2Y(molefrac, massfrac);

    amrex::Real rho = 0.0;
    eos.PYT2R(prob_parm.P_mean, massfrac, T, rho);

    amrex::Real h = 0.0;
    eos.TY2H(T, massfrac, h);

    state(i, j, k, URHO) = rho;
    state(i, j, k, UMX) = 0.0;
    state(i, j, k, UMY) = 0.0;
#if (AMREX_SPACEDIM == 3)
    state(i, j, k, UMZ) = 0.0;
#endif
    state(i, j, k, UEDEN) = rho * h;
    state(i, j, k, UTEMP) = T;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      state(i, j, k, UFS + n) = rho * massfrac[n];
    }
  }
};

#endif
