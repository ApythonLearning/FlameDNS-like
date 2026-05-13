#include "PeleLMeX.H"

#include <algorithm>
#include <cstring>
#include <string>

void
PeleLM::readProbParm()
{
  amrex::ParmParse pp("prob");

  pp.query("P_mean", PeleLM::prob_parm->P_mean);
  pp.query("T0", PeleLM::prob_parm->T0);
  pp.query("phi", PeleLM::prob_parm->phi);
  pp.query("h2_volume_fraction", PeleLM::prob_parm->h2_volume_fraction);
  pp.query("use_phi", PeleLM::prob_parm->use_phi);
  pp.query("ignition_radius", PeleLM::prob_parm->ignition_radius);
  pp.query("initial_flame_radius", PeleLM::prob_parm->initial_flame_radius);
  pp.query("ignition_temperature", PeleLM::prob_parm->ignition_temperature);
  pp.query("gravity_magnitude", PeleLM::prob_parm->gravity_magnitude);
  pp.queryarr("gravity_direction", PeleLM::prob_parm->gravity_direction, 0, 2);
  pp.queryarr("domain_size", PeleLM::prob_parm->domain_size, 0, 2);
  pp.queryarr("grid_resolution", PeleLM::prob_parm->grid_resolution, 0, 2);
  pp.query("max_step", PeleLM::prob_parm->max_step);
  pp.query("cfl", PeleLM::prob_parm->cfl);
  pp.query("plot_int", PeleLM::prob_parm->plot_int);
  std::string profile_csv(PeleLM::prob_parm->profile_csv);
  pp.query("profile_csv", profile_csv);
  std::strncpy(PeleLM::prob_parm->profile_csv, profile_csv.c_str(), sizeof(PeleLM::prob_parm->profile_csv) - 1);
  PeleLM::prob_parm->profile_csv[sizeof(PeleLM::prob_parm->profile_csv) - 1] = '\0';
}
