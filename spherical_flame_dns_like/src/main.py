from __future__ import annotations

from pathlib import Path

import pandas as pd

from .chemistry import compute_free_flame, equivalence_ratio_h2_o2
from .ignition import apply_hot_kernel, field_to_dataframe, map_profile_to_axisymmetric_mesh
from .mesh import create_axisymmetric_mesh, points_per_flame_thickness
from .postprocess import add_dimensionless_groups, check_cfl, fit_markstein, synthetic_expanding_history
from .utils import case_output_dir


def run_case(case: dict, config: dict) -> dict:
    case_name = case["name"]
    h2_fraction = float(case["h2_volume_fraction"])
    out_dir = case_output_dir(config, case_name)

    flame = compute_free_flame(h2_fraction, config, out_dir)
    mesh = create_axisymmetric_mesh(config["mesh"])
    ppft = points_per_flame_thickness(mesh, flame.delta_f_m)

    fields = map_profile_to_axisymmetric_mesh(mesh, flame.profile, float(config["ignition"]["initial_radius_m"]))
    fields = apply_hot_kernel(fields, mesh, config["ignition"])
    field_to_dataframe(mesh, fields).to_csv(out_dir / "initial_field.csv", index=False)

    history = synthetic_expanding_history(flame.S_L_m_s, flame.delta_f_m, config, h2_fraction)
    history.to_csv(out_dir / "flame_history.csv", index=False)

    fit = fit_markstein(history, config["postprocess"]["fit"])
    dt = float(history["t_s"].iloc[1] - history["t_s"].iloc[0]) if len(history) > 1 else 0.0
    dx = min(mesh.dr, mesh.dz)
    cfl = check_cfl(float(history["S_b_m_s"].abs().max()), dt, dx, float(config["time_history"]["cfl_limit"]))

    min_points = float(config["mesh"]["min_points_per_flame_thickness"])
    summary = {
        "case": case_name,
        "h2_volume_fraction": h2_fraction,
        "phi": equivalence_ratio_h2_o2(h2_fraction),
        "used_cantera": flame.used_cantera,
        "status": flame.status,
        "S_L_m_s": flame.S_L_m_s,
        "delta_f_m": flame.delta_f_m,
        "T_unburned_K": flame.T_unburned_K,
        "T_burned_K": flame.T_burned_K,
        "initial_radius_m": float(config["ignition"]["initial_radius_m"]),
        "points_per_flame_thickness": ppft,
        "grid_resolution_ok": bool(ppft >= min_points),
        "S_b0_model_m_s": float(history["S_b0_model_m_s"].iloc[0]),
        "L_b_model_m": float(history["L_b_model_m"].iloc[0]),
        "Ma_model": float(history["Ma_model"].iloc[0]),
        "Markstein_number_type": "real_scalar_from_Sb_kappa_fit",
        **fit,
        **cfl,
    }
    summary = add_dimensionless_groups(summary)
    pd.DataFrame([summary]).to_csv(out_dir / "case_summary.csv", index=False)
    return summary


def run_cases(cases: list[dict], config: dict, output_dir: Path) -> pd.DataFrame:
    rows = [run_case(case, config) for case in cases]
    summary = pd.DataFrame(rows)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_dir / "sweep_summary.csv", index=False)
    return summary
