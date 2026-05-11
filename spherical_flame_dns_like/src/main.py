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
        "rho_b_over_rho_u": flame.rho_b_over_rho_u,
        "S_b0_unstretched_input_m_s": flame.S_b0_unstretched_m_s,
        "initial_radius_m": float(config["ignition"]["initial_radius_m"]),
        "points_per_flame_thickness": ppft,
        "grid_resolution_ok": bool(ppft >= min_points),
        "S_b0_model_m_s": float(history["S_b0_model_m_s"].iloc[0]),
        "rho_b_over_rho_u_model": float(history["rho_b_over_rho_u_model"].iloc[0]),
        "S_L_model_m_s": float(history["S_L_model_m_s"].iloc[0]),
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
    summary = add_flame_speed_trend_diagnostics(summary, config)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_dir / "sweep_summary.csv", index=False)
    return summary


def add_flame_speed_trend_diagnostics(summary: pd.DataFrame, config: dict) -> pd.DataFrame:
    out = summary.sort_values("h2_volume_fraction").copy()
    speeds = out["S_L_m_s"].to_numpy()
    if len(speeds) > 1:
        increments = pd.Series(speeds).diff().fillna(0.0).to_numpy()
        out["S_L_increment_from_previous_m_s"] = increments
        out["S_L_monotonic_with_h2"] = bool((increments[1:] >= -1.0e-10).all())
    else:
        out["S_L_increment_from_previous_m_s"] = 0.0
        out["S_L_monotonic_with_h2"] = True
    speed_model = config.get("free_flame", {}).get("fallback_speed_model", {})
    out["S_L_reference_peak_phi"] = float(speed_model.get("phi_peak", 1.10))
    out["S_L_reference_peak_note"] = "fallback model peaks near phi=1.1; current 6-14% H2/O2 cases are far lean"
    return out
