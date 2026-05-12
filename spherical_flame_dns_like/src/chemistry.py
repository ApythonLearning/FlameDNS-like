from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .transport import placeholder_diffusivity_m2_s


@dataclass
class Flame1DResult:
    profile: pd.DataFrame
    S_L_m_s: float
    delta_f_m: float
    delta_f_gradient_m: float
    T_unburned_K: float
    T_burned_K: float
    used_cantera: bool
    status: str
    flame_thickness_method: str = "thermal_gradient"
    inner_layer_temperature_K: float | None = None
    thermal_diffusivity_T0_m2_s: float | None = None
    rho_unburned_kg_m3: float | None = None
    rho_b_over_rho_u: float | None = None
    S_b0_unstretched_m_s: float | None = None


def h2_o2_mole_fraction(h2_volume_fraction: float) -> dict[str, float]:
    x_h2 = float(h2_volume_fraction)
    if not 0.0 < x_h2 < 1.0:
        raise ValueError("h2_volume_fraction must be between 0 and 1.")
    return {"H2": x_h2, "O2": 1.0 - x_h2}


def equivalence_ratio_h2_o2(h2_volume_fraction: float) -> float:
    # Stoichiometric H2/O2 mole ratio is 2.0.
    return float((h2_volume_fraction / (1.0 - h2_volume_fraction)) / 2.0)


def reference_laminar_speed_h2_o2(phi: float, config: dict | None = None) -> float:
    """Smooth reference S_L(phi) used only for fallback and sanity checks.

    For H2/O2 mixtures the laminar speed increases on the lean branch and is
    represented with a broad peak near phi ~= 1.1. This is not a substitute for
    Cantera; it keeps non-converged very-lean cases on a physically ordered
    scale instead of mixing arbitrary placeholder values with converged cases.
    """
    model_config = {}
    if config is not None:
        model_config = config.get("free_flame", {}).get("fallback_speed_model", {})
    phi_peak = float(model_config.get("phi_peak", 1.10))
    S_peak = float(model_config.get("S_peak_m_s", 12.0))
    lean_width = float(model_config.get("lean_log_width", 0.72))
    rich_width = float(model_config.get("rich_log_width", 1.10))
    S_floor = float(model_config.get("S_floor_m_s", 1.0e-4))

    phi_safe = max(float(phi), 1.0e-8)
    width = lean_width if phi_safe <= phi_peak else rich_width
    log_distance = np.log(phi_safe / phi_peak)
    speed = S_peak * np.exp(-0.5 * (log_distance / width) ** 2)
    return float(max(speed, S_floor))


def reference_flame_thickness_h2_o2(h2_volume_fraction: float, config: dict | None = None) -> float:
    """Smooth effective thermal flame thickness used for calibrated Stage 1 histories.

    Cantera can return a numerically converged but overly diffuse profile for
    the very lean small-domain cases. The 2D synthetic expansion model should
    not let that outlier dominate ``L_b = Ma * delta_f`` and invert the speed
    trend. This surrogate follows the fallback profile width and is only used
    when the calibrated Stage 1 table is active and the extracted thickness is
    clearly inconsistent.
    """
    thickness_config = {}
    if config is not None:
        thickness_config = config.get("free_flame", {}).get("fallback_thickness_model", {})
    base_h2 = float(thickness_config.get("base_h2_volume_fraction", 0.10))
    base_delta = float(thickness_config.get("base_delta_f_m", 0.0020))
    min_ratio = float(thickness_config.get("min_h2_ratio", 0.7))
    ratio = max(float(h2_volume_fraction) / max(base_h2, 1.0e-12), min_ratio)
    return float(base_delta / ratio)


def preheat_zone_flame_thickness_h2_o2(
    h2_volume_fraction: float, S_L_m_s: float, config: dict
) -> dict[str, float | str]:
    """Evaluate delta_f = (k/cp)_{T0} / (rho_u S_L) = alpha(T0) / S_L.

    ``T0`` is a configurable characteristic inner-layer temperature. Replace
    the default in ``config/base.yaml`` with the gas-specific value from the
    reference used for a production study.
    """
    if S_L_m_s <= 0.0:
        raise ValueError("S_L_m_s must be positive for preheat-zone flame thickness.")
    thickness_config = config.get("free_flame", {}).get("flame_thickness", {})
    T0 = float(thickness_config.get("inner_layer_temperature_K", 650.0))
    pressure = float(config["gas"]["pressure_Pa"])
    Tu = float(config["gas"]["temperature_K"])
    composition = h2_o2_mole_fraction(h2_volume_fraction)
    try:
        import cantera as ct

        gas = ct.Solution(config["gas"].get("mechanism", "h2o2.yaml"))
        gas.TPX = T0, pressure, composition
        gas.transport_model = config["transport"].get("cantera_transport_model", "mixture-averaged")
        k_over_cp = float(gas.thermal_conductivity / gas.cp_mass)
        gas.TPX = Tu, pressure, composition
        rho_u = float(gas.density)
        alpha = k_over_cp / rho_u
        source = "cantera_transport_at_T0"
    except Exception:
        gas_constant = 8.31446261815324
        molecular_weights = {"H2": 2.01588e-3, "O2": 31.998e-3}
        mean_mw = sum(composition[name] * molecular_weights[name] for name in composition)
        rho_u = pressure * mean_mw / (gas_constant * Tu)
        alpha = placeholder_diffusivity_m2_s(T0)
        k_over_cp = alpha * rho_u
        source = "fallback_diffusivity_at_T0"
    return {
        "delta_f_m": float(alpha / S_L_m_s),
        "inner_layer_temperature_K": T0,
        "thermal_diffusivity_T0_m2_s": float(alpha),
        "k_over_cp_T0_kg_m_s": float(k_over_cp),
        "rho_unburned_kg_m3": float(rho_u),
        "source": source,
    }


def speed_calibration_from_h2(h2_volume_fraction: float, config: dict) -> dict[str, float] | None:
    """Return configured Stage 1 speed calibration, or ``None``.

    The input table stores the plotted no-stretch intercept ``S_b0`` and the
    density ratio ``rho_b/rho_u``. The laminar burning velocity is calculated as
    ``S_L = (rho_b/rho_u) S_b0``.
    """
    model = config.get("free_flame", {}).get("speed_calibration", {})
    if not bool(model.get("enabled", False)):
        return None
    h2_points = model.get("h2_volume_fraction_points", [])
    sb0_points = model.get("S_b0_unstretched_m_s_points", [])
    density_ratio_points = model.get("rho_b_over_rho_u_points", [])
    if len(h2_points) != len(sb0_points) or len(h2_points) != len(density_ratio_points) or len(h2_points) < 2:
        raise ValueError("speed_calibration requires matching h2, S_b0, and density-ratio point lists.")
    h2 = np.asarray(h2_points, dtype=float)
    sb0 = np.asarray(sb0_points, dtype=float)
    density_ratio = np.asarray(density_ratio_points, dtype=float)
    order = np.argsort(h2)
    h2_value = float(h2_volume_fraction)
    S_b0_value = float(np.interp(h2_value, h2[order], sb0[order]))
    density_ratio_value = float(np.interp(h2_value, h2[order], density_ratio[order]))
    return {
        "rho_b_over_rho_u": density_ratio_value,
        "S_b0_unstretched_m_s": S_b0_value,
        "S_L_m_s": density_ratio_value * S_b0_value,
    }


def compute_free_flame(h2_volume_fraction: float, config: dict, output_dir: Path) -> Flame1DResult:
    try:
        import cantera as ct
    except Exception as exc:
        return _fallback_profile(h2_volume_fraction, config, f"Cantera unavailable: {exc}")

    gas_config = config["gas"]
    flame_config = config["free_flame"]
    gas = ct.Solution(gas_config.get("mechanism", "h2o2.yaml"))
    gas.TPX = (
        float(gas_config["temperature_K"]),
        float(gas_config["pressure_Pa"]),
        h2_o2_mole_fraction(h2_volume_fraction),
    )

    flame = ct.FreeFlame(gas, width=float(flame_config["width_m"]))
    flame.transport_model = config["transport"].get("cantera_transport_model", "mixture-averaged")
    flame.set_refine_criteria(
        ratio=float(flame_config.get("ratio", 3.0)),
        slope=float(flame_config.get("slope", 0.06)),
        curve=float(flame_config.get("curve", 0.12)),
    )

    try:
        flame.solve(loglevel=int(flame_config.get("loglevel", 0)), auto=True)
        x = np.asarray(flame.grid, dtype=float)
        T = np.asarray(flame.T, dtype=float)
        velocity = np.asarray(flame.velocity, dtype=float)
        species = {}
        for name in gas.species_names:
            species[f"Y_{name}"] = np.asarray(flame.Y[gas.species_index(name), :], dtype=float)
        profile = pd.DataFrame({"x_m": x, "T_K": T, "u_m_s": velocity, **species})
        result = _make_result(profile, abs(float(velocity[0])), True, "cantera_converged", h2_volume_fraction, config)
    except Exception as exc:
        result = _fallback_profile(h2_volume_fraction, config, f"Cantera failed: {exc}")

    result = _apply_speed_calibration(result, h2_volume_fraction, config)
    result.profile.to_csv(output_dir / "free_flame_profile.csv", index=False)
    return result


def _thermal_gradient_flame_thickness(profile: pd.DataFrame) -> float:
    x = profile["x_m"].to_numpy()
    T = profile["T_K"].to_numpy()
    dTdx = np.gradient(T, x)
    max_grad = float(np.max(np.abs(dTdx)))
    return float((T.max() - T.min()) / max_grad) if max_grad > 0.0 else float("nan")


def _make_result(
    profile: pd.DataFrame,
    S_L_m_s: float,
    used_cantera: bool,
    status: str,
    h2_volume_fraction: float,
    config: dict,
) -> Flame1DResult:
    T = profile["T_K"].to_numpy()
    delta_f_gradient = _thermal_gradient_flame_thickness(profile)
    thickness_config = config.get("free_flame", {}).get("flame_thickness", {})
    method = str(thickness_config.get("method", "preheat_zone"))
    if method == "preheat_zone":
        thickness = preheat_zone_flame_thickness_h2_o2(h2_volume_fraction, float(S_L_m_s), config)
        delta_f = float(thickness["delta_f_m"])
        status = f"{status}; delta_f_preheat_zone_{thickness['source']}"
        inner_layer_temperature = float(thickness["inner_layer_temperature_K"])
        thermal_diffusivity = float(thickness["thermal_diffusivity_T0_m2_s"])
        rho_unburned = float(thickness["rho_unburned_kg_m3"])
    else:
        delta_f = delta_f_gradient
        inner_layer_temperature = None
        thermal_diffusivity = None
        rho_unburned = None
    return Flame1DResult(
        profile=profile,
        S_L_m_s=float(S_L_m_s),
        delta_f_m=delta_f,
        delta_f_gradient_m=delta_f_gradient,
        T_unburned_K=float(T[0]),
        T_burned_K=float(T[-1]),
        used_cantera=used_cantera,
        status=status,
        flame_thickness_method=method,
        inner_layer_temperature_K=inner_layer_temperature,
        thermal_diffusivity_T0_m2_s=thermal_diffusivity,
        rho_unburned_kg_m3=rho_unburned,
    )


def _fallback_profile(h2_volume_fraction: float, config: dict, status: str) -> Flame1DResult:
    Tu = float(config["gas"]["temperature_K"])
    phi = equivalence_ratio_h2_o2(h2_volume_fraction)
    Tb = Tu + 1700.0 * np.clip(phi / 0.08, 0.35, 1.4)
    delta = reference_flame_thickness_h2_o2(h2_volume_fraction, config)
    calibration = speed_calibration_from_h2(h2_volume_fraction, config)
    if calibration is None:
        S_L = reference_laminar_speed_h2_o2(phi, config)
        rho_b_over_rho_u = None
        S_b0_unstretched = None
    else:
        S_L = calibration["S_L_m_s"]
        rho_b_over_rho_u = calibration["rho_b_over_rho_u"]
        S_b0_unstretched = calibration["S_b0_unstretched_m_s"]
    x = np.linspace(0.0, float(config["free_flame"]["width_m"]), 400)
    x0 = 0.5 * x[-1]
    T = Tu + 0.5 * (Tb - Tu) * (1.0 + np.tanh((x - x0) / (0.25 * delta)))
    y_h2 = h2_volume_fraction * 0.02 * (1.0 - (T - Tu) / (Tb - Tu + 1e-12))
    y_o2 = (1.0 - h2_volume_fraction) * 0.98 * (1.0 - 0.4 * (T - Tu) / (Tb - Tu + 1e-12))
    y_h2o = np.clip(1.0 - y_h2 - y_o2, 0.0, 1.0)
    profile = pd.DataFrame({"x_m": x, "T_K": T, "u_m_s": S_L, "Y_H2": y_h2, "Y_O2": y_o2, "Y_H2O": y_h2o})
    result = _make_result(profile, float(S_L), False, status, h2_volume_fraction, config)
    result.rho_b_over_rho_u = rho_b_over_rho_u
    result.S_b0_unstretched_m_s = S_b0_unstretched
    return result


def _apply_speed_calibration(result: Flame1DResult, h2_volume_fraction: float, config: dict) -> Flame1DResult:
    calibration = speed_calibration_from_h2(h2_volume_fraction, config)
    if calibration is None:
        return result
    profile = result.profile.copy()
    profile["u_m_s"] = float(calibration["S_L_m_s"])
    status_parts = [result.status, "S_b0_to_S_L_calibrated_for_stage1"]
    thickness_config = config.get("free_flame", {}).get("flame_thickness", {})
    method = str(thickness_config.get("method", result.flame_thickness_method))
    if method == "preheat_zone":
        thickness = preheat_zone_flame_thickness_h2_o2(h2_volume_fraction, float(calibration["S_L_m_s"]), config)
        delta_f = float(thickness["delta_f_m"])
        inner_layer_temperature = float(thickness["inner_layer_temperature_K"])
        thermal_diffusivity = float(thickness["thermal_diffusivity_T0_m2_s"])
        rho_unburned = float(thickness["rho_unburned_kg_m3"])
        status_parts.append(f"delta_f_preheat_zone_recomputed_{thickness['source']}")
    else:
        delta_f = float(result.delta_f_m)
        inner_layer_temperature = result.inner_layer_temperature_K
        thermal_diffusivity = result.thermal_diffusivity_T0_m2_s
        rho_unburned = result.rho_unburned_kg_m3
    status = "; ".join(status_parts)
    return Flame1DResult(
        profile=profile,
        S_L_m_s=float(calibration["S_L_m_s"]),
        delta_f_m=delta_f,
        delta_f_gradient_m=result.delta_f_gradient_m,
        T_unburned_K=result.T_unburned_K,
        T_burned_K=result.T_burned_K,
        used_cantera=result.used_cantera,
        status=status,
        flame_thickness_method=method,
        inner_layer_temperature_K=inner_layer_temperature,
        thermal_diffusivity_T0_m2_s=thermal_diffusivity,
        rho_unburned_kg_m3=rho_unburned,
        rho_b_over_rho_u=float(calibration["rho_b_over_rho_u"]),
        S_b0_unstretched_m_s=float(calibration["S_b0_unstretched_m_s"]),
    )
