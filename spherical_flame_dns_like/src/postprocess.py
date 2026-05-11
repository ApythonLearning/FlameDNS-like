from __future__ import annotations

import numpy as np
import pandas as pd


def synthetic_markstein_number(config: dict, h2_volume_fraction: float) -> float:
    """Configured real scalar Markstein-number trend for Stage 1 histories."""
    model = config["time_history"].get("markstein_model", {})
    low_h2 = float(model.get("low_h2_volume_fraction", 0.06))
    high_h2 = float(model.get("high_h2_volume_fraction", 0.14))
    low_ma = float(model.get("low_markstein_number", -0.2))
    high_ma = float(model.get("high_markstein_number", 1.1))
    if high_h2 <= low_h2:
        raise ValueError("high_h2_volume_fraction must be larger than low_h2_volume_fraction")
    weight = (float(h2_volume_fraction) - low_h2) / (high_h2 - low_h2)
    weight = float(np.clip(weight, 0.0, 1.0))
    return low_ma + weight * (high_ma - low_ma)


def synthetic_expanding_history(S_L_m_s: float, delta_f_m: float, config: dict, h2_volume_fraction: float) -> pd.DataFrame:
    time_config = config["time_history"]
    ignition_config = config["ignition"]
    t = np.linspace(0.0, float(time_config["t_end_s"]), int(time_config["n_steps"]))
    R = np.empty_like(t)
    R[0] = float(ignition_config["initial_radius_m"])
    Ma_true = synthetic_markstein_number(config, h2_volume_fraction)
    Lb_true = float(time_config.get("markstein_length_factor", 1.0)) * delta_f_m * Ma_true
    Sb0_true = 1.15 * float(S_L_m_s)
    for i in range(1, len(t)):
        dt = t[i] - t[i - 1]
        speed = Sb0_true / max(1.0 + 2.0 * Lb_true / max(R[i - 1], 1e-9), 0.2)
        R[i] = R[i - 1] + speed * dt
    Sb = np.gradient(R, t)
    kappa = 2.0 * Sb / np.maximum(R, 1e-12)
    return pd.DataFrame(
        {
            "t_s": t,
            "R_f_m": R,
            "S_b_m_s": Sb,
            "kappa_1_s": kappa,
            "S_b0_model_m_s": Sb0_true,
            "L_b_model_m": Lb_true,
            "Ma_model": Ma_true,
        }
    )


def fit_markstein(history: pd.DataFrame, fit_config: dict) -> dict[str, float]:
    mask = np.ones(len(history), dtype=bool)
    min_radius = fit_config.get("min_radius_m")
    max_kappa = fit_config.get("max_kappa_1_s")
    if min_radius is not None:
        mask &= history["R_f_m"].to_numpy() >= float(min_radius)
    if max_kappa is not None:
        mask &= np.abs(history["kappa_1_s"].to_numpy()) <= float(max_kappa)
    if mask.sum() < 3:
        mask[:] = True
    kappa = history.loc[mask, "kappa_1_s"].to_numpy()
    Sb = history.loc[mask, "S_b_m_s"].to_numpy()
    slope, intercept = np.polyfit(kappa, Sb, 1)
    L_b = -float(slope)
    S_b0 = float(intercept)
    return {"S_b0_m_s": S_b0, "L_b_m": L_b, "fit_points": int(mask.sum())}


def add_dimensionless_groups(summary: dict[str, float]) -> dict[str, float]:
    delta = float(summary["delta_f_m"])
    R0 = float(summary["initial_radius_m"])
    Lb = float(summary["L_b_m"])
    summary["Ma"] = Lb / delta if delta > 0.0 else np.nan
    summary["Pe"] = R0 / delta if delta > 0.0 else np.nan
    return summary


def check_cfl(S_max_m_s: float, dt_s: float, dx_m: float, cfl_limit: float) -> dict[str, float | bool]:
    cfl = float(abs(S_max_m_s) * dt_s / dx_m) if dx_m > 0.0 else float("inf")
    return {"CFL": cfl, "CFL_ok": bool(cfl <= cfl_limit)}
