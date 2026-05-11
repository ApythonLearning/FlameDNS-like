from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass
class Flame1DResult:
    profile: pd.DataFrame
    S_L_m_s: float
    delta_f_m: float
    T_unburned_K: float
    T_burned_K: float
    used_cantera: bool
    status: str


def h2_o2_mole_fraction(h2_volume_fraction: float) -> dict[str, float]:
    x_h2 = float(h2_volume_fraction)
    if not 0.0 < x_h2 < 1.0:
        raise ValueError("h2_volume_fraction must be between 0 and 1.")
    return {"H2": x_h2, "O2": 1.0 - x_h2}


def equivalence_ratio_h2_o2(h2_volume_fraction: float) -> float:
    # Stoichiometric H2/O2 mole ratio is 2.0.
    return float((h2_volume_fraction / (1.0 - h2_volume_fraction)) / 2.0)


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
        result = _make_result(profile, abs(float(velocity[0])), True, "cantera_converged")
    except Exception as exc:
        result = _fallback_profile(h2_volume_fraction, config, f"Cantera failed: {exc}")

    result.profile.to_csv(output_dir / "free_flame_profile.csv", index=False)
    return result


def _make_result(profile: pd.DataFrame, S_L_m_s: float, used_cantera: bool, status: str) -> Flame1DResult:
    x = profile["x_m"].to_numpy()
    T = profile["T_K"].to_numpy()
    dTdx = np.gradient(T, x)
    max_grad = float(np.max(np.abs(dTdx)))
    delta_f = float((T.max() - T.min()) / max_grad) if max_grad > 0.0 else float("nan")
    return Flame1DResult(
        profile=profile,
        S_L_m_s=float(S_L_m_s),
        delta_f_m=delta_f,
        T_unburned_K=float(T[0]),
        T_burned_K=float(T[-1]),
        used_cantera=used_cantera,
        status=status,
    )


def _fallback_profile(h2_volume_fraction: float, config: dict, status: str) -> Flame1DResult:
    Tu = float(config["gas"]["temperature_K"])
    phi = equivalence_ratio_h2_o2(h2_volume_fraction)
    Tb = Tu + 1700.0 * np.clip(phi / 0.08, 0.35, 1.4)
    delta = 0.0020 / max(h2_volume_fraction / 0.10, 0.7)
    S_L = 0.08 + 1.2 * max(h2_volume_fraction - 0.05, 0.0) ** 0.8
    x = np.linspace(0.0, float(config["free_flame"]["width_m"]), 400)
    x0 = 0.5 * x[-1]
    T = Tu + 0.5 * (Tb - Tu) * (1.0 + np.tanh((x - x0) / (0.25 * delta)))
    y_h2 = h2_volume_fraction * 0.02 * (1.0 - (T - Tu) / (Tb - Tu + 1e-12))
    y_o2 = (1.0 - h2_volume_fraction) * 0.98 * (1.0 - 0.4 * (T - Tu) / (Tb - Tu + 1e-12))
    y_h2o = np.clip(1.0 - y_h2 - y_o2, 0.0, 1.0)
    profile = pd.DataFrame({"x_m": x, "T_K": T, "u_m_s": S_L, "Y_H2": y_h2, "Y_O2": y_o2, "Y_H2O": y_h2o})
    return _make_result(profile, float(S_L), False, status)
