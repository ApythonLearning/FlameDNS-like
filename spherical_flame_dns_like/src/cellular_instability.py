from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class LawModelInputs:
    """Inputs for the planar DL/TD/RT linear dispersion model.

    The model follows the Law/Sung-style discontinuity dispersion relation
    discussed by Zheng et al. (2024). It is used here as a DNS-like diagnostic
    for small flame-front perturbations, not as a replacement for a resolved
    two-dimensional cellular instability calculation.
    """

    S_L_m_s: float
    delta_f_m: float
    rho_b_over_rho_u: float
    markstein_number: float
    gravity_m_s2: float


def rt_modes(config: dict) -> list[str]:
    instability_config = config.get("cellular_instability", {})
    configured = instability_config.get("rt_modes")
    if configured is None:
        configured = [instability_config.get("rt_mode", "unstable")]
    modes = [str(mode).lower() for mode in configured]
    allowed = {"unstable", "neutral", "stable"}
    invalid = sorted(set(modes) - allowed)
    if invalid:
        raise ValueError(f"Unsupported RT modes: {invalid}. Use unstable, neutral, or stable.")
    return modes


def effective_gravity(config: dict, rt_mode: str | None = None) -> float:
    """Return signed gravity for the RT term in the Law dispersion relation.

    In the convention used here and in the referenced Law-model expression,
    ``g < 0`` is RT-unstable and ``g > 0`` is RT-stable. ``neutral`` suppresses
    the RT term even when gravity is globally enabled.
    """

    instability_config = config.get("cellular_instability", {})
    gravity_config = config.get("gravity", {})
    if not gravity_config.get("enabled", False):
        return 0.0
    magnitude = float(gravity_config.get("acceleration_m_s2", 9.80665))
    mode = str(rt_mode or instability_config.get("rt_mode", "unstable")).lower()
    if mode == "unstable":
        return -magnitude
    if mode == "stable":
        return magnitude
    if mode == "neutral":
        return 0.0
    raise ValueError("cellular_instability.rt_mode must be unstable, stable, or neutral")


def effective_lewis_number(phi: float, config: dict) -> float:
    """Return the deficient-reactant weighted effective Lewis number.

    Zheng et al. use the Matalon-Cui-Bechtold weighting
    ``Le_eff = (Le_O + A Le_F) / (1 + A)`` for lean mixtures, with
    ``A = 1 + beta(phi^-1 - 1)``. The constants are configurable because this
    project targets H2/O2 while the cited paper reports H2/air values.
    """

    if config.get("transport", {}).get("unity_lewis", False):
        return 1.0
    le_config = config.get("cellular_instability", {}).get("effective_lewis", {})
    le_fuel = float(le_config.get("fuel_lewis_number", 0.32))
    le_oxidizer = float(le_config.get("oxidizer_lewis_number", 1.15))
    beta = float(le_config.get("zeldovich_number", 6.0))
    phi = max(float(phi), 1.0e-12)
    if phi < 1.0:
        asymmetry = 1.0 + beta * (1.0 / phi - 1.0)
        return float((le_oxidizer + asymmetry * le_fuel) / (1.0 + asymmetry))
    asymmetry = 1.0 + beta * (phi - 1.0)
    return float((le_fuel + asymmetry * le_oxidizer) / (1.0 + asymmetry))


def law_growth_rate_complex(k_1_m: np.ndarray, inputs: LawModelInputs) -> np.ndarray:
    """Dimensional complex growth rate [1/s] from the Law DL/TD/RT model.

    Equation form:
        omega = sigma*S_L/(1+sigma) * (-k + sqrt(A*k^2 - B*k^3 - C*g*k/S_L^2))

    where ``sigma = rho_u/rho_b`` and ``Ma = L_b/delta_f``. A complex value is
    retained when the square-root argument is negative; the real part controls
    exponential amplitude growth/decay and the imaginary part represents
    oscillatory behavior in this low-order diagnostic.
    """

    k = np.asarray(k_1_m, dtype=float)
    S_L = max(float(inputs.S_L_m_s), 1.0e-12)
    delta = max(float(inputs.delta_f_m), 1.0e-12)
    rho_ratio = float(inputs.rho_b_over_rho_u)
    if not 0.0 < rho_ratio < 1.0:
        raise ValueError("rho_b_over_rho_u must be between 0 and 1")
    sigma = 1.0 / rho_ratio
    ma = float(inputs.markstein_number)
    g = float(inputs.gravity_m_s2)

    hydrodynamic = (sigma**2 + sigma - 1.0) / sigma * k**2
    stretch = 2.0 * (sigma + 1.0) * (1.0 - ma) * delta * k**3
    buoyancy = ((sigma**2 - 1.0) / sigma**2) * (g / S_L**2) * k
    radicand = hydrodynamic - stretch - buoyancy
    return sigma * S_L / (1.0 + sigma) * (-k + np.emath.sqrt(radicand))


def dl_growth_rate(k_1_m: np.ndarray, inputs: LawModelInputs) -> np.ndarray:
    """DL-only limit of the Law model, returned as a real growth rate [1/s]."""

    neutral_surface = LawModelInputs(
        S_L_m_s=inputs.S_L_m_s,
        delta_f_m=0.0,
        rho_b_over_rho_u=inputs.rho_b_over_rho_u,
        markstein_number=inputs.markstein_number,
        gravity_m_s2=0.0,
    )
    return np.real(law_growth_rate_complex(k_1_m, neutral_surface))


def build_dispersion_relation(
    history: pd.DataFrame,
    summary: dict | pd.Series,
    config: dict,
) -> pd.DataFrame:
    """Evaluate DL, TD/Markstein, and RT contributions over sampled times."""

    instability_config = config.get("cellular_instability", {})
    if not instability_config.get("enabled", True):
        return pd.DataFrame()

    n_k = int(instability_config.get("n_wavenumbers", 160))
    k_delta_min = float(instability_config.get("k_delta_min", 0.02))
    k_delta_max = float(instability_config.get("k_delta_max", 1.2))
    time_sample_count = int(instability_config.get("time_sample_count", 9))
    if n_k < 2:
        raise ValueError("cellular_instability.n_wavenumbers must be at least 2")

    h = history.sort_values("t_s").reset_index(drop=True)
    if len(h) == 0:
        return pd.DataFrame()
    sample_idx = np.unique(np.linspace(0, len(h) - 1, max(time_sample_count, 1)).round().astype(int))

    delta = float(summary["delta_f_m"])
    k_delta = np.linspace(k_delta_min, k_delta_max, n_k)
    k = k_delta / max(delta, 1.0e-12)
    base_inputs = LawModelInputs(
        S_L_m_s=float(summary["S_L_m_s"]),
        delta_f_m=delta,
        rho_b_over_rho_u=float(summary.get("rho_b_over_rho_u_model", summary.get("rho_b_over_rho_u", np.nan))),
        markstein_number=float(summary.get("Ma_model", summary.get("Ma", 0.0))),
        gravity_m_s2=0.0,
    )
    le_eff = effective_lewis_number(float(summary.get("phi", np.nan)), config)
    omega_no_rt = law_growth_rate_complex(k, base_inputs)
    omega_dl = dl_growth_rate(k, base_inputs)
    omega_td = np.real(omega_no_rt) - omega_dl
    tau_F = delta / max(base_inputs.S_L_m_s, 1.0e-12)

    rows = []
    for mode in rt_modes(config):
        inputs = replace(base_inputs, gravity_m_s2=effective_gravity(config, mode))
        omega = law_growth_rate_complex(k, inputs)
        omega_rt = np.real(omega) - np.real(omega_no_rt)
        for idx in sample_idx:
            row = h.iloc[int(idx)]
            R = float(row["R_f_m"])
            mode_n = k * max(R, 1.0e-12)
            for j in range(len(k)):
                wavelength = 2.0 * np.pi / k[j]
                rows.append(
                    {
                        "case": summary["case"],
                        "h2_volume_fraction": float(summary.get("h2_volume_fraction", np.nan)),
                        "phi": float(summary.get("phi", np.nan)),
                        "Le_eff": float(le_eff),
                        "t_s": float(row["t_s"]),
                        "R_f_m": R,
                        "k_delta": float(k_delta[j]),
                        "k_1_m": float(k[j]),
                        "mode_n_estimate": float(mode_n[j]),
                        "lambda_m": float(wavelength),
                        "lambda_over_delta": float(wavelength / delta),
                        "omega_real_1_s": float(np.real(omega[j])),
                        "omega_imag_1_s": float(np.imag(omega[j])),
                        "omega_dimensionless": float(np.real(omega[j]) * tau_F),
                        "omega_DL_1_s": float(omega_dl[j]),
                        "omega_TD_Markstein_1_s": float(omega_td[j]),
                        "omega_RT_1_s": float(omega_rt[j]),
                        "rt_mode": mode,
                        "gravity_effective_m_s2": float(inputs.gravity_m_s2),
                    }
                )
    return pd.DataFrame(rows)


def summarize_dispersion(dispersion: pd.DataFrame) -> dict[str, float | bool | str]:
    """Return compact cellular-instability metrics for a case summary row."""

    if dispersion.empty:
        return {
            "cellular_enabled": False,
            "cellular_unstable": False,
            "cellular_max_growth_1_s": np.nan,
            "cellular_k_delta_at_max": np.nan,
            "cellular_lambda_over_delta_at_max": np.nan,
            "cellular_mode_n_at_max": np.nan,
            "cellular_cutoff_k_delta": np.nan,
            "cellular_rt_mode": "disabled",
        }

    result: dict[str, float | bool | str] = {
        "cellular_enabled": True,
        "cellular_Le_eff": float(dispersion["Le_eff"].iloc[0]) if "Le_eff" in dispersion else np.nan,
    }
    final_t = float(dispersion["t_s"].max())
    final_all = dispersion[dispersion["t_s"] == final_t].copy()
    for mode, final in final_all.groupby("rt_mode"):
        idx = final["omega_real_1_s"].idxmax()
        peak = final.loc[idx]
        positive = final[final["omega_real_1_s"] > 0.0]
        cutoff = float(positive["k_delta"].max()) if not positive.empty else np.nan
        prefix = f"cellular_rt_{mode}"
        result.update(
            {
                f"{prefix}_unstable": bool(float(peak["omega_real_1_s"]) > 0.0),
                f"{prefix}_max_growth_1_s": float(peak["omega_real_1_s"]),
                f"{prefix}_k_delta_at_max": float(peak["k_delta"]),
                f"{prefix}_lambda_over_delta_at_max": float(peak["lambda_over_delta"]),
                f"{prefix}_mode_n_at_max": float(peak["mode_n_estimate"]),
                f"{prefix}_cutoff_k_delta": cutoff,
                f"{prefix}_gravity_effective_m_s2": float(peak["gravity_effective_m_s2"]),
            }
        )

    summary_mode = "unstable" if "unstable" in final_all["rt_mode"].unique() else str(final_all["rt_mode"].iloc[0])
    legacy_prefix = f"cellular_rt_{summary_mode}"
    result.update(
        {
            "cellular_unstable": bool(result.get(f"{legacy_prefix}_unstable", False)),
            "cellular_max_growth_1_s": float(result.get(f"{legacy_prefix}_max_growth_1_s", np.nan)),
            "cellular_k_delta_at_max": float(result.get(f"{legacy_prefix}_k_delta_at_max", np.nan)),
            "cellular_lambda_over_delta_at_max": float(result.get(f"{legacy_prefix}_lambda_over_delta_at_max", np.nan)),
            "cellular_mode_n_at_max": float(result.get(f"{legacy_prefix}_mode_n_at_max", np.nan)),
            "cellular_cutoff_k_delta": float(result.get(f"{legacy_prefix}_cutoff_k_delta", np.nan)),
            "cellular_rt_mode": summary_mode,
            "cellular_gravity_effective_m_s2": float(result.get(f"{legacy_prefix}_gravity_effective_m_s2", np.nan)),
        }
    )
    return result
