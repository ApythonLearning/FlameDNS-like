from __future__ import annotations


def lewis_number_mode(config: dict) -> str:
    return "unity" if config["transport"].get("unity_lewis", True) else "non-unity"


def placeholder_diffusivity_m2_s(temperature_K: float) -> float:
    return 2.0e-5 * (float(temperature_K) / 298.0) ** 1.7
