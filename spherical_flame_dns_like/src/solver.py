from __future__ import annotations


class AxisymmetricLowMachSolver:
    """Reserved interface for the Stage 2 low-Mach reacting-flow solver."""

    def __init__(self, config: dict):
        self.config = config

    def step(self, state: dict, dt_s: float) -> dict:
        raise NotImplementedError(
            "Stage 2 solver is not implemented yet. Current workflow uses "
            "Cantera 1D flames plus an expanding-flame post-processing model."
        )


def gravity_source_terms_enabled(config: dict) -> bool:
    return bool(config.get("gravity", {}).get("enabled", False))
