from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.main import run_cases
from src.utils import ensure_dir, load_yaml


def main() -> None:
    config = load_yaml(ROOT / "config" / "base.yaml")
    cases = load_yaml(ROOT / "config" / "cases.yaml")["cases"]
    output_dir = ensure_dir(ROOT / config["project"]["output_dir"])
    summary = run_cases(cases, config, output_dir)
    columns = [
        "case",
        "h2_volume_fraction",
        "used_cantera",
        "S_L_m_s",
        "S_b0_m_s",
        "L_b_m",
        "Ma",
        "Pe",
        "cellular_Le_eff",
        "cellular_rt_unstable_max_growth_1_s",
        "cellular_rt_neutral_max_growth_1_s",
        "cellular_rt_stable_max_growth_1_s",
        "nonlinear_rt_unstable_phase_final",
        "nonlinear_rt_unstable_onset_time_s",
        "nonlinear_rt_unstable_cellular_speed_factor_final",
    ]
    print(summary[[column for column in columns if column in summary.columns]])


if __name__ == "__main__":
    main()
