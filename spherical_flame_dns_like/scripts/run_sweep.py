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
    print(summary[["case", "h2_volume_fraction", "used_cantera", "S_L_m_s", "S_b0_m_s", "L_b_m", "Ma", "Pe"]])


if __name__ == "__main__":
    main()
