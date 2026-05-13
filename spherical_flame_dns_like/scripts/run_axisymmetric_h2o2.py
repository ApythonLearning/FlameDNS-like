from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.main import run_case
from src.utils import deep_update, ensure_dir, load_yaml


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the configured 2D axisymmetric premixed H2/O2 spherical-flame case."
    )
    parser.add_argument(
        "--config",
        default=str(ROOT / "config" / "axisymmetric_h2o2_spherical.yaml"),
        help="Case override YAML. Defaults to config/axisymmetric_h2o2_spherical.yaml.",
    )
    parser.add_argument(
        "--case",
        default=None,
        help="Optional case name inside the override YAML. Defaults to the first listed case.",
    )
    args = parser.parse_args()

    base_config = load_yaml(ROOT / "config" / "base.yaml")
    override = load_yaml(args.config)
    cases = override.pop("cases", None)
    config = deep_update(base_config, override)
    cases = cases or load_yaml(ROOT / "config" / "cases.yaml")["cases"]
    if not cases:
        raise SystemExit("No cases configured.")

    selected = cases[0] if args.case is None else next((case for case in cases if case["name"] == args.case), None)
    if selected is None:
        raise SystemExit(f"Unknown case {args.case!r}.")

    output_dir = ensure_dir(ROOT / config["project"]["output_dir"])
    summary = run_case(selected, config)
    (output_dir / "run_summary.txt").write_text(_format_summary(summary), encoding="utf-8")
    print(_format_summary(summary))


def _format_summary(summary: dict) -> str:
    keys = [
        "case",
        "h2_volume_fraction",
        "phi",
        "used_cantera",
        "status",
        "S_L_m_s",
        "S_b0_m_s",
        "delta_f_m",
        "points_per_flame_thickness",
        "grid_resolution_ok",
        "CFL",
        "CFL_ok",
        "Ma",
        "Pe",
    ]
    return "\n".join(f"{key}: {summary[key]}" for key in keys if key in summary)


if __name__ == "__main__":
    main()
