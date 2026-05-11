from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.main import run_case
from src.utils import load_yaml


def main() -> None:
    parser = argparse.ArgumentParser(description="Run one spherical flame DNS-like case.")
    parser.add_argument("--case", default="H2_10", help="Case name from config/cases.yaml.")
    args = parser.parse_args()

    config = load_yaml(ROOT / "config" / "base.yaml")
    cases = load_yaml(ROOT / "config" / "cases.yaml")["cases"]
    selected = next((case for case in cases if case["name"] == args.case), None)
    if selected is None:
        raise SystemExit(f"Unknown case {args.case!r}.")
    summary = run_case(selected, config)
    print(summary)


if __name__ == "__main__":
    main()
