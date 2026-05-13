from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Check DNS mesh and CFL requirements.")
    parser.add_argument("--summary", required=True, help="Cantera 1D summary JSON containing delta_T_m.")
    parser.add_argument("--dx", type=float, required=True)
    parser.add_argument("--cfl", type=float, required=True)
    args = parser.parse_args()

    summary = json.loads(Path(args.summary).read_text(encoding="utf-8"))
    delta = float(summary["delta_T_m"])
    required = delta / 10.0
    recommended = delta / 20.0
    print(f"delta_T_m: {delta:.8e}")
    print(f"dx_m: {args.dx:.8e}")
    print(f"dx <= delta_T/10: {args.dx <= required}")
    print(f"dx <= delta_T/20 recommended: {args.dx <= recommended}")
    print(f"CFL < 0.5: {args.cfl < 0.5}")
    if args.dx > required:
        print("WARNING: grid is under-resolved for DNS criterion dx <= delta_T/10.")
    elif args.dx > recommended:
        print("WARNING: grid passes dx <= delta_T/10 but not recommended dx <= delta_T/20.")
    if args.cfl >= 0.5:
        print("WARNING: CFL should be below 0.5.")


if __name__ == "__main__":
    main()
