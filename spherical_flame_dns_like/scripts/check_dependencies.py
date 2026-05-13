from __future__ import annotations

import importlib.util
import os
import shutil


def main() -> None:
    checks = [
        ("PeleLMeX executable", _which_any(["PeleLMeX2d.gnu.MPI.ex", "PeleLMeX2d.gnu.ex", "PeleLMeX2d.ex"])),
        ("PELE_HOME", os.environ.get("PELE_HOME")),
        ("AMREX_HOME", os.environ.get("AMREX_HOME")),
        ("Cantera Python", _python_import("cantera")),
        ("h5py Python", _python_import("h5py")),
        ("yt Python", _python_import("yt")),
    ]
    for name, value in checks:
        status = "FOUND" if value and not value.startswith("BROKEN") else "MISSING"
        detail = f" - {value}" if value else ""
        print(f"{status}: {name}{detail}")


def _which_any(names: list[str]) -> str | None:
    for name in names:
        path = shutil.which(name)
        if path:
            return path
    return None


def _python_import(name: str) -> str | None:
    spec = importlib.util.find_spec(name)
    if spec is None:
        return None
    try:
        module = __import__(name)
    except Exception as exc:
        return f"BROKEN import at {spec.origin}: {exc}"
    return getattr(module, "__file__", spec.origin)


if __name__ == "__main__":
    main()
