from __future__ import annotations

from pathlib import Path
from typing import Any

try:
    import yaml
except Exception:  # pragma: no cover - exercised only when PyYAML is absent.
    yaml = None


def project_root() -> Path:
    return Path(__file__).resolve().parents[1]


def load_yaml(path: str | Path) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as handle:
        text = handle.read()
    if yaml is not None:
        data = yaml.safe_load(text)
        return data or {}
    data = _load_project_yaml_subset(text)
    return data or {}


def _load_project_yaml_subset(text: str) -> dict[str, Any]:
    """Small YAML subset parser for this project's config files.

    The real dependency is PyYAML. This fallback keeps the framework runnable in
    bare Python environments for the simple nested dictionaries and case lists
    used in config/base.yaml and config/cases.yaml.
    """
    root: dict[str, Any] = {}
    stack: list[tuple[int, Any]] = [(-1, root)]
    current_list_item: dict[str, Any] | None = None

    for raw_line in text.splitlines():
        if not raw_line.strip() or raw_line.lstrip().startswith("#"):
            continue
        indent = len(raw_line) - len(raw_line.lstrip(" "))
        line = raw_line.strip()

        while stack and indent <= stack[-1][0]:
            stack.pop()
        parent = stack[-1][1]

        if line.startswith("- "):
            if not isinstance(parent, list):
                raise ValueError("Unsupported YAML list placement.")
            item_text = line[2:]
            item: dict[str, Any] = {}
            parent.append(item)
            current_list_item = item
            if item_text:
                key, value = item_text.split(":", 1)
                item[key.strip()] = _parse_scalar(value.strip())
            stack.append((indent, item))
            continue

        key, value = line.split(":", 1)
        key = key.strip()
        value = value.strip()
        if value == "":
            next_container: Any = [] if key == "cases" else {}
            parent[key] = next_container
            stack.append((indent, next_container))
            current_list_item = None
        else:
            target = current_list_item if isinstance(parent, list) and current_list_item is not None else parent
            target[key] = _parse_scalar(value)
    return root


def _parse_scalar(value: str) -> Any:
    if value.startswith("[") and value.endswith("]"):
        inner = value[1:-1].strip()
        if not inner:
            return []
        return [_parse_scalar(item.strip()) for item in inner.split(",")]
    lowered = value.lower()
    if lowered in {"true", "false"}:
        return lowered == "true"
    try:
        if any(char in value for char in [".", "e", "E"]):
            return float(value)
        return int(value)
    except ValueError:
        return value.strip("\"'")


def ensure_dir(path: str | Path) -> Path:
    out = Path(path)
    out.mkdir(parents=True, exist_ok=True)
    return out


def deep_update(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    result = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = deep_update(result[key], value)
        else:
            result[key] = value
    return result


def case_output_dir(config: dict[str, Any], case_name: str) -> Path:
    root = project_root()
    output_dir = root / config["project"]["output_dir"] / case_name
    return ensure_dir(output_dir)


def mechanism_file(config: dict[str, Any]) -> str:
    return str(config["gas"].get("mechanism", "h2o2.yaml"))
