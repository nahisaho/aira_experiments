from __future__ import annotations

import json
import random
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
FIGURES_DIR = ROOT / "figures"
RESULTS_DIR = ROOT / "results"
DATA_DIR = ROOT / "data"
LOGS_DIR = ROOT / "logs"
PROCESS_LOG = LOGS_DIR / "process-log.jsonl"
DEFAULT_SEED = 42


def ensure_output_dirs() -> None:
    for path in (FIGURES_DIR, RESULTS_DIR, DATA_DIR, LOGS_DIR):
        path.mkdir(parents=True, exist_ok=True)
    PROCESS_LOG.touch(exist_ok=True)


def seed_everything(seed: int = DEFAULT_SEED) -> None:
    random.seed(seed)
    np.random.seed(seed)


def _json_safe(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_json_safe(v) for v in value]
    return value


def log_event(
    phase: str,
    event_type: str,
    skill_or_tool: str,
    handoff_in: dict[str, Any] | None = None,
    handoff_out: dict[str, Any] | None = None,
    files_written: list[str] | None = None,
    status: str = "ok",
    actor: str = "co-scientist",
) -> None:
    ensure_output_dirs()
    record = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "phase": phase,
        "event_type": event_type,
        "actor": actor,
        "skill_or_tool": skill_or_tool,
        "handoff_in": _json_safe(handoff_in or {}),
        "handoff_out": _json_safe(handoff_out or {}),
        "files_written": _json_safe(files_written or []),
        "status": status,
    }
    with PROCESS_LOG.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, ensure_ascii=False) + "\n")


def save_json(path: Path, payload: dict[str, Any] | list[Any]) -> None:
    ensure_output_dirs()
    with path.open("w", encoding="utf-8") as handle:
        json.dump(_json_safe(payload), handle, indent=2, ensure_ascii=False)
    log_event(
        phase="io",
        event_type="file_written",
        skill_or_tool="save_json",
        files_written=[str(path.relative_to(ROOT))],
        handoff_out={"path": str(path.relative_to(ROOT))},
    )


ensure_output_dirs()
seed_everything(DEFAULT_SEED)
