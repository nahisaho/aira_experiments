from __future__ import annotations

import sys
sys.path.insert(0, '/home/nahisaho/GitHub/aira/projects/f78a410a-e891-4bfe-88c5-74583c82bef3/workspace')

import json
import logging
from datetime import datetime
from pathlib import Path


def setup_logger(name: str, log_file: str | None = None, level: int = logging.INFO) -> logging.Logger:
    """Create a console/file logger without duplicate handlers."""
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.propagate = False
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    if not any(isinstance(handler, logging.StreamHandler) and not isinstance(handler, logging.FileHandler) for handler in logger.handlers):
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

    if log_file:
        file_path = Path(log_file).resolve()
        if not any(isinstance(handler, logging.FileHandler) and Path(handler.baseFilename).resolve() == file_path for handler in logger.handlers):
            file_path.parent.mkdir(parents=True, exist_ok=True)
            file_handler = logging.FileHandler(file_path)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)

    return logger


class ProcessLogger:
    """Structured JSONL process logger."""

    def __init__(self, log_path: str = 'logs/process-log.jsonl') -> None:
        self.log_path = Path(log_path)
        self.log_path.parent.mkdir(parents=True, exist_ok=True)

    def log(
        self,
        phase: str,
        event_type: str,
        skill_or_tool: str = '',
        handoff_in: dict | None = None,
        handoff_out: dict | None = None,
        files_written: list[str] | None = None,
        status: str = 'ok',
    ) -> None:
        """Append a JSONL process log entry."""
        entry = {
            'timestamp': datetime.now().isoformat(),
            'phase': phase,
            'event_type': event_type,
            'actor': 'co-scientist',
            'skill_or_tool': skill_or_tool,
            'handoff_in': handoff_in or {},
            'handoff_out': handoff_out or {},
            'files_written': files_written or [],
            'status': status,
        }
        with self.log_path.open('a', encoding='utf-8') as file:
            file.write(json.dumps(entry, ensure_ascii=False) + '\\n')
