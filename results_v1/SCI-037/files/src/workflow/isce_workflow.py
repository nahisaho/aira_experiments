"""ISCE-based automated workflow orchestration."""

from __future__ import annotations

import argparse
import json
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any



def _simple_yaml_load(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    data: dict[str, Any] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if ":" not in stripped:
            continue
        key, value = stripped.split(":", 1)
        value = value.strip()
        if value.lower() in {"true", "false"}:
            parsed: Any = value.lower() == "true"
        else:
            try:
                parsed = json.loads(value)
            except Exception:
                parsed = value.strip('"\'')
        data[key.strip()] = parsed
    return data


@dataclass
class ISCEWorkflow:
    """Simple step-by-step ISCE workflow wrapper."""

    config: dict[str, Any]
    dry_run: bool = True

    @classmethod
    def from_yaml(cls, path: str | Path, dry_run: bool = True) -> "ISCEWorkflow":
        return cls(config=_simple_yaml_load(path), dry_run=dry_run)

    def _run_command(self, command: list[str]) -> dict[str, Any]:
        if self.dry_run:
            return {"command": command, "returncode": 0, "stdout": "DRY RUN", "stderr": ""}
        result = subprocess.run(command, check=False, capture_output=True, text=True)
        return {"command": command, "returncode": result.returncode, "stdout": result.stdout, "stderr": result.stderr}

    def preprocess_slc(self) -> dict[str, Any]:
        return self._run_command(["isce2_slc_preprocess", str(self.config.get("slc_dir", "./slc")), str(self.config.get("multilook", 4))])

    def coregister_esd(self) -> dict[str, Any]:
        return self._run_command(["isce2_coregister", "--method", "ESD", str(self.config.get("master", 0))])

    def generate_interferograms(self) -> dict[str, Any]:
        return self._run_command(["isce2_interferogram", str(self.config.get("pairing", "small-baseline"))])

    def goldstein_filter(self) -> dict[str, Any]:
        return self._run_command(["isce2_filter", "--type", "goldstein", "--alpha", str(self.config.get("goldstein_alpha", 0.8))])

    def unwrap_phase(self) -> dict[str, Any]:
        return self._run_command(["isce2_unwrap", "--engine", str(self.config.get("unwrap_engine", "snaphu"))])

    def geocode(self) -> dict[str, Any]:
        return self._run_command(["isce2_geocode", str(self.config.get("dem", "./dem.tif"))])

    def run_all(self) -> list[dict[str, Any]]:
        steps = [self.preprocess_slc, self.coregister_esd, self.generate_interferograms, self.goldstein_filter, self.unwrap_phase, self.geocode]
        return [step() for step in steps]


def main() -> None:
    parser = argparse.ArgumentParser(description="Run an ISCE-based InSAR workflow")
    parser.add_argument("config", help="Path to workflow YAML configuration")
    parser.add_argument("--execute", action="store_true", help="Run commands instead of dry-run mode")
    args = parser.parse_args()
    workflow = ISCEWorkflow.from_yaml(args.config, dry_run=not args.execute)
    for item in workflow.run_all():
        print(json.dumps(item))


if __name__ == "__main__":
    main()
