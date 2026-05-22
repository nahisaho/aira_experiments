"""StaMPS workflow helpers and MATLAB command wrappers."""

from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass
class StaMPSWorkflow:
    """Wrap StaMPS processing stages and SNAP-to-StaMPS preparation."""

    workdir: str
    matlab_executable: str = "matlab"
    dry_run: bool = True

    def _run_matlab(self, expression: str) -> dict[str, Any]:
        command = [self.matlab_executable, "-batch", expression]
        if self.dry_run:
            return {"command": command, "returncode": 0, "stdout": "DRY RUN", "stderr": ""}
        result = subprocess.run(command, check=False, capture_output=True, text=True, cwd=self.workdir)
        return {"command": command, "returncode": result.returncode, "stdout": result.stdout, "stderr": result.stderr}

    def mt_prep_snap(self, snap_export_dir: str, amplitude_dispersion_threshold: float = 0.4) -> dict[str, Any]:
        expression = f"mt_prep_snap('{snap_export_dir}', {amplitude_dispersion_threshold});"
        return self._run_matlab(expression)

    def export_to_stamps(self) -> dict[str, Any]:
        return self._run_matlab("stamps_export;")

    def run_stage(self, stage: int) -> dict[str, Any]:
        if not 1 <= stage <= 8:
            raise ValueError("StaMPS stage must be between 1 and 8")
        return self._run_matlab(f"stamps({stage},1);")

    def postprocess(self) -> list[dict[str, Any]]:
        return [self._run_matlab("ps_plot('v-do',1);") , self._run_matlab("ps_plot('ts',1);")]

    def generate_velocity_map(self, output_path: str | Path) -> dict[str, Any]:
        expression = f"save('{Path(output_path)}','ph_disp','ph_mm');"
        return self._run_matlab(expression)

    def run_all(self, snap_export_dir: str) -> list[dict[str, Any]]:
        results = [self.mt_prep_snap(snap_export_dir), self.export_to_stamps()]
        results.extend(self.run_stage(stage) for stage in range(1, 9))
        results.extend(self.postprocess())
        return results


if __name__ == "__main__":
    workflow = StaMPSWorkflow(workdir=".", dry_run=True)
    print(json.dumps(workflow.run_all("./snap_export"), indent=2))
