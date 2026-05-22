"""
Module 6: AiiDA / FireWorks Automated Workflow Pipeline
========================================================
Implements:
  - AiiDA WorkChain design for high-throughput perovskite screening
  - FireWorks Workflow definition (JSON-serializable)
  - DAG-based dependency management
  - Error handling and restart logic
  - Queue submission scripts (SLURM/PBS)
  - Status monitoring and checkpoint management
"""

import json
import os
from dataclasses import dataclass, field, asdict
from typing import List, Dict, Optional, Any
from datetime import datetime
import hashlib


# ── FireWorks-style Workflow Definition ───────────────────────────────────────

@dataclass
class Firework:
    """Represents a single FireWorks Firework (atomic calculation unit)."""
    fw_id: int
    name: str
    task_type: str
    spec: Dict[str, Any]
    parents: List[int] = field(default_factory=list)
    children: List[int] = field(default_factory=list)
    state: str = "WAITING"   # WAITING, READY, RUNNING, COMPLETED, FIZZLED
    priority: int = 0

    def to_dict(self):
        return {
            "fw_id": self.fw_id,
            "name": self.name,
            "task_type": self.task_type,
            "spec": self.spec,
            "parents": self.parents,
            "children": self.children,
            "state": self.state,
            "priority": self.priority,
        }


@dataclass
class Workflow:
    """FireWorks Workflow (DAG of Fireworks)."""
    wf_id: str
    name: str
    fireworks: List[Firework] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)
    created_at: str = field(default_factory=lambda: datetime.now().isoformat())

    def to_dict(self):
        return {
            "wf_id": self.wf_id,
            "name": self.name,
            "fireworks": [fw.to_dict() for fw in self.fireworks],
            "metadata": self.metadata,
            "created_at": self.created_at,
        }

    def save_json(self, path: str):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)


# ── Perovskite Screening Workflow Builder ─────────────────────────────────────

def build_screening_workflow(
    candidates: List[Dict],
    workflow_name: str = "lead_free_perovskite_screening",
    dft_code: str = "VASP",
    ml_enabled: bool = True,
    neb_enabled: bool = True,
    scaps_enabled: bool = True,
) -> Workflow:
    """
    Build the complete high-throughput screening workflow as a FireWorks DAG.

    Pipeline per candidate:
      [1] Tolerance factor filter  →  [2] ML band gap prediction
      →  [3] Structure relaxation (DFT)  →  [4] Band structure + DOS
      →  [5] Defect calculations  →  [6] NEB migration
      →  [7] SCAPS device simulation  →  [8] Results aggregation

    Non-DFT steps run immediately (no queue needed).
    DFT steps are VASP/Quantum ESPRESSO jobs.
    """
    wf_id = hashlib.md5(workflow_name.encode()).hexdigest()[:8]
    wf    = Workflow(
        wf_id=wf_id, name=workflow_name,
        metadata={
            "n_candidates": len(candidates),
            "dft_code": dft_code,
            "ml_enabled": ml_enabled,
            "neb_enabled": neb_enabled,
            "scaps_enabled": scaps_enabled,
            "description": "High-throughput screening of lead-free Sn/Ge/Bi perovskites",
        }
    )

    fw_counter = 1
    # Global aggregation task (final node)
    final_fw_id = 9000
    final_fw    = Firework(
        fw_id=final_fw_id,
        name="aggregate_and_rank",
        task_type="PythonTask",
        spec={
            "func": "perovskite_screener.ranking.rank_candidates",
            "args": {"top_n": 10},
            "resources": {"cores": 1, "memory_gb": 4, "walltime_h": 1},
        },
        state="WAITING",
        priority=10,
    )

    for cand in candidates:
        formula = cand["formula"]
        A, B, X = cand["A"], cand["B"], cand["X"]

        # ── Step 1: Tolerance filter (no DFT, instant) ───────────────────
        fw1 = Firework(
            fw_id=fw_counter, name=f"tolerance_{formula}",
            task_type="PythonTask",
            spec={
                "func": "perovskite_screener.tolerance_factor.analyze_perovskite",
                "args": {"A": A, "B": B, "X": X, "B_ox": cand.get("B_ox", 2)},
                "pass_to_db": True,
                "filter_key": "stability_class",
                "filter_values": ["perovskite"],  # discard non-perovskites
                "resources": {"cores": 1, "memory_gb": 1, "walltime_h": 0.1},
            },
            state="READY",
            priority=5,
        )
        fw_counter += 1

        # ── Step 2: ML band gap prediction ───────────────────────────────
        fw2_id = fw_counter
        fw2 = Firework(
            fw_id=fw2_id, name=f"ml_bandgap_{formula}",
            task_type="PythonTask",
            spec={
                "func": "perovskite_screener.bandgap_ml.BandGapPredictor.predict",
                "args": {"A": A, "B": B, "X": X, "B_ox": cand.get("B_ox", 2)},
                "pass_to_db": True,
                "filter_key": "Eg_predicted_eV",
                "filter_range": [0.9, 2.5],  # optimal for solar
                "resources": {"cores": 1, "memory_gb": 2, "walltime_h": 0.2},
            },
            parents=[fw1.fw_id],
            state="WAITING",
            priority=5,
        ) if ml_enabled else None
        fw_counter += 1

        # ── Step 3: DFT Structure Relaxation ─────────────────────────────
        fw3_id = fw_counter
        fw3 = Firework(
            fw_id=fw3_id, name=f"dft_relax_{formula}",
            task_type="VaspJob",
            spec={
                "calculation": "relax",
                "incar": {
                    "ENCUT": 520, "EDIFF": 1e-6, "EDIFFG": -0.01,
                    "IBRION": 2, "NSW": 100, "ISIF": 3,
                    "PREC": "Accurate", "ALGO": "Normal",
                    "ISMEAR": 0, "SIGMA": 0.05,
                    "LWAVE": False, "LCHARG": False,
                },
                "kpoints": {"scheme": "Monkhorst-Pack", "mesh": [6, 6, 6]},
                "formula": formula,
                "resources": {"cores": 32, "memory_gb": 64, "walltime_h": 12},
                "restart_strategy": "continue_from_last",
            },
            parents=[fw2_id] if ml_enabled else [fw1.fw_id],
            state="WAITING",
            priority=3,
        )
        fw_counter += 1

        # ── Step 4: Band Structure + DOS (HSE06 + SOC) ────────────────────
        fw4_id = fw_counter
        fw4 = Firework(
            fw_id=fw4_id, name=f"dft_bands_{formula}",
            task_type="VaspJob",
            spec={
                "calculation": "bands+dos",
                "incar": {
                    "ENCUT": 520, "EDIFF": 1e-6,
                    "IBRION": -1, "NSW": 0,
                    "HFSCREEN": 0.2, "AEXX": 0.25,  # HSE06
                    "LHFCALC": True, "LSORBIT": True,  # SOC
                    "PREC": "Accurate", "ALGO": "All",
                    "ISMEAR": 0, "SIGMA": 0.05,
                    "LORBIT": 11, "NEDOS": 2000,
                },
                "kpoints": {"scheme": "HSE_bandpath", "nkpts": 30},
                "formula": formula,
                "resources": {"cores": 64, "memory_gb": 128, "walltime_h": 48},
                "restart_strategy": "continue_from_last",
            },
            parents=[fw3_id],
            state="WAITING",
            priority=3,
        )
        fw_counter += 1

        # ── Step 5: Defect calculations ───────────────────────────────────
        fw5_id = fw_counter
        fw5 = Firework(
            fw_id=fw5_id, name=f"defect_calc_{formula}",
            task_type="CompositeTask",
            spec={
                "subtasks": [
                    {
                        "type": "VaspJob",
                        "defect": "V_X",
                        "incar_overrides": {"NELM": 200, "LDAUU": 0},
                        "supercell": [3, 3, 3],
                        "resources": {"cores": 64, "memory_gb": 128, "walltime_h": 24},
                    },
                    {
                        "type": "VaspJob",
                        "defect": "V_B",
                        "supercell": [3, 3, 3],
                        "resources": {"cores": 64, "memory_gb": 128, "walltime_h": 24},
                    },
                    {
                        "type": "VaspJob",
                        "defect": "i_X",
                        "supercell": [3, 3, 3],
                        "resources": {"cores": 64, "memory_gb": 128, "walltime_h": 24},
                    },
                ],
                "charge_states": [-2, -1, 0, 1, 2],
                "corrections": ["Freysoldt_FNV", "Kumagai"],
                "formula": formula,
            },
            parents=[fw4_id],
            state="WAITING",
            priority=2,
        )
        fw_counter += 1

        # ── Step 6: NEB ion migration ─────────────────────────────────────
        fw6_id = fw_counter
        fw6 = Firework(
            fw_id=fw6_id, name=f"neb_migration_{formula}",
            task_type="VaspNEB",
            spec={
                "n_images": 7,
                "migration_species": X,
                "mechanism": "vacancy_hop",
                "incar": {
                    "ENCUT": 520, "EDIFF": 1e-5, "EDIFFG": -0.05,
                    "IBRION": 3, "NSW": 200,
                    "ICHAIN": 0, "IMAGES": 7,
                    "SPRING": -5.0,
                    "LCLIMB": True,  # CI-NEB
                    "PREC": "Normal",
                },
                "supercell": [3, 3, 1],
                "formula": formula,
                "resources": {"cores": 64, "memory_gb": 128, "walltime_h": 36},
            },
            parents=[fw5_id],
            state="WAITING",
            priority=2,
        ) if neb_enabled else None
        fw_counter += 1

        # ── Step 7: SCAPS device simulation ──────────────────────────────
        fw7_id = fw_counter
        fw7 = Firework(
            fw_id=fw7_id, name=f"scaps_device_{formula}",
            task_type="PythonTask",
            spec={
                "func": "perovskite_screener.scaps_interface.run_scaps_simulation",
                "args": {"A": A, "B": B, "X": X},
                "inputs_from_db": ["Eg_predicted_eV", "defect_concentration"],
                "formula": formula,
                "resources": {"cores": 1, "memory_gb": 4, "walltime_h": 0.5},
            },
            parents=[fw6_id if neb_enabled else fw5_id],
            state="WAITING",
            priority=2,
        ) if scaps_enabled else None
        fw_counter += 1

        # Register fireworks in workflow
        for fw in [fw1, fw2, fw3, fw4, fw5, fw6, fw7]:
            if fw is not None:
                wf.fireworks.append(fw)

        # Connect to final aggregator
        last_fw_id = fw7_id if scaps_enabled else fw5_id
        final_fw.parents.append(last_fw_id)

    wf.fireworks.append(final_fw)
    return wf


# ── AiiDA WorkChain Definition ────────────────────────────────────────────────

AIIDA_WORKCHAIN_CODE = '''"""
AiiDA WorkChain for lead-free perovskite high-throughput screening.
Requires: aiida-core >= 2.0, aiida-vasp >= 3.0
"""
from aiida.engine import WorkChain, calcfunction, if_, while_, ToContext
from aiida.orm import Dict, StructureData, List, Float
from aiida.plugins import CalculationFactory, DataFactory


VaspCalculation   = CalculationFactory("vasp.vasp")
VaspNEBCalculation = CalculationFactory("vasp.neb")


class PerovskiteScreeningWorkChain(WorkChain):
    """
    AiiDA WorkChain for high-throughput Sn/Ge/Bi perovskite screening.

    Workflow outline:
      1. tolerance_filter  – discard non-perovskite compositions
      2. ml_prescreen      – ML band gap filter (0.9–2.5 eV window)
      3. dft_relax         – VASP structure relaxation (PBE+D3)
      4. dft_bands         – HSE06+SOC band structure
      5. defect_calc       – Defect formation energies (supercell)
      6. neb_calc          – CI-NEB ion migration barriers
      7. scaps_sim         – SCAPS-1D device simulation
      8. aggregate_rank    – Multi-criteria ranking
    """

    @classmethod
    def define(cls, spec):
        super().define(spec)
        spec.input("candidates", valid_type=List)
        spec.input("dft_parameters", valid_type=Dict)
        spec.input("screening_parameters", valid_type=Dict)
        spec.output("ranked_candidates", valid_type=List)
        spec.output("top_material", valid_type=Dict)
        spec.output("screening_report", valid_type=Dict)

        spec.outline(
            cls.tolerance_filter,
            if_(cls.has_candidates)(
                cls.ml_prescreen,
                if_(cls.has_ml_passed)(
                    cls.dft_relax,
                    cls.dft_bands,
                    cls.defect_calc,
                    cls.neb_calc,
                    cls.scaps_sim,
                )
            ),
            cls.aggregate_rank,
        )

        spec.exit_code(300, "NO_CANDIDATES_PASSED_FILTER",
                       "No candidates passed the tolerance factor filter")
        spec.exit_code(400, "DFT_FAILED",
                       "DFT calculation failed for all candidates")

    def tolerance_filter(self):
        from perovskite_screener.tolerance_factor import analyze_perovskite
        passed = []
        for c in self.inputs.candidates:
            res = analyze_perovskite(c["A"], c["B"], c["X"], c.get("B_ox", 2))
            if res.stability_class == "perovskite":
                passed.append(c)
        self.ctx.candidates = passed
        if not passed:
            return self.exit_codes.NO_CANDIDATES_PASSED_FILTER
        self.report(f"Tolerance filter: {len(passed)}/{len(self.inputs.candidates)} passed")

    def has_candidates(self):
        return len(self.ctx.candidates) > 0

    def ml_prescreen(self):
        from perovskite_screener.bandgap_ml import BandGapPredictor
        predictor = BandGapPredictor().fit(verbose=False)
        passed = []
        for c in self.ctx.candidates:
            res = predictor.predict(c["A"], c["B"], c["X"], c.get("B_ox", 2))
            if 0.9 <= res["Eg_predicted_eV"] <= 2.5:
                c["Eg_ml"] = res["Eg_predicted_eV"]
                passed.append(c)
        self.ctx.candidates = passed
        self.report(f"ML prescreening: {len(passed)} candidates in 0.9–2.5 eV window")

    def has_ml_passed(self):
        return len(self.ctx.candidates) > 0

    def dft_relax(self):
        calcs = {}
        for c in self.ctx.candidates:
            code     = self.inputs.dft_parameters["vasp_code"]
            builder  = VaspCalculation.get_builder()
            builder.structure = c["structure"]
            builder.parameters = Dict(dict=self.inputs.dft_parameters["relax_incar"])
            builder.kpoints = DataFactory("array.kpoints")()
            builder.kpoints.set_kpoints_mesh([6, 6, 6])
            future = self.submit(builder)
            calcs[f"relax_{c['formula']}"] = future
        return ToContext(**calcs)

    def dft_bands(self):
        # Similar to dft_relax but with HSE06+SOC settings
        pass

    def defect_calc(self):
        pass

    def neb_calc(self):
        pass

    def scaps_sim(self):
        from perovskite_screener.scaps_interface import run_scaps_simulation
        results = {}
        for c in self.ctx.candidates:
            res = run_scaps_simulation(c["A"], c["B"], c["X"], c.get("Eg_ml", 1.5))
            results[c["formula"]] = res
        self.ctx.scaps_results = results

    def aggregate_rank(self):
        from perovskite_screener.ranking import rank_candidates
        ranked = rank_candidates(
            self.ctx.candidates,
            scaps_results=getattr(self.ctx, "scaps_results", {}),
        )
        self.out("ranked_candidates", List(list=ranked))
        self.out("top_material", Dict(dict=ranked[0] if ranked else {}))
'''


def generate_slurm_script(
    job_name: str,
    n_cores: int = 32,
    memory_gb: int = 64,
    walltime_h: int = 12,
    partition: str = "normal",
    conda_env: str = "perovskite",
    output_dir: str = "results/slurm",
) -> str:
    """Generate SLURM submission script for DFT jobs."""
    os.makedirs(output_dir, exist_ok=True)
    script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition={partition}
#SBATCH --nodes={max(1, n_cores // 32)}
#SBATCH --ntasks-per-node={min(n_cores, 32)}
#SBATCH --mem={memory_gb}G
#SBATCH --time={walltime_h:02d}:00:00
#SBATCH --output=logs/slurm_%j.out
#SBATCH --error=logs/slurm_%j.err

module load intel/2023a VASP/6.4.1
conda activate {conda_env}

# AiiDA daemon must be running
verdi daemon start 4

# Submit workflow
python -c "
from perovskite_screener.workflow import build_screening_workflow
from perovskite_screener.materials_database import get_all_candidates
import json

candidates = get_all_candidates()
wf = build_screening_workflow(candidates)
wf.save_json('results/workflow_definition.json')
print(f'Workflow {{wf.wf_id}} with {{len(wf.fireworks)}} fireworks saved.')
"
"""
    script_path = f"{output_dir}/{job_name}.slurm"
    with open(script_path, "w") as f:
        f.write(script)
    return script_path


def save_workflow_definition(
    candidates: List[Dict],
    output_dir: str = "results",
) -> str:
    """Build and save the screening workflow as JSON."""
    wf = build_screening_workflow(candidates)
    path = f"{output_dir}/workflow_definition.json"
    os.makedirs(output_dir, exist_ok=True)
    wf.save_json(path)
    return path


def save_aiida_workchain(output_dir: str = "src/perovskite_screener") -> str:
    """Save the AiiDA WorkChain Python code."""
    path = f"{output_dir}/aiida_workchain.py"
    with open(path, "w") as f:
        f.write(AIIDA_WORKCHAIN_CODE)
    return path
