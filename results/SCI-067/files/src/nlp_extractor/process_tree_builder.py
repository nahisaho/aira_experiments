"""
NLP-based Process Tree Constructor for LCA Automation.

Extracts material flows, energy inputs, and process steps from
unstructured text (patents, technical reports, BOM documents)
and constructs a structured process tree compatible with Brightway2.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass, field, asdict
from typing import Optional


# ---------------------------------------------------------------------------
# Data Models
# ---------------------------------------------------------------------------
@dataclass
class Flow:
    """A single material or energy flow."""
    name: str
    amount: float
    unit: str
    flow_type: str  # "input" | "output" | "emission"
    category: str = ""
    uncertainty_type: str = "lognormal"
    uncertainty_params: dict = field(default_factory=lambda: {"loc": 0, "scale": 0.1})


@dataclass
class ProcessNode:
    """A unit process in the LCA process tree."""
    id: str
    name: str
    location: str = "GLO"
    inputs: list[Flow] = field(default_factory=list)
    outputs: list[Flow] = field(default_factory=list)
    emissions: list[Flow] = field(default_factory=list)
    children: list[str] = field(default_factory=list)  # child process IDs
    metadata: dict = field(default_factory=dict)


@dataclass
class ProcessTree:
    """Complete process tree for an LCA study."""
    product_name: str
    functional_unit: str
    root_process_id: str
    processes: dict[str, ProcessNode] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "product_name": self.product_name,
            "functional_unit": self.functional_unit,
            "root_process_id": self.root_process_id,
            "processes": {k: asdict(v) for k, v in self.processes.items()},
        }


# ---------------------------------------------------------------------------
# NLP Extraction Patterns
# ---------------------------------------------------------------------------
# Regex patterns for common material/energy mentions in technical documents
QUANTITY_PATTERN = re.compile(
    r"(\d+\.?\d*)\s*(kg|g|mg|t|kWh|MJ|GJ|L|m3|m²|m³|kW|MW)\b"
    r"\s+(?:of\s+)?([A-Za-z][A-Za-z0-9\s\-/]{2,40})",
    re.IGNORECASE,
)

PROCESS_STEP_PATTERN = re.compile(
    r"(?:step|stage|phase|process)\s*\d*[:\-\s]+([A-Za-z][A-Za-z\s\-]{3,60})",
    re.IGNORECASE,
)

EMISSION_PATTERN = re.compile(
    r"(?:emit|release|discharge|emission)[s]?\s+(?:of\s+)?"
    r"(\d+\.?\d*)\s*(kg|g|mg|t)\s+(?:of\s+)?([A-Za-z][A-Za-z0-9\s\-]{2,40})",
    re.IGNORECASE,
)


class NLPProcessTreeBuilder:
    """
    Builds LCA process trees from unstructured text using NLP.

    Architecture:
    1. Named Entity Recognition (NER) for materials, chemicals, energy carriers
    2. Relation Extraction for process-flow connections
    3. Quantity Extraction for amounts and units
    4. Graph Construction for process tree assembly

    In production, integrates with:
    - spaCy + custom NER model trained on LCA corpus
    - Sentence-BERT for semantic similarity matching
    - LLM (GPT-4 / Claude) for complex extraction via structured prompts
    """

    def __init__(self, use_llm: bool = False, model_name: str = "en_core_web_trf"):
        self.use_llm = use_llm
        self.model_name = model_name
        self._process_counter = 0

    def _next_id(self, prefix: str = "P") -> str:
        self._process_counter += 1
        return f"{prefix}_{self._process_counter:04d}"

    def extract_quantities(self, text: str) -> list[dict]:
        """Extract quantity-material pairs from text."""
        results = []
        for match in QUANTITY_PATTERN.finditer(text):
            results.append({
                "amount": float(match.group(1)),
                "unit": match.group(2),
                "material": match.group(3).strip(),
            })
        return results

    def extract_process_steps(self, text: str) -> list[str]:
        """Extract named process steps from text."""
        return [m.group(1).strip() for m in PROCESS_STEP_PATTERN.finditer(text)]

    def extract_emissions(self, text: str) -> list[dict]:
        """Extract emission flows from text."""
        results = []
        for match in EMISSION_PATTERN.finditer(text):
            results.append({
                "amount": float(match.group(1)),
                "unit": match.group(2),
                "substance": match.group(3).strip(),
            })
        return results

    def build_from_text(
        self, text: str, product_name: str, functional_unit: str
    ) -> ProcessTree:
        """
        Build a process tree from unstructured text description.

        Pipeline:
          text → tokenize → NER → relation extraction → graph assembly → ProcessTree
        """
        steps = self.extract_process_steps(text)
        quantities = self.extract_quantities(text)
        emissions = self.extract_emissions(text)

        tree = ProcessTree(
            product_name=product_name,
            functional_unit=functional_unit,
            root_process_id="",
        )

        # Create process nodes for each identified step
        prev_id: Optional[str] = None
        for step_name in steps:
            pid = self._next_id()
            node = ProcessNode(id=pid, name=step_name)
            tree.processes[pid] = node
            if prev_id:
                tree.processes[prev_id].children.append(pid)
            if tree.root_process_id == "":
                tree.root_process_id = pid
            prev_id = pid

        # Distribute extracted quantities across processes
        if tree.processes:
            process_list = list(tree.processes.values())
            for i, q in enumerate(quantities):
                target = process_list[i % len(process_list)]
                target.inputs.append(
                    Flow(
                        name=q["material"],
                        amount=q["amount"],
                        unit=q["unit"],
                        flow_type="input",
                    )
                )

            # Attach emissions
            for em in emissions:
                process_list[-1].emissions.append(
                    Flow(
                        name=em["substance"],
                        amount=em["amount"],
                        unit=em["unit"],
                        flow_type="emission",
                    )
                )

        return tree

    def build_from_bom(self, bom: list[dict], product_name: str) -> ProcessTree:
        """
        Build process tree from a Bill of Materials (structured input).

        Each BOM entry: {"component": str, "material": str, "mass_kg": float,
                         "process": str, "supplier_location": str}
        """
        tree = ProcessTree(
            product_name=product_name,
            functional_unit=f"1 unit of {product_name}",
            root_process_id="",
        )

        assembly_id = self._next_id("ASM")
        assembly = ProcessNode(id=assembly_id, name=f"{product_name} Assembly")
        tree.root_process_id = assembly_id

        for entry in bom:
            pid = self._next_id("CMP")
            node = ProcessNode(
                id=pid,
                name=entry.get("process", f"Produce {entry['component']}"),
                location=entry.get("supplier_location", "GLO"),
                inputs=[
                    Flow(
                        name=entry["material"],
                        amount=entry["mass_kg"],
                        unit="kg",
                        flow_type="input",
                    )
                ],
                outputs=[
                    Flow(
                        name=entry["component"],
                        amount=entry["mass_kg"],
                        unit="kg",
                        flow_type="output",
                    )
                ],
            )
            tree.processes[pid] = node
            assembly.children.append(pid)

        tree.processes[assembly_id] = assembly
        return tree


# ---------------------------------------------------------------------------
# LLM-Assisted Extraction Prompt Template
# ---------------------------------------------------------------------------
LLM_EXTRACTION_PROMPT = """
You are an LCA (Life Cycle Assessment) expert. Extract structured process information
from the following technical description.

Return a JSON object with:
{{
  "processes": [
    {{
      "name": "process name",
      "inputs": [{{"material": "...", "amount": 0.0, "unit": "kg"}}],
      "outputs": [{{"product": "...", "amount": 0.0, "unit": "kg"}}],
      "emissions": [{{"substance": "...", "amount": 0.0, "unit": "kg"}}],
      "energy": [{{"type": "electricity|heat|fuel", "amount": 0.0, "unit": "kWh"}}]
    }}
  ],
  "process_connections": [
    {{"from": "process A", "to": "process B", "flow": "intermediate product"}}
  ]
}}

Text:
{text}
"""


def generate_llm_prompt(text: str) -> str:
    """Generate a structured extraction prompt for LLM-assisted process tree building."""
    return LLM_EXTRACTION_PROMPT.format(text=text)
