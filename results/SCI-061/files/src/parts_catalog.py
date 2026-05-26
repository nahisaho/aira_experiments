"""
Parts Catalog Module — Promoter, RBS, Terminator, and CDS definitions
for the synthetic gene circuit design pipeline (SBOL-inspired).
"""

import json
from dataclasses import dataclass, field, asdict
from typing import List, Dict, Optional
from enum import Enum


class PartType(Enum):
    PROMOTER = "promoter"
    RBS = "rbs"
    CDS = "cds"
    TERMINATOR = "terminator"
    INSULATOR = "insulator"


@dataclass
class GeneticPart:
    name: str
    part_type: PartType
    sequence: str
    parameters: Dict[str, float] = field(default_factory=dict)
    sbol_uri: str = ""
    description: str = ""


# ---- Built-in parts library (curated from literature) ----

PROMOTER_CATALOG = [
    GeneticPart(
        name="pTac", part_type=PartType.PROMOTER,
        sequence="AATTGTGAGCGGATAACAATTGACATTGTGAGCGGATAACAAGATACTGAGCACA",
        parameters={"k_max": 2.8, "K_d": 50.0, "n": 2.0, "leak": 0.01},
        sbol_uri="https://synbiohub.org/public/igem/BBa_R0010",
        description="IPTG-inducible promoter, strong"
    ),
    GeneticPart(
        name="pTet", part_type=PartType.PROMOTER,
        sequence="TCCCTATCAGTGATAGAGATTGACATCCCTATCAGTGATAGAGATACTGAGCAC",
        parameters={"k_max": 3.2, "K_d": 40.0, "n": 2.5, "leak": 0.008},
        sbol_uri="https://synbiohub.org/public/igem/BBa_R0040",
        description="aTc-inducible promoter"
    ),
    GeneticPart(
        name="pLac", part_type=PartType.PROMOTER,
        sequence="AATTGTGAGCGGATAACAATTGACATTGTGAGCGGATAACAAGATACTGAGCACA",
        parameters={"k_max": 2.2, "K_d": 60.0, "n": 1.8, "leak": 0.015},
        sbol_uri="https://synbiohub.org/public/igem/BBa_R0011",
        description="Lac promoter, moderate strength"
    ),
    GeneticPart(
        name="pBAD", part_type=PartType.PROMOTER,
        sequence="ACATTGATTATTTGCACGGCGTCACACTTTGCTATGCCATAGCATTTTTATCCATAAG",
        parameters={"k_max": 3.5, "K_d": 35.0, "n": 2.2, "leak": 0.005},
        sbol_uri="https://synbiohub.org/public/igem/BBa_I0500",
        description="Arabinose-inducible promoter"
    ),
    GeneticPart(
        name="pLambda", part_type=PartType.PROMOTER,
        sequence="TAATACGACTCACTATAGGGAGACCACAACGGTTTCCCTCTAGAAATAATTTTG",
        parameters={"k_max": 4.0, "K_d": 25.0, "n": 3.0, "leak": 0.003},
        sbol_uri="https://synbiohub.org/public/igem/BBa_R0051",
        description="Lambda cI-regulated promoter, strong repression"
    ),
]

RBS_CATALOG = [
    GeneticPart(
        name="B0034", part_type=PartType.RBS,
        sequence="AAAGAGGAGAAA",
        parameters={"translation_rate": 1.0, "efficiency": 0.85},
        sbol_uri="https://synbiohub.org/public/igem/BBa_B0034",
        description="Strong RBS, community standard"
    ),
    GeneticPart(
        name="B0032", part_type=PartType.RBS,
        sequence="TCACACAGGAAAG",
        parameters={"translation_rate": 0.3, "efficiency": 0.45},
        sbol_uri="https://synbiohub.org/public/igem/BBa_B0032",
        description="Medium RBS"
    ),
    GeneticPart(
        name="B0031", part_type=PartType.RBS,
        sequence="TCACACAGGAAAGTACTAG",
        parameters={"translation_rate": 0.07, "efficiency": 0.15},
        sbol_uri="https://synbiohub.org/public/igem/BBa_B0031",
        description="Weak RBS"
    ),
]

TERMINATOR_CATALOG = [
    GeneticPart(
        name="B0015", part_type=PartType.TERMINATOR,
        sequence="CCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTT",
        parameters={"termination_efficiency": 0.98},
        sbol_uri="https://synbiohub.org/public/igem/BBa_B0015",
        description="Double terminator (B0010+B0012)"
    ),
    GeneticPart(
        name="B0010", part_type=PartType.TERMINATOR,
        sequence="CCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCG",
        parameters={"termination_efficiency": 0.92},
        sbol_uri="https://synbiohub.org/public/igem/BBa_B0010",
        description="T1 terminator from E. coli rrnB"
    ),
]

CDS_CATALOG = [
    GeneticPart(
        name="LacI", part_type=PartType.CDS,
        sequence="ATGAAACCAGTAACGTTATACGATGTCGCAGAGTATGCCGGTGTCTCTT",
        parameters={"degradation_rate": 0.0023, "maturation_time": 6.0},
        sbol_uri="https://synbiohub.org/public/igem/BBa_C0012",
        description="LacI repressor"
    ),
    GeneticPart(
        name="TetR", part_type=PartType.CDS,
        sequence="ATGTCTAGATTAGATAAAAGTAAAGTGATTAACAGCGCATTAGAGCTGCTT",
        parameters={"degradation_rate": 0.0023, "maturation_time": 6.0},
        sbol_uri="https://synbiohub.org/public/igem/BBa_C0040",
        description="TetR repressor"
    ),
    GeneticPart(
        name="cI", part_type=PartType.CDS,
        sequence="ATGAGCACAAAAAAGAAACCATTAACACAAGAGCAGCTTGAGGACGCACGT",
        parameters={"degradation_rate": 0.0023, "maturation_time": 6.0},
        sbol_uri="https://synbiohub.org/public/igem/BBa_C0051",
        description="Lambda cI repressor"
    ),
    GeneticPart(
        name="GFP", part_type=PartType.CDS,
        sequence="ATGCGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATG",
        parameters={"degradation_rate": 0.0023, "maturation_time": 12.0},
        sbol_uri="https://synbiohub.org/public/igem/BBa_E0040",
        description="Green fluorescent protein reporter"
    ),
]


class PartsCatalog:
    """Central registry for genetic parts with query capabilities."""

    def __init__(self):
        self.parts: Dict[str, GeneticPart] = {}
        self._load_defaults()

    def _load_defaults(self):
        for part in (PROMOTER_CATALOG + RBS_CATALOG +
                     TERMINATOR_CATALOG + CDS_CATALOG):
            self.parts[part.name] = part

    def get(self, name: str) -> GeneticPart:
        return self.parts[name]

    def query(self, part_type: Optional[PartType] = None) -> List[GeneticPart]:
        results = list(self.parts.values())
        if part_type:
            results = [p for p in results if p.part_type == part_type]
        return results

    def add(self, part: GeneticPart):
        self.parts[part.name] = part

    def to_sbol_json(self) -> str:
        """Export catalog in simplified SBOL-like JSON format."""
        data = []
        for p in self.parts.values():
            data.append({
                "displayId": p.name,
                "type": p.part_type.value,
                "sequence": p.sequence,
                "parameters": p.parameters,
                "sbol_uri": p.sbol_uri,
                "description": p.description,
            })
        return json.dumps(data, indent=2)
