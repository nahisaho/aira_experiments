#!/usr/bin/env python3
"""Automated FAIR compliance assessment simulation and analysis."""

from __future__ import annotations

import json
import math
import random
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)
random.seed(RANDOM_SEED)
RNG = np.random.default_rng(RANDOM_SEED)

BASE_DIR = Path(__file__).resolve().parents[1]
FIGURES_DIR = BASE_DIR / "figures"
RESULTS_DIR = BASE_DIR / "results"
DATA_DIR = BASE_DIR / "data"
LOGS_DIR = BASE_DIR / "logs"
REPORT_PATH = BASE_DIR / "report.md"
RESULTS_PATH = RESULTS_DIR / "fair_results.json"
DATASET_PATH = DATA_DIR / "simulated_fair_datasets.csv"
PREPROCESSING_LOG_PATH = DATA_DIR / "preprocessing-log.md"
STAT_SUMMARY_PATH = RESULTS_DIR / "statistical-summary.md"
PROCESS_LOG_PATH = LOGS_DIR / "process-log.jsonl"

OPEN_PROTOCOLS = {"https", "http", "ftp", "ftps", "rsync"}
STANDARD_PROTOCOLS = OPEN_PROTOCOLS | {"s3", "globus", "api"}
OPEN_LICENSES = {
    "CC-BY-4.0",
    "CC0-1.0",
    "MIT",
    "Apache-2.0",
    "ODC-BY-1.0",
    "ODbL-1.0",
    "BSD-3-Clause",
}
CLEAR_BUT_RESTRICTED_LICENSES = {
    "Custom-Restricted",
    "Data-Use-Agreement",
    "CC-BY-NC-4.0",
    "Institutional-Access",
}
RECOGNIZED_VOCABS = {
    "Gene Ontology",
    "NCBI Taxonomy",
    "ChEBI",
    "EDAM",
    "Dublin Core",
    "Schema.org",
    "EnvO",
    "MeSH",
    "SNOMED CT",
    "UMLS",
    "ORCID",
    "DataCite",
    "OBO",
}
FORMAL_SCHEMAS = {"DataCite", "Dublin Core", "JSON-LD", "RDF", "Schema.org", "XML"}
OPEN_INTEROPERABLE_FORMATS = {
    "CSV",
    "TSV",
    "JSON",
    "JSON-LD",
    "XML",
    "RDF",
    "TTL",
    "Parquet",
    "NetCDF",
    "HDF5",
    "GeoTIFF",
    "OME-TIFF",
    "TIFF",
    "FASTA",
    "FASTQ",
    "BAM",
    "NWB",
    "SDF",
    "mzML",
}
PROPRIETARY_OR_LIMITED_FORMATS = {"XLSX", "SAV", "DTA", "MAT", "CZI", "LIF", "CDX", "DOCX"}
REQUIRED_METADATA_FIELDS = [
    "title",
    "description",
    "creator",
    "keywords",
    "publication_year",
    "version",
    "license",
    "identifier",
    "repository",
    "field",
    "methodology",
    "funding",
]
CORE_METADATA_FIELDS = {"title", "description", "creator", "publication_year", "identifier"}

REPOSITORIES = {
    "Zenodo": {"quality": 0.06, "pid": 0.08, "search": 0.08, "auth": 0.22, "curation": 0.05, "reuse": 0.12},
    "Dryad": {"quality": 0.10, "pid": 0.12, "search": 0.10, "auth": 0.18, "curation": 0.09, "reuse": 0.10},
    "Figshare": {"quality": 0.03, "pid": 0.06, "search": 0.07, "auth": 0.20, "curation": 0.03, "reuse": 0.08},
    "Dataverse": {"quality": 0.09, "pid": 0.10, "search": 0.09, "auth": 0.35, "curation": 0.10, "reuse": 0.09},
    "GEO": {"quality": 0.12, "pid": 0.10, "search": 0.12, "auth": 0.28, "curation": 0.12, "reuse": 0.13},
    "PANGAEA": {"quality": 0.11, "pid": 0.11, "search": 0.11, "auth": 0.16, "curation": 0.12, "reuse": 0.11},
}
REPOSITORY_WEIGHTS = np.array([0.20, 0.13, 0.20, 0.17, 0.15, 0.15])

FIELDS = {
    "Genomics": {"quality": 0.10, "reuse": 0.14, "standard_bonus": 0.16},
    "Climate Science": {"quality": 0.12, "reuse": 0.11, "standard_bonus": 0.17},
    "Social Science": {"quality": -0.02, "reuse": 0.08, "standard_bonus": 0.05},
    "Imaging": {"quality": 0.04, "reuse": 0.10, "standard_bonus": 0.12},
    "Chemistry": {"quality": 0.02, "reuse": 0.09, "standard_bonus": 0.10},
    "Neuroscience": {"quality": 0.05, "reuse": 0.12, "standard_bonus": 0.14},
}
FIELD_WEIGHTS = np.array([0.18, 0.15, 0.17, 0.17, 0.15, 0.18])

FIELD_FORMATS = {
    "Genomics": {
        "open": ["FASTQ", "FASTA", "BAM", "TSV", "CSV", "HDF5"],
        "limited": ["XLSX", "DOCX"],
        "vocabs": ["Gene Ontology", "NCBI Taxonomy", "EDAM", "ORCID"],
        "standard": "MIxS",
    },
    "Climate Science": {
        "open": ["NetCDF", "CSV", "GeoTIFF", "JSON", "HDF5"],
        "limited": ["XLSX", "MAT"],
        "vocabs": ["EnvO", "Dublin Core", "Schema.org", "ORCID"],
        "standard": "CF Conventions",
    },
    "Social Science": {
        "open": ["CSV", "TSV", "JSON", "Parquet"],
        "limited": ["SAV", "DTA", "XLSX"],
        "vocabs": ["Dublin Core", "Schema.org", "ORCID", "DataCite"],
        "standard": "DDI",
    },
    "Imaging": {
        "open": ["OME-TIFF", "TIFF", "CSV", "JSON", "HDF5"],
        "limited": ["CZI", "LIF", "JPEG"],
        "vocabs": ["MeSH", "OBO", "Schema.org", "ORCID"],
        "standard": "OME",
    },
    "Chemistry": {
        "open": ["SDF", "CSV", "JSON", "mzML"],
        "limited": ["CDX", "XLSX"],
        "vocabs": ["ChEBI", "MeSH", "Schema.org", "ORCID"],
        "standard": "JCAMP-DX",
    },
    "Neuroscience": {
        "open": ["NWB", "CSV", "JSON", "HDF5"],
        "limited": ["MAT", "XLSX"],
        "vocabs": ["SNOMED CT", "UMLS", "OBO", "ORCID"],
        "standard": "BIDS",
    },
}

RELATION_TYPES = ["IsDerivedFrom", "References", "IsSupplementTo", "IsVersionOf"]
PROVENANCE_FIELDS = ["collection_protocol", "processing_pipeline", "instrument", "version_history"]


@dataclass(frozen=True)
class FAIRSubPrinciple:
    code: str
    principle: str
    description: str


SUBPRINCIPLES = [
    FAIRSubPrinciple("F1", "Findable", "Globally unique persistent identifier"),
    FAIRSubPrinciple("F2", "Findable", "Rich metadata are described"),
    FAIRSubPrinciple("F3", "Findable", "Metadata explicitly include the data identifier"),
    FAIRSubPrinciple("F4", "Findable", "Indexed in a searchable resource"),
    FAIRSubPrinciple("A1", "Accessible", "Retrievable by identifier using a standard protocol"),
    FAIRSubPrinciple("A1.1", "Accessible", "Protocol is open, free, and universally implementable"),
    FAIRSubPrinciple("A1.2", "Accessible", "Protocol allows authentication and authorization where needed"),
    FAIRSubPrinciple("A2", "Accessible", "Metadata remain accessible when the data are unavailable"),
    FAIRSubPrinciple("I1", "Interoperable", "Uses a formal, accessible knowledge representation language"),
    FAIRSubPrinciple("I2", "Interoperable", "Uses FAIR-compliant vocabularies and ontologies"),
    FAIRSubPrinciple("I3", "Interoperable", "Includes qualified references to related data"),
    FAIRSubPrinciple("R1", "Reusable", "Richly described with relevant and accurate attributes"),
    FAIRSubPrinciple("R1.1", "Reusable", "Clear and accessible data usage license"),
    FAIRSubPrinciple("R1.2", "Reusable", "Detailed provenance is associated with the data"),
    FAIRSubPrinciple("R1.3", "Reusable", "Conforms to community standards"),
]
SUBPRINCIPLE_CODES = [item.code for item in SUBPRINCIPLES]
PRINCIPLE_MAP = {
    "Findable": ["F1", "F2", "F3", "F4"],
    "Accessible": ["A1", "A1.1", "A1.2", "A2"],
    "Interoperable": ["I1", "I2", "I3"],
    "Reusable": ["R1", "R1.1", "R1.2", "R1.3"],
}


def ensure_output_dirs() -> None:
    for path in [FIGURES_DIR, RESULTS_DIR, DATA_DIR, LOGS_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def log_event(event_type: str, phase: str, skill_or_tool: str, handoff_in: Mapping[str, Any], handoff_out: Mapping[str, Any], files_written: Sequence[str] | None = None, status: str = "ok") -> None:
    entry = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "phase": phase,
        "event_type": event_type,
        "actor": "co-scientist",
        "skill_or_tool": skill_or_tool,
        "handoff_in": handoff_in,
        "handoff_out": handoff_out,
        "files_written": list(files_written or []),
        "status": status,
    }
    with PROCESS_LOG_PATH.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(entry) + "\n")


class FAIRAssessmentTool:
    def __init__(self) -> None:
        self.framework = {item.code: item for item in SUBPRINCIPLES}

    @staticmethod
    def check_persistent_identifier(identifier: str | None) -> bool:
        if not identifier:
            return False
        normalized = identifier.lower()
        return normalized.startswith("10.") or normalized.startswith("doi:") or "doi.org/10." in normalized or normalized.startswith("hdl:") or normalized.startswith("ark:") or normalized.startswith("urn:uuid:")

    @staticmethod
    def evaluate_metadata_completeness(metadata: Mapping[str, Any]) -> Dict[str, Any]:
        present = [field for field in REQUIRED_METADATA_FIELDS if metadata.get(field) not in (None, "", [], {})]
        completeness = len(present) / len(REQUIRED_METADATA_FIELDS)
        missing = [field for field in REQUIRED_METADATA_FIELDS if field not in present]
        return {"present": present, "missing": missing, "completeness": completeness}

    @staticmethod
    def check_license(license_name: str | None) -> Dict[str, Any]:
        if not license_name:
            return {"present": False, "clear": False, "open": False, "category": "missing"}
        if license_name in OPEN_LICENSES:
            return {"present": True, "clear": True, "open": True, "category": "open"}
        if license_name in CLEAR_BUT_RESTRICTED_LICENSES:
            return {"present": True, "clear": True, "open": False, "category": "restricted"}
        return {"present": True, "clear": False, "open": False, "category": "unknown"}

    @staticmethod
    def assess_format_interoperability(file_formats: Sequence[str]) -> Dict[str, Any]:
        if not file_formats:
            return {"open_ratio": 0.0, "interoperable": False, "formal_representation": False}
        open_count = sum(fmt in OPEN_INTEROPERABLE_FORMATS for fmt in file_formats)
        formal = any(fmt in OPEN_INTEROPERABLE_FORMATS for fmt in file_formats)
        open_ratio = open_count / len(file_formats)
        return {
            "open_ratio": open_ratio,
            "interoperable": open_ratio >= 0.5,
            "formal_representation": formal,
        }

    @staticmethod
    def check_readme(documentation: str | None) -> bool:
        return bool(documentation and documentation.strip())

    @staticmethod
    def evaluate_vocabulary_usage(vocabularies: Sequence[str]) -> Dict[str, Any]:
        recognized = sorted(set(vocabularies).intersection(RECOGNIZED_VOCABS))
        return {
            "present": bool(vocabularies),
            "recognized": recognized,
            "fair_compliant": bool(recognized),
        }

    @staticmethod
    def evaluate_provenance(provenance: Mapping[str, Any]) -> Dict[str, Any]:
        completed = [field for field in PROVENANCE_FIELDS if provenance.get(field)]
        completeness = len(completed) / len(PROVENANCE_FIELDS)
        return {"present": bool(completed), "completed": completed, "completeness": completeness}

    @staticmethod
    def has_qualified_references(related_resources: Sequence[Mapping[str, Any]]) -> bool:
        return any(item.get("relation_type") in RELATION_TYPES and item.get("target") for item in related_resources)

    def assess_dataset(self, record: Mapping[str, Any]) -> Dict[str, Any]:
        metadata_eval = self.evaluate_metadata_completeness(record["metadata"])
        license_eval = self.check_license(record.get("license"))
        format_eval = self.assess_format_interoperability(record.get("file_formats", []))
        readme_present = self.check_readme(record.get("documentation"))
        vocab_eval = self.evaluate_vocabulary_usage(record.get("vocabularies", []))
        provenance_eval = self.evaluate_provenance(record.get("provenance", {}))
        id_in_metadata = record.get("identifier") and record["metadata"].get("identifier") == record.get("identifier")

        scores = {
            "F1": int(self.check_persistent_identifier(record.get("identifier"))),
            "F2": int(metadata_eval["completeness"] >= 0.70),
            "F3": int(bool(id_in_metadata)),
            "F4": int(bool(record.get("indexed_in_searchable_resource"))),
            "A1": int(bool(record.get("identifier")) and record.get("access_protocol") in STANDARD_PROTOCOLS),
            "A1.1": int(record.get("access_protocol") in OPEN_PROTOCOLS),
            "A1.2": int(bool(record.get("protocol_supports_auth"))),
            "A2": int(bool(record.get("metadata_persists"))),
            "I1": int(record.get("metadata_schema") in FORMAL_SCHEMAS or format_eval["formal_representation"]),
            "I2": int(vocab_eval["fair_compliant"]),
            "I3": int(self.has_qualified_references(record.get("related_resources", []))),
            "R1": int(metadata_eval["completeness"] >= 0.75 and readme_present),
            "R1.1": int(license_eval["clear"]),
            "R1.2": int(provenance_eval["completeness"] >= 0.60),
            "R1.3": int(bool(record.get("community_standard")) and format_eval["interoperable"]),
        }

        principle_scores = {
            principle: 100.0 * float(np.mean([scores[code] for code in codes]))
            for principle, codes in PRINCIPLE_MAP.items()
        }
        total_score = 100.0 * float(np.mean(list(scores.values())))

        return {
            "subprinciple_scores": scores,
            "principle_scores": principle_scores,
            "total_score": total_score,
            "metadata_completeness": metadata_eval["completeness"],
            "provenance_completeness": provenance_eval["completeness"],
            "format_open_ratio": format_eval["open_ratio"],
            "recognized_vocabularies": vocab_eval["recognized"],
        }


def clamp(value: float, lower: float = 0.02, upper: float = 0.98) -> float:
    return float(max(lower, min(upper, value)))


def choose_identifier(repository: str, year: int, has_pid: bool, index: int) -> str:
    if has_pid:
        suffix = f"{year}{index:04d}"
        prefix = random.choice(["10.5281", "10.5061", "10.6084", "10.1594", "10.18112"])
        return f"{prefix}/dataset.{suffix}"
    return f"local-{repository.lower()}-{year}-{index:04d}"


def choose_access_protocol(repository: str, quality: float) -> str:
    standard_protocols = ["https", "https", "api", "ftp", "s3", "globus"]
    restricted_protocols = ["email-request", "web-form", "custom-client"]
    if repository in {"GEO", "PANGAEA"} and quality > 0.55:
        standard_protocols.append("ftps")
    if quality >= 0.72:
        pool = standard_protocols + ["https"]
    elif quality >= 0.45:
        pool = standard_protocols + restricted_protocols[:2]
    else:
        pool = ["https", "ftp", "email-request", "web-form", "custom-client"]
    return random.choice(pool)


def build_metadata(record_id: str, repository: str, field: str, year: int, quality: float, license_name: str | None) -> Dict[str, Any]:
    metadata: Dict[str, Any] = {}
    for name in REQUIRED_METADATA_FIELDS:
        base_probability = 0.34 + 0.62 * quality + (0.12 if name in CORE_METADATA_FIELDS else 0.0)
        include = RNG.random() < clamp(base_probability, 0.18, 0.99)
        if name == "identifier":
            include = include and record_id.startswith(("10.", "doi:", "ark:", "hdl:", "urn:"))
        if not include:
            continue
        if name == "title":
            metadata[name] = f"{field} benchmark dataset {record_id.split('.')[-1]}"
        elif name == "description":
            metadata[name] = f"Curated {field.lower()} dataset deposited in {repository}."
        elif name == "creator":
            metadata[name] = "Consortium A"
        elif name == "keywords":
            metadata[name] = [field.lower(), repository.lower(), "fair"]
        elif name == "publication_year":
            metadata[name] = year
        elif name == "version":
            metadata[name] = f"v{1 + int(quality * 3)}.0"
        elif name == "license":
            metadata[name] = license_name
        elif name == "identifier":
            metadata[name] = record_id
        elif name == "repository":
            metadata[name] = repository
        elif name == "field":
            metadata[name] = field
        elif name == "methodology":
            metadata[name] = "Documented acquisition and preprocessing workflow"
        elif name == "funding":
            metadata[name] = "Grant-supported"
    return metadata


def choose_license(quality: float) -> str | None:
    roll = RNG.random()
    if roll < clamp(0.58 * quality + 0.18):
        return random.choice(sorted(OPEN_LICENSES))
    if roll < clamp(0.78 * quality + 0.30):
        return random.choice(sorted(CLEAR_BUT_RESTRICTED_LICENSES))
    return None


def choose_file_formats(field: str, quality: float) -> List[str]:
    pools = FIELD_FORMATS[field]
    format_count = int(RNG.integers(1, 4))
    formats: List[str] = []
    for _ in range(format_count):
        source = pools["open"] if RNG.random() < clamp(0.38 + 0.62 * quality) else pools["limited"]
        formats.append(random.choice(source))
    return sorted(set(formats))


def choose_vocabularies(field: str, quality: float) -> List[str]:
    vocabularies: List[str] = []
    for vocab in FIELD_FORMATS[field]["vocabs"]:
        if RNG.random() < clamp(0.20 + 0.75 * quality):
            vocabularies.append(vocab)
    return vocabularies


def build_related_resources(quality: float) -> List[Dict[str, str]]:
    count = int(RNG.integers(0, 4))
    resources = []
    for idx in range(count):
        if RNG.random() < clamp(0.25 + 0.70 * quality):
            resources.append(
                {
                    "relation_type": random.choice(RELATION_TYPES),
                    "target": f"10.1234/related.{RNG.integers(1000, 9999)}.{idx}",
                }
            )
    return resources


def build_provenance(quality: float) -> Dict[str, str]:
    provenance: Dict[str, str] = {}
    for field_name in PROVENANCE_FIELDS:
        if RNG.random() < clamp(0.18 + 0.82 * quality):
            provenance[field_name] = "documented"
    return provenance


def simulate_datasets(n_datasets: int = 500) -> pd.DataFrame:
    assessment_tool = FAIRAssessmentTool()
    rows: List[Dict[str, Any]] = []
    repository_names = list(REPOSITORIES.keys())
    field_names = list(FIELDS.keys())

    for index in range(n_datasets):
        repository = str(RNG.choice(repository_names, p=REPOSITORY_WEIGHTS))
        field = str(RNG.choice(field_names, p=FIELD_WEIGHTS))
        year = int(RNG.integers(2018, 2025))
        year_effect = (year - 2018) / 6.0 * 0.18
        quality = clamp(0.42 + REPOSITORIES[repository]["quality"] + FIELDS[field]["quality"] + year_effect + RNG.normal(0.0, 0.11), 0.05, 0.95)

        has_pid = RNG.random() < clamp(0.30 + 0.62 * quality + REPOSITORIES[repository]["pid"])
        identifier = choose_identifier(repository, year, has_pid, index)
        license_name = choose_license(quality)
        metadata = build_metadata(identifier, repository, field, year, quality, license_name)
        file_formats = choose_file_formats(field, quality)
        vocabularies = choose_vocabularies(field, quality)
        metadata_schema = random.choice(sorted(FORMAL_SCHEMAS)) if RNG.random() < clamp(0.22 + 0.70 * quality) else "free-text"
        documentation = "README available" if RNG.random() < clamp(0.22 + 0.72 * quality) else ""
        related_resources = build_related_resources(quality)
        provenance = build_provenance(quality)
        indexed = RNG.random() < clamp(0.48 + 0.38 * quality + REPOSITORIES[repository]["search"])
        metadata_persists = RNG.random() < clamp(0.50 + 0.38 * quality + REPOSITORIES[repository]["curation"])
        protocol = choose_access_protocol(repository, quality)
        protocol_supports_auth = RNG.random() < clamp(0.06 + 0.35 * quality + REPOSITORIES[repository]["auth"])
        community_standard = RNG.random() < clamp(0.12 + 0.58 * quality + FIELDS[field]["standard_bonus"])

        record = {
            "dataset_id": f"DS-{index + 1:04d}",
            "repository": repository,
            "field": field,
            "year": year,
            "identifier": identifier,
            "metadata": metadata,
            "license": license_name,
            "file_formats": file_formats,
            "documentation": documentation,
            "vocabularies": vocabularies,
            "metadata_schema": metadata_schema,
            "indexed_in_searchable_resource": indexed,
            "access_protocol": protocol,
            "protocol_supports_auth": protocol_supports_auth,
            "metadata_persists": metadata_persists,
            "related_resources": related_resources,
            "provenance": provenance,
            "community_standard": community_standard,
            "community_standard_name": FIELD_FORMATS[field]["standard"],
        }

        assessment = assessment_tool.assess_dataset(record)
        total_score = assessment["total_score"]
        reuse_mean = math.exp(
            0.65
            + 0.020 * total_score
            + REPOSITORIES[repository]["reuse"]
            + FIELDS[field]["reuse"]
            + 0.020 * (year - 2018)
        )
        reuse_count = int(RNG.negative_binomial(6, 6 / (6 + reuse_mean)))

        row: Dict[str, Any] = {
            "dataset_id": record["dataset_id"],
            "repository": repository,
            "field": field,
            "year": year,
            "identifier": identifier,
            "access_protocol": protocol,
            "license": license_name or "Missing",
            "metadata_schema": metadata_schema,
            "documentation_present": bool(documentation),
            "indexed_in_searchable_resource": indexed,
            "protocol_supports_auth": protocol_supports_auth,
            "metadata_persists": metadata_persists,
            "community_standard": community_standard,
            "community_standard_name": FIELD_FORMATS[field]["standard"],
            "metadata_field_count": len(metadata),
            "metadata_completeness": assessment["metadata_completeness"],
            "provenance_completeness": assessment["provenance_completeness"],
            "format_open_ratio": assessment["format_open_ratio"],
            "file_formats": "|".join(file_formats),
            "vocabularies": "|".join(vocabularies),
            "related_resource_count": len(related_resources),
            "reuse_count": reuse_count,
            "fair_total_score": total_score,
            **assessment["principle_scores"],
            **assessment["subprinciple_scores"],
        }
        rows.append(row)

    return pd.DataFrame(rows)


def confidence_interval(values: Sequence[float], alpha: float = 0.05) -> tuple[float, float]:
    arr = np.asarray(values, dtype=float)
    mean = float(np.mean(arr))
    if arr.size < 2:
        return mean, mean
    sem = stats.sem(arr)
    margin = float(stats.t.ppf(1 - alpha / 2, arr.size - 1) * sem)
    return mean - margin, mean + margin


def pearson_confidence_interval(r_value: float, n_obs: int, alpha: float = 0.05) -> tuple[float, float]:
    if n_obs <= 3:
        return r_value, r_value
    bounded_r = float(np.clip(r_value, -0.999999, 0.999999))
    z_value = np.arctanh(bounded_r)
    z_error = 1 / math.sqrt(n_obs - 3)
    z_critical = stats.norm.ppf(1 - alpha / 2)
    lower, upper = np.tanh([z_value - z_critical * z_error, z_value + z_critical * z_error])
    return float(lower), float(upper)


def to_native(value: Any) -> Any:
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, pd.Timestamp):
        return value.isoformat()
    if isinstance(value, dict):
        return {str(key): to_native(inner) for key, inner in value.items()}
    if isinstance(value, list):
        return [to_native(item) for item in value]
    return value


def build_repository_summary(data: pd.DataFrame) -> pd.DataFrame:
    return (
        data.groupby("repository")
        .agg(
            dataset_count=("dataset_id", "count"),
            mean_fair_score=("fair_total_score", "mean"),
            median_fair_score=("fair_total_score", "median"),
            std_fair_score=("fair_total_score", "std"),
            mean_reuse_count=("reuse_count", "mean"),
        )
        .round(3)
        .reset_index()
        .sort_values("mean_fair_score", ascending=False)
    )


def build_compliance_summary(data: pd.DataFrame) -> pd.DataFrame:
    compliance = data[SUBPRINCIPLE_CODES].mean().mul(100).round(2)
    principle_lookup = {item.code: item.principle for item in SUBPRINCIPLES}
    description_lookup = {item.code: item.description for item in SUBPRINCIPLES}
    return pd.DataFrame(
        {
            "subprinciple": compliance.index,
            "principle": [principle_lookup[code] for code in compliance.index],
            "description": [description_lookup[code] for code in compliance.index],
            "compliance_rate": compliance.values,
        }
    )


def build_temporal_summary(data: pd.DataFrame) -> pd.DataFrame:
    records: List[Dict[str, Any]] = []
    for year, subset in data.groupby("year"):
        lower, upper = confidence_interval(subset["fair_total_score"])
        records.append(
            {
                "year": int(year),
                "mean_fair_score": round(float(subset["fair_total_score"].mean()), 3),
                "ci_lower": round(lower, 3),
                "ci_upper": round(upper, 3),
                "dataset_count": int(len(subset)),
            }
        )
    return pd.DataFrame(records).sort_values("year")


def build_field_summary(data: pd.DataFrame) -> pd.DataFrame:
    return (
        data.groupby("field")
        .agg(
            dataset_count=("dataset_id", "count"),
            mean_fair_score=("fair_total_score", "mean"),
            mean_findable=("Findable", "mean"),
            mean_accessible=("Accessible", "mean"),
            mean_interoperable=("Interoperable", "mean"),
            mean_reusable=("Reusable", "mean"),
            mean_reuse_count=("reuse_count", "mean"),
        )
        .round(3)
        .reset_index()
        .sort_values("mean_fair_score", ascending=False)
    )


def build_correlation_summary(data: pd.DataFrame) -> Dict[str, Any]:
    fair_scores = data["fair_total_score"].to_numpy(dtype=float)
    reuse_counts = data["reuse_count"].to_numpy(dtype=float)
    pearson_r, pearson_p = stats.pearsonr(fair_scores, reuse_counts)
    spearman_rho, spearman_p = stats.spearmanr(fair_scores, reuse_counts)
    regression = stats.linregress(fair_scores, reuse_counts)
    lower, upper = pearson_confidence_interval(float(pearson_r), len(data))
    return {
        "pearson_r": round(float(pearson_r), 4),
        "pearson_p_value": round(float(pearson_p), 8),
        "pearson_ci95": [round(lower, 4), round(upper, 4)],
        "spearman_rho": round(float(spearman_rho), 4),
        "spearman_p_value": round(float(spearman_p), 8),
        "regression_slope": round(float(regression.slope), 4),
        "regression_intercept": round(float(regression.intercept), 4),
        "regression_r_squared": round(float(regression.rvalue ** 2), 4),
    }


def configure_plot_style() -> None:
    plt.style.use("default")
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "grid.linestyle": "--",
            "font.size": 10,
            "axes.titlesize": 13,
            "axes.labelsize": 11,
            "legend.frameon": False,
        }
    )


def save_radar_chart(repository_scores: pd.DataFrame) -> Path:
    labels = ["Findable", "Accessible", "Interoperable", "Reusable"]
    angles = np.linspace(0, 2 * np.pi, len(labels), endpoint=False)
    angles = np.concatenate([angles, [angles[0]]])
    fig, ax = plt.subplots(figsize=(8.5, 8.5), subplot_kw={"projection": "polar"})
    colors = plt.cm.tab10(np.linspace(0, 1, len(repository_scores)))

    for color, (_, row) in zip(colors, repository_scores.iterrows()):
        values = row[labels].to_numpy(dtype=float)
        values = np.concatenate([values, [values[0]]])
        ax.plot(angles, values, linewidth=2, color=color, label=row["repository"])
        ax.fill(angles, values, color=color, alpha=0.08)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels)
    ax.set_ylim(0, 100)
    ax.set_yticks([20, 40, 60, 80, 100])
    ax.set_title("FAIR principle scores by repository", pad=20)
    ax.legend(loc="upper right", bbox_to_anchor=(1.25, 1.12))
    output_path = FIGURES_DIR / "fair_radar_by_repo.png"
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path


def save_heatmap(data: pd.DataFrame) -> Path:
    heatmap_data = data.groupby("field")[SUBPRINCIPLE_CODES].mean().mul(100).loc[list(FIELDS.keys())]
    fig, ax = plt.subplots(figsize=(14, 6.5))
    im = ax.imshow(heatmap_data.to_numpy(), cmap="cividis", aspect="auto", vmin=0, vmax=100)
    ax.set_xticks(np.arange(len(SUBPRINCIPLE_CODES)))
    ax.set_xticklabels(SUBPRINCIPLE_CODES, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(heatmap_data.index)))
    ax.set_yticklabels(heatmap_data.index)
    ax.set_title("FAIR sub-principle compliance by field")
    ax.set_xlabel("Sub-principle")
    ax.set_ylabel("Field")
    for row in range(heatmap_data.shape[0]):
        for col in range(heatmap_data.shape[1]):
            value = heatmap_data.iloc[row, col]
            ax.text(col, row, f"{value:.0f}", ha="center", va="center", color="white" if value < 55 else "black", fontsize=8)
    colorbar = fig.colorbar(im, ax=ax)
    colorbar.set_label("Compliance rate (%)")
    output_path = FIGURES_DIR / "fair_compliance_heatmap.png"
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path


def save_temporal_trends(temporal_summary: pd.DataFrame) -> Path:
    fig, ax = plt.subplots(figsize=(10, 5.8))
    years = temporal_summary["year"].to_numpy(dtype=int)
    mean_scores = temporal_summary["mean_fair_score"].to_numpy(dtype=float)
    lower = temporal_summary["ci_lower"].to_numpy(dtype=float)
    upper = temporal_summary["ci_upper"].to_numpy(dtype=float)
    ax.plot(years, mean_scores, marker="o", linewidth=2.5, color=plt.cm.viridis(0.68), label="Mean FAIR score")
    ax.fill_between(years, lower, upper, color=plt.cm.viridis(0.45), alpha=0.25, label="95% CI")
    ax.set_title("Temporal FAIR compliance trends (2018-2024)")
    ax.set_xlabel("Year")
    ax.set_ylabel("Mean FAIR score")
    ax.set_xticks(years)
    ax.set_ylim(max(0, lower.min() - 5), min(100, upper.max() + 5))
    ax.legend()
    output_path = FIGURES_DIR / "fair_temporal_trends.png"
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path


def save_reuse_scatter(data: pd.DataFrame, correlation_summary: Mapping[str, Any]) -> Path:
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    fields = list(FIELDS.keys())
    colors = plt.cm.tab10(np.linspace(0, 1, len(fields)))
    for color, field in zip(colors, fields):
        subset = data[data["field"] == field]
        ax.scatter(subset["fair_total_score"], subset["reuse_count"], s=36, alpha=0.72, color=color, label=field, edgecolor="none")

    slope = correlation_summary["regression_slope"]
    intercept = correlation_summary["regression_intercept"]
    x_values = np.linspace(data["fair_total_score"].min(), data["fair_total_score"].max(), 200)
    y_values = slope * x_values + intercept
    ax.plot(x_values, y_values, color=plt.cm.viridis(0.15), linewidth=2.5, label="Linear fit")
    annotation = (
        f"Pearson r = {correlation_summary['pearson_r']:.2f}\n"
        f"95% CI [{correlation_summary['pearson_ci95'][0]:.2f}, {correlation_summary['pearson_ci95'][1]:.2f}]\n"
        f"R² = {correlation_summary['regression_r_squared']:.2f}"
    )
    ax.text(0.02, 0.98, annotation, transform=ax.transAxes, va="top", ha="left", bbox={"boxstyle": "round,pad=0.4", "facecolor": "white", "edgecolor": "#808080", "alpha": 0.95})
    ax.set_title("FAIR score and data reuse correlation")
    ax.set_xlabel("FAIR score")
    ax.set_ylabel("Reuse count")
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0))
    output_path = FIGURES_DIR / "fair_reuse_correlation.png"
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return output_path


def write_preprocessing_log(data: pd.DataFrame) -> None:
    content = f"""# Preprocessing Log

- Timestamp: {datetime.now(timezone.utc).isoformat()}
- Random seed: {RANDOM_SEED}
- Dataset count: {len(data)}
- Simulation scope: FAIR assessment records across repositories, fields, and years 2018-2024.
- Generated variables: repository, field, year, identifier, metadata completeness, provenance completeness, interoperability indicators, FAIR sub-principles, total FAIR score, reuse count.
- Normalization: all FAIR sub-principles encoded as binary 0/1 values and aggregated to percentage scores.
- File formats and vocabularies were generated from field-specific pools to emulate realistic interoperability profiles.
- No raw external data were ingested; all records are simulated for reproducible benchmarking.
"""
    PREPROCESSING_LOG_PATH.write_text(content, encoding="utf-8")


def write_statistical_summary(data: pd.DataFrame, correlation_summary: Mapping[str, Any], temporal_summary: pd.DataFrame) -> None:
    year_groups = [group["fair_total_score"].to_numpy(dtype=float) for _, group in data.groupby("year")]
    anova = stats.f_oneway(*year_groups)
    eta_squared = float((anova.statistic * (len(year_groups) - 1)) / (anova.statistic * (len(year_groups) - 1) + (len(data) - len(year_groups))))
    summary = f"""# Statistical Summary

## Correlation between FAIR score and reuse
- Pearson correlation: r = {correlation_summary['pearson_r']}, 95% CI {tuple(correlation_summary['pearson_ci95'])}, p = {correlation_summary['pearson_p_value']}
- Spearman correlation: rho = {correlation_summary['spearman_rho']}, p = {correlation_summary['spearman_p_value']}
- Linear trend: slope = {correlation_summary['regression_slope']} reuse-count units per FAIR-score point, R² = {correlation_summary['regression_r_squared']}

## Temporal FAIR trend
- One-way ANOVA across years: F = {anova.statistic:.3f}, p = {anova.pvalue:.6f}, eta² = {eta_squared:.3f}
- Mean FAIR score in 2018: {temporal_summary.iloc[0]['mean_fair_score']:.2f}
- Mean FAIR score in 2024: {temporal_summary.iloc[-1]['mean_fair_score']:.2f}

## Notes
- Pearson and linear regression were used for continuous score-to-reuse analysis; Spearman correlation was added as a rank-based sensitivity check.
- Temporal comparisons were summarized with ANOVA effect size (eta²) to complement p-values.
"""
    STAT_SUMMARY_PATH.write_text(summary, encoding="utf-8")


def write_report(repository_summary: pd.DataFrame, compliance_summary: pd.DataFrame, field_summary: pd.DataFrame, temporal_summary: pd.DataFrame, correlation_summary: Mapping[str, Any]) -> None:
    best_repo = repository_summary.iloc[0]
    top_field = field_summary.iloc[0]
    weakest_subprinciple = compliance_summary.sort_values("compliance_rate").iloc[0]
    strongest_subprinciple = compliance_summary.sort_values("compliance_rate", ascending=False).iloc[0]
    report = f"""# DRAFT — NOT FOR DISTRIBUTION

## FAIR Assessment Tool Report

- Timestamp: {datetime.now(timezone.utc).isoformat()}
- Scope: Simulated FAIR compliance assessment for 500 datasets across 6 repositories, 6 research fields, and publication years 2018-2024.
- Random seed: {RANDOM_SEED}

## Methods
This workflow implemented automated checkers for persistent identifiers, metadata completeness, license clarity, interoperability of file formats, README/documentation presence, vocabulary usage, qualified cross-references, provenance coverage, and community-standard alignment. Each dataset was assessed against the 15 FAIR sub-principles (F1-F4, A1-A2, I1-I3, R1-R1.3). Total FAIR score was calculated as the mean percentage compliance across sub-principles.

## Key Results
- Highest-scoring repository: {best_repo['repository']} (mean FAIR score {best_repo['mean_fair_score']:.2f}).
- Highest-scoring field: {top_field['field']} (mean FAIR score {top_field['mean_fair_score']:.2f}).
- Strongest sub-principle: {strongest_subprinciple['subprinciple']} ({strongest_subprinciple['compliance_rate']:.1f}% compliance).
- Weakest sub-principle: {weakest_subprinciple['subprinciple']} ({weakest_subprinciple['compliance_rate']:.1f}% compliance).
- FAIR score improved from {temporal_summary.iloc[0]['mean_fair_score']:.2f} in 2018 to {temporal_summary.iloc[-1]['mean_fair_score']:.2f} in 2024.
- FAIR score vs reuse correlation: Pearson r = {correlation_summary['pearson_r']}, 95% CI {tuple(correlation_summary['pearson_ci95'])}.

## Figure Inventory
- `figures/fair_radar_by_repo.png`: Repository-level FAIR radar chart.
- `figures/fair_compliance_heatmap.png`: Field-by-sub-principle compliance heatmap.
- `figures/fair_temporal_trends.png`: Annual FAIR trend with 95% confidence interval.
- `figures/fair_reuse_correlation.png`: FAIR score vs reuse count scatter plot with regression line.

## File Inventory
- `src/fair_assessment_tool.py`
- `results/fair_results.json`
- `results/statistical-summary.md`
- `data/simulated_fair_datasets.csv`
- `data/preprocessing-log.md`
- `logs/process-log.jsonl`
"""
    REPORT_PATH.write_text(report, encoding="utf-8")


def generate_outputs(data: pd.DataFrame) -> Dict[str, Any]:
    repository_summary = build_repository_summary(data)
    compliance_summary = build_compliance_summary(data)
    temporal_summary = build_temporal_summary(data)
    field_summary = build_field_summary(data)
    correlation_summary = build_correlation_summary(data)

    repository_principles = (
        data.groupby("repository")[["Findable", "Accessible", "Interoperable", "Reusable"]]
        .mean()
        .round(3)
        .reset_index()
        .sort_values("repository")
    )

    figures = {
        "fair_radar_by_repo": str(save_radar_chart(repository_principles).relative_to(BASE_DIR)),
        "fair_compliance_heatmap": str(save_heatmap(data).relative_to(BASE_DIR)),
        "fair_temporal_trends": str(save_temporal_trends(temporal_summary).relative_to(BASE_DIR)),
        "fair_reuse_correlation": str(save_reuse_scatter(data, correlation_summary).relative_to(BASE_DIR)),
    }

    write_preprocessing_log(data)
    write_statistical_summary(data, correlation_summary, temporal_summary)
    write_report(repository_summary, compliance_summary, field_summary, temporal_summary, correlation_summary)

    results_payload = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "random_seed": RANDOM_SEED,
        "simulation": {
            "dataset_count": int(len(data)),
            "repositories": list(REPOSITORIES.keys()),
            "fields": list(FIELDS.keys()),
            "year_range": [2018, 2024],
        },
        "framework": [item.__dict__ for item in SUBPRINCIPLES],
        "repository_summary": repository_summary.to_dict(orient="records"),
        "principle_compliance": compliance_summary.to_dict(orient="records"),
        "temporal_trends": temporal_summary.to_dict(orient="records"),
        "field_patterns": field_summary.to_dict(orient="records"),
        "fair_reuse_correlation": correlation_summary,
        "figures": figures,
        "data_file": str(DATASET_PATH.relative_to(BASE_DIR)),
        "report_file": str(REPORT_PATH.relative_to(BASE_DIR)),
        "statistical_summary_file": str(STAT_SUMMARY_PATH.relative_to(BASE_DIR)),
    }

    RESULTS_PATH.write_text(json.dumps(to_native(results_payload), indent=2), encoding="utf-8")
    data.to_csv(DATASET_PATH, index=False)
    return to_native(results_payload)


def main() -> None:
    ensure_output_dirs()
    log_event(
        event_type="run_started",
        phase="plan",
        skill_or_tool="co-scientist-data-analysis",
        handoff_in={"request": "Create automated FAIR assessment tool and outputs"},
        handoff_out={"status": "initialized"},
    )
    log_event(
        event_type="prompt_received",
        phase="plan",
        skill_or_tool="co-scientist-data-analysis",
        handoff_in={"dataset_target": 500, "year_range": [2018, 2024]},
        handoff_out={"framework": SUBPRINCIPLE_CODES},
    )
    log_event(
        event_type="skill_selected",
        phase="plan",
        skill_or_tool="co-scientist-data-analysis",
        handoff_in={"analysis_type": "simulated FAIR benchmarking"},
        handoff_out={"outputs": [str(RESULTS_PATH), str(FIGURES_DIR)]},
    )

    configure_plot_style()
    log_event(
        event_type="handoff_started",
        phase="execute",
        skill_or_tool="fair_assessment_tool.py",
        handoff_in={"random_seed": RANDOM_SEED},
        handoff_out={"status": "simulating datasets"},
    )
    dataset = simulate_datasets(500)
    results_payload = generate_outputs(dataset)
    log_event(
        event_type="handoff_completed",
        phase="verify",
        skill_or_tool="fair_assessment_tool.py",
        handoff_in={"dataset_count": len(dataset)},
        handoff_out={"mean_fair_score": round(float(dataset['fair_total_score'].mean()), 3)},
        files_written=[
            str(RESULTS_PATH.relative_to(BASE_DIR)),
            str(DATASET_PATH.relative_to(BASE_DIR)),
            str(REPORT_PATH.relative_to(BASE_DIR)),
            str(STAT_SUMMARY_PATH.relative_to(BASE_DIR)),
        ],
    )
    for file_name in [
        REPORT_PATH,
        RESULTS_PATH,
        DATASET_PATH,
        PREPROCESSING_LOG_PATH,
        STAT_SUMMARY_PATH,
        FIGURES_DIR / "fair_radar_by_repo.png",
        FIGURES_DIR / "fair_compliance_heatmap.png",
        FIGURES_DIR / "fair_temporal_trends.png",
        FIGURES_DIR / "fair_reuse_correlation.png",
    ]:
        log_event(
            event_type="file_written",
            phase="report",
            skill_or_tool="fair_assessment_tool.py",
            handoff_in={"file": str(file_name.relative_to(BASE_DIR))},
            handoff_out={"exists": file_name.exists()},
            files_written=[str(file_name.relative_to(BASE_DIR))],
        )
    log_event(
        event_type="report_finalized",
        phase="report",
        skill_or_tool="fair_assessment_tool.py",
        handoff_in={"report": str(REPORT_PATH.relative_to(BASE_DIR))},
        handoff_out={"summary_file": str(STAT_SUMMARY_PATH.relative_to(BASE_DIR))},
        files_written=[str(REPORT_PATH.relative_to(BASE_DIR)), str(STAT_SUMMARY_PATH.relative_to(BASE_DIR))],
    )
    log_event(
        event_type="run_completed",
        phase="log",
        skill_or_tool="fair_assessment_tool.py",
        handoff_in={"result_file": str(RESULTS_PATH.relative_to(BASE_DIR))},
        handoff_out={"figures": results_payload["figures"]},
        files_written=[str(PROCESS_LOG_PATH.relative_to(BASE_DIR))],
    )
    print(json.dumps({"status": "ok", "results": str(RESULTS_PATH), "figures": results_payload["figures"]}, indent=2))


if __name__ == "__main__":
    main()
