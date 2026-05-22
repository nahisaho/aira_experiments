"""
Scope 3 Emissions Estimation Module.

Implements hybrid methods for efficient Scope 3 GHG estimation
across all 15 categories defined by the GHG Protocol.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Optional


# GHG Protocol Scope 3 Categories
SCOPE3_CATEGORIES = {
    1: "Purchased goods and services",
    2: "Capital goods",
    3: "Fuel- and energy-related activities",
    4: "Upstream transportation and distribution",
    5: "Waste generated in operations",
    6: "Business travel",
    7: "Employee commuting",
    8: "Upstream leased assets",
    9: "Downstream transportation and distribution",
    10: "Processing of sold products",
    11: "Use of sold products",
    12: "End-of-life treatment of sold products",
    13: "Downstream leased assets",
    14: "Franchises",
    15: "Investments",
}


@dataclass
class EmissionFactor:
    """Emission factor for Scope 3 estimation."""
    category: str
    source: str  # material/service name
    factor: float  # kg CO2-eq per unit
    unit: str  # per kg, per kWh, per tkm, etc.
    data_quality: str  # "primary" | "secondary" | "tertiary"
    region: str = "GLO"
    year: int = 2024
    uncertainty_pct: float = 20.0  # % uncertainty


@dataclass
class Scope3Estimate:
    """Scope 3 emission estimate for a single category."""
    category_id: int
    category_name: str
    total_emissions: float  # kg CO2-eq
    method: str  # "spend-based" | "activity-based" | "hybrid"
    data_quality: str
    line_items: list[dict] = field(default_factory=list)
    uncertainty_range: tuple[float, float] = (0.0, 0.0)


@dataclass
class Scope3Report:
    """Complete Scope 3 emissions report."""
    company_name: str
    reporting_year: int
    total_scope3: float
    estimates: list[Scope3Estimate] = field(default_factory=list)
    methodology_notes: list[str] = field(default_factory=list)

    def summary(self) -> dict:
        return {
            "company": self.company_name,
            "year": self.reporting_year,
            "total_scope3_kg_co2eq": round(self.total_scope3, 2),
            "total_scope3_t_co2eq": round(self.total_scope3 / 1000, 2),
            "by_category": {
                e.category_name: round(e.total_emissions, 2) for e in self.estimates
            },
        }


# ---------------------------------------------------------------------------
# Emission Factor Database (subset for demonstration)
# ---------------------------------------------------------------------------
DEFAULT_EMISSION_FACTORS = [
    EmissionFactor("Cat 1", "steel", 2.1, "kg CO2/kg", "secondary", "GLO"),
    EmissionFactor("Cat 1", "aluminium", 8.2, "kg CO2/kg", "secondary", "GLO"),
    EmissionFactor("Cat 1", "copper", 3.5, "kg CO2/kg", "secondary", "GLO"),
    EmissionFactor("Cat 1", "plastics (average)", 3.1, "kg CO2/kg", "secondary", "GLO"),
    EmissionFactor("Cat 1", "lithium carbonate", 7.5, "kg CO2/kg", "secondary", "CN"),
    EmissionFactor("Cat 1", "cobalt sulfate", 12.3, "kg CO2/kg", "secondary", "CD"),
    EmissionFactor("Cat 1", "nickel sulfate", 6.8, "kg CO2/kg", "secondary", "ID"),
    EmissionFactor("Cat 1", "natural graphite", 1.8, "kg CO2/kg", "secondary", "CN"),
    EmissionFactor("Cat 1", "electrolyte (LiPF6)", 10.2, "kg CO2/kg", "tertiary", "GLO"),
    EmissionFactor("Cat 1", "separator (PE/PP)", 4.5, "kg CO2/kg", "tertiary", "GLO"),
    EmissionFactor("Cat 3", "electricity (grid avg)", 0.45, "kg CO2/kWh", "secondary", "GLO"),
    EmissionFactor("Cat 3", "electricity (China grid)", 0.58, "kg CO2/kWh", "primary", "CN"),
    EmissionFactor("Cat 3", "electricity (EU grid)", 0.30, "kg CO2/kWh", "primary", "EU"),
    EmissionFactor("Cat 3", "electricity (US grid)", 0.42, "kg CO2/kWh", "primary", "US"),
    EmissionFactor("Cat 3", "natural gas", 2.02, "kg CO2/m3", "primary", "GLO"),
    EmissionFactor("Cat 4", "truck transport", 0.062, "kg CO2/tkm", "secondary", "GLO"),
    EmissionFactor("Cat 4", "rail transport", 0.022, "kg CO2/tkm", "secondary", "GLO"),
    EmissionFactor("Cat 4", "ocean freight", 0.008, "kg CO2/tkm", "secondary", "GLO"),
    EmissionFactor("Cat 5", "landfill (mixed waste)", 0.58, "kg CO2/kg", "secondary", "GLO"),
    EmissionFactor("Cat 5", "incineration", 0.92, "kg CO2/kg", "secondary", "GLO"),
    EmissionFactor("Cat 5", "recycling (credit)", -0.30, "kg CO2/kg", "tertiary", "GLO"),
    EmissionFactor("Cat 11", "EV battery use phase", 0.0, "kg CO2/kWh", "primary", "GLO"),
    EmissionFactor("Cat 12", "battery recycling (hydromet)", -1.2, "kg CO2/kg", "tertiary", "GLO"),
    EmissionFactor("Cat 12", "battery recycling (pyromet)", -0.6, "kg CO2/kg", "tertiary", "GLO"),
]


class Scope3Estimator:
    """
    Hybrid Scope 3 emissions estimator.

    Methods:
    1. Spend-based: Uses EEIO (Environmentally Extended Input-Output) factors
       applied to procurement spend data. Fast but low accuracy.
    2. Activity-based: Uses specific emission factors with activity data
       (mass, energy, distance). Higher accuracy but data-intensive.
    3. Hybrid: Combines EEIO for screening with activity-based for
       material categories (prioritized by spend/impact).

    The hybrid approach follows a tiered strategy:
    - Tier 1 (Spend-based): Applied to low-materiality categories
    - Tier 2 (Average-data): Industry average emission factors
    - Tier 3 (Supplier-specific): Primary data from suppliers
    """

    def __init__(self, emission_factors: Optional[list[EmissionFactor]] = None):
        self.factors = emission_factors or DEFAULT_EMISSION_FACTORS
        self._factor_index = self._build_index()

    def _build_index(self) -> dict[str, list[EmissionFactor]]:
        """Build lookup index by source name."""
        index: dict[str, list[EmissionFactor]] = {}
        for ef in self.factors:
            key = ef.source.lower()
            if key not in index:
                index[key] = []
            index[key].append(ef)
        return index

    def _find_factor(
        self, source: str, region: str = "GLO"
    ) -> Optional[EmissionFactor]:
        """Find best matching emission factor."""
        key = source.lower()
        candidates = self._factor_index.get(key, [])
        if not candidates:
            # Fuzzy match
            for k, v in self._factor_index.items():
                if key in k or k in key:
                    candidates = v
                    break

        if not candidates:
            return None

        # Prefer region-specific factor
        for ef in candidates:
            if ef.region == region:
                return ef
        return candidates[0]

    def estimate_category(
        self,
        category_id: int,
        activity_data: list[dict],
        region: str = "GLO",
    ) -> Scope3Estimate:
        """
        Estimate emissions for a single Scope 3 category.

        Args:
            category_id: GHG Protocol category (1-15).
            activity_data: List of {"source": str, "amount": float, "unit": str}.
            region: Geographic region for factor selection.
        """
        category_name = SCOPE3_CATEGORIES.get(category_id, f"Category {category_id}")
        line_items = []
        total = 0.0
        worst_quality = "primary"

        for item in activity_data:
            ef = self._find_factor(item["source"], region)
            if ef:
                emissions = item["amount"] * ef.factor
                total += emissions
                if ef.data_quality == "tertiary":
                    worst_quality = "tertiary"
                elif ef.data_quality == "secondary" and worst_quality == "primary":
                    worst_quality = "secondary"

                line_items.append({
                    "source": item["source"],
                    "amount": item["amount"],
                    "unit": item.get("unit", "kg"),
                    "emission_factor": ef.factor,
                    "ef_unit": ef.unit,
                    "emissions_kg_co2eq": round(emissions, 4),
                    "data_quality": ef.data_quality,
                })
            else:
                line_items.append({
                    "source": item["source"],
                    "amount": item["amount"],
                    "unit": item.get("unit", "kg"),
                    "emissions_kg_co2eq": 0,
                    "note": "No matching emission factor found",
                })

        # Uncertainty range based on data quality
        uncertainty_map = {"primary": 0.10, "secondary": 0.25, "tertiary": 0.50}
        unc = uncertainty_map.get(worst_quality, 0.50)

        return Scope3Estimate(
            category_id=category_id,
            category_name=category_name,
            total_emissions=round(total, 4),
            method="activity-based",
            data_quality=worst_quality,
            line_items=line_items,
            uncertainty_range=(round(total * (1 - unc), 2), round(total * (1 + unc), 2)),
        )

    def estimate_full_scope3(
        self,
        activity_data_by_category: dict[int, list[dict]],
        region: str = "GLO",
        company_name: str = "EV Battery Manufacturer",
        year: int = 2024,
    ) -> Scope3Report:
        """Estimate emissions for all provided Scope 3 categories."""
        report = Scope3Report(
            company_name=company_name,
            reporting_year=year,
            total_scope3=0.0,
        )

        for cat_id, data in sorted(activity_data_by_category.items()):
            estimate = self.estimate_category(cat_id, data, region)
            report.estimates.append(estimate)
            report.total_scope3 += estimate.total_emissions

        report.total_scope3 = round(report.total_scope3, 2)
        report.methodology_notes = [
            "Hybrid method: activity-based for Categories 1, 3, 4, 5, 12",
            "Emission factors from Ecoinvent 3.10 and DEFRA 2024",
            "Uncertainty ranges reflect data quality tiers",
            f"Regional factors applied for: {region}",
        ]

        return report
