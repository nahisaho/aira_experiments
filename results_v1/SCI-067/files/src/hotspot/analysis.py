"""
Hotspot Analysis and Scenario Comparison Module.

Identifies environmental impact hotspots and generates
comparative scenario analyses for LCA decision support.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class ImpactContribution:
    """Contribution of a process/flow to total impact."""
    process_name: str
    process_id: str
    impact_value: float
    impact_unit: str
    percentage: float
    category: str  # "material" | "energy" | "transport" | "emission" | "waste"
    rank: int = 0


@dataclass
class Hotspot:
    """An identified environmental hotspot."""
    process_name: str
    process_id: str
    impact_category: str
    contribution_pct: float
    absolute_value: float
    unit: str
    improvement_potential: str
    suggested_alternatives: list[str] = field(default_factory=list)


@dataclass
class ScenarioDefinition:
    """Definition of an LCA scenario for comparison."""
    name: str
    description: str
    parameter_changes: dict[str, float] = field(default_factory=dict)
    process_substitutions: dict[str, str] = field(default_factory=dict)
    metadata: dict = field(default_factory=dict)


@dataclass
class ScenarioResult:
    """Results for a single scenario."""
    scenario_name: str
    total_gwp: float  # kg CO2-eq
    total_ap: float   # kg SO2-eq (Acidification)
    total_ep: float   # kg PO4-eq (Eutrophication)
    total_ced: float  # MJ (Cumulative Energy Demand)
    process_contributions: list[ImpactContribution] = field(default_factory=list)
    hotspots: list[Hotspot] = field(default_factory=list)


@dataclass
class ScenarioComparison:
    """Comparison across multiple scenarios."""
    baseline_name: str
    scenarios: list[ScenarioResult] = field(default_factory=list)
    relative_changes: dict[str, dict[str, float]] = field(default_factory=dict)

    def compute_relative_changes(self) -> None:
        """Compute % change relative to baseline for each impact category."""
        baseline = next(
            (s for s in self.scenarios if s.scenario_name == self.baseline_name), None
        )
        if not baseline:
            return

        for scenario in self.scenarios:
            if scenario.scenario_name == self.baseline_name:
                continue
            changes = {}
            if baseline.total_gwp != 0:
                changes["GWP"] = ((scenario.total_gwp - baseline.total_gwp)
                                  / baseline.total_gwp * 100)
            if baseline.total_ap != 0:
                changes["AP"] = ((scenario.total_ap - baseline.total_ap)
                                 / baseline.total_ap * 100)
            if baseline.total_ep != 0:
                changes["EP"] = ((scenario.total_ep - baseline.total_ep)
                                 / baseline.total_ep * 100)
            if baseline.total_ced != 0:
                changes["CED"] = ((scenario.total_ced - baseline.total_ced)
                                  / baseline.total_ced * 100)
            self.relative_changes[scenario.scenario_name] = changes


class HotspotAnalyzer:
    """
    Performs contribution analysis and hotspot identification.

    Methods:
    1. Process contribution analysis (Pareto-based)
    2. Flow contribution analysis
    3. Upstream tracing (supply chain hotspots)
    4. Multi-criteria hotspot ranking
    """

    HOTSPOT_THRESHOLD = 0.10  # 10% contribution = hotspot

    def analyze(
        self,
        process_impacts: dict[str, dict],
        impact_category: str = "GWP",
    ) -> list[Hotspot]:
        """
        Identify hotspots from process-level impact contributions.

        Args:
            process_impacts: {process_id: {"name": str, "gwp": float, ...}}
            impact_category: Impact category to analyze.
        """
        cat_key = impact_category.lower()
        total = sum(
            p.get(cat_key, 0) for p in process_impacts.values()
        )
        if total == 0:
            return []

        hotspots = []
        for pid, data in process_impacts.items():
            value = data.get(cat_key, 0)
            pct = value / total
            if pct >= self.HOTSPOT_THRESHOLD:
                hotspots.append(
                    Hotspot(
                        process_name=data.get("name", pid),
                        process_id=pid,
                        impact_category=impact_category,
                        contribution_pct=round(pct * 100, 1),
                        absolute_value=round(value, 4),
                        unit=data.get("unit", "kg CO2-eq"),
                        improvement_potential=self._assess_improvement(pct),
                        suggested_alternatives=self._suggest_alternatives(
                            data.get("name", ""), data.get("category", "")
                        ),
                    )
                )

        hotspots.sort(key=lambda h: h.contribution_pct, reverse=True)
        return hotspots

    def _assess_improvement(self, pct: float) -> str:
        if pct >= 0.30:
            return "High — primary target for reduction"
        elif pct >= 0.15:
            return "Medium — significant reduction potential"
        else:
            return "Moderate — incremental improvement possible"

    def _suggest_alternatives(self, process_name: str, category: str) -> list[str]:
        """Suggest alternatives based on process type."""
        suggestions_db = {
            "electricity": [
                "Switch to renewable energy (solar/wind PPA)",
                "Improve energy efficiency",
                "On-site generation with battery storage",
            ],
            "aluminium": [
                "Use recycled aluminium (secondary)",
                "Substitute with lightweight composites",
                "Optimize material usage (lightweighting)",
            ],
            "steel": [
                "Use EAF steel (recycled scrap)",
                "Source from low-carbon steelmakers",
                "Material substitution where feasible",
            ],
            "transport": [
                "Nearshore/local sourcing",
                "Modal shift (road → rail/sea)",
                "Optimize logistics (load factor improvement)",
            ],
            "cathode": [
                "Transition NMC811 → LFP chemistry",
                "Increase recycled content",
                "Optimize synthesis process efficiency",
            ],
        }
        name_lower = process_name.lower()
        for key, suggestions in suggestions_db.items():
            if key in name_lower:
                return suggestions
        return ["Conduct detailed sub-process analysis", "Explore material substitution"]


class ScenarioComparator:
    """
    Generates and compares LCA scenarios automatically.

    Supports:
    - Parameter sensitivity scenarios
    - Technology substitution scenarios
    - Geographic variation scenarios
    - Temporal projection scenarios
    """

    def create_baseline(
        self,
        process_impacts: dict[str, dict],
        name: str = "Baseline",
    ) -> ScenarioResult:
        """Create baseline scenario result from process impacts."""
        total_gwp = sum(p.get("gwp", 0) for p in process_impacts.values())
        total_ap = sum(p.get("ap", 0) for p in process_impacts.values())
        total_ep = sum(p.get("ep", 0) for p in process_impacts.values())
        total_ced = sum(p.get("ced", 0) for p in process_impacts.values())

        contributions = []
        for pid, data in process_impacts.items():
            contributions.append(
                ImpactContribution(
                    process_name=data.get("name", pid),
                    process_id=pid,
                    impact_value=data.get("gwp", 0),
                    impact_unit="kg CO2-eq",
                    percentage=round(
                        data.get("gwp", 0) / total_gwp * 100, 1
                    ) if total_gwp > 0 else 0,
                    category=data.get("category", "material"),
                )
            )
        contributions.sort(key=lambda c: c.percentage, reverse=True)
        for i, c in enumerate(contributions):
            c.rank = i + 1

        analyzer = HotspotAnalyzer()
        hotspots = analyzer.analyze(process_impacts)

        return ScenarioResult(
            scenario_name=name,
            total_gwp=round(total_gwp, 4),
            total_ap=round(total_ap, 4),
            total_ep=round(total_ep, 4),
            total_ced=round(total_ced, 2),
            process_contributions=contributions,
            hotspots=hotspots,
        )

    def apply_scenario(
        self,
        baseline_impacts: dict[str, dict],
        scenario: ScenarioDefinition,
    ) -> ScenarioResult:
        """Apply a scenario definition to baseline impacts."""
        modified = {}
        for pid, data in baseline_impacts.items():
            modified[pid] = dict(data)
            for param, multiplier in scenario.parameter_changes.items():
                if param in pid or param in data.get("name", ""):
                    for key in ["gwp", "ap", "ep", "ced"]:
                        if key in modified[pid]:
                            modified[pid][key] *= multiplier

        return self.create_baseline(modified, name=scenario.name)

    def compare(
        self,
        baseline_impacts: dict[str, dict],
        scenarios: list[ScenarioDefinition],
    ) -> ScenarioComparison:
        """Run full scenario comparison."""
        baseline = self.create_baseline(baseline_impacts, "Baseline")
        comparison = ScenarioComparison(
            baseline_name="Baseline",
            scenarios=[baseline],
        )

        for scenario in scenarios:
            result = self.apply_scenario(baseline_impacts, scenario)
            comparison.scenarios.append(result)

        comparison.compute_relative_changes()
        return comparison
