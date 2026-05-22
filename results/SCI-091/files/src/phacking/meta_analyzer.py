"""
P-hacking/HARKing メタ分析モジュール

複数の論文にわたるP-hacking指標を集約し、
研究分野全体の信頼性を評価する。
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Optional
from .phacking_detector import PHackingDetector, PHackingResult
from .harking_detector import HARKingDetector, HARKingResult


@dataclass
class MetaAnalysisResult:
    """メタ分析結果"""
    num_papers: int
    phacking_results: List[PHackingResult]
    harking_results: List[HARKingResult]
    aggregate_p_curve_power: float
    aggregate_excess_significance: float
    field_risk_score: float  # 分野全体のリスク
    high_risk_papers: List[str]
    summary_statistics: Dict
    recommendations: List[str]


class PHackingMetaAnalyzer:
    """
    複数論文にわたるP-hacking/HARKingのメタ分析。

    Parameters
    ----------
    significance_level : float
        有意水準
    """

    def __init__(self, significance_level: float = 0.05):
        self.phacking_detector = PHackingDetector(significance_level)
        self.harking_detector = HARKingDetector()

    def analyze_corpus(
        self,
        papers: List[Dict],
    ) -> MetaAnalysisResult:
        """
        論文コーパスのメタ分析を実行。

        Parameters
        ----------
        papers : list of dict
            各dictは以下のキーを持つ:
            - 'id': 論文ID
            - 'p_values': p値リスト
            - 'introduction': 序論テキスト (optional)
            - 'results': 結果テキスト (optional)
            - 'discussion': 考察テキスト (optional)
        """
        phacking_results = []
        harking_results = []
        high_risk_papers = []
        all_p_values = []

        for paper in papers:
            # P-hacking分析
            p_vals = paper.get("p_values", [])
            if p_vals:
                all_p_values.extend(p_vals)
                ph_result = self.phacking_detector.analyze(p_vals)
                phacking_results.append(ph_result)

                if ph_result.overall_risk in ("high", "very_high"):
                    high_risk_papers.append(paper.get("id", "unknown"))

            # HARKing分析
            intro = paper.get("introduction", "")
            results = paper.get("results", "")
            discussion = paper.get("discussion", "")
            if intro and results:
                hk_result = self.harking_detector.analyze(
                    intro, results, discussion
                )
                harking_results.append(hk_result)

                if hk_result.risk_level == "high":
                    pid = paper.get("id", "unknown")
                    if pid not in high_risk_papers:
                        high_risk_papers.append(pid)

        # 集約統計
        if all_p_values:
            aggregate_ph = self.phacking_detector.analyze(all_p_values)
            agg_power = (
                aggregate_ph.p_curve.power_estimate
                if aggregate_ph.p_curve else 0.0
            )
            agg_excess = 1.0 if aggregate_ph.excess_significance else 0.0
        else:
            agg_power = 0.0
            agg_excess = 0.0

        # 分野リスクスコア
        if phacking_results:
            avg_ph_score = np.mean(
                [r.overall_score for r in phacking_results]
            )
        else:
            avg_ph_score = 0.0

        if harking_results:
            avg_hk_score = np.mean(
                [r.risk_score for r in harking_results]
            )
        else:
            avg_hk_score = 0.0

        field_risk = float(0.6 * avg_ph_score + 0.4 * avg_hk_score)

        # 推奨事項
        recommendations = self._generate_recommendations(
            field_risk, phacking_results, harking_results,
            len(papers), len(high_risk_papers)
        )

        summary_stats = {
            "total_p_values": len(all_p_values),
            "significant_ratio": (
                sum(1 for p in all_p_values if p < 0.05) / len(all_p_values)
                if all_p_values else 0.0
            ),
            "mean_phacking_score": float(avg_ph_score),
            "mean_harking_score": float(avg_hk_score),
            "high_risk_ratio": (
                len(high_risk_papers) / len(papers) if papers else 0.0
            ),
        }

        return MetaAnalysisResult(
            num_papers=len(papers),
            phacking_results=phacking_results,
            harking_results=harking_results,
            aggregate_p_curve_power=agg_power,
            aggregate_excess_significance=agg_excess,
            field_risk_score=field_risk,
            high_risk_papers=high_risk_papers,
            summary_statistics=summary_stats,
            recommendations=recommendations,
        )

    def _generate_recommendations(
        self, risk: float,
        ph_results: List[PHackingResult],
        hk_results: List[HARKingResult],
        n_papers: int, n_high_risk: int,
    ) -> List[str]:
        """分析結果に基づく推奨事項を生成"""
        recs = []

        if risk > 0.5:
            recs.append(
                "CRITICAL: Field-level risk score is high. "
                "Systematic replication efforts are strongly recommended."
            )
        elif risk > 0.3:
            recs.append(
                "WARNING: Moderate field-level risk detected. "
                "Consider targeted replication of key findings."
            )

        if n_high_risk > n_papers * 0.3:
            recs.append(
                f"{n_high_risk}/{n_papers} papers flagged as high-risk. "
                "Journal-level policy review may be warranted."
            )

        # P-curve検出力
        powers = [
            r.p_curve.power_estimate
            for r in ph_results
            if r.p_curve is not None
        ]
        if powers and np.mean(powers) < 0.5:
            recs.append(
                f"Average estimated power is low ({np.mean(powers):.1%}). "
                "Studies may be underpowered — larger samples needed."
            )

        # 事前登録の推奨
        harking_high = sum(1 for r in hk_results if r.risk_level == "high")
        if harking_high > 0:
            recs.append(
                f"{harking_high} papers show HARKing indicators. "
                "Preregistration should be encouraged."
            )

        if not recs:
            recs.append("No major concerns detected at the meta-analysis level.")

        return recs
