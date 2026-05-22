"""
再現性予測スコア統合モジュール

複数の指標を統合して論文の再現可能性を予測する。

入力指標:
1. 方法論の詳細度スコア
2. 統計的不整合スコア
3. P-hacking リスクスコア
4. サンプルサイズの妥当性
5. 効果量の妥当性
6. 事前登録の有無
7. データ・コードの公開状況
8. ジャーナルのインパクトファクター（逆相関の知見あり）

出力:
- 再現性予測スコア (0.0 - 1.0)
- 信頼区間
- リスク要因の内訳
"""

import math
from dataclasses import dataclass, field
from typing import List, Dict, Optional


@dataclass
class ReproducibilityScore:
    """再現性予測スコア"""
    overall_score: float  # 0.0 - 1.0
    confidence_interval: tuple  # (lower, upper)
    risk_factors: Dict[str, float]  # 各リスク要因のスコア
    protective_factors: Dict[str, float]  # 保護因子のスコア
    prediction_class: str  # likely_reproducible, uncertain, unlikely_reproducible
    component_scores: Dict[str, float]
    recommendations: List[str]
    details: str


class ReproducibilityScorer:
    """
    再現性予測スコアの統合計算。

    重み付きロジスティック回帰モデルで再現確率を推定する。
    重みは Open Science Collaboration (2015) の再現プロジェクトおよび
    Camerer et al. (2018) の知見に基づいて設定。
    """

    # 各コンポーネントの重みと方向
    COMPONENTS = {
        "methodology_quality": {
            "weight": 0.20,
            "direction": "positive",  # 高い = 再現性高い
            "description": "方法論の詳細度と品質",
        },
        "statistical_consistency": {
            "weight": 0.18,
            "direction": "positive",
            "description": "統計的整合性（GRIM/SPRITE通過率）",
        },
        "phacking_risk": {
            "weight": 0.15,
            "direction": "negative",  # 高い = 再現性低い
            "description": "P-hacking リスクスコア",
        },
        "effect_size_plausibility": {
            "weight": 0.12,
            "direction": "positive",
            "description": "効果量の妥当性",
        },
        "sample_size_adequacy": {
            "weight": 0.10,
            "direction": "positive",
            "description": "サンプルサイズの妥当性",
        },
        "preregistration": {
            "weight": 0.08,
            "direction": "positive",
            "description": "事前登録の有無",
        },
        "data_code_availability": {
            "weight": 0.08,
            "direction": "positive",
            "description": "データ・コードの公開状況",
        },
        "harking_risk": {
            "weight": 0.05,
            "direction": "negative",
            "description": "HARKing リスクスコア",
        },
        "journal_rigor": {
            "weight": 0.04,
            "direction": "positive",
            "description": "ジャーナルの方法論的厳密性",
        },
    }

    def score(
        self,
        methodology_quality: float = 0.5,
        statistical_consistency: float = 0.5,
        phacking_risk: float = 0.0,
        effect_size_plausibility: float = 0.5,
        sample_size_adequacy: float = 0.5,
        preregistration: float = 0.0,
        data_code_availability: float = 0.0,
        harking_risk: float = 0.0,
        journal_rigor: float = 0.5,
    ) -> ReproducibilityScore:
        """
        再現性予測スコアを計算。

        すべてのパラメータは 0.0 - 1.0 の範囲。
        """
        inputs = {
            "methodology_quality": methodology_quality,
            "statistical_consistency": statistical_consistency,
            "phacking_risk": phacking_risk,
            "effect_size_plausibility": effect_size_plausibility,
            "sample_size_adequacy": sample_size_adequacy,
            "preregistration": preregistration,
            "data_code_availability": data_code_availability,
            "harking_risk": harking_risk,
            "journal_rigor": journal_rigor,
        }

        # ロジスティック回帰の線形結合
        logit = -1.0  # 基本バイアス（ベースレート ≈ 0.27: Camerer et al.）

        for name, config in self.COMPONENTS.items():
            value = inputs[name]
            weight = config["weight"]

            if config["direction"] == "positive":
                contribution = value * weight * 8.0
            else:
                contribution = -(value * weight * 8.0)

            logit += contribution

        # シグモイド変換
        overall = 1.0 / (1.0 + math.exp(-logit))

        # 信頼区間（ブートストラップ近似）
        se = math.sqrt(overall * (1 - overall) / 10)  # 仮想サンプルN=10
        ci_lower = max(0.0, overall - 1.96 * se)
        ci_upper = min(1.0, overall + 1.96 * se)

        # リスク・保護因子の特定
        risk_factors = {}
        protective_factors = {}

        for name, config in self.COMPONENTS.items():
            value = inputs[name]
            if config["direction"] == "negative" and value > 0.3:
                risk_factors[config["description"]] = value
            elif config["direction"] == "positive" and value < 0.3:
                risk_factors[config["description"]] = 1.0 - value
            elif config["direction"] == "positive" and value > 0.7:
                protective_factors[config["description"]] = value

        # 予測クラス
        if overall >= 0.6:
            pred_class = "likely_reproducible"
        elif overall >= 0.35:
            pred_class = "uncertain"
        else:
            pred_class = "unlikely_reproducible"

        # 推奨事項
        recommendations = self._generate_recommendations(inputs, overall)

        return ReproducibilityScore(
            overall_score=overall,
            confidence_interval=(ci_lower, ci_upper),
            risk_factors=risk_factors,
            protective_factors=protective_factors,
            prediction_class=pred_class,
            component_scores=inputs,
            recommendations=recommendations,
            details=(
                f"Reproducibility Score: {overall:.2f} "
                f"[{ci_lower:.2f}, {ci_upper:.2f}]\n"
                f"Class: {pred_class}\n"
                f"Risk factors: {len(risk_factors)}\n"
                f"Protective factors: {len(protective_factors)}"
            ),
        )

    def _generate_recommendations(self, inputs: Dict[str, float],
                                  overall: float) -> List[str]:
        """改善推奨事項を生成"""
        recs = []

        if inputs["preregistration"] < 0.5:
            recs.append(
                "Consider preregistering hypotheses and analysis plans "
                "(e.g., OSF, AsPredicted)"
            )

        if inputs["data_code_availability"] < 0.5:
            recs.append(
                "Share data and analysis code in a public repository "
                "(e.g., Zenodo, GitHub, Dryad)"
            )

        if inputs["sample_size_adequacy"] < 0.5:
            recs.append(
                "Conduct a priori power analysis to justify sample sizes"
            )

        if inputs["methodology_quality"] < 0.5:
            recs.append(
                "Improve methodological reporting following ARRIVE/"
                "CONSORT/STROBE guidelines"
            )

        if inputs["phacking_risk"] > 0.5:
            recs.append(
                "Address potential p-hacking concerns: "
                "report all conducted analyses, including null results"
            )

        if inputs["statistical_consistency"] < 0.5:
            recs.append(
                "Verify reported statistics — inconsistencies detected "
                "in means, SDs, or p-values"
            )

        if overall < 0.4:
            recs.append(
                "PRIORITY: Multiple reproducibility risk factors present. "
                "Independent replication is strongly recommended."
            )

        return recs
