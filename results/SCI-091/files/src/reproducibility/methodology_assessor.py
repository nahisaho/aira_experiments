"""
方法論の詳細度評価モジュール

論文の方法セクションの記述充実度を定量的に評価し、
再現性への寄与を推定する。

評価軸:
1. プロトコル詳細度（試薬、手順、パラメータの記載）
2. 統計手法の報告品質
3. データ可用性・コード共有
4. サンプルサイズの正当化
5. ブラインド化・ランダム化の記載
6. 除外基準の事前定義
"""

import re
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple


@dataclass
class MethodologyScore:
    """方法論評価スコア"""
    overall_score: float  # 0.0 - 1.0
    dimension_scores: Dict[str, float]
    missing_elements: List[str]
    present_elements: List[str]
    quality_level: str  # insufficient, basic, good, excellent
    reproducibility_prediction: float  # 再現確率の推定
    details: str


class MethodologyAssessor:
    """
    方法論の詳細度を評価するアセッサー。

    チェック項目はNIH Rigor and Reproducibility Guidelines、
    ARRIVE/CONSORT/STROBE等のレポーティングガイドラインに基づく。
    """

    CHECKLIST = {
        "sample_size_justification": {
            "weight": 0.12,
            "patterns": [
                r'(?:power analysis|sample size.*calculation|a priori.*power)',
                r'(?:effect size.*d\s*=|Cohen.*d|η²|partial eta)',
                r'(?:G\*Power|required sample|minimum.*sample)',
            ],
            "description": "Sample size justification / power analysis",
        },
        "randomization": {
            "weight": 0.08,
            "patterns": [
                r'(?:random(?:ly|ized|ization)|pseudo-random)',
                r'(?:block randomiz|stratified random|permuted block)',
                r'(?:random number generator|randomization sequence)',
            ],
            "description": "Randomization procedure",
        },
        "blinding": {
            "weight": 0.08,
            "patterns": [
                r'(?:blind(?:ed|ing)|mask(?:ed|ing)|double.blind)',
                r'(?:single.blind|assessor.*blind|participant.*blind)',
            ],
            "description": "Blinding/masking",
        },
        "inclusion_exclusion": {
            "weight": 0.07,
            "patterns": [
                r'(?:inclusion criteria|exclusion criteria)',
                r'(?:eligibility|eligible participants)',
                r'(?:excluded.*because|inclusion.*defined)',
            ],
            "description": "Inclusion/exclusion criteria",
        },
        "statistical_methods": {
            "weight": 0.10,
            "patterns": [
                r'(?:ANOVA|t-test|chi-square|regression|mixed.model)',
                r'(?:Bonferroni|Tukey|FDR|multiple comparison)',
                r'(?:confidence interval|effect size|Bayesian)',
                r'(?:non-parametric|Mann-Whitney|Kruskal-Wallis)',
            ],
            "description": "Statistical methods specification",
        },
        "software_tools": {
            "weight": 0.06,
            "patterns": [
                r'(?:SPSS|R\s+(?:version|v\d)|Python|MATLAB|Stata|SAS)',
                r'(?:GraphPad|ImageJ|FIJI|CellProfiler)',
                r'(?:version\s+\d|v\d+\.\d+)',
            ],
            "description": "Software/tools with versions",
        },
        "reagents_materials": {
            "weight": 0.08,
            "patterns": [
                r'(?:catalog\s*(?:number|#|no)|Cat\.?\s*#)',
                r'(?:RRID|manufacturer|supplier|vendor)',
                r'(?:antibody|primer|plasmid|cell line)',
                r'(?:concentration|dilution|dose)',
            ],
            "description": "Reagents/materials identification",
        },
        "data_availability": {
            "weight": 0.10,
            "patterns": [
                r'(?:data.*available|deposited.*(?:GEO|SRA|Zenodo|Dryad|figshare))',
                r'(?:accession.*number|repository|open.*data)',
                r'(?:supplementary.*data|raw.*data.*available)',
            ],
            "description": "Data availability statement",
        },
        "code_availability": {
            "weight": 0.08,
            "patterns": [
                r'(?:code.*available|github\.com|gitlab|bitbucket)',
                r'(?:source.*code|script.*available|notebook.*available)',
                r'(?:reproducible.*code|analysis.*code)',
            ],
            "description": "Code availability",
        },
        "protocol_detail": {
            "weight": 0.08,
            "patterns": [
                r'(?:protocol.*(?:available|described|published))',
                r'(?:step.by.step|detailed.*procedure|following.*protocol)',
                r'(?:incubat(?:ed|ion).*(?:min|hour|°C)|centrifug)',
            ],
            "description": "Detailed experimental protocol",
        },
        "replication_info": {
            "weight": 0.07,
            "patterns": [
                r'(?:biological.*replicate|technical.*replicate)',
                r'(?:independent.*experiment|repeated.*(?:three|3|twice))',
                r'(?:n\s*=\s*\d+.*(?:per|each|group))',
            ],
            "description": "Replication information",
        },
        "ethics_approval": {
            "weight": 0.04,
            "patterns": [
                r'(?:IRB|ethics.*committee|institutional.*review)',
                r'(?:informed.*consent|ethical.*approval)',
                r'(?:IACUC|animal.*protocol|Helsinki)',
            ],
            "description": "Ethics approval",
        },
        "preregistration": {
            "weight": 0.04,
            "patterns": [
                r'(?:preregist|pre.regist|OSF|AsPredicted|ClinicalTrials)',
                r'(?:registered.*(?:protocol|plan|analysis))',
            ],
            "description": "Preregistration",
        },
    }

    def assess(self, methods_text: str,
               full_text: Optional[str] = None) -> MethodologyScore:
        """
        方法セクションの詳細度を評価する。

        Parameters
        ----------
        methods_text : str
            方法セクションのテキスト
        full_text : str, optional
            論文全文（データ可用性等は方法以外にもある場合）
        """
        search_text = full_text if full_text else methods_text
        dimension_scores = {}
        present = []
        missing = []

        weighted_total = 0.0
        weight_sum = 0.0

        for dim_name, dim_info in self.CHECKLIST.items():
            score = self._check_dimension(search_text, dim_info["patterns"])
            dimension_scores[dim_name] = score
            weight = dim_info["weight"]
            weighted_total += score * weight
            weight_sum += weight

            if score >= 0.5:
                present.append(dim_info["description"])
            else:
                missing.append(dim_info["description"])

        overall = weighted_total / weight_sum if weight_sum > 0 else 0.0

        # テキストの量的評価（短すぎる方法セクションはペナルティ）
        word_count = len(methods_text.split())
        length_factor = min(word_count / 500, 1.0)  # 500語以上で満点
        overall *= (0.7 + 0.3 * length_factor)

        # 品質レベル
        if overall >= 0.75:
            quality = "excellent"
        elif overall >= 0.50:
            quality = "good"
        elif overall >= 0.25:
            quality = "basic"
        else:
            quality = "insufficient"

        # 再現性予測（ロジスティック回帰風の変換）
        # Ioannidis (2005) の知見に基づく概算
        import math
        logit = -2.0 + 5.0 * overall  # 基本ロジット
        reproducibility_pred = 1.0 / (1.0 + math.exp(-logit))

        return MethodologyScore(
            overall_score=overall,
            dimension_scores=dimension_scores,
            missing_elements=missing,
            present_elements=present,
            quality_level=quality,
            reproducibility_prediction=reproducibility_pred,
            details=(
                f"Methodology Assessment: {quality}\n"
                f"Overall score: {overall:.2f}\n"
                f"Word count: {word_count}\n"
                f"Present: {len(present)}/{len(self.CHECKLIST)}\n"
                f"Missing: {', '.join(missing[:5])}"
                f"{'...' if len(missing) > 5 else ''}\n"
                f"Reproducibility prediction: {reproducibility_pred:.1%}"
            ),
        )

    def _check_dimension(self, text: str,
                         patterns: List[str]) -> float:
        """パターンマッチングで各次元のスコアを計算"""
        matches = 0
        for pattern in patterns:
            if re.search(pattern, text, re.IGNORECASE):
                matches += 1
        return min(matches / max(len(patterns) * 0.5, 1), 1.0)
