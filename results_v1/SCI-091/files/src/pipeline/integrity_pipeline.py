"""
Research Integrity Assessment System (RIAS) — 統合パイプライン

NLPとコンピュータビジョンを統合した研究公正性評価パイプライン。

全6モジュールを統合し、論文単位の包括的な公正性評価を行う。

パイプライン:
  Paper Input
       │
       ├──→ [1] Image Forensics Module
       │         ├── Duplicate Detection (pHash/dHash + Block Matching)
       │         ├── Manipulation Detection (ELA + Noise Analysis)
       │         └── CNN Forensics (ManTraNet-style)
       │
       ├──→ [2] Statistical Checks Module
       │         ├── GRIM Test (Mean Granularity)
       │         ├── SPRITE Test (Descriptive Stats Consistency)
       │         └── DF Consistency / P-value Distribution
       │
       ├──→ [3] Plagiarism Detection Module
       │         ├── Winnowing Fingerprint (Fast Screening)
       │         ├── Citation-Aware Similarity (Detailed Comparison)
       │         └── Section-level Risk Assessment
       │
       ├──→ [4] P-hacking / HARKing Module
       │         ├── P-curve Analysis
       │         ├── Z-curve Analysis
       │         ├── Caliper Test / Excess Significance
       │         └── HARKing Linguistic Markers
       │
       ├──→ [5] Reproducibility Score Module
       │         ├── Methodology Detail Assessment
       │         ├── Reporting Quality (ARRIVE/CONSORT/STROBE)
       │         └── Integrated Reproducibility Prediction
       │
       └──→ [6] Validation Module
                 ├── PubPeer External Signals
                 ├── Retraction Watch Ground Truth
                 └── System Performance Evaluation
       │
       ▼
  Integrated Integrity Report
"""

import json
from dataclasses import dataclass, field, asdict
from typing import List, Dict, Optional, Any
from datetime import datetime
import numpy as np

from ..statistical_checks import StatisticalAnalyzer
from ..plagiarism import CitationAwareSimilarity, TextFingerprinter
from ..phacking import PHackingDetector, HARKingDetector
from ..reproducibility import ReproducibilityScorer, MethodologyAssessor
from ..validation import IntegrityValidator


@dataclass
class PaperInput:
    """論文入力データ"""
    doi: str
    title: str
    authors: List[str]
    sections: Dict[str, str]  # introduction, methods, results, discussion
    images: Optional[Dict[str, np.ndarray]] = None
    p_values: Optional[List[float]] = None
    metadata: Optional[Dict] = None


@dataclass
class ModuleResult:
    """個別モジュール結果"""
    module_name: str
    risk_score: float  # 0.0 - 1.0
    risk_level: str  # low, medium, high, critical
    flags: List[str]
    details: Dict[str, Any]


@dataclass
class IntegrityReport:
    """統合公正性レポート"""
    doi: str
    title: str
    timestamp: str
    overall_integrity_score: float  # 0.0 (不正) - 1.0 (公正)
    overall_risk_level: str
    module_results: Dict[str, ModuleResult]
    reproducibility_prediction: float
    top_concerns: List[str]
    recommendations: List[str]
    confidence: float  # レポートの信頼度


class IntegrityPipeline:
    """
    研究公正性評価の統合パイプライン。

    全6モジュールを順次実行し、統合スコアを算出する。

    Module Weights (検出寄与度):
        - image_forensics: 0.20
        - statistical_checks: 0.20
        - plagiarism: 0.15
        - phacking: 0.20
        - reproducibility: 0.15
        - external_signals: 0.10
    """

    MODULE_WEIGHTS = {
        "image_forensics": 0.20,
        "statistical_checks": 0.20,
        "plagiarism": 0.15,
        "phacking": 0.20,
        "reproducibility": 0.15,
        "external_signals": 0.10,
    }

    def __init__(self):
        self.stat_analyzer = StatisticalAnalyzer()
        self.similarity_engine = CitationAwareSimilarity()
        self.fingerprinter = TextFingerprinter()
        self.phacking_detector = PHackingDetector()
        self.harking_detector = HARKingDetector()
        self.repro_scorer = ReproducibilityScorer()
        self.method_assessor = MethodologyAssessor()
        self.validator = IntegrityValidator()

    def analyze(self, paper: PaperInput) -> IntegrityReport:
        """
        論文の包括的な公正性分析を実行。

        Parameters
        ----------
        paper : PaperInput
            分析対象の論文データ
        """
        module_results = {}

        # Module 1: 画像フォレンジクス
        if paper.images:
            img_result = self._run_image_forensics(paper.images)
            module_results["image_forensics"] = img_result
        else:
            module_results["image_forensics"] = ModuleResult(
                module_name="image_forensics",
                risk_score=0.0,
                risk_level="not_applicable",
                flags=["No images provided for analysis"],
                details={},
            )

        # Module 2: 統計的不整合
        stat_result = self._run_statistical_checks(paper)
        module_results["statistical_checks"] = stat_result

        # Module 3: 盗作検出（コーパスが必要なためプレースホルダー）
        plag_result = self._run_plagiarism_check(paper)
        module_results["plagiarism"] = plag_result

        # Module 4: P-hacking / HARKing
        phack_result = self._run_phacking_analysis(paper)
        module_results["phacking"] = phack_result

        # Module 5: 再現性予測
        repro_result = self._run_reproducibility_assessment(
            paper, stat_result, phack_result
        )
        module_results["reproducibility"] = repro_result

        # Module 6: 外部シグナル
        ext_result = ModuleResult(
            module_name="external_signals",
            risk_score=0.0,
            risk_level="pending",
            flags=["External validation requires API access"],
            details={"note": "PubPeer/Retraction Watch lookup pending"},
        )
        module_results["external_signals"] = ext_result

        # 統合スコア計算
        overall_score, risk_level, concerns = self._integrate_results(
            module_results
        )

        # 推奨事項
        recommendations = self._generate_recommendations(
            module_results, overall_score
        )

        # 再現性予測値の取得
        repro_pred = repro_result.details.get("reproducibility_prediction", 0.5)

        # 信頼度（利用可能モジュール数に基づく）
        active_modules = sum(
            1 for r in module_results.values()
            if r.risk_level not in ("not_applicable", "pending")
        )
        confidence = active_modules / len(self.MODULE_WEIGHTS)

        return IntegrityReport(
            doi=paper.doi,
            title=paper.title,
            timestamp=datetime.utcnow().isoformat(),
            overall_integrity_score=overall_score,
            overall_risk_level=risk_level,
            module_results=module_results,
            reproducibility_prediction=repro_pred,
            top_concerns=concerns[:5],
            recommendations=recommendations,
            confidence=confidence,
        )

    def _run_image_forensics(
        self, images: Dict[str, np.ndarray]
    ) -> ModuleResult:
        """画像フォレンジクスモジュールの実行"""
        from ..image_forensics import (
            DuplicateDetector, ManipulationDetector
        )

        dup_detector = DuplicateDetector()
        manip_detector = ManipulationDetector()

        flags = []
        max_risk = 0.0

        # 重複検出
        dup_results = dup_detector.find_duplicates(images)
        for dr in dup_results:
            if dr.similarity_score > 0.8:
                flags.append(
                    f"Image duplicate: {dr.image_id_a} ↔ {dr.image_id_b} "
                    f"({dr.duplication_type.value}, sim={dr.similarity_score:.2f})"
                )
                max_risk = max(max_risk, dr.similarity_score)

        # 加工検出
        manip_results = {}
        for img_id, img in images.items():
            result = manip_detector.detect(img)
            if result.is_manipulated:
                flags.append(
                    f"Image manipulation detected in {img_id}: "
                    f"{result.manipulation_type} "
                    f"(confidence={result.confidence:.2f})"
                )
                max_risk = max(max_risk, result.confidence)
            manip_results[img_id] = {
                "is_manipulated": result.is_manipulated,
                "confidence": result.confidence,
                "type": result.manipulation_type,
            }

        risk_level = self._score_to_level(max_risk)

        return ModuleResult(
            module_name="image_forensics",
            risk_score=max_risk,
            risk_level=risk_level,
            flags=flags,
            details={
                "num_images": len(images),
                "duplicates_found": len(dup_results),
                "manipulation_results": manip_results,
            },
        )

    def _run_statistical_checks(self, paper: PaperInput) -> ModuleResult:
        """統計的不整合チェック"""
        results_text = paper.sections.get("results", "")
        full_text = " ".join(paper.sections.values())

        analysis = self.stat_analyzer.analyze_text(
            results_text if results_text else full_text
        )

        risk_score = analysis.overall_inconsistency_score

        return ModuleResult(
            module_name="statistical_checks",
            risk_score=min(risk_score, 1.0),
            risk_level=self._score_to_level(risk_score),
            flags=analysis.flags,
            details={
                "grim_inconsistencies": sum(
                    1 for r in analysis.grim_results if not r.is_consistent
                ),
                "grim_total": len(analysis.grim_results),
                "sprite_inconsistencies": sum(
                    1 for r in analysis.sprite_results if not r.is_consistent
                ),
                "sprite_total": len(analysis.sprite_results),
                "p_value_suspicious": (
                    analysis.p_value_analysis.suspicious_clustering
                    if analysis.p_value_analysis else False
                ),
                "num_tests_checked": analysis.num_tests_checked,
            },
        )

    def _run_plagiarism_check(self, paper: PaperInput) -> ModuleResult:
        """盗作チェック（コーパスなしの自己分析）"""
        sections = paper.sections
        flags = []

        # セクション間の自己類似度チェック
        fp_results = {}
        for name, text in sections.items():
            if len(text.strip()) > 50:
                fp = self.fingerprinter.fingerprint(text, name)
                fp_results[name] = fp

        # セクション間の不自然な類似度
        names = list(fp_results.keys())
        for i in range(len(names)):
            for j in range(i + 1, len(names)):
                sim = self.fingerprinter.similarity(
                    fp_results[names[i]], fp_results[names[j]]
                )
                if sim > 0.5:
                    flags.append(
                        f"High intra-paper similarity between "
                        f"{names[i]} and {names[j]}: {sim:.2f}"
                    )

        risk_score = min(len(flags) * 0.2, 1.0)

        return ModuleResult(
            module_name="plagiarism",
            risk_score=risk_score,
            risk_level=self._score_to_level(risk_score),
            flags=flags if flags else ["No plagiarism indicators (corpus comparison pending)"],
            details={
                "sections_analyzed": len(fp_results),
                "note": "Full plagiarism check requires reference corpus",
            },
        )

    def _run_phacking_analysis(self, paper: PaperInput) -> ModuleResult:
        """P-hacking/HARKing分析"""
        flags = []
        risk_score = 0.0

        # P-hacking
        p_values = paper.p_values or []
        if not p_values:
            # テキストからp値を抽出
            full_text = " ".join(paper.sections.values())
            extractions = self.stat_analyzer.extract_statistics(full_text)
            p_values = [
                e.p_value for e in extractions
                if e.p_value is not None and 0 < e.p_value < 1
            ]

        ph_result = None
        if len(p_values) >= 3:
            ph_result = self.phacking_detector.analyze(p_values)
            flags.extend(ph_result.flags)
            risk_score = max(risk_score, ph_result.overall_score)

        # HARKing
        intro = paper.sections.get("introduction", "")
        results = paper.sections.get("results", "")
        discussion = paper.sections.get("discussion", "")

        hk_result = None
        if intro and results:
            hk_result = self.harking_detector.analyze(
                intro, results, discussion
            )
            flags.extend(hk_result.flags)
            risk_score = max(risk_score, hk_result.risk_score)

        return ModuleResult(
            module_name="phacking",
            risk_score=min(risk_score, 1.0),
            risk_level=self._score_to_level(risk_score),
            flags=flags if flags else ["No P-hacking/HARKing indicators"],
            details={
                "num_p_values": len(p_values),
                "phacking_risk": (
                    ph_result.overall_risk if ph_result else "insufficient_data"
                ),
                "harking_risk": (
                    hk_result.risk_level if hk_result else "insufficient_data"
                ),
            },
        )

    def _run_reproducibility_assessment(
        self, paper: PaperInput,
        stat_result: ModuleResult,
        phack_result: ModuleResult,
    ) -> ModuleResult:
        """再現性予測"""
        methods = paper.sections.get("methods", "")
        full_text = " ".join(paper.sections.values())

        # 方法論評価
        method_score = self.method_assessor.assess(methods, full_text)

        # 統合再現性スコア
        repro = self.repro_scorer.score(
            methodology_quality=method_score.overall_score,
            statistical_consistency=1.0 - stat_result.risk_score,
            phacking_risk=phack_result.risk_score,
        )

        flags = repro.recommendations[:3]

        return ModuleResult(
            module_name="reproducibility",
            risk_score=1.0 - repro.overall_score,
            risk_level=self._score_to_level(1.0 - repro.overall_score),
            flags=flags,
            details={
                "reproducibility_prediction": repro.overall_score,
                "confidence_interval": repro.confidence_interval,
                "prediction_class": repro.prediction_class,
                "methodology_score": method_score.overall_score,
                "methodology_quality": method_score.quality_level,
                "missing_elements": method_score.missing_elements[:5],
            },
        )

    def _integrate_results(
        self, results: Dict[str, ModuleResult]
    ) -> tuple:
        """モジュール結果を統合"""
        weighted_risk = 0.0
        weight_sum = 0.0
        all_flags = []

        for name, weight in self.MODULE_WEIGHTS.items():
            if name in results:
                r = results[name]
                if r.risk_level not in ("not_applicable", "pending"):
                    weighted_risk += r.risk_score * weight
                    weight_sum += weight
                    all_flags.extend(r.flags)

        if weight_sum > 0:
            avg_risk = weighted_risk / weight_sum
        else:
            avg_risk = 0.0

        integrity_score = 1.0 - avg_risk

        if integrity_score < 0.3:
            level = "critical"
        elif integrity_score < 0.5:
            level = "high_risk"
        elif integrity_score < 0.7:
            level = "moderate_risk"
        else:
            level = "low_risk"

        # 重要な懸念事項をフィルタ
        concerns = [f for f in all_flags if not f.startswith("No ")]

        return integrity_score, level, concerns

    def _generate_recommendations(
        self, results: Dict[str, ModuleResult],
        overall_score: float,
    ) -> List[str]:
        """統合推奨事項"""
        recs = []

        if overall_score < 0.5:
            recs.append(
                "CRITICAL: Multiple integrity concerns detected. "
                "Detailed manual review is strongly recommended."
            )

        for name, result in results.items():
            if result.risk_level in ("high", "critical"):
                recs.append(
                    f"[{name}] High risk detected — "
                    f"investigate: {result.flags[0] if result.flags else 'see details'}"
                )

        if not recs:
            recs.append(
                "No major concerns detected. "
                "Standard quality checks passed."
            )

        return recs

    @staticmethod
    def _score_to_level(score: float) -> str:
        """スコアをリスクレベルに変換"""
        if score >= 0.7:
            return "critical"
        elif score >= 0.4:
            return "high"
        elif score >= 0.2:
            return "medium"
        else:
            return "low"
