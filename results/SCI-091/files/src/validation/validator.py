"""
研究公正性検証パイプライン統合バリデータ

全モジュールの結果を統合し、PubPeer/Retraction Watch データで
システム全体の検出性能を評価する。
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional
from .pubpeer_client import PubPeerClient
from .retraction_watch import RetractionWatchAnalyzer


@dataclass
class ValidationResult:
    """統合検証結果"""
    num_papers_validated: int
    detection_performance: Dict[str, float]
    module_performance: Dict[str, Dict[str, float]]
    confusion_matrix: Dict[str, int]
    false_positive_analysis: List[Dict]
    false_negative_analysis: List[Dict]
    external_signal_correlation: float
    recommendations: List[str]


class IntegrityValidator:
    """
    研究公正性検出システムの統合バリデータ。

    検証プロセス:
    1. PubPeerで既知の問題論文を収集
    2. Retraction Watch で撤回論文のグラウンドトゥルースを構築
    3. 各検出モジュールの予測と比較
    4. 性能指標の計算と分析
    """

    def __init__(self):
        self.pubpeer = PubPeerClient()
        self.retraction_watch = RetractionWatchAnalyzer()

    def validate_system(
        self,
        predictions: Dict[str, Dict[str, float]],
        ground_truth: Dict[str, bool],
    ) -> ValidationResult:
        """
        システム全体の検証を実行。

        Parameters
        ----------
        predictions : dict
            DOI -> { module_name: score } の辞書
        ground_truth : dict
            DOI -> 撤回/問題あり (True/False)
        """
        # 統合スコアの計算
        integrated_scores = {}
        for doi, module_scores in predictions.items():
            integrated_scores[doi] = sum(module_scores.values()) / max(len(module_scores), 1)

        # 全体性能
        overall_perf = self.retraction_watch.evaluate_detector(
            integrated_scores, ground_truth
        )

        # モジュール別性能
        module_perf = {}
        module_names = set()
        for scores in predictions.values():
            module_names.update(scores.keys())

        for module in module_names:
            module_scores = {
                doi: scores.get(module, 0.0)
                for doi, scores in predictions.items()
            }
            module_perf[module] = self.retraction_watch.evaluate_detector(
                module_scores, ground_truth
            )

        # 偽陽性・偽陰性分析
        fp_analysis = []
        fn_analysis = []
        threshold = 0.5

        for doi, score in integrated_scores.items():
            actual = ground_truth.get(doi, False)
            if score >= threshold and not actual:
                fp_analysis.append({
                    "doi": doi,
                    "score": score,
                    "module_scores": predictions.get(doi, {}),
                    "note": "False positive — flagged but not retracted",
                })
            elif score < threshold and actual:
                fn_analysis.append({
                    "doi": doi,
                    "score": score,
                    "module_scores": predictions.get(doi, {}),
                    "note": "False negative — retracted but not detected",
                })

        # 外部シグナルとの相関
        ext_corr = self._external_signal_correlation(
            integrated_scores, ground_truth
        )

        recs = self._generate_validation_recommendations(
            overall_perf, module_perf, len(fp_analysis), len(fn_analysis)
        )

        return ValidationResult(
            num_papers_validated=len(predictions),
            detection_performance=overall_perf,
            module_performance=module_perf,
            confusion_matrix={
                "TP": overall_perf.get("true_positives", 0),
                "FP": overall_perf.get("false_positives", 0),
                "TN": overall_perf.get("true_negatives", 0),
                "FN": overall_perf.get("false_negatives", 0),
            },
            false_positive_analysis=fp_analysis[:20],
            false_negative_analysis=fn_analysis[:20],
            external_signal_correlation=ext_corr,
            recommendations=recs,
        )

    def _external_signal_correlation(
        self, scores: Dict[str, float], truth: Dict[str, bool]
    ) -> float:
        """検出スコアと外部シグナルの相関"""
        import numpy as np
        common_dois = set(scores.keys()) & set(truth.keys())
        if len(common_dois) < 5:
            return 0.0

        x = [scores[d] for d in common_dois]
        y = [1.0 if truth[d] else 0.0 for d in common_dois]

        x_arr = np.array(x)
        y_arr = np.array(y)

        if np.std(x_arr) == 0 or np.std(y_arr) == 0:
            return 0.0

        corr = float(np.corrcoef(x_arr, y_arr)[0, 1])
        return corr

    def _generate_validation_recommendations(
        self, overall: Dict, module: Dict,
        n_fp: int, n_fn: int
    ) -> List[str]:
        """検証結果に基づく推奨事項"""
        recs = []

        f1 = overall.get("f1_score", 0)
        if f1 < 0.5:
            recs.append(
                "Overall F1 score is below 0.5 — "
                "model requires significant improvement"
            )

        precision = overall.get("precision", 0)
        recall = overall.get("recall", 0)

        if precision < 0.5:
            recs.append(
                f"High false positive rate ({n_fp} FPs). "
                "Consider raising detection thresholds."
            )

        if recall < 0.5:
            recs.append(
                f"High false negative rate ({n_fn} FNs). "
                "Detection sensitivity needs improvement."
            )

        # 最良・最悪のモジュール
        if module:
            best = max(module.items(), key=lambda x: x[1].get("f1_score", 0))
            worst = min(module.items(), key=lambda x: x[1].get("f1_score", 0))
            recs.append(
                f"Best module: {best[0]} (F1={best[1].get('f1_score', 0):.2f})"
            )
            recs.append(
                f"Weakest module: {worst[0]} "
                f"(F1={worst[1].get('f1_score', 0):.2f}) — prioritize improvement"
            )

        return recs
