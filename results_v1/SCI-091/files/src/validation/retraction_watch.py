"""
Retraction Watch データベース分析モジュール

Retraction Watch Database から撤回論文データを取得・分析し、
検出システムの検証用データセットとして活用する。

主要な撤回理由カテゴリ:
- 画像の問題 (Image Issues)
- データの捏造・改竄 (Fabrication/Falsification)
- 盗作 (Plagiarism)
- 統計エラー (Statistical Errors)
- 著者の問題 (Authorship Issues)
- 倫理的問題 (Ethical Issues)
"""

import re
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple
from collections import Counter
from datetime import datetime


@dataclass
class RetractionRecord:
    """撤回論文レコード"""
    doi: str
    title: str
    authors: List[str]
    journal: str
    retraction_date: str
    original_pub_date: str
    retraction_reasons: List[str]
    retraction_type: str  # retraction, correction, expression_of_concern
    is_self_retraction: bool
    has_investigation: bool
    country: str = ""
    subject_area: str = ""
    citations_at_retraction: int = 0
    time_to_retraction_days: int = 0


@dataclass
class RetractionStatistics:
    """撤回統計"""
    total_retractions: int
    by_reason: Dict[str, int]
    by_year: Dict[int, int]
    by_journal: Dict[str, int]
    by_country: Dict[str, int]
    by_subject: Dict[str, int]
    mean_time_to_retraction: float
    median_citations: float
    self_retraction_rate: float


class RetractionWatchAnalyzer:
    """
    Retraction Watch データの分析と検証用データセット構築。

    機能:
    1. 撤回データの統計分析
    2. 撤回理由によるフィルタリング
    3. 検出システムの検証用正例/負例データセット構築
    4. 検出性能の評価指標計算
    """

    REASON_CATEGORIES = {
        "image_issues": [
            "image manipulation", "figure manipulation", "duplicate image",
            "image fabrication", "doctored images",
        ],
        "data_fabrication": [
            "fabrication", "falsification", "data fabrication",
            "fake data", "manipulated data",
        ],
        "plagiarism": [
            "plagiarism", "duplication of text", "text overlap",
            "copied text", "self-plagiarism",
        ],
        "statistical_errors": [
            "errors in data", "statistical errors", "calculation errors",
            "incorrect results", "errors in analysis",
        ],
        "ethical_issues": [
            "ethical violations", "no IRB approval", "informed consent",
            "conflict of interest",
        ],
        "authorship": [
            "authorship dispute", "fake authorship",
            "undisclosed conflict",
        ],
        "reproducibility": [
            "unreproducible", "cannot be reproduced",
            "failed replication",
        ],
    }

    def __init__(self):
        self._records: List[RetractionRecord] = []

    def load_records(self, records: List[Dict]):
        """レコードをロード"""
        for r in records:
            reasons = r.get("reasons", [])
            categorized = self._categorize_reasons(reasons)

            self._records.append(RetractionRecord(
                doi=r.get("doi", ""),
                title=r.get("title", ""),
                authors=r.get("authors", []),
                journal=r.get("journal", ""),
                retraction_date=r.get("retraction_date", ""),
                original_pub_date=r.get("original_pub_date", ""),
                retraction_reasons=categorized,
                retraction_type=r.get("type", "retraction"),
                is_self_retraction=r.get("self_retraction", False),
                has_investigation=r.get("investigation", False),
                country=r.get("country", ""),
                subject_area=r.get("subject", ""),
                citations_at_retraction=r.get("citations", 0),
                time_to_retraction_days=r.get("time_to_retraction", 0),
            ))

    def compute_statistics(self) -> RetractionStatistics:
        """撤回データの統計を計算"""
        if not self._records:
            return RetractionStatistics(
                total_retractions=0, by_reason={}, by_year={},
                by_journal={}, by_country={}, by_subject={},
                mean_time_to_retraction=0.0, median_citations=0.0,
                self_retraction_rate=0.0,
            )

        by_reason = Counter()
        by_year = Counter()
        by_journal = Counter()
        by_country = Counter()
        by_subject = Counter()

        times = []
        citations = []
        self_count = 0

        for r in self._records:
            for reason in r.retraction_reasons:
                by_reason[reason] += 1
            if r.retraction_date:
                try:
                    year = int(r.retraction_date[:4])
                    by_year[year] += 1
                except (ValueError, IndexError):
                    pass
            by_journal[r.journal] += 1
            if r.country:
                by_country[r.country] += 1
            if r.subject_area:
                by_subject[r.subject_area] += 1
            if r.time_to_retraction_days > 0:
                times.append(r.time_to_retraction_days)
            citations.append(r.citations_at_retraction)
            if r.is_self_retraction:
                self_count += 1

        import numpy as np
        mean_time = float(np.mean(times)) if times else 0.0
        median_cite = float(np.median(citations)) if citations else 0.0

        return RetractionStatistics(
            total_retractions=len(self._records),
            by_reason=dict(by_reason.most_common()),
            by_year=dict(sorted(by_year.items())),
            by_journal=dict(by_journal.most_common(20)),
            by_country=dict(by_country.most_common(20)),
            by_subject=dict(by_subject.most_common(20)),
            mean_time_to_retraction=mean_time,
            median_citations=median_cite,
            self_retraction_rate=self_count / len(self._records),
        )

    def build_validation_set(
        self,
        reason_filter: Optional[List[str]] = None,
        max_samples: int = 1000,
    ) -> Dict[str, List[str]]:
        """
        検証用データセットを構築。

        Returns
        -------
        dict
            'positive': 撤回論文のDOIリスト
            'negative': 対照群の指定（同じジャーナル・時期の非撤回論文）
        """
        positives = []
        for r in self._records:
            if reason_filter:
                if any(reason in r.retraction_reasons for reason in reason_filter):
                    positives.append(r.doi)
            else:
                positives.append(r.doi)

        positives = positives[:max_samples]

        return {
            "positive": positives,
            "negative_criteria": {
                "strategy": "matched_controls",
                "matching_variables": ["journal", "year", "subject_area"],
                "ratio": "1:3",
                "description": (
                    "For each retracted paper, select 3 non-retracted papers "
                    "from the same journal, year, and subject area"
                ),
            },
        }

    def evaluate_detector(
        self,
        predictions: Dict[str, float],
        ground_truth: Dict[str, bool],
        threshold: float = 0.5,
    ) -> Dict[str, float]:
        """
        検出器の性能を評価。

        Parameters
        ----------
        predictions : dict
            DOI -> 予測スコア (0.0-1.0)
        ground_truth : dict
            DOI -> 撤回された (True/False)
        threshold : float
            判定閾値
        """
        tp = fp = tn = fn = 0

        for doi, score in predictions.items():
            if doi not in ground_truth:
                continue
            predicted = score >= threshold
            actual = ground_truth[doi]

            if predicted and actual:
                tp += 1
            elif predicted and not actual:
                fp += 1
            elif not predicted and actual:
                fn += 1
            else:
                tn += 1

        total = tp + fp + tn + fn
        if total == 0:
            return {"error": "No overlapping DOIs"}

        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1 = (
            2 * precision * recall / (precision + recall)
            if (precision + recall) > 0 else 0.0
        )
        accuracy = (tp + tn) / total
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0

        # AUROCの概算（閾値を変化させて）
        thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        tpr_list = []
        fpr_list = []
        for t in thresholds:
            tp_t = sum(
                1 for d, s in predictions.items()
                if s >= t and ground_truth.get(d, False)
            )
            fn_t = sum(
                1 for d, s in predictions.items()
                if s < t and ground_truth.get(d, False)
            )
            fp_t = sum(
                1 for d, s in predictions.items()
                if s >= t and not ground_truth.get(d, True)
            )
            tn_t = sum(
                1 for d, s in predictions.items()
                if s < t and not ground_truth.get(d, True)
            )
            tpr_t = tp_t / (tp_t + fn_t) if (tp_t + fn_t) > 0 else 0.0
            fpr_t = fp_t / (fp_t + tn_t) if (fp_t + tn_t) > 0 else 0.0
            tpr_list.append(tpr_t)
            fpr_list.append(fpr_t)

        # 台形則でAUROC近似
        auroc = 0.5  # デフォルト
        if len(fpr_list) > 1:
            import numpy as np
            sorted_pairs = sorted(zip(fpr_list, tpr_list))
            fpr_sorted = [p[0] for p in sorted_pairs]
            tpr_sorted = [p[1] for p in sorted_pairs]
            auroc = float(np.trapz(tpr_sorted, fpr_sorted))
            auroc = max(0.0, min(1.0, abs(auroc)))

        return {
            "accuracy": accuracy,
            "precision": precision,
            "recall": recall,
            "f1_score": f1,
            "specificity": specificity,
            "auroc": auroc,
            "true_positives": tp,
            "false_positives": fp,
            "true_negatives": tn,
            "false_negatives": fn,
        }

    def _categorize_reasons(self, reasons: List[str]) -> List[str]:
        """撤回理由をカテゴリに分類"""
        categories = set()
        for reason in reasons:
            reason_lower = reason.lower()
            for cat, keywords in self.REASON_CATEGORIES.items():
                if any(kw in reason_lower for kw in keywords):
                    categories.add(cat)
        if not categories:
            categories.add("other")
        return list(categories)
