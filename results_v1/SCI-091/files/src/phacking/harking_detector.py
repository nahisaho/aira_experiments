"""
HARKing (Hypothesizing After Results are Known) 検出モジュール

検出指標:
1. 仮説の具体性と結果の一致度（過度な一致は疑わしい）
2. Introduction-Results-Discussion の論理的整合性
3. 予測的表現 vs 探索的表現の比率
4. 効果量の大きさと仮説の精度の不自然な相関
5. 事前登録との比較（利用可能な場合）
"""

import re
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple
from collections import Counter


@dataclass
class HARKingResult:
    """HARKing検出結果"""
    risk_level: str  # low, moderate, high
    risk_score: float  # 0.0 - 1.0
    predictive_ratio: float  # 予測的表現の割合
    exploratory_ratio: float  # 探索的表現の割合
    hypothesis_result_alignment: float  # 仮説と結果の一致度
    narrative_consistency: float  # ストーリーの一貫性
    flags: List[str] = field(default_factory=list)
    linguistic_markers: Dict[str, int] = field(default_factory=dict)
    preregistration_deviations: List[str] = field(default_factory=list)
    details: str = ""


class HARKingDetector:
    """
    HARKing検出器。

    テキスト分析により事後的仮説構築の兆候を検出する。

    検出手法:
    1. 言語マーカー分析（予測的 vs 探索的表現）
    2. Introduction-Results整合性分析
    3. 仮説-結果アラインメント評価
    4. 事前登録比較（OSF/AsPredicted）
    """

    # 予測的表現パターン（HARKingの兆候となりうる過度に確信的な言語）
    PREDICTIVE_MARKERS = [
        r'\b(?:we predicted|we hypothesized|as expected|consistent with our hypothesis)\b',
        r'\b(?:as anticipated|we expected|in line with our prediction)\b',
        r'\b(?:confirming our hypothesis|supporting our prediction)\b',
        r'\b(?:we proposed that|our model predicts)\b',
    ]

    # 探索的表現パターン（正直な探索的研究のマーカー）
    EXPLORATORY_MARKERS = [
        r'\b(?:unexpectedly|surprisingly|we did not expect)\b',
        r'\b(?:contrary to our expectations|an unexpected finding)\b',
        r'\b(?:exploratory analysis|post-hoc analysis)\b',
        r'\b(?:we explored|we examined whether|to our surprise)\b',
        r'\b(?:interestingly|notably|remarkably)\b',
    ]

    # 曖昧化表現（結果に合わせた曖昧な仮説）
    HEDGING_MARKERS = [
        r'\b(?:may|might|could|possibly|potentially)\b',
        r'\b(?:it is possible that|one explanation is)\b',
        r'\b(?:we speculated|we reasoned)\b',
        r'\b(?:it seems|appears to|tends to)\b',
    ]

    # 結果の強調パターン
    EMPHASIS_MARKERS = [
        r'\b(?:importantly|critically|crucially|strikingly)\b',
        r'\b(?:robust|strong|clear|compelling|definitive)\b',
        r'\b(?:unambiguous|decisive|conclusive)\b',
    ]

    def __init__(self):
        self.predictive_re = [re.compile(p, re.IGNORECASE) for p in self.PREDICTIVE_MARKERS]
        self.exploratory_re = [re.compile(p, re.IGNORECASE) for p in self.EXPLORATORY_MARKERS]
        self.hedging_re = [re.compile(p, re.IGNORECASE) for p in self.HEDGING_MARKERS]
        self.emphasis_re = [re.compile(p, re.IGNORECASE) for p in self.EMPHASIS_MARKERS]

    def analyze(self, introduction: str, results: str,
                discussion: str = "",
                preregistration: Optional[str] = None) -> HARKingResult:
        """
        HARKing分析を実行。

        Parameters
        ----------
        introduction : str
            序論セクション
        results : str
            結果セクション
        discussion : str
            考察セクション
        preregistration : str, optional
            事前登録テキスト
        """
        flags = []
        markers = {}

        # 1. 言語マーカー分析
        intro_markers = self._count_markers(introduction)
        results_markers = self._count_markers(results)
        disc_markers = self._count_markers(discussion)

        total_text = introduction + " " + results + " " + discussion
        total_markers = self._count_markers(total_text)
        markers = total_markers

        # 予測的・探索的表現の比率
        pred_count = total_markers.get("predictive", 0)
        expl_count = total_markers.get("exploratory", 0)
        total_assertions = pred_count + expl_count + 1

        predictive_ratio = pred_count / total_assertions
        exploratory_ratio = expl_count / total_assertions

        # 2. 過度に予測的（HARKingの兆候）
        if predictive_ratio > 0.8 and pred_count >= 3:
            flags.append(
                f"Excessively predictive language: {pred_count} predictive "
                f"vs {expl_count} exploratory markers"
            )

        # 探索的表現が皆無（不自然）
        if expl_count == 0 and pred_count >= 3:
            flags.append(
                "No exploratory language detected despite multiple predictions"
            )

        # 3. 仮説-結果アラインメント
        alignment = self._hypothesis_result_alignment(introduction, results)
        if alignment > 0.9:
            flags.append(
                f"Suspiciously high hypothesis-result alignment: {alignment:.2f}"
            )

        # 4. Introduction内の曖昧な仮説
        hedging_in_intro = sum(
            len(p.findall(introduction)) for p in self.hedging_re
        )
        emphasis_in_results = sum(
            len(p.findall(results)) for p in self.emphasis_re
        )

        if hedging_in_intro > 3 and emphasis_in_results > 3:
            flags.append(
                "Pattern: vague hypotheses + emphatic results "
                "(potential narrative shaping)"
            )

        # 5. 事前登録との比較
        prereg_deviations = []
        if preregistration:
            prereg_deviations = self._compare_preregistration(
                preregistration, introduction, results
            )
            if prereg_deviations:
                flags.extend(
                    f"Preregistration deviation: {d}"
                    for d in prereg_deviations
                )

        # 6. ストーリーの一貫性分析
        narrative_consistency = self._narrative_consistency(
            introduction, results, discussion
        )

        # リスクスコアの計算
        risk_score = 0.0
        risk_score += 0.2 * min(predictive_ratio * 1.5, 1.0)
        risk_score += 0.2 * max(0, alignment - 0.7) / 0.3
        risk_score += 0.1 * (1.0 if expl_count == 0 and pred_count >= 2 else 0.0)
        risk_score += 0.2 * min(hedging_in_intro / 5, 1.0)
        risk_score += 0.15 * min(len(prereg_deviations) / 3, 1.0)
        risk_score += 0.15 * max(0, narrative_consistency - 0.8) / 0.2
        risk_score = min(risk_score, 1.0)

        if risk_score > 0.6:
            risk_level = "high"
        elif risk_score > 0.3:
            risk_level = "moderate"
        else:
            risk_level = "low"

        return HARKingResult(
            risk_level=risk_level,
            risk_score=risk_score,
            predictive_ratio=predictive_ratio,
            exploratory_ratio=exploratory_ratio,
            hypothesis_result_alignment=alignment,
            narrative_consistency=narrative_consistency,
            flags=flags,
            linguistic_markers=markers,
            preregistration_deviations=prereg_deviations,
            details=(
                f"HARKing Analysis: risk={risk_level} ({risk_score:.2f})\n"
                f"Predictive ratio: {predictive_ratio:.1%}\n"
                f"Exploratory ratio: {exploratory_ratio:.1%}\n"
                f"Hypothesis-result alignment: {alignment:.2f}\n"
                f"Narrative consistency: {narrative_consistency:.2f}\n"
                f"Flags: {len(flags)}"
            ),
        )

    def _count_markers(self, text: str) -> Dict[str, int]:
        """テキスト中の各種マーカーをカウント"""
        counts = {
            "predictive": sum(
                len(p.findall(text)) for p in self.predictive_re
            ),
            "exploratory": sum(
                len(p.findall(text)) for p in self.exploratory_re
            ),
            "hedging": sum(
                len(p.findall(text)) for p in self.hedging_re
            ),
            "emphasis": sum(
                len(p.findall(text)) for p in self.emphasis_re
            ),
        }
        return counts

    def _hypothesis_result_alignment(self, introduction: str,
                                     results: str) -> float:
        """仮説と結果の一致度を評価"""
        # Introduction内の主要な名詞句を抽出
        intro_keywords = self._extract_keywords(introduction)
        result_keywords = self._extract_keywords(results)

        if not intro_keywords or not result_keywords:
            return 0.5

        # キーワードの重複度
        overlap = len(intro_keywords & result_keywords)
        total = len(intro_keywords | result_keywords)

        return overlap / total if total > 0 else 0.5

    def _extract_keywords(self, text: str, min_len: int = 4) -> set:
        """テキストからキーワードを抽出"""
        words = re.findall(r'\b[a-zA-Z]{' + str(min_len) + r',}\b', text.lower())
        # ストップワード除去
        stopwords = {
            'that', 'this', 'with', 'from', 'were', 'have', 'been',
            'their', 'which', 'would', 'than', 'they', 'more',
            'also', 'between', 'other', 'these', 'into', 'some',
            'when', 'will', 'each', 'about', 'such',
        }
        return set(w for w in words if w not in stopwords)

    def _compare_preregistration(self, prereg: str, intro: str,
                                 results: str) -> List[str]:
        """事前登録との比較"""
        deviations = []

        prereg_keywords = self._extract_keywords(prereg)
        intro_keywords = self._extract_keywords(intro)
        result_keywords = self._extract_keywords(results)

        # 事前登録にないが結果にある分析
        new_in_results = result_keywords - prereg_keywords
        if len(new_in_results) > len(prereg_keywords) * 0.5:
            deviations.append(
                "Substantial new analyses not in preregistration"
            )

        # 事前登録にあるが結果にない分析
        missing = prereg_keywords - result_keywords - intro_keywords
        if len(missing) > len(prereg_keywords) * 0.3:
            deviations.append(
                "Preregistered analyses appear to be unreported"
            )

        return deviations

    def _narrative_consistency(self, intro: str, results: str,
                               discussion: str) -> float:
        """ストーリーの一貫性を評価（過度に一貫的＝疑わしい）"""
        if not discussion:
            return 0.5

        intro_kw = self._extract_keywords(intro)
        result_kw = self._extract_keywords(results)
        disc_kw = self._extract_keywords(discussion)

        if not intro_kw or not disc_kw:
            return 0.5

        # 3セクション間のキーワード重複
        all_overlap = intro_kw & result_kw & disc_kw
        all_union = intro_kw | result_kw | disc_kw

        if not all_union:
            return 0.5

        return len(all_overlap) / len(all_union)
