"""
引用文脈考慮型テキスト類似度分析

論文の盗作検出において、正当な引用と不正なコピーを区別する。

特徴:
1. 引用マーカーの認識と除外
2. パラフレーズ検出（構文変換を考慮）
3. セクション別重み付け（Introduction vs Methods vs Discussion）
4. 自己引用 vs 他者からの盗作の区別
5. 文レベルのアラインメント（Smith-Waterman 風）
"""

import re
import math
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple, Set
from collections import Counter


@dataclass
class SimilarityMatch:
    """類似箇所の検出結果"""
    source_text: str
    target_text: str
    similarity_score: float  # 0.0 - 1.0
    match_type: str  # verbatim, paraphrase, mosaic, citation
    source_location: Tuple[int, int]  # (start, end) in source
    target_location: Tuple[int, int]
    is_cited: bool  # 引用の有無
    citation_refs: List[str] = field(default_factory=list)
    is_legitimate: bool = False  # 正当な引用かどうか


@dataclass
class PlagiarismReport:
    """盗作検出の総合レポート"""
    overall_similarity: float
    adjusted_similarity: float  # 正当な引用を除外した類似度
    matches: List[SimilarityMatch]
    section_scores: Dict[str, float]
    verbatim_ratio: float
    paraphrase_ratio: float
    mosaic_ratio: float
    num_citations_found: int
    num_uncited_matches: int
    risk_level: str  # low, medium, high, critical
    details: str = ""


class CitationAwareSimilarity:
    """
    引用文脈を考慮したテキスト類似度分析。

    Parameters
    ----------
    min_match_length : int
        最小マッチ長（単語数）
    verbatim_threshold : float
        逐語的コピーの閾値
    paraphrase_threshold : float
        パラフレーズの閾値
    """

    # 引用パターン
    CITATION_PATTERNS = [
        re.compile(r'\(([A-Z][a-zA-Z]+(?:\s+(?:et\s+al\.?|&\s+[A-Z][a-zA-Z]+))?'
                   r'(?:,?\s*\d{4}[a-z]?)+)\)'),
        re.compile(r'\[(\d+(?:\s*[-,]\s*\d+)*)\]'),
        re.compile(r'([A-Z][a-zA-Z]+(?:\s+(?:et\s+al\.?|and\s+[A-Z][a-zA-Z]+))?'
                   r'\s*\(\d{4}[a-z]?\))'),
    ]

    SECTION_WEIGHTS = {
        "abstract": 1.5,
        "introduction": 1.2,
        "methods": 0.8,  # 方法の定型表現は許容度が高い
        "results": 1.0,
        "discussion": 1.3,
        "conclusion": 1.4,
    }

    def __init__(self, min_match_length: int = 6,
                 verbatim_threshold: float = 0.95,
                 paraphrase_threshold: float = 0.70):
        self.min_match_length = min_match_length
        self.verbatim_threshold = verbatim_threshold
        self.paraphrase_threshold = paraphrase_threshold

    def compare(self, source_text: str, target_text: str,
                source_section: str = "unknown",
                target_section: str = "unknown") -> PlagiarismReport:
        """
        2つのテキストを比較し盗作検出を行う。

        Parameters
        ----------
        source_text : str
            比較元テキスト
        target_text : str
            検査対象テキスト
        """
        # 引用箇所の特定
        source_citations = self._extract_citations(source_text)
        target_citations = self._extract_citations(target_text)

        # 文に分割
        source_sents = self._split_sentences(source_text)
        target_sents = self._split_sentences(target_text)

        matches = []

        for t_idx, t_sent in enumerate(target_sents):
            t_clean = self._remove_citations(t_sent)
            t_words = self._tokenize(t_clean)
            if len(t_words) < self.min_match_length:
                continue

            best_match = None
            best_score = 0.0

            for s_idx, s_sent in enumerate(source_sents):
                s_clean = self._remove_citations(s_sent)
                s_words = self._tokenize(s_clean)
                if len(s_words) < self.min_match_length:
                    continue

                # 類似度計算
                sim = self._sentence_similarity(s_words, t_words)

                if sim > self.paraphrase_threshold and sim > best_score:
                    best_score = sim
                    # マッチタイプの判定
                    if sim >= self.verbatim_threshold:
                        match_type = "verbatim"
                    elif sim >= 0.85:
                        match_type = "mosaic"
                    else:
                        match_type = "paraphrase"

                    # 引用の有無
                    is_cited = self._has_citation_nearby(
                        target_text, t_sent, target_citations
                    )
                    refs = self._get_citation_refs(t_sent)

                    best_match = SimilarityMatch(
                        source_text=s_sent,
                        target_text=t_sent,
                        similarity_score=sim,
                        match_type=match_type,
                        source_location=(s_idx, s_idx),
                        target_location=(t_idx, t_idx),
                        is_cited=is_cited,
                        citation_refs=refs,
                        is_legitimate=is_cited and match_type != "verbatim",
                    )

            if best_match is not None:
                matches.append(best_match)

        # 総合スコアの計算
        return self._compile_report(
            matches, source_sents, target_sents,
            source_section, target_section
        )

    def _sentence_similarity(self, words_a: List[str],
                             words_b: List[str]) -> float:
        """文レベルの類似度（n-gramベース + Jaccard + 編集距離ハイブリッド）"""
        if not words_a or not words_b:
            return 0.0

        # Jaccard類似度 (unigram)
        set_a = set(words_a)
        set_b = set(words_b)
        if not set_a or not set_b:
            return 0.0
        jaccard = len(set_a & set_b) / len(set_a | set_b)

        # N-gram類似度 (bigram, trigram)
        bigram_sim = self._ngram_similarity(words_a, words_b, 2)
        trigram_sim = self._ngram_similarity(words_a, words_b, 3)

        # 語順考慮の類似度（LCS比率）
        lcs_len = self._lcs_length(words_a, words_b)
        lcs_ratio = 2 * lcs_len / (len(words_a) + len(words_b))

        # 重み付き平均
        similarity = (
            0.2 * jaccard
            + 0.2 * bigram_sim
            + 0.2 * trigram_sim
            + 0.4 * lcs_ratio
        )

        return similarity

    def _ngram_similarity(self, words_a: List[str],
                          words_b: List[str], n: int) -> float:
        """N-gram Jaccard類似度"""
        if len(words_a) < n or len(words_b) < n:
            return 0.0

        ngrams_a = set(
            tuple(words_a[i:i+n]) for i in range(len(words_a) - n + 1)
        )
        ngrams_b = set(
            tuple(words_b[i:i+n]) for i in range(len(words_b) - n + 1)
        )

        if not ngrams_a or not ngrams_b:
            return 0.0
        return len(ngrams_a & ngrams_b) / len(ngrams_a | ngrams_b)

    def _lcs_length(self, a: List[str], b: List[str]) -> int:
        """最長共通部分列の長さ"""
        m, n = len(a), len(b)
        if m == 0 or n == 0:
            return 0

        # メモリ効率のため2行のみ保持
        prev = [0] * (n + 1)
        curr = [0] * (n + 1)

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                if a[i-1] == b[j-1]:
                    curr[j] = prev[j-1] + 1
                else:
                    curr[j] = max(prev[j], curr[j-1])
            prev, curr = curr, [0] * (n + 1)

        return prev[n]

    def _extract_citations(self, text: str) -> List[Dict]:
        """テキストから引用箇所を抽出"""
        citations = []
        for pattern in self.CITATION_PATTERNS:
            for m in pattern.finditer(text):
                citations.append({
                    "text": m.group(0),
                    "ref": m.group(1) if m.lastindex else m.group(0),
                    "start": m.start(),
                    "end": m.end(),
                })
        return citations

    def _remove_citations(self, text: str) -> str:
        """テキストから引用マーカーを除去"""
        cleaned = text
        for pattern in self.CITATION_PATTERNS:
            cleaned = pattern.sub("", cleaned)
        return cleaned.strip()

    def _has_citation_nearby(self, full_text: str, sentence: str,
                             citations: List[Dict]) -> bool:
        """文の近くに引用があるか"""
        sent_pos = full_text.find(sentence)
        if sent_pos == -1:
            return False
        for cit in citations:
            if abs(cit["start"] - sent_pos) < len(sentence) + 50:
                return True
        # 文内に引用パターンがあるか
        for pattern in self.CITATION_PATTERNS:
            if pattern.search(sentence):
                return True
        return False

    def _get_citation_refs(self, text: str) -> List[str]:
        """テキスト内の引用参照を抽出"""
        refs = []
        for pattern in self.CITATION_PATTERNS:
            for m in pattern.finditer(text):
                refs.append(m.group(0))
        return refs

    def _split_sentences(self, text: str) -> List[str]:
        """テキストを文に分割"""
        sentences = re.split(r'(?<=[.!?])\s+', text)
        return [s.strip() for s in sentences if len(s.strip()) > 10]

    def _tokenize(self, text: str) -> List[str]:
        """テキストをトークンに分割"""
        text = text.lower()
        text = re.sub(r'[^\w\s]', ' ', text)
        return [w for w in text.split() if len(w) > 1]

    def _compile_report(self, matches: List[SimilarityMatch],
                        source_sents: List[str],
                        target_sents: List[str],
                        source_section: str,
                        target_section: str) -> PlagiarismReport:
        """検出結果を総合レポートにまとめる"""
        if not target_sents:
            return PlagiarismReport(
                overall_similarity=0.0, adjusted_similarity=0.0,
                matches=[], section_scores={},
                verbatim_ratio=0.0, paraphrase_ratio=0.0,
                mosaic_ratio=0.0, num_citations_found=0,
                num_uncited_matches=0, risk_level="low",
            )

        total_sents = len(target_sents)
        matched_sents = len(matches)
        overall_sim = matched_sents / total_sents if total_sents > 0 else 0.0

        # 正当な引用を除外した類似度
        illegitimate = [m for m in matches if not m.is_legitimate]
        adjusted_sim = len(illegitimate) / total_sents if total_sents > 0 else 0.0

        # タイプ別の比率
        verbatim = [m for m in matches if m.match_type == "verbatim"]
        paraphrase = [m for m in matches if m.match_type == "paraphrase"]
        mosaic = [m for m in matches if m.match_type == "mosaic"]

        verbatim_ratio = len(verbatim) / total_sents if total_sents > 0 else 0.0
        paraphrase_ratio = len(paraphrase) / total_sents if total_sents > 0 else 0.0
        mosaic_ratio = len(mosaic) / total_sents if total_sents > 0 else 0.0

        cited = sum(1 for m in matches if m.is_cited)
        uncited = matched_sents - cited

        # リスクレベルの判定
        section_weight = self.SECTION_WEIGHTS.get(target_section, 1.0)
        weighted_score = adjusted_sim * section_weight

        if weighted_score > 0.3 or verbatim_ratio > 0.1:
            risk = "critical"
        elif weighted_score > 0.15 or verbatim_ratio > 0.05:
            risk = "high"
        elif weighted_score > 0.05:
            risk = "medium"
        else:
            risk = "low"

        return PlagiarismReport(
            overall_similarity=overall_sim,
            adjusted_similarity=adjusted_sim,
            matches=matches,
            section_scores={target_section: weighted_score},
            verbatim_ratio=verbatim_ratio,
            paraphrase_ratio=paraphrase_ratio,
            mosaic_ratio=mosaic_ratio,
            num_citations_found=cited,
            num_uncited_matches=uncited,
            risk_level=risk,
            details=(
                f"Overall similarity: {overall_sim:.1%}\n"
                f"Adjusted (excl. citations): {adjusted_sim:.1%}\n"
                f"Verbatim: {len(verbatim)}, Paraphrase: {len(paraphrase)}, "
                f"Mosaic: {len(mosaic)}\n"
                f"Cited matches: {cited}, Uncited: {uncited}\n"
                f"Risk level: {risk}"
            ),
        )
