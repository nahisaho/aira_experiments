"""
盗作検出統合モジュール

パイプライン:
1. テキストフィンガープリントで候補文書を高速検索
2. 引用文脈考慮型の詳細比較
3. セクション別・文レベルの類似度マッピング
4. リスクスコアの統合判定
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple
from .citation_aware_similarity import (
    CitationAwareSimilarity, PlagiarismReport, SimilarityMatch
)
from .text_fingerprint import TextFingerprinter, FingerprintResult


@dataclass
class DocumentProfile:
    """文書のプロファイル"""
    doc_id: str
    title: str
    sections: Dict[str, str]  # section_name -> text
    fingerprint: Optional[FingerprintResult] = None
    authors: List[str] = field(default_factory=list)
    references: List[str] = field(default_factory=list)


@dataclass
class ComprehensivePlagiarismResult:
    """包括的盗作検出結果"""
    target_doc_id: str
    overall_risk: str  # low, medium, high, critical
    overall_score: float  # 0.0 - 1.0
    section_reports: Dict[str, PlagiarismReport]
    top_matches: List[Dict]  # 最も類似度の高い文書
    self_plagiarism_flag: bool
    verbatim_count: int
    paraphrase_count: int
    total_sentences_checked: int
    summary: str


class PlagiarismDetector:
    """
    包括的な盗作検出システム。

    ステージ1: フィンガープリントによる高速スクリーニング
    ステージ2: 候補文書との詳細比較（引用文脈考慮）
    ステージ3: セクション別リスク評価

    Parameters
    ----------
    fingerprint_k : int
        フィンガープリントのk-gramサイズ
    min_similarity : float
        詳細比較の閾値
    max_candidates : int
        詳細比較する候補文書の最大数
    """

    def __init__(self, fingerprint_k: int = 5,
                 min_similarity: float = 0.1,
                 max_candidates: int = 20):
        self.fingerprinter = TextFingerprinter(k=fingerprint_k)
        self.similarity_engine = CitationAwareSimilarity()
        self.min_similarity = min_similarity
        self.max_candidates = max_candidates
        self._corpus: Dict[str, DocumentProfile] = {}

    def add_to_corpus(self, doc: DocumentProfile):
        """比較用コーパスに文書を追加"""
        full_text = " ".join(doc.sections.values())
        doc.fingerprint = self.fingerprinter.fingerprint(full_text, doc.doc_id)
        self.fingerprinter.index_document(doc.fingerprint)
        self._corpus[doc.doc_id] = doc

    def check_document(self, target: DocumentProfile) -> ComprehensivePlagiarismResult:
        """文書の盗作チェックを実行"""
        # ステージ1: フィンガープリントスクリーニング
        full_text = " ".join(target.sections.values())
        target_fp = self.fingerprinter.fingerprint(full_text, target.doc_id)
        candidates = self.fingerprinter.query(target_fp, min_shared=2)

        # 候補を制限
        candidates = candidates[:self.max_candidates]

        # ステージ2: セクション別詳細比較
        section_reports: Dict[str, PlagiarismReport] = {}
        all_matches: List[Dict] = []

        for section_name, section_text in target.sections.items():
            if len(section_text.strip()) < 50:
                continue

            best_report = None
            best_score = 0.0

            for candidate_id, _ in candidates:
                if candidate_id not in self._corpus:
                    continue
                candidate = self._corpus[candidate_id]

                for c_section, c_text in candidate.sections.items():
                    if len(c_text.strip()) < 50:
                        continue

                    report = self.similarity_engine.compare(
                        c_text, section_text, c_section, section_name
                    )

                    if report.adjusted_similarity > best_score:
                        best_score = report.adjusted_similarity
                        best_report = report
                        all_matches.append({
                            "source_doc": candidate_id,
                            "source_section": c_section,
                            "target_section": section_name,
                            "similarity": report.adjusted_similarity,
                            "risk": report.risk_level,
                        })

            if best_report:
                section_reports[section_name] = best_report

        # ステージ3: 統合評価
        self_plagiarism = self._check_self_plagiarism(target, candidates)

        # 統計集計
        verbatim_total = sum(
            r.verbatim_ratio * 100 for r in section_reports.values()
        )
        paraphrase_total = sum(
            r.paraphrase_ratio * 100 for r in section_reports.values()
        )
        total_sents = sum(
            len(r.matches) for r in section_reports.values()
        )

        # 全体リスク
        if section_reports:
            avg_adjusted = sum(
                r.adjusted_similarity for r in section_reports.values()
            ) / len(section_reports)
        else:
            avg_adjusted = 0.0

        if avg_adjusted > 0.25:
            overall_risk = "critical"
        elif avg_adjusted > 0.15:
            overall_risk = "high"
        elif avg_adjusted > 0.05:
            overall_risk = "medium"
        else:
            overall_risk = "low"

        # トップマッチのソート
        all_matches.sort(key=lambda x: x["similarity"], reverse=True)

        return ComprehensivePlagiarismResult(
            target_doc_id=target.doc_id,
            overall_risk=overall_risk,
            overall_score=avg_adjusted,
            section_reports=section_reports,
            top_matches=all_matches[:10],
            self_plagiarism_flag=self_plagiarism,
            verbatim_count=int(verbatim_total),
            paraphrase_count=int(paraphrase_total),
            total_sentences_checked=total_sents,
            summary=(
                f"Document: {target.doc_id}\n"
                f"Overall risk: {overall_risk} "
                f"(adjusted similarity: {avg_adjusted:.1%})\n"
                f"Sections checked: {len(section_reports)}\n"
                f"Top matches: {len(all_matches)}\n"
                f"Self-plagiarism: {'Yes' if self_plagiarism else 'No'}"
            ),
        )

    def _check_self_plagiarism(self, target: DocumentProfile,
                               candidates: List[Tuple[str, int]]) -> bool:
        """自己盗作の検出"""
        for candidate_id, _ in candidates:
            if candidate_id not in self._corpus:
                continue
            candidate = self._corpus[candidate_id]
            # 著者の重複チェック
            if set(target.authors) & set(candidate.authors):
                return True
        return False
