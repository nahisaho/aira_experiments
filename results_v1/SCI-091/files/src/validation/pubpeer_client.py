"""
PubPeer APIクライアント

PubPeer.comから論文に対するポストパブリケーションレビューデータを取得し、
研究公正性の外部シグナルとして活用する。
"""

import json
from dataclasses import dataclass, field
from typing import List, Dict, Optional
from datetime import datetime


@dataclass
class PubPeerComment:
    """PubPeerコメント"""
    comment_id: str
    content: str
    date: str
    is_anonymous: bool
    category: str  # concern, question, praise, other
    sentiment: float  # -1.0 (negative) to 1.0 (positive)


@dataclass
class PubPeerRecord:
    """PubPeer論文レコード"""
    doi: str
    title: str
    num_comments: int
    comments: List[PubPeerComment]
    has_integrity_concerns: bool
    concern_categories: List[str]
    first_comment_date: Optional[str] = None
    latest_comment_date: Optional[str] = None


class PubPeerClient:
    """
    PubPeer API クライアント。

    PubPeer の公開APIを使用して論文に対するコメントを取得し、
    研究公正性に関する懸念の有無を評価する。

    API endpoint: https://pubpeer.com/api/v1/

    Parameters
    ----------
    api_base : str
        API のベースURL
    timeout : int
        リクエストタイムアウト（秒）
    """

    CONCERN_KEYWORDS = [
        "manipulation", "duplicate", "fabricat", "falsif",
        "plagiaris", "error", "incorrect", "retract",
        "concern", "questionable", "suspicious", "overlap",
        "inconsistent", "impossible", "copy", "splice",
    ]

    CATEGORY_KEYWORDS = {
        "image_concern": [
            "image", "figure", "western blot", "gel", "microscopy",
            "duplicate", "manipulation", "splice",
        ],
        "statistical_concern": [
            "statistics", "p-value", "mean", "standard deviation",
            "GRIM", "SPRITE", "inconsistent", "impossible",
        ],
        "plagiarism_concern": [
            "plagiarism", "copied", "overlap", "verbatim",
            "text recycling", "self-plagiarism",
        ],
        "data_concern": [
            "data", "fabricat", "falsif", "raw data",
            "made up", "generated",
        ],
        "methodological_concern": [
            "methodology", "protocol", "reproducib",
            "cannot replicate", "failed to reproduce",
        ],
    }

    def __init__(self, api_base: str = "https://pubpeer.com/api/v1",
                 timeout: int = 30):
        self.api_base = api_base
        self.timeout = timeout

    def search_by_doi(self, doi: str) -> Optional[PubPeerRecord]:
        """
        DOIでPubPeerレコードを検索する。

        注: 実際のAPI呼び出しはネットワーク接続が必要。
        ここではデータ構造とロジックを定義する。

        Parameters
        ----------
        doi : str
            論文のDOI

        Returns
        -------
        PubPeerRecord or None
        """
        # API呼び出しのテンプレート
        # 実際にはrequestsで取得
        url = f"{self.api_base}/publications?doi={doi}"
        # response = requests.get(url, timeout=self.timeout)
        # data = response.json()

        # シミュレーションデータ（実運用時はAPI結果を使用）
        return None

    def analyze_comments(self, comments: List[Dict]) -> PubPeerRecord:
        """
        PubPeerコメントを分析し構造化する。

        Parameters
        ----------
        comments : list of dict
            APIから取得したコメントデータ
        """
        parsed_comments = []
        concern_cats = set()

        for c in comments:
            content = c.get("content", "")
            category = self._classify_comment(content)
            sentiment = self._estimate_sentiment(content)
            is_concern = any(
                kw in content.lower() for kw in self.CONCERN_KEYWORDS
            )

            if is_concern:
                for cat, keywords in self.CATEGORY_KEYWORDS.items():
                    if any(kw in content.lower() for kw in keywords):
                        concern_cats.add(cat)

            parsed_comments.append(PubPeerComment(
                comment_id=c.get("id", ""),
                content=content,
                date=c.get("created_at", ""),
                is_anonymous=c.get("is_anonymous", True),
                category=category,
                sentiment=sentiment,
            ))

        has_concerns = len(concern_cats) > 0

        dates = [c.date for c in parsed_comments if c.date]
        dates.sort()

        return PubPeerRecord(
            doi=c.get("doi", "") if comments else "",
            title=c.get("title", "") if comments else "",
            num_comments=len(parsed_comments),
            comments=parsed_comments,
            has_integrity_concerns=has_concerns,
            concern_categories=list(concern_cats),
            first_comment_date=dates[0] if dates else None,
            latest_comment_date=dates[-1] if dates else None,
        )

    def _classify_comment(self, text: str) -> str:
        """コメントをカテゴリに分類"""
        text_lower = text.lower()
        for cat, keywords in self.CATEGORY_KEYWORDS.items():
            if any(kw in text_lower for kw in keywords):
                return cat
        if any(kw in text_lower for kw in self.CONCERN_KEYWORDS):
            return "general_concern"
        return "other"

    def _estimate_sentiment(self, text: str) -> float:
        """簡易センチメント推定"""
        negative_words = [
            "concern", "error", "incorrect", "suspicious",
            "questionable", "fraud", "fabricat", "manipulat",
        ]
        positive_words = [
            "excellent", "impressive", "well-done", "robust",
            "thorough", "confirm",
        ]

        text_lower = text.lower()
        neg = sum(1 for w in negative_words if w in text_lower)
        pos = sum(1 for w in positive_words if w in text_lower)

        if neg + pos == 0:
            return 0.0
        return (pos - neg) / (pos + neg)
