"""
テキストフィンガープリント — 高速な文書類似度スクリーニング

手法:
1. Winnowing (Schleimer, Wilkerson & Aiken, 2003) — ローカルフィンガープリント
2. SimHash — 近似最近傍探索用のハッシュ
3. MinHash + LSH — 大規模コーパスでの高速検索
"""

import hashlib
from dataclasses import dataclass, field
from typing import List, Set, Dict, Tuple, Optional
from collections import defaultdict


@dataclass
class FingerprintResult:
    """フィンガープリント結果"""
    document_id: str
    fingerprints: Set[int]
    simhash: int
    num_tokens: int
    num_fingerprints: int


class TextFingerprinter:
    """
    テキストフィンガープリント生成器。

    Winnowing アルゴリズムを使用してローカルフィンガープリントを生成し、
    文書間の類似度を高速に判定する。

    Parameters
    ----------
    k : int
        k-gramのサイズ（デフォルト: 5）
    window_size : int
        Winnowingのウィンドウサイズ（デフォルト: 4）
    """

    def __init__(self, k: int = 5, window_size: int = 4):
        self.k = k
        self.window_size = window_size
        self._index: Dict[int, List[str]] = defaultdict(list)

    def fingerprint(self, text: str,
                    document_id: str = "unknown") -> FingerprintResult:
        """テキストのフィンガープリントを生成"""
        tokens = self._tokenize(text)
        kgrams = self._generate_kgrams(tokens)
        hashes = [self._hash_kgram(kg) for kg in kgrams]

        # Winnowingでフィンガープリントを選択
        fingerprints = self._winnow(hashes)

        # SimHash
        simhash = self._compute_simhash(tokens)

        return FingerprintResult(
            document_id=document_id,
            fingerprints=fingerprints,
            simhash=simhash,
            num_tokens=len(tokens),
            num_fingerprints=len(fingerprints),
        )

    def similarity(self, fp_a: FingerprintResult,
                   fp_b: FingerprintResult) -> float:
        """2つのフィンガープリント間のJaccard類似度"""
        if not fp_a.fingerprints or not fp_b.fingerprints:
            return 0.0
        intersection = fp_a.fingerprints & fp_b.fingerprints
        union = fp_a.fingerprints | fp_b.fingerprints
        return len(intersection) / len(union) if union else 0.0

    def simhash_distance(self, fp_a: FingerprintResult,
                         fp_b: FingerprintResult) -> int:
        """SimHashのハミング距離"""
        xor = fp_a.simhash ^ fp_b.simhash
        return bin(xor).count("1")

    def index_document(self, fp: FingerprintResult):
        """フィンガープリントをインデックスに追加"""
        for h in fp.fingerprints:
            self._index[h].append(fp.document_id)

    def query(self, fp: FingerprintResult,
              min_shared: int = 3) -> List[Tuple[str, int]]:
        """インデックスから類似文書を検索"""
        candidates = defaultdict(int)
        for h in fp.fingerprints:
            for doc_id in self._index.get(h, []):
                if doc_id != fp.document_id:
                    candidates[doc_id] += 1

        results = [
            (doc_id, count)
            for doc_id, count in candidates.items()
            if count >= min_shared
        ]
        results.sort(key=lambda x: x[1], reverse=True)
        return results

    def _tokenize(self, text: str) -> List[str]:
        """正規化トークン化"""
        import re
        text = text.lower()
        text = re.sub(r'[^\w\s]', '', text)
        text = re.sub(r'\s+', ' ', text)
        return text.split()

    def _generate_kgrams(self, tokens: List[str]) -> List[Tuple[str, ...]]:
        """k-gramの生成"""
        return [
            tuple(tokens[i:i+self.k])
            for i in range(len(tokens) - self.k + 1)
        ]

    def _hash_kgram(self, kgram: Tuple[str, ...]) -> int:
        """k-gramのハッシュ値"""
        text = " ".join(kgram)
        return int(hashlib.md5(text.encode()).hexdigest()[:8], 16)

    def _winnow(self, hashes: List[int]) -> Set[int]:
        """Winnowingアルゴリズム"""
        if len(hashes) < self.window_size:
            return set(hashes)

        fingerprints = set()
        prev_min_idx = -1

        for i in range(len(hashes) - self.window_size + 1):
            window = hashes[i:i + self.window_size]
            # ウィンドウ内の最小値の最右端を選択
            min_val = min(window)
            min_idx = i + len(window) - 1 - window[::-1].index(min_val)

            if min_idx != prev_min_idx:
                fingerprints.add(min_val)
                prev_min_idx = min_idx

        return fingerprints

    def _compute_simhash(self, tokens: List[str], bits: int = 64) -> int:
        """SimHashの計算"""
        v = [0] * bits

        for token in tokens:
            h = int(hashlib.md5(token.encode()).hexdigest()[:16], 16)
            for i in range(bits):
                if h & (1 << i):
                    v[i] += 1
                else:
                    v[i] -= 1

        simhash = 0
        for i in range(bits):
            if v[i] > 0:
                simhash |= (1 << i)
        return simhash
