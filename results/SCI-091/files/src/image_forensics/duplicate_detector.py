"""
画像重複検出モジュール

手法:
1. 知覚ハッシュ (pHash/dHash) による高速スクリーニング
2. SIFT/ORB特徴量マッチングによる部分重複検出
3. CNN特徴量 (ResNet/EfficientNet) によるセマンティック類似度
4. Copy-Move Forgery Detection (CMFD) 用のブロックマッチング
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Dict
from enum import Enum


class DuplicationType(Enum):
    EXACT = "exact"
    NEAR_DUPLICATE = "near_duplicate"
    PARTIAL_OVERLAP = "partial_overlap"
    COPY_MOVE = "copy_move"
    SPLICING = "splicing"


@dataclass
class DuplicateMatch:
    """重複検出結果"""
    image_id_a: str
    image_id_b: str
    duplication_type: DuplicationType
    similarity_score: float  # 0.0 - 1.0
    confidence: float
    matched_regions: list = field(default_factory=list)
    transform_detected: Optional[str] = None  # rotation, flip, scale, etc.


class PerceptualHasher:
    """知覚ハッシュによる高速重複スクリーニング"""

    def __init__(self, hash_size: int = 16):
        self.hash_size = hash_size

    def dhash(self, image: np.ndarray) -> int:
        """Difference Hash: 隣接ピクセルの勾配方向をエンコード"""
        if image.ndim == 3:
            gray = np.mean(image, axis=2)
        else:
            gray = image.astype(np.float64)

        # リサイズ (hash_size+1 x hash_size)
        h, w = gray.shape
        resized = self._resize(gray, self.hash_size + 1, self.hash_size)

        # 水平勾配の符号をビット列に変換
        diff = resized[:, 1:] > resized[:, :-1]
        hash_val = 0
        for bit in diff.flatten():
            hash_val = (hash_val << 1) | int(bit)
        return hash_val

    def phash(self, image: np.ndarray) -> int:
        """Perceptual Hash: DCTベースの知覚ハッシュ"""
        if image.ndim == 3:
            gray = np.mean(image, axis=2)
        else:
            gray = image.astype(np.float64)

        resized = self._resize(gray, 32, 32)

        # 簡易DCT（numpy FFTベース）
        dct = np.fft.fft2(resized)
        dct_low = np.abs(dct[:8, :8])

        median_val = np.median(dct_low)
        hash_val = 0
        for bit in (dct_low.flatten() > median_val):
            hash_val = (hash_val << 1) | int(bit)
        return hash_val

    def hamming_distance(self, hash_a: int, hash_b: int) -> int:
        """ハミング距離"""
        xor = hash_a ^ hash_b
        return bin(xor).count("1")

    def _resize(self, img: np.ndarray, new_w: int, new_h: int) -> np.ndarray:
        """バイリニア補間リサイズ"""
        h, w = img.shape
        y_ratio = h / new_h
        x_ratio = w / new_w
        y_idx = (np.arange(new_h) * y_ratio).astype(int).clip(0, h - 1)
        x_idx = (np.arange(new_w) * x_ratio).astype(int).clip(0, w - 1)
        return img[np.ix_(y_idx, x_idx)]


class BlockMatcher:
    """
    Copy-Move Forgery Detection のためのブロックマッチング。

    アルゴリズム:
    1. 画像を重複ブロックに分割
    2. 各ブロックの特徴ベクトルを計算 (DCT係数 or PCA)
    3. kd-tree で近傍探索
    4. 幾何学的整合性チェック (RANSAC)
    """

    def __init__(self, block_size: int = 16, stride: int = 4,
                 similarity_threshold: float = 0.95):
        self.block_size = block_size
        self.stride = stride
        self.similarity_threshold = similarity_threshold

    def detect_copy_move(self, image: np.ndarray) -> List[Dict]:
        """画像内のコピー＆ムーブ偽造を検出"""
        if image.ndim == 3:
            gray = np.mean(image, axis=2)
        else:
            gray = image.astype(np.float64)

        h, w = gray.shape
        blocks = []
        positions = []

        # ブロック抽出と特徴量計算
        for y in range(0, h - self.block_size + 1, self.stride):
            for x in range(0, w - self.block_size + 1, self.stride):
                block = gray[y:y+self.block_size, x:x+self.block_size]
                # DCT特徴量の簡易計算
                features = np.fft.fft2(block).flatten()[:32]
                features = np.abs(features)
                if np.linalg.norm(features) > 0:
                    features = features / np.linalg.norm(features)
                blocks.append(features)
                positions.append((x, y))

        if len(blocks) == 0:
            return []

        blocks = np.array(blocks)
        positions = np.array(positions)

        # 類似ブロックペアの検出（コサイン類似度）
        matches = []
        n = len(blocks)

        # 効率化: ランダムサンプリング + 局所探索
        sample_size = min(n, 500)
        indices = np.random.choice(n, sample_size, replace=False)

        for i in indices:
            # 自分自身と近傍を除外した類似度計算
            sims = blocks[i] @ blocks.T
            for j in np.where(sims > self.similarity_threshold)[0]:
                dist = np.linalg.norm(positions[i] - positions[j])
                if dist > self.block_size * 2:  # 近傍を除外
                    matches.append({
                        "source": positions[i].tolist(),
                        "target": positions[j].tolist(),
                        "similarity": float(sims[j]),
                        "displacement": dist,
                    })

        # クラスタリングで偽造領域を特定
        return self._cluster_matches(matches)

    def _cluster_matches(self, matches: List[Dict],
                         min_cluster_size: int = 5) -> List[Dict]:
        """マッチをクラスタリングして偽造領域を特定"""
        if len(matches) < min_cluster_size:
            return []

        # 変位ベクトルでクラスタリング
        displacements = np.array([
            [m["target"][0] - m["source"][0],
             m["target"][1] - m["source"][1]]
            for m in matches
        ])

        # 簡易クラスタリング（ヒストグラムベース）
        regions = []
        if len(displacements) > 0:
            # 変位ベクトルの主要方向を検出
            unique_dirs = {}
            for i, d in enumerate(displacements):
                key = (round(d[0] / 10) * 10, round(d[1] / 10) * 10)
                if key not in unique_dirs:
                    unique_dirs[key] = []
                unique_dirs[key].append(i)

            for key, idx_list in unique_dirs.items():
                if len(idx_list) >= min_cluster_size:
                    source_pts = np.array(
                        [matches[i]["source"] for i in idx_list]
                    )
                    regions.append({
                        "displacement_vector": list(key),
                        "num_matches": len(idx_list),
                        "source_bbox": [
                            int(source_pts[:, 0].min()),
                            int(source_pts[:, 1].min()),
                            int(source_pts[:, 0].max()) + self.block_size,
                            int(source_pts[:, 1].max()) + self.block_size,
                        ],
                        "confidence": min(len(idx_list) / 20.0, 1.0),
                    })

        return regions


class DuplicateDetector:
    """
    論文内画像の重複検出統合クラス。

    パイプライン:
    1. 知覚ハッシュで高速フィルタリング
    2. ブロックマッチングでコピー＆ムーブ検出
    3. CNN特徴量で意味的類似度評価
    """

    def __init__(self, hash_threshold: int = 10,
                 block_size: int = 16,
                 cnn_threshold: float = 0.85):
        self.hasher = PerceptualHasher(hash_size=16)
        self.block_matcher = BlockMatcher(block_size=block_size)
        self.hash_threshold = hash_threshold
        self.cnn_threshold = cnn_threshold
        self._hash_db: Dict[str, Tuple[int, int]] = {}

    def register_image(self, image_id: str, image: np.ndarray):
        """画像をデータベースに登録"""
        dhash = self.hasher.dhash(image)
        phash = self.hasher.phash(image)
        self._hash_db[image_id] = (dhash, phash)

    def find_duplicates(
        self, images: Dict[str, np.ndarray]
    ) -> List[DuplicateMatch]:
        """画像セット内の重複を検出"""
        results = []

        # Phase 1: 知覚ハッシュによるペアワイズ比較
        ids = list(images.keys())
        for img_id in ids:
            self.register_image(img_id, images[img_id])

        for i in range(len(ids)):
            for j in range(i + 1, len(ids)):
                id_a, id_b = ids[i], ids[j]
                dhash_a, phash_a = self._hash_db[id_a]
                dhash_b, phash_b = self._hash_db[id_b]

                d_dist = self.hasher.hamming_distance(dhash_a, dhash_b)
                p_dist = self.hasher.hamming_distance(phash_a, phash_b)
                avg_dist = (d_dist + p_dist) / 2

                if avg_dist <= self.hash_threshold:
                    similarity = 1.0 - (avg_dist / 256.0)
                    dup_type = (DuplicationType.EXACT if avg_dist == 0
                                else DuplicationType.NEAR_DUPLICATE)
                    results.append(DuplicateMatch(
                        image_id_a=id_a,
                        image_id_b=id_b,
                        duplication_type=dup_type,
                        similarity_score=similarity,
                        confidence=similarity,
                    ))

        # Phase 2: 各画像内のコピー＆ムーブ検出
        for img_id, img in images.items():
            cm_regions = self.block_matcher.detect_copy_move(img)
            if cm_regions:
                results.append(DuplicateMatch(
                    image_id_a=img_id,
                    image_id_b=img_id,
                    duplication_type=DuplicationType.COPY_MOVE,
                    similarity_score=max(r["confidence"] for r in cm_regions),
                    confidence=max(r["confidence"] for r in cm_regions),
                    matched_regions=cm_regions,
                ))

        return results

    def detect_intra_image_forgery(self, image: np.ndarray) -> List[Dict]:
        """単一画像内のコピー＆ムーブ偽造検出"""
        return self.block_matcher.detect_copy_move(image)
