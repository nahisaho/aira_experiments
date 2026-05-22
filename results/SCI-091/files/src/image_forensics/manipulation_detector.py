"""
画像加工検出モジュール

検出手法:
1. ELA (Error Level Analysis)
2. ノイズ不整合分析
3. EXIF メタデータ整合性チェック
4. JPEG ゴースト検出
5. CNN ベースの加工検出 (ManTraNet 風アーキテクチャ)
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple
from .ela_analyzer import ELAAnalyzer, ELAResult


@dataclass
class ManipulationResult:
    """加工検出の総合結果"""
    is_manipulated: bool
    confidence: float  # 0.0 - 1.0
    manipulation_type: str  # splicing, copy-move, retouching, generation
    ela_result: Optional[ELAResult] = None
    noise_inconsistency: float = 0.0
    jpeg_ghost_score: float = 0.0
    metadata_flags: list = field(default_factory=list)
    region_mask: Optional[np.ndarray] = None
    details: dict = field(default_factory=dict)


class NoiseAnalyzer:
    """
    ノイズレベル不整合分析。

    原理: 自然画像は一様なノイズ特性を持つ。加工された領域は
    異なるノイズレベルを示すため、局所的なノイズ推定の不整合から
    加工を検出できる。
    """

    def __init__(self, block_size: int = 32):
        self.block_size = block_size

    def estimate_local_noise(self, image: np.ndarray) -> np.ndarray:
        """局所ノイズレベルマップを推定"""
        if image.ndim == 3:
            gray = np.mean(image, axis=2)
        else:
            gray = image.astype(np.float64)

        h, w = gray.shape
        noise_map = np.zeros(
            (h // self.block_size, w // self.block_size)
        )

        for i in range(noise_map.shape[0]):
            for j in range(noise_map.shape[1]):
                y_start = i * self.block_size
                x_start = j * self.block_size
                block = gray[
                    y_start:y_start + self.block_size,
                    x_start:x_start + self.block_size
                ]
                # Median Absolute Deviation ベースのノイズ推定
                # (Donoho & Johnstone, 1994)
                # 高周波成分からノイズレベルを推定
                hf = self._high_frequency(block)
                noise_map[i, j] = np.median(np.abs(hf)) / 0.6745

        return noise_map

    def detect_inconsistency(self, image: np.ndarray) -> Tuple[float, np.ndarray]:
        """ノイズ不整合を検出しスコアとマスクを返す"""
        noise_map = self.estimate_local_noise(image)
        global_noise = np.median(noise_map)

        if global_noise == 0:
            return 0.0, np.zeros_like(noise_map, dtype=bool)

        # ノイズレベルの局所偏差
        deviation = np.abs(noise_map - global_noise) / (global_noise + 1e-8)
        inconsistency_score = float(np.mean(deviation > 0.5))

        suspicious_mask = deviation > 0.5
        return inconsistency_score, suspicious_mask

    def _high_frequency(self, block: np.ndarray) -> np.ndarray:
        """ラプラシアンフィルタで高周波成分を抽出"""
        kernel = np.array([[0, 1, 0], [1, -4, 1], [0, 1, 0]], dtype=np.float64)
        h, w = block.shape
        result = np.zeros_like(block)
        for i in range(1, h - 1):
            for j in range(1, w - 1):
                result[i, j] = np.sum(
                    block[i-1:i+2, j-1:j+2] * kernel
                )
        return result


class JPEGGhostDetector:
    """
    JPEGゴースト検出。

    原理: 二重圧縮された画像は特定の品質レベルで再圧縮した際に
    最小のエラーを示す。加工された領域はこのパターンが崩れる。
    """

    def __init__(self, quality_range: Tuple[int, int] = (50, 99), step: int = 5):
        self.quality_range = quality_range
        self.step = step

    def detect(self, image: np.ndarray) -> Tuple[float, Dict]:
        """JPEGゴーストを検出"""
        try:
            from PIL import Image as PILImage
            from io import BytesIO
        except ImportError:
            return 0.0, {"error": "PIL not available"}

        img = PILImage.fromarray(image.astype(np.uint8))
        errors = {}

        for q in range(self.quality_range[0],
                       self.quality_range[1] + 1, self.step):
            buf = BytesIO()
            img.save(buf, format="JPEG", quality=q)
            buf.seek(0)
            recomp = np.array(PILImage.open(buf)).astype(np.float32)
            error = np.mean((image.astype(np.float32) - recomp) ** 2)
            errors[q] = float(error)

        # 最小エラーの品質レベルが元の圧縮品質を示唆
        if errors:
            min_q = min(errors, key=errors.get)
            min_error = errors[min_q]
            mean_error = np.mean(list(errors.values()))
            ghost_score = 1.0 - (min_error / (mean_error + 1e-8))
            ghost_score = max(0.0, min(1.0, ghost_score))
        else:
            ghost_score = 0.0
            min_q = -1

        return ghost_score, {
            "estimated_original_quality": min_q,
            "quality_error_curve": errors,
            "ghost_score": ghost_score,
        }


class ManipulationDetector:
    """
    画像加工検出の統合クラス。

    複数の検出手法を統合し、総合的な加工判定を行う。

    重み付けスコア:
    - ELA疑わしさ: 0.30
    - ノイズ不整合: 0.25
    - JPEGゴースト: 0.20
    - メタデータ異常: 0.15
    - CNNスコア: 0.10 (モデルがある場合)
    """

    WEIGHTS = {
        "ela": 0.30,
        "noise": 0.25,
        "jpeg_ghost": 0.20,
        "metadata": 0.15,
        "cnn": 0.10,
    }

    THRESHOLD = 0.45  # この閾値を超えると加工と判定

    def __init__(self):
        self.ela_analyzer = ELAAnalyzer(quality=95, threshold_factor=3.0)
        self.noise_analyzer = NoiseAnalyzer(block_size=32)
        self.jpeg_ghost = JPEGGhostDetector()

    def detect(self, image: np.ndarray,
               metadata: Optional[Dict] = None) -> ManipulationResult:
        """
        総合的な画像加工検出を実行。

        Parameters
        ----------
        image : np.ndarray
            入力画像 (H, W, 3) uint8
        metadata : dict, optional
            EXIF等のメタデータ

        Returns
        -------
        ManipulationResult
        """
        scores = {}

        # 1. ELA解析
        ela_result = self.ela_analyzer.analyze(image)
        scores["ela"] = min(ela_result.suspicious_ratio * 5, 1.0)

        # 2. ノイズ不整合
        noise_score, noise_mask = self.noise_analyzer.detect_inconsistency(image)
        scores["noise"] = min(noise_score * 3, 1.0)

        # 3. JPEGゴースト
        ghost_score, ghost_details = self.jpeg_ghost.detect(image)
        scores["jpeg_ghost"] = ghost_score

        # 4. メタデータチェック
        metadata_flags = []
        metadata_score = 0.0
        if metadata:
            metadata_flags, metadata_score = self._check_metadata(metadata)
        scores["metadata"] = metadata_score

        # 5. CNN (プレースホルダー)
        scores["cnn"] = 0.0

        # 重み付き総合スコア
        total = sum(
            scores[k] * self.WEIGHTS[k] for k in self.WEIGHTS
        )

        # 加工タイプの推定
        manipulation_type = self._infer_type(scores, ela_result)

        return ManipulationResult(
            is_manipulated=total > self.THRESHOLD,
            confidence=min(total, 1.0),
            manipulation_type=manipulation_type,
            ela_result=ela_result,
            noise_inconsistency=scores["noise"],
            jpeg_ghost_score=scores["jpeg_ghost"],
            metadata_flags=metadata_flags,
            details={
                "individual_scores": scores,
                "weighted_total": total,
                "threshold": self.THRESHOLD,
                "ghost_details": ghost_details,
            },
        )

    def _check_metadata(self, metadata: Dict) -> Tuple[List[str], float]:
        """EXIFメタデータの整合性チェック"""
        flags = []
        score = 0.0

        # ソフトウェア編集痕跡
        software = metadata.get("Software", "")
        editing_tools = [
            "photoshop", "gimp", "paint", "affinity",
            "pixelmator", "illustrator",
        ]
        if any(tool in software.lower() for tool in editing_tools):
            flags.append(f"Editing software detected: {software}")
            score += 0.3

        # 日付整合性
        date_original = metadata.get("DateTimeOriginal")
        date_modified = metadata.get("DateTime")
        if date_original and date_modified and date_modified < date_original:
            flags.append("Modified date precedes original date")
            score += 0.4

        # サムネイルと本体の不一致
        if metadata.get("thumbnail_mismatch"):
            flags.append("Thumbnail does not match main image")
            score += 0.5

        # GPS情報の不整合
        if metadata.get("gps_inconsistency"):
            flags.append("GPS data inconsistent with other metadata")
            score += 0.2

        return flags, min(score, 1.0)

    def _infer_type(self, scores: Dict, ela_result: ELAResult) -> str:
        """検出パターンから加工タイプを推定"""
        if scores["ela"] > 0.5 and ela_result.suspicious_ratio > 0.05:
            if len(ela_result.suspicious_regions) > 1:
                return "splicing"
            return "retouching"
        if scores["noise"] > 0.5:
            return "splicing"
        if scores["jpeg_ghost"] > 0.5:
            return "double_compression"
        if scores["metadata"] > 0.3:
            return "metadata_tampering"
        return "unknown"
