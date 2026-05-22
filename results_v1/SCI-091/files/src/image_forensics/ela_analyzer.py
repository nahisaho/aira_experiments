"""
Error Level Analysis (ELA) — JPEG再圧縮差分による画像加工領域検出。

原理: JPEG画像を既知の品質で再圧縮し、元画像との差分を計算。
加工された領域は周囲と異なるエラーレベルを示す。
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional
from io import BytesIO

try:
    from PIL import Image
except ImportError:
    Image = None


@dataclass
class ELAResult:
    """ELA解析結果"""
    ela_map: np.ndarray
    mean_error: float
    max_error: float
    suspicious_ratio: float  # 閾値超過ピクセルの割合
    suspicious_regions: list = field(default_factory=list)
    quality_used: int = 95


class ELAAnalyzer:
    """
    Error Level Analysis を実行し、画像加工の痕跡を検出する。

    Parameters
    ----------
    quality : int
        再圧縮のJPEG品質 (デフォルト: 95)
    threshold_factor : float
        平均エラーに対する閾値倍率。この倍率を超えるピクセルを疑わしいとする。
    """

    def __init__(self, quality: int = 95, threshold_factor: float = 3.0):
        self.quality = quality
        self.threshold_factor = threshold_factor

    def analyze(self, image_array: np.ndarray) -> ELAResult:
        """
        ELA解析を実行する。

        Parameters
        ----------
        image_array : np.ndarray
            入力画像 (H, W, 3) uint8

        Returns
        -------
        ELAResult
            解析結果
        """
        if Image is None:
            return self._analyze_without_pil(image_array)

        img = Image.fromarray(image_array.astype(np.uint8))

        # JPEG再圧縮
        buffer = BytesIO()
        img.save(buffer, format="JPEG", quality=self.quality)
        buffer.seek(0)
        recompressed = np.array(Image.open(buffer)).astype(np.float32)

        original = image_array.astype(np.float32)

        # 差分計算（絶対値）
        ela_map = np.abs(original - recompressed)

        # チャンネル平均
        if ela_map.ndim == 3:
            ela_single = np.mean(ela_map, axis=2)
        else:
            ela_single = ela_map

        mean_error = float(np.mean(ela_single))
        max_error = float(np.max(ela_single))

        # 疑わしい領域の検出
        threshold = mean_error * self.threshold_factor
        suspicious_mask = ela_single > threshold
        suspicious_ratio = float(np.mean(suspicious_mask))

        # 連結成分で疑わしい領域を特定
        regions = self._find_suspicious_regions(suspicious_mask)

        return ELAResult(
            ela_map=ela_map,
            mean_error=mean_error,
            max_error=max_error,
            suspicious_ratio=suspicious_ratio,
            suspicious_regions=regions,
            quality_used=self.quality,
        )

    def _analyze_without_pil(self, image_array: np.ndarray) -> ELAResult:
        """PIL不在時のフォールバック: DCT近似によるELA"""
        h, w = image_array.shape[:2]

        # ブロック単位のDCT近似（8x8ブロック）
        block_size = 8
        ela_map = np.zeros_like(image_array, dtype=np.float32)

        for i in range(0, h - block_size + 1, block_size):
            for j in range(0, w - block_size + 1, block_size):
                block = image_array[i:i+block_size, j:j+block_size].astype(np.float32)
                # 量子化シミュレーション
                quantized = np.round(block / 8.0) * 8.0
                ela_map[i:i+block_size, j:j+block_size] = np.abs(block - quantized)

        if ela_map.ndim == 3:
            ela_single = np.mean(ela_map, axis=2)
        else:
            ela_single = ela_map

        mean_error = float(np.mean(ela_single))
        max_error = float(np.max(ela_single))
        threshold = mean_error * self.threshold_factor
        suspicious_mask = ela_single > threshold
        suspicious_ratio = float(np.mean(suspicious_mask))
        regions = self._find_suspicious_regions(suspicious_mask)

        return ELAResult(
            ela_map=ela_map,
            mean_error=mean_error,
            max_error=max_error,
            suspicious_ratio=suspicious_ratio,
            suspicious_regions=regions,
            quality_used=self.quality,
        )

    def _find_suspicious_regions(
        self, mask: np.ndarray, min_area: int = 100
    ) -> list:
        """連結成分ラベリングで疑わしい領域を抽出"""
        try:
            from scipy import ndimage
            labeled, num_features = ndimage.label(mask)
            regions = []
            for i in range(1, num_features + 1):
                region_mask = labeled == i
                area = int(np.sum(region_mask))
                if area >= min_area:
                    ys, xs = np.where(region_mask)
                    regions.append({
                        "bbox": [int(xs.min()), int(ys.min()),
                                 int(xs.max()), int(ys.max())],
                        "area": area,
                        "centroid": [int(np.mean(xs)), int(np.mean(ys))],
                    })
            return regions
        except ImportError:
            return []
