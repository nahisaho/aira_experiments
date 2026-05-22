"""
GRIM Test (Granularity-Related Inconsistency of Means)

Brown & Heathers (2017) によって提案された手法。
報告された平均値がサンプルサイズと整合するかを検証する。

原理: N人のサンプルの平均は 1/N の精度でしか取りえない。
例: N=25 の場合、平均値は 0.04 の倍数でなければならない。
報告された平均値がこの制約を満たさない場合、不整合と判定。
"""

import math
from dataclasses import dataclass
from typing import List, Optional, Tuple


@dataclass
class GRIMResult:
    """GRIM Test の結果"""
    reported_mean: float
    sample_size: int
    decimal_places: int
    is_consistent: bool
    closest_possible_mean: float
    granularity: float
    deviation: float
    details: str


class GRIMTest:
    """
    GRIM Test の実装。

    整数スケール（例: Likert尺度 1-7）の平均値に対する
    粒度整合性を検証する。

    Parameters
    ----------
    tolerance : float
        丸め誤差の許容範囲（デフォルト: 1e-6）
    scale_items : int
        尺度の項目数（複数項目の合計の場合）
    """

    def __init__(self, tolerance: float = 1e-6, scale_items: int = 1):
        self.tolerance = tolerance
        self.scale_items = scale_items

    def test(self, mean: float, n: int,
             decimal_places: Optional[int] = None) -> GRIMResult:
        """
        GRIM Testを実行する。

        Parameters
        ----------
        mean : float
            報告された平均値
        n : int
            サンプルサイズ
        decimal_places : int, optional
            報告された小数桁数。未指定の場合は自動検出。

        Returns
        -------
        GRIMResult
        """
        if decimal_places is None:
            decimal_places = self._detect_decimal_places(mean)

        # 粒度: 1/(N * items)
        granularity = 1.0 / (n * self.scale_items)

        # 報告された平均値から推定される合計値
        total = mean * n * self.scale_items
        rounded_total = round(total)
        reconstructed_mean = rounded_total / (n * self.scale_items)

        # 小数桁数に合わせて丸め
        reconstructed_rounded = round(reconstructed_mean, decimal_places)
        reported_rounded = round(mean, decimal_places)

        deviation = abs(reconstructed_rounded - reported_rounded)
        is_consistent = deviation <= self.tolerance

        # 最も近い可能な平均値を探索
        closest = self._find_closest_possible(mean, n, decimal_places)

        return GRIMResult(
            reported_mean=mean,
            sample_size=n,
            decimal_places=decimal_places,
            is_consistent=is_consistent,
            closest_possible_mean=closest,
            granularity=granularity,
            deviation=deviation,
            details=self._generate_details(
                mean, n, decimal_places, is_consistent,
                reconstructed_mean, closest
            ),
        )

    def batch_test(self, means_and_ns: List[Tuple[float, int]],
                   decimal_places: Optional[int] = None) -> List[GRIMResult]:
        """複数の平均値を一括テスト"""
        return [
            self.test(mean, n, decimal_places)
            for mean, n in means_and_ns
        ]

    def inconsistency_rate(self, results: List[GRIMResult]) -> float:
        """不整合率を計算"""
        if not results:
            return 0.0
        inconsistent = sum(1 for r in results if not r.is_consistent)
        return inconsistent / len(results)

    def _find_closest_possible(self, mean: float, n: int,
                               decimal_places: int) -> float:
        """最も近い可能な平均値を探索"""
        total = mean * n * self.scale_items
        lower = math.floor(total)
        upper = math.ceil(total)

        candidates = []
        for t in [lower, upper]:
            possible_mean = t / (n * self.scale_items)
            rounded = round(possible_mean, decimal_places)
            candidates.append((abs(rounded - round(mean, decimal_places)), rounded))

        candidates.sort()
        return candidates[0][1]

    def _detect_decimal_places(self, value: float) -> int:
        """小数桁数を自動検出"""
        s = f"{value:.10f}".rstrip("0")
        if "." not in s:
            return 0
        return len(s.split(".")[1])

    def _generate_details(self, mean, n, dp, consistent,
                          reconstructed, closest) -> str:
        """詳細レポートを生成"""
        status = "CONSISTENT ✓" if consistent else "INCONSISTENT ✗"
        return (
            f"GRIM Test Result: {status}\n"
            f"  Reported mean: {mean}\n"
            f"  Sample size (N): {n}\n"
            f"  Decimal places: {dp}\n"
            f"  Reconstructed mean: {reconstructed:.{dp}f}\n"
            f"  Closest possible mean: {closest:.{dp}f}\n"
            f"  Granularity: {1.0/(n*self.scale_items):.6f}"
        )
