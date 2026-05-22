"""
SPRITE Test (Sample Parameter Reconstruction via Iterative TEchniques)

Heathers et al. (2018) による手法。
報告された記述統計量（平均、標準偏差、範囲）が
互いに整合するかを検証する。

原理: 与えられた制約（N, 範囲, 平均, SD）を満たすデータセットを
反復的に構築し、そのようなデータセットが存在し得るかを確認する。
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional, Tuple


@dataclass
class SPRITEResult:
    """SPRITE Test の結果"""
    reported_mean: float
    reported_sd: float
    sample_size: int
    scale_range: Tuple[int, int]
    is_consistent: bool
    best_achieved_mean: float
    best_achieved_sd: float
    mean_deviation: float
    sd_deviation: float
    num_solutions_found: int
    iterations_used: int
    example_distribution: Optional[np.ndarray] = None
    details: str = ""


class SPRITETest:
    """
    SPRITE Test の実装。

    反復的なデータ再構成により、報告された記述統計量の整合性を検証する。

    Parameters
    ----------
    max_iterations : int
        最大反復回数
    num_attempts : int
        独立な試行回数
    tolerance_mean : float
        平均値の許容誤差
    tolerance_sd : float
        標準偏差の許容誤差
    """

    def __init__(self, max_iterations: int = 10000,
                 num_attempts: int = 100,
                 tolerance_mean: float = 0.01,
                 tolerance_sd: float = 0.01):
        self.max_iterations = max_iterations
        self.num_attempts = num_attempts
        self.tolerance_mean = tolerance_mean
        self.tolerance_sd = tolerance_sd

    def test(self, mean: float, sd: float, n: int,
             scale_min: int = 1, scale_max: int = 7,
             seed: Optional[int] = None) -> SPRITEResult:
        """
        SPRITE Testを実行する。

        Parameters
        ----------
        mean : float
            報告された平均値
        sd : float
            報告された標準偏差
        n : int
            サンプルサイズ
        scale_min, scale_max : int
            尺度の最小値・最大値
        seed : int, optional
            乱数シード

        Returns
        -------
        SPRITEResult
        """
        rng = np.random.RandomState(seed)

        # 基本的な整合性チェック
        if mean < scale_min or mean > scale_max:
            return SPRITEResult(
                reported_mean=mean, reported_sd=sd, sample_size=n,
                scale_range=(scale_min, scale_max),
                is_consistent=False,
                best_achieved_mean=mean, best_achieved_sd=0.0,
                mean_deviation=abs(mean - (scale_min + scale_max) / 2),
                sd_deviation=sd,
                num_solutions_found=0, iterations_used=0,
                details="Mean outside scale range — trivially inconsistent",
            )

        # 理論的SD上限チェック
        max_possible_sd = self._max_possible_sd(mean, n, scale_min, scale_max)
        if sd > max_possible_sd + self.tolerance_sd:
            return SPRITEResult(
                reported_mean=mean, reported_sd=sd, sample_size=n,
                scale_range=(scale_min, scale_max),
                is_consistent=False,
                best_achieved_mean=mean, best_achieved_sd=max_possible_sd,
                mean_deviation=0.0,
                sd_deviation=sd - max_possible_sd,
                num_solutions_found=0, iterations_used=0,
                details=(
                    f"SD ({sd:.2f}) exceeds theoretical maximum "
                    f"({max_possible_sd:.2f}) for given mean and range"
                ),
            )

        # SPRITE反復
        best_solution = None
        best_distance = float("inf")
        solutions_found = 0
        total_iterations = 0

        for attempt in range(self.num_attempts):
            data, achieved_mean, achieved_sd, iters = self._sprite_iterate(
                mean, sd, n, scale_min, scale_max, rng
            )
            total_iterations += iters

            distance = (
                abs(achieved_mean - mean) / (self.tolerance_mean + 1e-10)
                + abs(achieved_sd - sd) / (self.tolerance_sd + 1e-10)
            )

            if (abs(achieved_mean - mean) <= self.tolerance_mean
                    and abs(achieved_sd - sd) <= self.tolerance_sd):
                solutions_found += 1
                if distance < best_distance:
                    best_distance = distance
                    best_solution = data.copy()

            if distance < best_distance:
                best_distance = distance
                best_solution = data.copy()

        if best_solution is not None:
            best_mean = float(np.mean(best_solution))
            best_sd = float(np.std(best_solution, ddof=1))
        else:
            best_mean = mean
            best_sd = 0.0

        is_consistent = solutions_found > 0

        return SPRITEResult(
            reported_mean=mean,
            reported_sd=sd,
            sample_size=n,
            scale_range=(scale_min, scale_max),
            is_consistent=is_consistent,
            best_achieved_mean=best_mean,
            best_achieved_sd=best_sd,
            mean_deviation=abs(best_mean - mean),
            sd_deviation=abs(best_sd - sd),
            num_solutions_found=solutions_found,
            iterations_used=total_iterations,
            example_distribution=best_solution,
            details=self._generate_details(
                mean, sd, n, scale_min, scale_max,
                is_consistent, solutions_found, best_mean, best_sd
            ),
        )

    def _sprite_iterate(self, target_mean, target_sd, n,
                        scale_min, scale_max, rng):
        """SPRITE反復アルゴリズム"""
        # 初期データ: ターゲット平均に近いデータを生成
        target_sum = round(target_mean * n)
        data = np.full(n, scale_min, dtype=np.float64)

        # 合計値をターゲットに合わせる
        remaining = target_sum - int(np.sum(data))
        for i in range(n):
            add = min(remaining, scale_max - scale_min)
            add = max(add, 0)
            data[i] += add
            remaining -= add
            if remaining <= 0:
                break

        # 反復的にSDを調整
        for iteration in range(self.max_iterations):
            current_mean = np.mean(data)
            current_sd = np.std(data, ddof=1) if n > 1 else 0.0

            if (abs(current_mean - target_mean) <= self.tolerance_mean
                    and abs(current_sd - target_sd) <= self.tolerance_sd):
                return data, float(current_mean), float(current_sd), iteration

            # SDが低すぎる場合: 分散を増やす
            if current_sd < target_sd - self.tolerance_sd:
                idx = rng.randint(0, n)
                if data[idx] < scale_max and data[idx] > scale_min:
                    # 中央値に近い値を極端な値に変更
                    if rng.random() < 0.5 and data[idx] < scale_max:
                        data[idx] += 1
                        # 合計を保持するため別の値を減らす
                        other = rng.randint(0, n)
                        while other == idx or data[other] <= scale_min:
                            other = rng.randint(0, n)
                            if data[other] > scale_min:
                                break
                        if data[other] > scale_min:
                            data[other] -= 1
                    elif data[idx] > scale_min:
                        data[idx] -= 1
                        other = rng.randint(0, n)
                        while other == idx or data[other] >= scale_max:
                            other = rng.randint(0, n)
                            if data[other] < scale_max:
                                break
                        if data[other] < scale_max:
                            data[other] += 1

            # SDが高すぎる場合: 分散を減らす
            elif current_sd > target_sd + self.tolerance_sd:
                idx_high = np.argmax(data)
                idx_low = np.argmin(data)
                if data[idx_high] > scale_min + 1:
                    data[idx_high] -= 1
                    data[idx_low] += 1

            # ランダム摂動
            else:
                i, j = rng.choice(n, 2, replace=False)
                if data[i] < scale_max and data[j] > scale_min:
                    data[i] += 1
                    data[j] -= 1

        final_mean = float(np.mean(data))
        final_sd = float(np.std(data, ddof=1)) if n > 1 else 0.0
        return data, final_mean, final_sd, self.max_iterations

    def _max_possible_sd(self, mean, n, scale_min, scale_max):
        """与えられた制約下での理論的最大SD"""
        target_sum = round(mean * n)

        # 極端な分布: できるだけ多くの値を端に配置
        data = np.full(n, scale_min, dtype=np.float64)
        remaining = target_sum - int(np.sum(data))

        # 上限値をできるだけ多く配置
        for i in range(n):
            if remaining >= (scale_max - scale_min):
                data[i] = scale_max
                remaining -= (scale_max - scale_min)
            elif remaining > 0:
                data[i] = scale_min + remaining
                remaining = 0
                break

        return float(np.std(data, ddof=1)) if n > 1 else 0.0

    def _generate_details(self, mean, sd, n, s_min, s_max,
                          consistent, n_solutions, best_mean, best_sd):
        status = "CONSISTENT ✓" if consistent else "INCONSISTENT ✗"
        return (
            f"SPRITE Test Result: {status}\n"
            f"  Reported: M={mean}, SD={sd}, N={n}\n"
            f"  Scale range: [{s_min}, {s_max}]\n"
            f"  Solutions found: {n_solutions}/{self.num_attempts}\n"
            f"  Best achieved: M={best_mean:.4f}, SD={best_sd:.4f}"
        )
