"""
P-hacking 検出モジュール

検出指標:
1. P-curve分析 (Simonsohn et al., 2014)
2. Z-curve分析 (Brunner & Schimmack, 2020)
3. Caliper test (Masicampo & Lalande, 2012)
4. Excess significance test (Ioannidis & Trikalinos, 2007)
5. 多重比較の報告パターン分析
"""

import math
import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple


@dataclass
class PCurveResult:
    """P-curve分析結果"""
    p_values: List[float]
    num_significant: int
    right_skew_z: float  # 右偏度のZ値
    flat_test_z: float   # 平坦性のZ値
    has_evidential_value: bool  # 真の効果がある証拠
    is_consistent_with_phacking: bool
    power_estimate: float  # 推定検出力
    binomial_test_p: float
    details: str


@dataclass
class ZCurveResult:
    """Z-curve分析結果"""
    z_values: List[float]
    expected_discovery_rate: float  # EDR
    expected_replication_rate: float  # ERR
    observed_discovery_rate: float  # ODR
    file_drawer_ratio: float
    is_inflated: bool
    details: str


@dataclass
class PHackingResult:
    """P-hacking検出の総合結果"""
    p_curve: Optional[PCurveResult] = None
    z_curve: Optional[ZCurveResult] = None
    caliper_suspicious: bool = False
    excess_significance: bool = False
    multiple_testing_issues: List[str] = field(default_factory=list)
    overall_risk: str = "low"  # low, moderate, high, very_high
    overall_score: float = 0.0  # 0.0 - 1.0
    flags: List[str] = field(default_factory=list)


class PHackingDetector:
    """
    P-hacking の検出と分析。

    Parameters
    ----------
    significance_level : float
        有意水準（デフォルト: 0.05）
    caliper_width : float
        Caliper testの幅（デフォルト: 0.005）
    """

    def __init__(self, significance_level: float = 0.05,
                 caliper_width: float = 0.005):
        self.alpha = significance_level
        self.caliper_width = caliper_width

    def analyze(self, p_values: List[float],
                test_statistics: Optional[List[Tuple[float, float]]] = None
                ) -> PHackingResult:
        """
        P-hacking分析を実行。

        Parameters
        ----------
        p_values : list of float
            論文から抽出されたp値のリスト
        test_statistics : list of (statistic, df), optional
            検定統計量と自由度のペア
        """
        result = PHackingResult()

        if len(p_values) < 3:
            result.flags.append("Insufficient p-values for analysis (need ≥ 3)")
            return result

        # 1. P-curve分析
        sig_p = [p for p in p_values if 0 < p < self.alpha]
        if len(sig_p) >= 3:
            result.p_curve = self._p_curve_analysis(sig_p)
            if result.p_curve.is_consistent_with_phacking:
                result.flags.append("P-curve suggests p-hacking")
                result.overall_score += 0.3

        # 2. Z-curve分析
        if test_statistics:
            z_values = [self._to_z(stat, df) for stat, df in test_statistics]
            z_values = [z for z in z_values if z is not None]
            if len(z_values) >= 3:
                result.z_curve = self._z_curve_analysis(z_values)
                if result.z_curve.is_inflated:
                    result.flags.append("Z-curve indicates inflation")
                    result.overall_score += 0.2
        else:
            z_values = [self._p_to_z(p) for p in p_values if 0 < p < 1]
            z_values = [z for z in z_values if z is not None]
            if len(z_values) >= 3:
                result.z_curve = self._z_curve_analysis(z_values)
                if result.z_curve.is_inflated:
                    result.flags.append("Z-curve indicates inflation")
                    result.overall_score += 0.2

        # 3. Caliper test
        result.caliper_suspicious = self._caliper_test(p_values)
        if result.caliper_suspicious:
            result.flags.append(
                f"Caliper test: suspicious clustering around p={self.alpha}"
            )
            result.overall_score += 0.25

        # 4. Excess significance
        result.excess_significance = self._excess_significance_test(p_values)
        if result.excess_significance:
            result.flags.append("Excess significance detected")
            result.overall_score += 0.15

        # 5. 多重比較
        result.multiple_testing_issues = self._check_multiple_testing(p_values)
        if result.multiple_testing_issues:
            result.overall_score += 0.1 * len(result.multiple_testing_issues)

        # 総合リスク
        result.overall_score = min(result.overall_score, 1.0)
        if result.overall_score > 0.6:
            result.overall_risk = "very_high"
        elif result.overall_score > 0.4:
            result.overall_risk = "high"
        elif result.overall_score > 0.2:
            result.overall_risk = "moderate"
        else:
            result.overall_risk = "low"

        return result

    def _p_curve_analysis(self, sig_p_values: List[float]) -> PCurveResult:
        """P-curve分析（Simonsohn et al., 2014）"""
        n = len(sig_p_values)

        # 有意なp値を (0, alpha) の一様分布と比較
        # pp値: p値のp値（一様分布下での累積確率）
        pp_values = [p / self.alpha for p in sig_p_values]

        # 右偏度テスト: 多くのpp値が0.5未満なら真の効果あり
        below_half = sum(1 for pp in pp_values if pp < 0.5)
        # 二項検定
        binomial_p = self._binomial_test(below_half, n, 0.5)

        # 検出力推定（pp値の中央値から）
        median_pp = float(np.median(pp_values))
        power_estimate = max(0.0, 1.0 - 2 * median_pp)

        # Z値の計算
        right_skew_z = (below_half - n * 0.5) / math.sqrt(n * 0.25)
        flat_z = abs(right_skew_z)

        has_evidential = right_skew_z > 1.96
        is_phacking = right_skew_z < -1.96  # 左偏はp-hackingの兆候

        return PCurveResult(
            p_values=sig_p_values,
            num_significant=n,
            right_skew_z=right_skew_z,
            flat_test_z=flat_z,
            has_evidential_value=has_evidential,
            is_consistent_with_phacking=is_phacking,
            power_estimate=power_estimate,
            binomial_test_p=binomial_p,
            details=(
                f"P-curve: {n} significant p-values\n"
                f"Right-skew Z = {right_skew_z:.2f}\n"
                f"Power estimate: {power_estimate:.1%}\n"
                f"Evidential value: {has_evidential}\n"
                f"P-hacking consistent: {is_phacking}"
            ),
        )

    def _z_curve_analysis(self, z_values: List[float]) -> ZCurveResult:
        """Z-curve分析"""
        z_arr = np.array([abs(z) for z in z_values])
        z_critical = abs(self._p_to_z(self.alpha))

        # 観察された発見率
        odr = float(np.mean(z_arr > z_critical))

        # 期待発見率（EM推定の簡略版）
        sig_z = z_arr[z_arr > z_critical]
        if len(sig_z) > 0:
            mean_sig_z = float(np.mean(sig_z))
            edr = float(1 - self._normal_cdf(z_critical - mean_sig_z + z_critical))
            edr = max(0.05, min(1.0, edr))
        else:
            edr = 0.05

        # 期待再現率
        err = edr * 0.8  # 簡略化

        # ファイルドロワー比率
        if edr > 0:
            file_drawer = max(0, (odr - edr) / edr)
        else:
            file_drawer = 0.0

        is_inflated = odr > edr * 1.5 and odr > 0.3

        return ZCurveResult(
            z_values=z_values,
            expected_discovery_rate=edr,
            expected_replication_rate=err,
            observed_discovery_rate=odr,
            file_drawer_ratio=file_drawer,
            is_inflated=is_inflated,
            details=(
                f"Z-curve: {len(z_values)} tests\n"
                f"ODR: {odr:.1%}, EDR: {edr:.1%}\n"
                f"ERR: {err:.1%}\n"
                f"File drawer ratio: {file_drawer:.2f}\n"
                f"Inflated: {is_inflated}"
            ),
        )

    def _caliper_test(self, p_values: List[float]) -> bool:
        """Caliper test: p=0.05周辺の不自然な集中"""
        just_below = sum(
            1 for p in p_values
            if self.alpha - self.caliper_width <= p < self.alpha
        )
        just_above = sum(
            1 for p in p_values
            if self.alpha <= p < self.alpha + self.caliper_width
        )
        # 比率が3:1以上なら疑わしい
        ratio = (just_below + 1) / (just_above + 1)
        return ratio > 3.0

    def _excess_significance_test(self, p_values: List[float]) -> bool:
        """過剰有意性テスト"""
        n = len(p_values)
        observed_sig = sum(1 for p in p_values if p < self.alpha)

        # 期待有意数（中程度の検出力を仮定: 50%）
        expected_power = 0.5
        expected_sig = n * expected_power

        if n < 5:
            return False

        # 二項検定
        p = self._binomial_test(observed_sig, n, expected_power)
        return p < 0.10 and observed_sig > expected_sig

    def _check_multiple_testing(self, p_values: List[float]) -> List[str]:
        """多重検定の問題を検出"""
        issues = []
        n_tests = len(p_values)

        if n_tests > 10:
            sig_count = sum(1 for p in p_values if p < self.alpha)
            # Bonferroni補正後
            bonferroni_sig = sum(
                1 for p in p_values if p < self.alpha / n_tests
            )
            if sig_count > 0 and bonferroni_sig == 0:
                issues.append(
                    f"{sig_count} tests significant at α={self.alpha}, "
                    f"but none survive Bonferroni correction "
                    f"(α={self.alpha/n_tests:.4f})"
                )

        # p値が.05直下に集中
        near_boundary = sum(
            1 for p in p_values if 0.04 <= p < 0.05
        )
        if near_boundary >= 3:
            issues.append(
                f"{near_boundary} p-values in [.04, .05) range"
            )

        return issues

    def _to_z(self, statistic: float, df: float) -> Optional[float]:
        """検定統計量をZ値に変換"""
        try:
            if df > 30:
                return statistic * (1 - 1/(4*df)) / math.sqrt(1 + statistic**2/(2*df))
            else:
                p = 2 * (1 - self._t_cdf(abs(statistic), df))
                return self._p_to_z(p)
        except (ValueError, OverflowError):
            return None

    def _p_to_z(self, p: float) -> Optional[float]:
        """p値をZ値に変換"""
        if p <= 0 or p >= 1:
            return None
        # Abramowitz & Stegun近似
        if p > 0.5:
            return -self._p_to_z(1 - p)
        t = math.sqrt(-2 * math.log(p))
        c0, c1, c2 = 2.515517, 0.802853, 0.010328
        d1, d2, d3 = 1.432788, 0.189269, 0.001308
        z = t - (c0 + c1*t + c2*t**2) / (1 + d1*t + d2*t**2 + d3*t**3)
        return z

    def _normal_cdf(self, x: float) -> float:
        return 0.5 * (1 + math.erf(x / math.sqrt(2)))

    def _t_cdf(self, t: float, df: float) -> float:
        """t分布の累積分布関数（正規近似）"""
        z = t * (1 - 1/(4*df)) / math.sqrt(1 + t**2/(2*df))
        return self._normal_cdf(z)

    def _binomial_test(self, k: int, n: int, p: float) -> float:
        """二項検定のp値（正規近似）"""
        if n == 0:
            return 1.0
        mean = n * p
        std = math.sqrt(n * p * (1 - p))
        if std == 0:
            return 1.0
        z = (k - mean) / std
        return 2 * (1 - self._normal_cdf(abs(z)))
