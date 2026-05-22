"""
統計的不整合の包括的分析モジュール

検出対象:
1. GRIM不整合（平均値の粒度チェック）
2. SPRITE不整合（記述統計の整合性）
3. p値分布の異常（Caliper test）
4. 効果量の異常パターン
5. 自由度と検定統計量の不整合
6. サンプルサイズの丸め/繰り返しパターン
"""

import re
import math
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple
from .grim_test import GRIMTest, GRIMResult
from .sprite_test import SPRITETest, SPRITEResult


@dataclass
class StatisticalExtraction:
    """テキストから抽出された統計値"""
    test_type: str  # t, F, chi2, r, etc.
    statistic_value: float
    df1: Optional[float] = None
    df2: Optional[float] = None
    p_value: Optional[float] = None
    effect_size: Optional[float] = None
    sample_size: Optional[int] = None
    mean: Optional[float] = None
    sd: Optional[float] = None
    source_text: str = ""


@dataclass
class PValueDistribution:
    """p値分布の分析結果"""
    p_values: List[float]
    num_significant: int
    num_total: int
    proportion_significant: float
    caliper_test_result: Optional[Dict] = None
    uniformity_test: Optional[Dict] = None
    suspicious_clustering: bool = False
    details: str = ""


@dataclass
class StatisticalAnalysisResult:
    """統計的不整合分析の総合結果"""
    grim_results: List[GRIMResult] = field(default_factory=list)
    sprite_results: List[SPRITEResult] = field(default_factory=list)
    p_value_analysis: Optional[PValueDistribution] = None
    df_consistency: List[Dict] = field(default_factory=list)
    overall_inconsistency_score: float = 0.0
    num_tests_checked: int = 0
    num_inconsistencies: int = 0
    flags: List[str] = field(default_factory=list)


class StatisticalAnalyzer:
    """
    論文の統計的不整合を包括的に分析する。

    テキストから統計値を自動抽出し、複数の整合性テストを実行する。
    """

    # 統計量抽出用の正規表現パターン
    PATTERNS = {
        "t_test": re.compile(
            r't\s*\((\d+(?:\.\d+)?)\)\s*=\s*(-?\d+\.?\d*)\s*,?\s*'
            r'p\s*[=<>]\s*\.?(\d+\.?\d*)',
            re.IGNORECASE
        ),
        "f_test": re.compile(
            r'F\s*\((\d+(?:\.\d+)?)\s*,\s*(\d+(?:\.\d+)?)\)\s*=\s*'
            r'(\d+\.?\d*)\s*,?\s*p\s*[=<>]\s*\.?(\d+\.?\d*)',
            re.IGNORECASE
        ),
        "chi2": re.compile(
            r'[χχ²X²]\s*\((\d+)\)\s*=\s*(\d+\.?\d*)\s*,?\s*'
            r'p\s*[=<>]\s*\.?(\d+\.?\d*)',
            re.IGNORECASE
        ),
        "correlation": re.compile(
            r'r\s*[=]\s*(-?\.?\d+\.?\d*)\s*,?\s*'
            r'p\s*[=<>]\s*\.?(\d+\.?\d*)',
            re.IGNORECASE
        ),
        "mean_sd": re.compile(
            r'M\s*=\s*(\d+\.?\d*)\s*,?\s*SD\s*=\s*(\d+\.?\d*)',
            re.IGNORECASE
        ),
        "n_size": re.compile(
            r'[Nn]\s*=\s*(\d+)',
        ),
        "p_value": re.compile(
            r'p\s*[=<>]\s*\.?(\d+\.?\d*(?:e[+-]?\d+)?)',
            re.IGNORECASE
        ),
    }

    def __init__(self):
        self.grim = GRIMTest()
        self.sprite = SPRITETest(max_iterations=5000, num_attempts=50)

    def analyze_text(self, text: str,
                     scale_range: Tuple[int, int] = (1, 7)) -> StatisticalAnalysisResult:
        """
        論文テキストから統計値を抽出し整合性を分析する。

        Parameters
        ----------
        text : str
            論文テキスト（結果セクション）
        scale_range : tuple
            尺度の範囲（SPRITE テスト用）
        """
        result = StatisticalAnalysisResult()

        # 1. 統計値の抽出
        extractions = self.extract_statistics(text)

        # 2. GRIM Test
        mean_n_pairs = [
            (e.mean, e.sample_size)
            for e in extractions
            if e.mean is not None and e.sample_size is not None
        ]
        if mean_n_pairs:
            result.grim_results = self.grim.batch_test(mean_n_pairs)
            grim_failures = sum(
                1 for r in result.grim_results if not r.is_consistent
            )
            if grim_failures > 0:
                result.flags.append(
                    f"GRIM: {grim_failures}/{len(result.grim_results)} "
                    f"means inconsistent with sample sizes"
                )
                result.num_inconsistencies += grim_failures

        # 3. SPRITE Test
        for e in extractions:
            if (e.mean is not None and e.sd is not None
                    and e.sample_size is not None):
                sprite_result = self.sprite.test(
                    e.mean, e.sd, e.sample_size,
                    scale_min=scale_range[0], scale_max=scale_range[1]
                )
                result.sprite_results.append(sprite_result)
                if not sprite_result.is_consistent:
                    result.num_inconsistencies += 1
                    result.flags.append(
                        f"SPRITE: M={e.mean}, SD={e.sd}, N={e.sample_size} "
                        f"inconsistent"
                    )

        # 4. p値分布分析
        p_values = [
            e.p_value for e in extractions if e.p_value is not None
        ]
        if len(p_values) >= 3:
            result.p_value_analysis = self._analyze_p_distribution(p_values)
            if result.p_value_analysis.suspicious_clustering:
                result.flags.append(
                    "P-value distribution shows suspicious clustering"
                )
                result.num_inconsistencies += 1

        # 5. 自由度の整合性チェック
        result.df_consistency = self._check_df_consistency(extractions)
        for df_issue in result.df_consistency:
            if not df_issue["consistent"]:
                result.flags.append(df_issue["message"])
                result.num_inconsistencies += 1

        # 総合スコア
        result.num_tests_checked = (
            len(result.grim_results)
            + len(result.sprite_results)
            + (1 if result.p_value_analysis else 0)
            + len(result.df_consistency)
        )

        if result.num_tests_checked > 0:
            result.overall_inconsistency_score = (
                result.num_inconsistencies / result.num_tests_checked
            )

        return result

    def extract_statistics(self, text: str) -> List[StatisticalExtraction]:
        """テキストから統計値を抽出"""
        extractions = []

        # サンプルサイズの抽出（コンテキスト全体から）
        n_matches = self.PATTERNS["n_size"].findall(text)
        default_n = int(n_matches[0]) if n_matches else None

        # t検定
        for m in self.PATTERNS["t_test"].finditer(text):
            df = float(m.group(1))
            t_val = float(m.group(2))
            p_val = float(m.group(3))
            if p_val > 1:
                p_val = p_val / (10 ** len(str(int(p_val))))
            extractions.append(StatisticalExtraction(
                test_type="t",
                statistic_value=t_val,
                df1=df,
                p_value=p_val,
                sample_size=int(df + 1) if default_n is None else default_n,
                source_text=m.group(0),
            ))

        # F検定
        for m in self.PATTERNS["f_test"].finditer(text):
            df1 = float(m.group(1))
            df2 = float(m.group(2))
            f_val = float(m.group(3))
            p_val = float(m.group(4))
            if p_val > 1:
                p_val = p_val / (10 ** len(str(int(p_val))))
            extractions.append(StatisticalExtraction(
                test_type="F",
                statistic_value=f_val,
                df1=df1,
                df2=df2,
                p_value=p_val,
                sample_size=int(df1 + df2 + 1) if default_n is None else default_n,
                source_text=m.group(0),
            ))

        # 平均値・標準偏差
        for m in self.PATTERNS["mean_sd"].finditer(text):
            mean_val = float(m.group(1))
            sd_val = float(m.group(2))
            extractions.append(StatisticalExtraction(
                test_type="descriptive",
                statistic_value=mean_val,
                mean=mean_val,
                sd=sd_val,
                sample_size=default_n,
                source_text=m.group(0),
            ))

        # p値のみ
        for m in self.PATTERNS["p_value"].finditer(text):
            p_val = float(m.group(1))
            if p_val > 1:
                p_val = p_val / (10 ** len(str(int(p_val))))
            # 既に抽出されたものと重複チェック
            existing_texts = [e.source_text for e in extractions]
            if not any(m.group(0) in t for t in existing_texts):
                extractions.append(StatisticalExtraction(
                    test_type="p_only",
                    statistic_value=0.0,
                    p_value=p_val,
                    source_text=m.group(0),
                ))

        return extractions

    def _analyze_p_distribution(self, p_values: List[float]) -> PValueDistribution:
        """p値分布の分析"""
        import numpy as np

        p_arr = np.array([p for p in p_values if 0 < p <= 1])
        if len(p_arr) == 0:
            return PValueDistribution(
                p_values=p_values,
                num_significant=0,
                num_total=len(p_values),
                proportion_significant=0.0,
            )

        num_sig = int(np.sum(p_arr < 0.05))
        prop_sig = num_sig / len(p_arr)

        # Caliper test: p=0.05付近の密度を検査
        caliper_width = 0.005
        just_below = np.sum((p_arr >= 0.045) & (p_arr < 0.05))
        just_above = np.sum((p_arr >= 0.05) & (p_arr < 0.055))
        caliper_ratio = (just_below + 1) / (just_above + 1)

        caliper_result = {
            "just_below_05": int(just_below),
            "just_above_05": int(just_above),
            "ratio": float(caliper_ratio),
            "suspicious": caliper_ratio > 3.0,
        }

        # p値分布の一様性テスト（非有意p値）
        non_sig = p_arr[p_arr >= 0.05]
        uniformity = None
        if len(non_sig) >= 5:
            # KS検定のシンプルな近似
            sorted_p = np.sort(non_sig)
            expected = np.linspace(0.05, 1.0, len(sorted_p))
            ks_stat = float(np.max(np.abs(sorted_p - expected)))
            uniformity = {
                "ks_statistic": ks_stat,
                "suspicious": ks_stat > 0.3,
            }

        suspicious = (
            caliper_result["suspicious"]
            or (uniformity is not None and uniformity["suspicious"])
            or prop_sig > 0.9
        )

        return PValueDistribution(
            p_values=p_values,
            num_significant=num_sig,
            num_total=len(p_arr),
            proportion_significant=prop_sig,
            caliper_test_result=caliper_result,
            uniformity_test=uniformity,
            suspicious_clustering=suspicious,
            details=(
                f"P-value analysis: {num_sig}/{len(p_arr)} significant\n"
                f"Caliper ratio (below/above .05): {caliper_ratio:.2f}\n"
                f"Suspicious: {suspicious}"
            ),
        )

    def _check_df_consistency(
        self, extractions: List[StatisticalExtraction]
    ) -> List[Dict]:
        """自由度と検定統計量の整合性チェック"""
        issues = []

        for e in extractions:
            if e.test_type == "t" and e.df1 is not None and e.p_value is not None:
                # t値からp値を概算して整合性を確認
                expected_p = self._approx_t_to_p(
                    abs(e.statistic_value), e.df1
                )
                if expected_p is not None and e.p_value is not None:
                    ratio = max(expected_p, e.p_value) / (min(expected_p, e.p_value) + 1e-15)
                    consistent = ratio < 10  # 1桁以内
                    issues.append({
                        "test": f"t({e.df1}) = {e.statistic_value}",
                        "reported_p": e.p_value,
                        "expected_p_approx": expected_p,
                        "consistent": consistent,
                        "message": (
                            f"df consistency: t({e.df1})={e.statistic_value}, "
                            f"reported p={e.p_value}, expected ≈{expected_p:.4f}"
                            if not consistent else ""
                        ),
                    })

        return issues

    def _approx_t_to_p(self, t_val: float, df: float) -> Optional[float]:
        """t値からp値を近似計算（Abramowitz & Stegun近似）"""
        if df <= 0:
            return None
        try:
            x = df / (df + t_val ** 2)
            # 不完全ベータ関数の近似
            if t_val == 0:
                return 1.0
            # 正規近似 (df > 30)
            if df > 30:
                z = t_val * (1 - 1/(4*df)) / math.sqrt(1 + t_val**2/(2*df))
                p = 2 * (1 - self._normal_cdf(abs(z)))
                return max(0.0, min(1.0, p))
            # 小さい自由度の場合の粗い近似
            z = t_val / math.sqrt(df)
            p = 2 * (1 - self._normal_cdf(abs(t_val) * math.sqrt(df/(df+2))))
            return max(0.0, min(1.0, p))
        except (ValueError, OverflowError):
            return None

    def _normal_cdf(self, x: float) -> float:
        """標準正規分布の累積分布関数（近似）"""
        return 0.5 * (1 + math.erf(x / math.sqrt(2)))
