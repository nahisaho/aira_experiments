"""
RIAS デモスクリプト — 統合パイプラインの動作確認

合成データを使用して全モジュールの動作を検証する。
"""

import sys
import json
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from src.image_forensics import ELAAnalyzer, DuplicateDetector, ManipulationDetector
from src.image_forensics.model import ImageForensicsNet, ModelConfig
from src.statistical_checks import GRIMTest, SPRITETest, StatisticalAnalyzer
from src.plagiarism import CitationAwareSimilarity, TextFingerprinter
from src.phacking import PHackingDetector, HARKingDetector, PHackingMetaAnalyzer
from src.reproducibility import ReproducibilityScorer, MethodologyAssessor
from src.pipeline import IntegrityPipeline
from src.pipeline.integrity_pipeline import PaperInput


def demo_image_forensics():
    """画像フォレンジクスのデモ"""
    print("=" * 60)
    print("Module 1: Image Forensics")
    print("=" * 60)

    # 合成テスト画像（小サイズでデモ高速化）
    np.random.seed(42)
    clean_image = np.random.randint(0, 255, (64, 64, 3), dtype=np.uint8)

    # コピー&ムーブ偽造のシミュレーション
    manipulated = clean_image.copy()
    manipulated[10:30, 10:30] = manipulated[35:55, 35:55]

    # ELA分析
    ela = ELAAnalyzer(quality=95)
    result = ela.analyze(clean_image)
    print(f"\n[ELA] Clean image:")
    print(f"  Mean error: {result.mean_error:.2f}")
    print(f"  Suspicious ratio: {result.suspicious_ratio:.4f}")

    result_manip = ela.analyze(manipulated)
    print(f"\n[ELA] Manipulated image:")
    print(f"  Mean error: {result_manip.mean_error:.2f}")
    print(f"  Suspicious ratio: {result_manip.suspicious_ratio:.4f}")

    # 重複検出
    detector = DuplicateDetector()
    images = {"fig1": clean_image, "fig2": clean_image.copy(), "fig3": manipulated}
    duplicates = detector.find_duplicates(images)
    print(f"\n[Duplicate Detection]")
    print(f"  Duplicates found: {len(duplicates)}")
    for d in duplicates:
        print(f"  - {d.image_id_a} ↔ {d.image_id_b}: "
              f"{d.duplication_type.value} (sim={d.similarity_score:.2f})")

    # モデルアーキテクチャ
    model = ImageForensicsNet(ModelConfig(backbone="resnet50"))
    print(f"\n[CNN Model Architecture]")
    print(model.get_architecture_summary())

    return True


def demo_statistical_checks():
    """統計チェックのデモ"""
    print("\n" + "=" * 60)
    print("Module 2: Statistical Checks (GRIM/SPRITE)")
    print("=" * 60)

    # GRIM Test
    grim = GRIMTest()
    print("\n[GRIM Test]")

    test_cases = [
        (3.48, 25),   # 整合
        (3.47, 25),   # 不整合 (3.47 × 25 = 86.75, not integer)
        (4.57, 30),   # チェック
        (2.33, 15),   # チェック
    ]

    for mean, n in test_cases:
        result = grim.test(mean, n)
        status = "✓" if result.is_consistent else "✗"
        print(f"  M={mean}, N={n}: {status} (closest: {result.closest_possible_mean})")

    # SPRITE Test
    sprite = SPRITETest(max_iterations=5000, num_attempts=50)
    print("\n[SPRITE Test]")

    sprite_cases = [
        (3.50, 1.20, 20, 1, 5),  # 整合的
        (3.50, 0.10, 20, 1, 5),  # 非常に小さいSD — 疑わしい
        (4.00, 2.50, 10, 1, 7),  # チェック
    ]

    for mean, sd, n, s_min, s_max in sprite_cases:
        result = sprite.test(mean, sd, n, s_min, s_max, seed=42)
        status = "✓" if result.is_consistent else "✗"
        print(f"  M={mean}, SD={sd}, N={n}, [{s_min}-{s_max}]: {status} "
              f"(solutions: {result.num_solutions_found})")

    # テキストからの統計抽出
    analyzer = StatisticalAnalyzer()
    sample_text = """
    The results showed a significant effect, t(48) = 2.45, p = .018.
    Group A (M = 3.47, SD = 1.23, N = 25) scored higher than
    Group B (M = 2.89, SD = 1.15, N = 25). An ANOVA revealed
    F(2, 87) = 4.56, p = .013. Post-hoc tests showed p = .042
    and p = .048.
    """
    result = analyzer.analyze_text(sample_text)
    print(f"\n[Statistical Text Analysis]")
    print(f"  Tests checked: {result.num_tests_checked}")
    print(f"  Inconsistencies: {result.num_inconsistencies}")
    print(f"  Flags: {result.flags}")

    return True


def demo_plagiarism():
    """盗作検出のデモ"""
    print("\n" + "=" * 60)
    print("Module 3: Plagiarism Detection")
    print("=" * 60)

    similarity = CitationAwareSimilarity()

    source = (
        "The experiment was conducted using a randomized controlled design. "
        "Participants were randomly assigned to one of three conditions. "
        "Data were collected over a period of six weeks. "
        "Statistical analyses were performed using SPSS version 25."
    )

    target_verbatim = (
        "The experiment was conducted using a randomized controlled design. "
        "Participants were randomly assigned to one of three conditions. "
        "Data were collected over a period of six weeks."
    )

    target_cited = (
        "As described by Smith et al. (2020), the experiment was conducted "
        "using a randomized controlled design. Participants were randomly "
        "assigned to one of three conditions (Johnson & Lee, 2019)."
    )

    target_paraphrase = (
        "A randomized controlled methodology was employed for this study. "
        "The subjects were placed into three different groups at random. "
        "The data collection phase lasted approximately six weeks."
    )

    print("\n[Verbatim Copy]")
    report = similarity.compare(source, target_verbatim)
    print(f"  Overall similarity: {report.overall_similarity:.1%}")
    print(f"  Risk level: {report.risk_level}")

    print("\n[Cited Copy]")
    report = similarity.compare(source, target_cited)
    print(f"  Overall similarity: {report.overall_similarity:.1%}")
    print(f"  Adjusted (excl. citations): {report.adjusted_similarity:.1%}")
    print(f"  Risk level: {report.risk_level}")

    print("\n[Paraphrase]")
    report = similarity.compare(source, target_paraphrase)
    print(f"  Overall similarity: {report.overall_similarity:.1%}")
    print(f"  Risk level: {report.risk_level}")

    # フィンガープリント
    fp = TextFingerprinter(k=3)
    fp1 = fp.fingerprint(source, "source")
    fp2 = fp.fingerprint(target_verbatim, "verbatim")
    fp3 = fp.fingerprint(target_paraphrase, "paraphrase")

    print(f"\n[Fingerprint Similarity]")
    print(f"  Source vs Verbatim: {fp.similarity(fp1, fp2):.2f}")
    print(f"  Source vs Paraphrase: {fp.similarity(fp1, fp3):.2f}")

    return True


def demo_phacking():
    """P-hacking検出のデモ"""
    print("\n" + "=" * 60)
    print("Module 4: P-hacking / HARKing Detection")
    print("=" * 60)

    detector = PHackingDetector()

    # 正常なp値分布
    np.random.seed(42)
    normal_p = np.random.uniform(0, 1, 20).tolist()
    normal_p[:5] = [0.001, 0.008, 0.012, 0.023, 0.031]  # 一部有意

    # P-hackingが疑われる分布（.05直下に集中）
    suspicious_p = [
        0.048, 0.049, 0.042, 0.045, 0.039, 0.044, 0.047,
        0.12, 0.34, 0.67, 0.89, 0.23, 0.56, 0.78,
    ]

    print("\n[Normal P-value Distribution]")
    result = detector.analyze(normal_p)
    print(f"  Risk: {result.overall_risk} (score={result.overall_score:.2f})")
    print(f"  Caliper suspicious: {result.caliper_suspicious}")

    print("\n[Suspicious P-value Distribution]")
    result = detector.analyze(suspicious_p)
    print(f"  Risk: {result.overall_risk} (score={result.overall_score:.2f})")
    print(f"  Caliper suspicious: {result.caliper_suspicious}")
    print(f"  Flags: {result.flags}")

    # HARKing
    harking = HARKingDetector()

    intro = (
        "We hypothesized that treatment A would significantly improve "
        "outcomes compared to control. As expected, based on prior "
        "literature, we predicted a strong positive effect. "
        "We anticipated that the effect would be moderated by age."
    )

    results_text = (
        "Consistent with our hypothesis, treatment A significantly "
        "improved outcomes (p < .001). As predicted, the effect was "
        "robust and compelling. Confirming our hypothesis, the "
        "moderation by age was significant and strong."
    )

    hk_result = harking.analyze(intro, results_text)
    print(f"\n[HARKing Detection]")
    print(f"  Risk: {hk_result.risk_level} (score={hk_result.risk_score:.2f})")
    print(f"  Predictive ratio: {hk_result.predictive_ratio:.1%}")
    print(f"  Flags: {hk_result.flags}")

    return True


def demo_reproducibility():
    """再現性予測のデモ"""
    print("\n" + "=" * 60)
    print("Module 5: Reproducibility Score")
    print("=" * 60)

    assessor = MethodologyAssessor()

    # 詳細な方法セクション
    good_methods = """
    Participants (N = 120) were recruited from the university campus.
    Sample size was determined by a priori power analysis using G*Power
    (effect size d = 0.5, alpha = 0.05, power = 0.80). Participants
    were randomly assigned to conditions using a computer-generated
    randomization sequence. The study was double-blind; neither
    participants nor experimenters knew the assignment. Inclusion
    criteria required age 18-65 and no prior diagnosis. Exclusion
    criteria included current medication use. Data were analyzed
    using R version 4.2.1 with the lme4 package. A mixed-effects
    model was fitted with Bonferroni correction for multiple comparisons.
    Effect sizes (Cohen's d) and 95% confidence intervals are reported.
    Three biological replicates were performed for each condition.
    Raw data are available on Zenodo (DOI: 10.5281/zenodo.XXXXXXX).
    Analysis code is available at github.com/example/repo.
    The study was preregistered on OSF (osf.io/XXXXX).
    IRB approval was obtained (Protocol #2023-001).
    """

    poor_methods = """
    We tested our hypothesis using standard methods. Data were
    collected and analyzed. Statistical tests were performed.
    Results were significant.
    """

    print("\n[Good Methods Section]")
    result = assessor.assess(good_methods)
    print(f"  Quality: {result.quality_level}")
    print(f"  Score: {result.overall_score:.2f}")
    print(f"  Reproducibility prediction: {result.reproducibility_prediction:.1%}")
    print(f"  Present: {len(result.present_elements)}/{len(assessor.CHECKLIST)}")

    print("\n[Poor Methods Section]")
    result = assessor.assess(poor_methods)
    print(f"  Quality: {result.quality_level}")
    print(f"  Score: {result.overall_score:.2f}")
    print(f"  Reproducibility prediction: {result.reproducibility_prediction:.1%}")
    print(f"  Missing: {result.missing_elements[:5]}")

    # 統合スコア
    scorer = ReproducibilityScorer()

    print("\n[Integrated Reproducibility Score]")
    # 高品質論文
    score = scorer.score(
        methodology_quality=0.85,
        statistical_consistency=0.90,
        phacking_risk=0.10,
        effect_size_plausibility=0.80,
        sample_size_adequacy=0.75,
        preregistration=1.0,
        data_code_availability=1.0,
    )
    print(f"  High-quality paper: {score.overall_score:.2f} "
          f"[{score.confidence_interval[0]:.2f}-{score.confidence_interval[1]:.2f}] "
          f"({score.prediction_class})")

    # 低品質論文
    score = scorer.score(
        methodology_quality=0.20,
        statistical_consistency=0.30,
        phacking_risk=0.70,
        effect_size_plausibility=0.30,
        sample_size_adequacy=0.20,
        preregistration=0.0,
        data_code_availability=0.0,
    )
    print(f"  Low-quality paper:  {score.overall_score:.2f} "
          f"[{score.confidence_interval[0]:.2f}-{score.confidence_interval[1]:.2f}] "
          f"({score.prediction_class})")

    return True


def demo_integrated_pipeline():
    """統合パイプラインのデモ"""
    print("\n" + "=" * 60)
    print("Module 6: Integrated Pipeline")
    print("=" * 60)

    pipeline = IntegrityPipeline()

    paper = PaperInput(
        doi="10.1234/example.2024",
        title="Effects of Treatment X on Outcome Y: A Randomized Trial",
        authors=["Smith, J.", "Doe, A."],
        sections={
            "introduction": (
                "We hypothesized that treatment X would improve outcome Y. "
                "Based on prior work, we predicted a moderate effect. "
                "The theoretical framework suggests positive outcomes."
            ),
            "methods": (
                "Participants (N = 60) were randomly assigned to two groups. "
                "The study used a double-blind design. Data were analyzed "
                "using t-tests with alpha = 0.05. Effect sizes are reported "
                "as Cohen's d. R version 4.1 was used for analysis."
            ),
            "results": (
                "The main effect was significant, t(58) = 2.12, p = .038. "
                "Group A (M = 4.23, SD = 1.45, N = 30) outperformed "
                "Group B (M = 3.67, SD = 1.38, N = 30). A secondary "
                "analysis showed F(2, 57) = 3.45, p = .039."
            ),
            "discussion": (
                "These findings support our hypothesis. The effect size "
                "was moderate, consistent with prior predictions. "
                "Limitations include the sample size."
            ),
        },
        p_values=[0.038, 0.039, 0.045, 0.23, 0.67],
    )

    report = pipeline.analyze(paper)

    print(f"\n[Integrity Report]")
    print(f"  DOI: {report.doi}")
    print(f"  Integrity Score: {report.overall_integrity_score:.2f}")
    print(f"  Risk Level: {report.overall_risk_level}")
    print(f"  Reproducibility: {report.reproducibility_prediction:.1%}")
    print(f"  Confidence: {report.confidence:.1%}")
    print(f"\n  Module Results:")
    for name, result in report.module_results.items():
        print(f"    {name}: {result.risk_level} (score={result.risk_score:.2f})")
    print(f"\n  Top Concerns:")
    for c in report.top_concerns[:3]:
        print(f"    - {c}")
    print(f"\n  Recommendations:")
    for r in report.recommendations[:3]:
        print(f"    - {r}")

    return report


def main():
    """全モジュールのデモを実行"""
    print("=" * 60)
    print("Research Integrity Assessment System (RIAS) v1.0.0")
    print("NLP + Computer Vision 統合型 研究公正性評価システム")
    print("=" * 60)

    results = {}

    results["image_forensics"] = demo_image_forensics()
    results["statistical_checks"] = demo_statistical_checks()
    results["plagiarism"] = demo_plagiarism()
    results["phacking"] = demo_phacking()
    results["reproducibility"] = demo_reproducibility()
    report = demo_integrated_pipeline()
    results["pipeline"] = report is not None

    print("\n" + "=" * 60)
    print("Demo Summary")
    print("=" * 60)
    for module, success in results.items():
        status = "✓ PASS" if success else "✗ FAIL"
        print(f"  {module}: {status}")

    # 結果をJSONとして保存
    output_dir = Path(__file__).parent / "results"
    output_dir.mkdir(exist_ok=True)

    summary = {
        "system": "RIAS v1.0.0",
        "modules_tested": len(results),
        "all_passed": all(results.values()),
        "demo_report": {
            "doi": report.doi if report else None,
            "integrity_score": report.overall_integrity_score if report else None,
            "risk_level": report.overall_risk_level if report else None,
            "reproducibility": report.reproducibility_prediction if report else None,
        } if report else None,
    }

    with open(output_dir / "demo_results.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n  Results saved to {output_dir / 'demo_results.json'}")
    return all(results.values())


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
