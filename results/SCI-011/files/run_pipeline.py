"""
run_pipeline.py
---------------
全脳コネクトーム解析パイプライン - メインランナー
全ステップを順次実行し、results/ と report.md を生成する。
"""

import sys
import os
import json
import time
from datetime import datetime

# プロジェクトルートを sys.path に追加
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

LOGS_DIR = "logs"
os.makedirs(LOGS_DIR, exist_ok=True)


def log_event(event: dict):
    with open(f"{LOGS_DIR}/process-log.jsonl", "a") as f:
        f.write(json.dumps(event, ensure_ascii=False) + "\n")


def main():
    start_time = time.time()
    print("=" * 60)
    print("  全脳コネクトーム解析パイプライン 開始")
    print(f"  {datetime.utcnow().isoformat()} UTC")
    print("=" * 60)

    log_event({
        "timestamp": datetime.utcnow().isoformat(),
        "phase": "pipeline",
        "event_type": "run_started",
        "actor": "co-scientist",
        "skill_or_tool": "run_pipeline.py",
        "status": "ok",
    })

    # ─── ステップ 1: 前処理 ─────────────────────────────────────────────────
    print("\n[STEP 1/7] 前処理パイプライン")
    from src._01_preprocessing import main as step1
    qc_summary, optim = step1()

    # ─── ステップ 2: 構造的コネクティビティ ─────────────────────────────────
    print("\n[STEP 2/7] 構造的コネクティビティ")
    from src._02_structural_connectivity import main as step2
    sc_data = step2()

    # ─── ステップ 3: 機能的コネクティビティ ─────────────────────────────────
    print("\n[STEP 3/7] 機能的コネクティビティ")
    from src._03_functional_connectivity import main as step3
    fc_data, dfc_results, ts_data = step3()

    # ─── ステップ 4: グラフ理論解析 ─────────────────────────────────────────
    print("\n[STEP 4/7] グラフ理論解析")
    from src._04_graph_analysis import main as step4
    graph_metrics = step4()

    # ─── ステップ 5: バイオマーカー ─────────────────────────────────────────
    print("\n[STEP 5/7] 疾患バイオマーカー同定")
    from src._05_biomarkers import main as step5
    biomarker_results = step5()

    # ─── ステップ 6: 信頼性評価 ──────────────────────────────────────────────
    print("\n[STEP 6/7] テスト-リテスト信頼性")
    from src._06_reliability import main as step6
    reliability_results = step6()

    # ─── ステップ 7: 可視化 ──────────────────────────────────────────────────
    print("\n[STEP 7/7] 可視化")
    from src._07_visualization import main as step7
    step7()

    elapsed = time.time() - start_time
    print(f"\n{'='*60}")
    print(f"  パイプライン完了 ({elapsed:.1f}秒)")
    print("=" * 60)

    log_event({
        "timestamp": datetime.utcnow().isoformat(),
        "phase": "pipeline",
        "event_type": "run_completed",
        "actor": "co-scientist",
        "elapsed_seconds": round(elapsed, 1),
        "status": "ok",
    })

    return {
        "qc_summary": qc_summary,
        "graph_metrics": graph_metrics,
        "biomarker_results": biomarker_results,
        "reliability_results": reliability_results,
    }


if __name__ == "__main__":
    results = main()
