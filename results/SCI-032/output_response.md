Good. Now let me build the full simulation framework.Now let me create the main simulation runner and visualization scripts:Now run the simulations:Simulations complete. Now generating figures:All 7 figures generated. Now creating the final report:全ファイルが揃いました。以下に実験結果をまとめます。

---

## 実験完了サマリー

**総実行時間**: ~303秒 | **フレームワーク**: Stim 1.16.0 / PyMatching 2.3.1

### 主要な数値結果

| 項目 | 結果 |
|------|------|
| **MWPM 閾値エラー率** | **0.72% ± 0.07%**（理論値 ~0.57%） |
| d=9, p=0.1% の論理エラー率 | < 0.05%（10,000ショット以内で0件） |
| PyMatching デコード速度 | **0.5〜26 μs/shot**（d=3〜9） |
| Python UF デコーダ速度 | 188μs〜466ms/shot（C++比 3,000〜18,000x 遅い） |
| 格子手術 CNOT (d=7, p=0.1%) | **~0.00%** エラー率 |
| T ゲート (15-to-1蒸留, d=7, p=0.1%) | 理論値 < 10⁻⁸ |

### 生成ファイル構造

```
workspace/
├── report.md                    ← 本レポート
├── src/                         ← 全ソースコード（6モジュール）
├── run_simulations.py           ← メインランナー
├── plot_results.py              ← 可視化スクリプト
├── figures/                     ← 7図 × PNG+SVG = 14ファイル
├── results/                     ← 7つの JSON 結果ファイル
└── logs/process-log.jsonl       ← 実行トレース
```

### 重要な知見

1. **閾値確認**：p < ~0.7% で大きな距離ほど論理エラー率が低下 → QEC利得を実証
2. **脱分極 vs 偏り雑音**：T1/T2*偏りの雑音はMWPMで純脱分極より若干低い論理エラー率
3. **非パウリ雑音**：リーケージのPauli近似では過少評価になる可能性あり（密度行列シミュレーションが必要）
4. **格子手術**：d=7, p=0.1% では理論的にT ゲートエラー率 < 10⁻⁸ を達成可能