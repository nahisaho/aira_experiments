Now let's create the main analysis runner and visualization script:All 6 analyses completed successfully. Now let's generate the report:全ての解析が正常に完了しました。以下に結果をまとめます。

---

## SOA 反応ネットワーク解析システム — 完了

### 主要な結果

| モジュール | 成果 |
|-----------|------|
| **反応ネットワーク** | 41 化学種 × 45 反応（3 世代、SOA 前駆体 27 種） |
| **分配熱力学** | ELVOC（Fpart ≈ 1.0）〜 IVOC（Fpart < 0.001）を 34 種で計算 |
| **ML 速度定数** | R² = 0.997, RMSE = 0.034 log 単位（Evans-Polanyi GPR） |
| **箱モデル (8h)** | α-ピネン: 6.66 μg/m³、リモネン: 3.15、イソプレン: 18.0 μg/m³ |
| **感度解析** | **温度**が最支配的（Sobol S₁ = 0.67）、VOC 濃度が 2 番目 |
| **SOA 収率** | リモネン最高（OH: 0.26、O₃: 0.28）、イソプレン最低（OH: 0.028） |

### 生成ファイル（計 18 ファイル）
- `report.md` — 本レポート（日本語、全手法・考察）
- `figures/fig01〜06_*.png` — 6 枚の出版品質図表
- `results/` — 統計サマリー、収率テーブル、感度解析 JSON
- `data/` — 反応ネットワーク JSON、分配 CSV、ML 予測 CSV、箱モデル時系列 5 種
- `logs/process-log.jsonl` — 完全実行トレース
- `src/` — 6 モジュールの完全 Python ソースコード