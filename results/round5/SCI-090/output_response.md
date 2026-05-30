Now I have enough foundation to proceed. Let me also search for a few more specific papers and then build the simulation system.Now I have sufficient literature. Let me implement the simulation system and create the required output files.The results have several physically unrealistic values that need fixing:
- ACH of 50-79 is unrealistically high (should be 2-15 for offices)
- UDI values >100% (double multiplication bug)  
- Mean DF 144% (impossible)

Let me fix these bugs:Now fix the daylighting bugs:Let me fix the remaining daylighting attenuation issue and then rerun:Good results now. CFD ACH is 1.9–5.2 (realistic), ASE reduced significantly. Let me now create the paper and report files.すべてのファイルが生成されました。以下にまとめます。

---

## 完了サマリー

### ステップ1: 先行研究調査（7件確認）

| # | 著者 | 年 | DOI | 主要知見 |
|---|------|-----|-----|---------|
| 1 | Alexandrou et al. | 2023 | 10.1080/09613218.2023.2222856 | IFC→gbXML変換で15-30%誤差 |
| 2 | Aydin & Koçlar Oral | 2025 | 10.2139/ssrn.5111348 | BIM双方向EnergyPlus連携 |
| 3 | Mehraban et al. | 2025 | 10.3390/buildings15213950 | BIM+AI熱性能予測 |
| 4 | Kharvari | 2020 | 10.1016/j.solener.2020.07.054 | Radiance昼光検証（RMSE 8-12%）|
| 5 | Li et al. | 2024 | 10.1007/s12273-024-1178-6 | 多目的最適化（採光-エネルギーTOT）|
| 6 | Mathew et al. | 2026 | 10.15627/jd.2026.9 | CFD換気：40% WWRが最適 |
| 7 | Sarkar & Solanki | 2025 | 10.1007/s12053-025-10398-y | ZEB設計で83%削減 |

### ステップ2〜3: 主要実験結果（5分割交差検証付き）

| 指標 | ZEB | ベースライン |
|------|-----|------------|
| Site EUI [kWh/(m²·yr)] | **217.2 ± 0.9** | 508.7 ± 3.2 |
| ZEB比 | 0.253 ± 0.002 | 0 |
| 年間ACH（平均） | 4.05 | 2.1 |
| sDA(300/50%) | 100% | 100% |
| ASE(1000,250h) | **40.3%** | 88.9% |
| 平均DF | 4.79% | 5.99% |

**⚠️ 自己批判**: ZEB比25.3%は目標（100%）未達。EUI推定値は簡略化モデルのため実際より40-60%高い可能性。sDA=100%は昼光モデルの過剰楽観性が疑われる。

### ステップ4: 成果物

- **`paper.md`** — 学術論文（英語、Abstract 250語以上、References 10件、DOI付き）
- **`report.md`** — 日本語実験レポート
- **7枚のFigure** (`figures/fig1〜fig7.png`)
- **ソースコード** (`src/` 4モジュール + `run_simulation.py`)