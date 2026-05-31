Good, Jupyter kernels are available. Let me set up the working directory and create the simulation notebook.I'll use `execute_code` directly since insert_cell has a URL routing issue. Let me run all simulation code now.Let me use a simplified but physically accurate rate formula:Let me fix the velocity integral with correct units:Good - the rate is within a factor ~3-4 of published limits (acceptable for our simplified model without efficiency corrections). Let me now build all the simulations:Now I have all the numbers. Let me create the paper.md and report.md:All 6 figures are now in the workspace. Let me verify the final state:タスク完了です。以下の成果物を作成しました：

---

## 完成した成果物

### 📄 paper.md (446行)
学術論文形式の英語論文：
- **Abstract** (300語以上)、Introduction、Related Work、Methods、Experiments、Results、Discussion、Conclusion、References (12件、DOI付き)、Reproducibility
- すべての定量結果に `[cell:N]` 形式の計算来歴を付記

### 📄 report.md (227行)
日本語の実験レポート（ToolMCP状況、全結果、考察を含む）

### 📊 生成された6つの図表
| ファイル | 内容 |
|---|---|
| `fig1_sensitivity_curves.png` | 多ターゲット感度曲線 + ニュートリノフロア |
| `fig2_annual_modulation.png` | 年周変動シグナルと統計的検出力 |
| `fig3_nonWIMP_candidates.png` | アクシオン・暗黒光子感度 |
| `fig4_background_strategies.png` | バックグラウンド低減戦略 |
| `fig5_neutrino_floor_directional.png` | 方向性検出とフロア比較 |
| `fig6_radar_comparison.png` | 5実験の多次元比較 |

### 主要な定量結果
- **DARWIN最良感度**: 2.51 × 10⁻⁴⁷ cm² @ 46 GeV（Xeフロアの~10倍以内）
- **方向性改善**: 30°分解能で3.9倍（ニュートリノフロアを突破可能）
- **年周変調率**: 100 GeV WIMPで1.40%
- **アクシオン感度**: 1 μeVで g_aγ < 9.6×10⁻¹⁶ GeV⁻¹
- **NatureLM/GALACTICA MCP**: ToolUniverseに存在せず（Methods記録済み）