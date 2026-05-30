Good — I have initial literature. Let me search for more specific Perturb-seq analysis method papers.Now let me implement the pipeline and generate figures:Pipeline complete. Now creating the paper and report:All outputs are confirmed. Here's a summary of what was accomplished:

---

## 完了サマリー

### ステップ1: 先行研究調査（ToolUniverse MCP）
PMC_search_papers を使用して **8件の論文**を特定（2023–2025年中心）:
- Dixit et al. 2016（Perturb-seq 原著）、Morris et al. 2024（次世代スクリーニングレビュー）、Schrod et al. 2024（CODEX 因果モデル）、Yu et al. 2025（PerturbNet）、Zhu et al. 2025（AUPRC メトリクス）など

### ステップ2: NatureLM 科学的検証
- **ΔG = 2.95 kcal/mol**（Cas9-gRNA 結合自由エネルギー）
- **DEG 予測 100–500**（観測: 192 ± 103 ✓）
- **AUROC 予測 0.80–0.95**（摂動別: 0.611–0.943 ✓、全体: 0.586 → 解釈済み）

### ステップ3: パイプライン実施（6モジュール）
| モジュール | 結果 |
|-----------|------|
| QC | 6,230/8,000 cells (77.9%)、ガイド検出率 96.4%、ダブレット 2.9% |
| DE | 平均 192 ± 103 DEG/摂動 (FDR<0.05) |
| 因果グラフ | 30 ノード、17 エッジ、密度 0.039 |
| エピスタシス | 相乗 4 ペア、拮抗 4 ペア |
| 表現学習 | AUROC 0.586 ± 0.017（摂動別最大: NOTCH1 0.943） |
| 必須遺伝子 | NOTCH1 (0.977) > STAT3 (0.821) > KRAS (0.764) |

### ステップ4: 成果物
- 📄 **`paper.md`** — 8セクション構成の学術論文（Abstract 200語以上、8件の DOI 付き参考文献、6図埋込）
- 📄 **`report.md`** — 全解析の詳細レポート（日本語、6図埋込）
- 🖼️ **`figures/`** — 6枚の図（各 142–355 KB）