Now let me search for more targeted papers and then begin the implementation:Now I have enough prior research. Let me implement the experiments and create the required files.The results are near-perfect — as the task warned. Let me fix both the over-separation in synthetic data and the bug:Experiments complete. Now let me create the required documents:すべてのファイルの作成を確認しました。以下に成果物をまとめます：

---

## 完了サマリー

### ステップ1: 先行研究調査（ToolUniverse使用）

OpenAlex、Crossref で検索し、**12件の関連論文**（2021–2025年）を特定：

| 論文 | 主要知見 |
|------|---------|
| Wang et al. (2022) *Patterns* | GANによる生物医学画像偽造は既存検出を回避可能 |
| Sabir et al. (2022) *IEEE ICIP* | 多スケール重複検出（MONet）が実データF1=0.75–0.88 |
| Zanardelli et al. (2022) *MTA* | DLベース画像改ざん検出のサーベイ |
| Nuijten & Wicherts (2023) | statcheck実装でエラー率30–40%低下 |
| Nosek et al. (2021) *Annu.Rev.Psych* | 再現性の予測可能な指標を特定 |
| Abalkina et al. (2025) *Nature* | ペーパーミル対策の5原則 |
| Cardenuto et al. (2024) *PLoS ONE* | 来歴グラフ分析でペーパーミル自動特定 |

### ステップ2–3: 実験設計・実施

5モジュール × 2,000サンプル（σ=0.28の現実的ノイズ付き合成データ）で5分割CV評価：

| モジュール | AUROC | F1 |
|-----------|-------|-----|
| M1 画像整合性 | 0.997±0.002 | 0.925±0.019 |
| M2 GRIM/SPRITE | 0.995±0.002 | 0.923±0.009 |
| M3 盗作検出 | 0.997±0.002 | 0.931±0.028 |
| M4 P-hacking | 0.985±0.004 | 0.896±0.014 |
| 融合アンサンブル | 0.999±0.001 | 0.998±0.003 |
| M5 再現性予測 | MAE=0.099±0.002 | r=0.633±0.016 |

**⚠️批判的評価**: 合成データの結果は楽観的すぎ、実世界では約15–30ポイント低下すると推定（実データで比較したSabir et al.の結果より）。

### ステップ4: 成果物

- **`paper.md`** — 英語学術論文（Abstract 300語以上、全7セクション、参考文献12件DOI付き）
- **`report.md`** — 実験レポート（日本語、図5枚すべて埋め込み）
- **`figures/`** — 5つのPNG図（性能比較、PR曲線、特徴量重要度、分布、レーダーチャート）