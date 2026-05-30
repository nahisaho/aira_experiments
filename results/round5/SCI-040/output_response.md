Now let me run parallel literature searches and start the experiment implementation:I have enough literature. Now let me implement the full experiment:Excellent! Now let me generate the figures:All 6 figures generated. Now let me create the paper.md and report.md:すべてのファイルが正常に生成されました。以下に結果をまとめます：

---

## 実験完了サマリー

### ステップ1: 先行研究調査（Crossref/Semantic Scholar使用）
6件の関連論文を特定：
- **Yamasaki et al. (2022)** — 桜島/姶良カルデラのMaxwell粘弾性モデル
- **Yamasaki et al. (2023)** — 低粘性ゾーンの影響
- **Bonaccorso & Aloisi (2021)** — エトナ火山40年の変形ソースモデリング
- **Xue et al. (2020)** — GPS+InSAR+カルマンフィルタ（Okmok）
- **Nishiyama (2022)** — 円錐形火山体への圧力源変形解
- **Saito et al. (2018)** — 阿蘇中岳の岩石学的特性

### ステップ2〜3: 実験結果（主要数値）

| 指標 | 桜島（Mogi） | 桜島（楕円体） | 阿蘇（Mogi） | 阿蘇（楕円体） |
|------|------------|-------------|------------|-------------|
| CV RMSE | 9.49±1.73 mm | 7.21±1.25 mm | 29.52±5.76 mm | **8.78±2.43 mm** |
| R̂ (max) | **1.011** | 3.70 | **1.010** | 7.72 |
| 深さ誤差 | +1.1% | — | −1.5% | — |
| ΔV誤差 | +0.7% | — | −2.9% | — |

**⚠️ 自己批判的評価**：弾性半空間仮定・単一ソース仮定・楽観的合成ノイズのため、実世界適用時のRMSEは報告値より大幅に大きくなると予想。楕円体モデルはR̂ > 3で収束不十分——HMCサンプラー必須。

### ステップ4: 成果物
- **`paper.md`** — 英語学術論文（Abstract 300語超、7セクション、図6枚埋め込み、文献12件）
- **`report.md`** — 日本語実験レポート（図6枚埋め込み）
- **`figures/`** — fig1〜fig6（PNG 6枚）
- **`src/volcanic_inversion.py`** — PyMC互換ベイズインバージョンフレ��ムワーク