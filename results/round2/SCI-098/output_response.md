Now let me query NatureLM for scientific knowledge and start building the simulation framework in parallel.Now let me build the comprehensive simulation framework:All 6 figures generated. Now let me create the paper.md and report.md:すべてのファイルが正常に作成されました。作業内容をまとめます：

---

## 完了サマリー

### ステップ1: 先行研究調査
ToolUniverse MCP（Crossref, Semantic Scholar）を使用し、**11件の関連論文**（2020–2026年）を特定：
- LZ 2023, XENONnT 2023, PandaX-4T 2021（最新液体キセノン結果）
- DAMIC-M（スキッパーCCD軽いDM）
- DAMA/COSINE年周変動論文
- プリモーディアルBH（LVK O4, 2026）
- 暗黒光子 × JWST（2026）

### ステップ2: NatureLM科学的検証
`ask_naturelm`を3回試行 → **接続成功（HTTP 200）だが定量パラメータが返されず**（カテゴリラベルのみ）。すべてのパラメータを査読論文から取得（Methodsに記録済み）。

### ステップ3: シミュレーション実行
`dm_simulation.py`（約440行）で6種類の図を生成：
- **fig1**: 多標的感度曲線 + ニュートリノフロア
- **fig2**: 核反跳スペクトル（Xe/Ar/Ge/NaI）
- **fig3**: 年周変動解析（変調分率~3.3%、10⁵ kg·yr で5σ検出）
- **fig4**: アクシオン（ADMX KSVZ/DFSZ）+ 暗黒光子（ε~10⁻¹⁷）
- **fig5**: バックグラウンド低減（−97.6%）+ CF₄方向性検出（19°@10keV）
- **fig6**: ニュートリノフロア（Ar: Xeより6倍低い）

### ステップ4: 成果物
- **`paper.md`**: 30 KB, 410行 — Abstract 300語以上、7セクション、参照文献15件（DOI付き）
- **`report.md`**: 14 KB, 254行 — 全実験結果・考察（日本語）、全図埋め込み済み