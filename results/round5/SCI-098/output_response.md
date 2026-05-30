I have enough initial papers. Now let me implement the full simulation framework and find more references.Now let me implement the comprehensive simulation framework:All 6 figures generated. Now let me create the paper.md and report.md:全ファイルが正常に作成されました。以下に結果をまとめます。

---

## 完了報告

### ステップ1: 先行研究調査（Crossref MCP使用）

7件の論文を特定：
1. **LZ Collaboration (2022)** — LZ初期結果 σ_SI < 9.2×10⁻⁴⁸ cm²
2. **O'Hare (2021, PRL)** — ニュートリノフロアの新定義（DOI: 10.1103/physrevlett.127.251802）
3. **Nikolic et al. (2022, EPJC)** — 修正ニュートリノフロア（DOI: 10.1140/epjc/s10052-022-10534-3）
4. **Miuchi (2024, JAIS)** — 方向感度検出の課題
5. **Akerib (2024, NPB)** — LZ開発史・現状
6. **Casali (2025)** — NUCLEUS実験でのCEνNS
7. **Ko et al. (2023)** — NaI検出器とDAMAクロスチェック

### ステップ2–3: シミュレーション実施

- **6枚の図** を生成（感度曲線・年周変調・方向感度・非WIMP候補・背景低減・多ターゲット相補性）
- **重要な自己批判**: 絶対感度値はLZ実験値から5–6桁劣る（速度積分の単純化、スペクトル最尤推定の未実装による）。CV分散がゼロという問題も検出し、エネルギー積分感度の設計上の限界として明記

### ステップ4: 成果物

| ファイル | 内容 |
|---|---|
| **`paper.md`** | 368行の学術論文（Abstract 250語以上、全セクション、参考文献10件） |
| **`report.md`** | 日本語実験レポート（目的・手法・結果・考察・ファイル一覧） |
| `dm_simulation.py` | コアフレームワーク（Helm形状因子・MB速度分布・感度計算） |
| `run_analysis.py` | 全解析・図生成スクリプト |
| `figures/fig1–6.png` | 6枚の定量的図表 |