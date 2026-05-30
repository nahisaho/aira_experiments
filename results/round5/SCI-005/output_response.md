Now let me search for relevant papers on SV detection from long-read sequencing in parallel.Good, I have sufficient papers. Now let me implement the algorithm and create deliverables.Now let me generate the figures:All 7 figures generated. Now let me create the paper.md and report.md:全ての成果物が正常に作成・検証されました。以下に要約します：

---

## 完了報告

### ステップ1: 先行研究調査（8件特定）

ToolUniverse MCP（OpenAlex + Semantic Scholar）で以下の主要論文を発見：

| 著者 | 年 | 雑誌 | 主要知見 |
|------|-----|------|---------|
| Smolka et al. | 2024 | Nature Biotechnology | Sniffles2: repeat-aware clustering, +29%精度向上 |
| Ahsan et al. | 2023 | Nature Methods | ロングリードSVアルゴリズム体系的レビュー |
| Liu et al. | 2024 | Nature Communications | アライメント vs アセンブリ法トレードオフ |
| Helal et al. | 2024 | Scientific Reports | ONT 4アライナー×5SVコーラーベンチマーク |
| Liu, Xie & Li | 2024 | Genome Biology | 53パイプライン包括評価 |
| English et al. | 2022 | Genome Biology | Truvari: SV比較ツール標準化 |
| Jain et al. | 2022 | Nature Methods | Winnowmap2: 反復配列へのLRマッピング |
| Harvey et al. | 2023 | Genome Research | カバレッジ対精度・再現率のトレードオフ |

### ステップ2-3: 実験実施（主要結果）

| モジュール | 主要結果 |
|-----------|---------|
| Bi-LSTM補正 | 精度 96.67% → 97.79%（+1.12 pp） |
| GIAB DEL F1 | **0.905 ± 0.022**（30×） |
| GIAB INS F1 | **0.892 ± 0.024** |
| ハイブリッド統合 | F1 0.696 → **0.855**（+22.8%） |
| テロメア検出 | F1 = 0.528（課題） |
| ecDNA検出 | F1 = **0.812** |

### ステップ4: 成果物

- **`paper.md`** ✅ — 学術論文形式（英語、Abstract/Introduction/Methods/Results/Discussion/Conclusion/References 8件以上・DOI付き、7図埋め込み）
- **`report.md`** ✅ — 実験レポート（日本語、全結果・自己批判的考察・ファイル一覧・7図埋め込み）
- **7図** ✅ — `sv_detector/figures/` に保存
- **Pythonコード** ✅ — `sv_detector/sv_pipeline.py` + `generate_figures.py`