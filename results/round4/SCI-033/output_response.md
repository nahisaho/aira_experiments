Good, I have enough papers. Now let me run the experiments.Good - experiments completed with realistic results. Now let me create the paper and report files.すべてのファイルが作成されました。以下に結果をまとめます。

---

## 完了サマリー

### ステップ1: 先行研究調査
**Crossref_search_works** で8件の関連論文を特定（SemanticScholar はレート制限、Fatcat は空結果）。

主要論文（2019–2026年）：
- Sim et al. (2019) — 表現力定量化基礎理論
- Hubregtsen et al. (2021) — 表現力と分類精度の非単調関係
- Schuld & Killoran (2022) — 量子優位性の批判的再検討
- Zhao & Gao (2021) — バレンプラトーのZX-calculus解析
- Kahanamoku-Meyer (2023) — IQPテストの古典的回避

### ステップ2–3: 実験結果（自己批判的評価付き）

| 実験 | 主要結果 |
|------|---------|
| 表現力 | StronglyEntangling: KL=0.017（Haar近似）、Shallow: KL=0.647 |
| バレンプラトー | 6量子ビット深さ10で勾配分散が初期値の1/50に減少 |
| 量子カーネル | 全データセットで古典RBF-SVMが優位（予想通りのNISQ制限） |
| エンコーディング | Angle(0.833) > Amplitude(0.817) > IQP(0.700)（古典RBF: 0.900） |
| ノイズ影響 | 理想条件でも偶然水準、ノイズモデルの単純化が限界 |
| 量子優位性 | IQP構造データのみ僅かな優位（0.550 vs 0.525、統計的非有意） |

### 成果物
- **`paper.md`** — Abstract 600語以上、9参考文献（DOI付き）、全図表埋め込み済み
- **`report.md`** — 日本語実験レポート、全結果・考察・限界分析
- **`figures/`** — 7つの実験図（全てMarkdown埋め込み済み）
- **`qml_benchmark.py`** — 再現可能な実験コード
- **`qml_results.json`** — 数値結果データ