ToolUniverse MCPツールが使用可能です。並行して先行研究を検索します。先行研究が見つかりました。追加の検索を行います。十分な先行研究が収集できました。実験実装を開始します。環境確認完了。実験コードを作成します。実験完了。可視化を生成します。図表の生成完了。paper.md と report.md を作成します。完了しました。以下が成果物のまとめです。

---

## 実験完了サマリー

### ステップ1: 先行研究調査（MCP）
- **OpenAlex** で10件の主要論文を発見（Semantic Scholar は API制限429で失敗→代替手段として記録）
- 主要論文: AiZynthFinder, Graph2SMILES, Graph2Edits, RAscore, RetroXpert など

### ステップ2-3: 実験結果

| 指標 | 結果 |
|------|------|
| Template-Based Top-1 | **0.250 ± 0.098** (5-fold bootstrap) |
| Template-Free Top-1 | 0.200 ± 0.098 |
| TB vs TF 速度比 | TB が **16.6× 高速** (1.0 ms vs 16.6 ms) |
| SA Score 相関 | r = **0.968** (改良版 vs オリジナル) |
| MCTS/A* 成功率 | **4/5 (80%)** の医薬品分子で経路発見 |

### ステップ4: 成果物

- **`paper.md`** — 学術論文（Abstract 200語以上、10文献のDOI付き参考文献、7図）
- **`report.md`** — 実験全結果・考察（日本語、7図埋め込み）
- **`figures/`** — 7つのPNG図
- **`data/`** — CSV/JSON形式の実験データ