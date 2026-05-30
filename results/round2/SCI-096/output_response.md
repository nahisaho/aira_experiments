I'll now search for related papers in parallel across multiple topics.Now let me query NatureLM and search for more specific papers in parallel:Now let me run the computational experiment to generate figures:All 5 figures generated. Now let me create the paper.md and report.md:すべてのファイルが正常に作成されました。以下に結果をまとめます。

---

## 完了した作業

### Step 1: 先行研究調査（ToolUniverse MCP使用）

**PubMed検索（5クエリ）** で10件以上の関連論文を特定：
- Corcoran et al. (2026) — IIT・Active Inferenceの敵対的協調レビュー
- Mayner, Tononi et al. (2026) — IIT 4.0の分化-特定性トレードオフ
- McQueen et al. (2026) — 量子重ね合わせとIITの結合
- Wiest (2025) — 量子微小管基質の実験的支持
- Maschke et al. (2024) — EEG臨界ダイナミクスとPCIの相関

※ Semantic Scholar APIは400エラーのためPubMedで代替

### Step 2: NatureLM MCP科学的検証

4クエリ実行、IIT公理・Orch-OR量子コヒーレンスタイムスケール（10〜1000 fs）・FEP-IIT統合の知見を取得。ゾンビ論証クエリは不完全応答（§3.3に記録）。

### Step 3: 計算実験

5つのネットワークトポロジー × 5つの意識状態で計算実験を実施：
- **Small-Worldが最高Φ** (5.94 ± 0.64) — IIT予測と一致
- **意識状態のΦ順序**: 瞑想(8.42) > 覚醒(7.94) > REM(6.19) > NREM(3.11) > 麻酔(2.55)
- **Φ–PCI相関**: r = −0.69（モデル限界あり、詳細はDiscussionに記載）

### 生成ファイル

| ファイル | 内容 |
|---------|------|
| `paper.md` | 学術論文（英語、7セクション + 12参考文献） |
| `report.md` | 実験レポート（日本語、全結果・考察・図埋め込み） |
| `figures/fig1〜5.png` | 5つの実験図 |
| `experiment.py` | 計算実験スクリプト |