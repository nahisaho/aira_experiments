I now have sufficient literature. Let me set up the experiments and create all deliverables.Excellent! All experiments succeeded. Now let me create the tests and write the deliverables.All 6 tests pass. Now let me write report.md and paper.md.すべての成果物が生成されました。結果をまとめます。

---

## Experiment Complete: 因果効果推定手法の体系的比較

### 主要科学的知見

1. **PSMが最小バイアス** — 5手法中、傾向スコアマッチングが最小絶対バイアス 0.0476 ± 0.0036（CV）を達成。5共変量・高重複率（99.8%）の設定でロジスティック回帰によるPS推定が有効だった
2. **全手法でポジティブバイアス** — 真のATE（−0.1469）に対し全手法が過小評価（IV最大バイアス 0.0742）。健康状態が悪い患者ほどスタチンを処方される「適応による交絡」が残存
3. **IV法の弱操作変数問題なし、DiDの平行トレンド仮定は違反** — 第1段階F統計量=27.62（閾値>10クリア）；DiD並行トレンド検定 t=8.74, p<0.001（仮定棄却）
4. **DMLは最低分散** — 交差検証バイアスのStd=0.0011で最も安定した推定（ただしバイアス自体は大きい）
5. **因果フォレストが異質性を可視化** — CATE標準偏差=0.036、高齢Q4群の治療効果がQ1群比約30%大（年齢・LDL依存の異質性を正確に捕捉）

### 最も重要な図
![ATE比較（95% CI）](figures/fig1_ate_comparison.png)

### 成果物一覧
- **ソースコード**: 4モジュール（828行）
- **テスト**: 6件全通過
- **report.md**: 1,061語（日本語）
- **paper.md**: 3,339語（英語、IMRaD形式、16文献・DOI付き）
- **図**: 5枚（fig1〜fig5）
- **結果**: CSV/JSON形式で `results/` に保存

### 限界
- 合成データ：実際のRWDが持つ欠測・測定誤差・時変交絡を再現していない
- N=2,000はDML・因果フォレストの真の優位性を発揮するには小規模
- DiD平行トレンド違反への感度分析（HonestDiD）は実施していない