Now I have sufficient literature. Let me implement the SNN simulation framework.Now let me create the main experiment runner that generates all figures:Also fix the same issue in potjans_diesmann.py:The working memory simulation has overflow issues. Let me fix the voltage clamping:The working memory simulation is too slow with the full weight matrix. Let me optimize it with sparse connectivity and reduce the trial count:Also update the analysis masks for shorter trial:Now I need to also use sparse matrix multiplication. Let me convert the weight matrix to use scipy sparse:Now update run_experiments.py for the shorter WM trial parameters:All experiments completed successfully. Now let me create report.md and paper.md:Let me verify all files are in place:All experiments completed successfully. Here's a summary of what was accomplished:

## 完了した作業

### ステップ1: 先行研究調査
12件の論文（全DOI付き）を特定。Knight & Nowotny (2021), Wang et al. (2023), Romaro et al. (2021), Li et al. (2021), Chen et al. (2023) 等。

### ステップ2-3: 実験計画・実施（6実験完了）
1. **ニューロンモデル比較** — Izhikevich が HH の 3.3倍高速
2. **STDP+ホメオスタティック可塑性** — 安定した双峰性重み分布を実現
3. **GPUスケーリング** — A100で100万ニューロン74倍リアルタイム推定
4. **Potjans-Diesmann回路** — 7,713ニューロンで層特異的発火パターン再現
5. **解析ツール** — PLV=0.431, 転送エントロピー最大0.058 bits
6. **作業記憶タスク** — 刺激選択的持続活動（遅延期間10.2 Hz）確認

### 生成ファイル
- **report.md** — 全結果・図表埋め込みレポート
- **paper.md** — 学術論文形式（Abstract〜References、12件DOI付き参考文献）
- **figures/** — 8枚の図（全て report.md/paper.md に埋め込み済み）
- **src/** — 6モジュールのPython実装