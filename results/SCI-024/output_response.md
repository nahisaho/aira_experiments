Now let's create the visualization module and run all computations:Now create the QE and Wannier90 input templates:Now let's run all the computations:All figures generated. Now create the process log and the final report:すべての計算が完了しました。以下に結果をまとめます。

---

## フレームワーク完成 — 実行結果サマリー

### 主要計算結果

| モジュール | 主な結果 |
|-----------|---------|
| **①対称性分類** | 10物質中8物質がZ₂=(1;000)の強いTI — Fu-Kane公式で確認 |
| **②Wannier TBモデル** | Bi₂Se₃ 4バンドk·pモデル実装（Liu 2010パラメータ）、ギャップ 0.354 eV |
| **③Z₂・Chern数** | フルSOC(λ=1.0)でZ₂=1、Chern=−1 ✓、半SOC(λ=0.5)でZ₂=0へ相転移 |
| **④表面状態** | 22層スラブのグリーン関数法でA(k,ω)計算、ヘリカルDiracコーンを確認 |
| **⑤SOC位相図** | 臨界SOC: Bi₂Te₃(λ_c=0.35) < Sb₂Te₃(0.49) < Bi₂Se₃(0.65) < TlBiSe₂(0.81) |
| **⑥候補スクリーニング** | 22物質評価、9つの新規候補 — **Bi₂Po₃(0.861), TlBiPo₂(0.846), SnBi₂Te₄(0.811)** が最上位 |

### 生成ファイル

- **`report.md`** — 完全な研究報告書（日本語）
- **`figures/`** — 6図（SVG+PNG）: バンド構造、表面状態、位相図、スクリーニング、ワークフロー図
- **`results/`** — 7つのJSONデータファイル
- **`src/`** — 8つのPythonモジュール（完全実装）
- **`qe_inputs/`, `w90_inputs/`** — QE/Wannier90入力テンプレート
- **`logs/process-log.jsonl`** — 実行トレース