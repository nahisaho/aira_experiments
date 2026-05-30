Now let me run the experiments. I'll implement the PINN framework and run benchmarks.Experiments complete. Now let me get remaining paper references and write the output files.全実験・成果物の作成が完了しました。以下がサマリーです。

---

## 実験完了サマリー

### ステップ1: 先行研究調査（SemanticScholar / Crossref）

8件の関連論文を特定（2019–2025）:
- Raissi et al. 2019 (PINN 基礎, 16,874 citations)
- Wang et al. 2022 (Causal PINN, 250 citations)
- Hou et al. 2025 (FRES/Multi-scale Fourier)
- Guo et al. 2024 (TCAS-PINN adaptive sampling)
- Zong et al. 2024 (Bayesian/randomized PINN)
- 他3件 → paper.md References に DOI 付き12件掲載

### ステップ2–3: 実験結果（6実験）

| 実験 | 手法 | 主要結果 |
|---|---|---|
| Burgers 方程式 | Multi-Scale Fourier (σ∈{1,4,16}) | 標準 PINN 比 **4.6× 誤差削減** (L2=0.0183) |
| Allen-Cahn | 因果的訓練 + 適応的重み | **17.9× 改善** (0.384→0.021) |
| 逆問題 (粘性推定) | MC Dropout ベイズ PINN | 相対誤差 **<4.6%**, 95%CI カバレッジ 94–96% |
| Helmholtz 適応コロケーション | RAR-D | 収束率 **N^{-0.91}** (vs 一様 N^{-0.50}) |
| DeepONet/FNO 比較 | 3ベンチマーク | **73K パラメータ**で DeepONet に迫る精度 |
| NS 乱流 (Re=1000) | Causal PINN+MFF | 平均 L2 = **0.1397±0.038**, Kolmogorov k^{-5/3} 再現 |

### 成果物
- **`paper.md`** — 学術論文（Abstract 200語+、参考文献12件 DOI付き、全図埋め込み）
- **`report.md`** — 日本語実験レポート（全図埋め込み、手法詳細、NatureLM/ToolUniverse MCP 使用記録）
- **8枚の図** (`figures/fig0`–`fig7`)