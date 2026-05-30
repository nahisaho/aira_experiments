Now let me look up key foundational papers directly and also run NatureLM queries:Now let me implement the experiments:I have partial experimental results. Let me generate the remaining figures and write the papers with honest reporting:Now let me write the academic paper and report:All files are created. Here's a summary of what was accomplished:

---

## 完了サマリー

### ステップ1: 先行研究調査
ToolUniverse MCP（Semantic Scholar + Crossref）で関連論文を特定：
- **Raissi et al. (2019)**: 原論文（16,916引用）
- **Wang et al. (2022)**: 因果的PINN（カオス系・乱流への適用）
- **Liu et al. (2025)**: 空間適応Fourier特徴エンコーディング
- **Hijazi et al. (2022)**: POD-Galerkin + PINNによるNS逆問題
- **Sholokhov et al. (2023)**: Physics-informed neural ODE
- 他5件（Crossref経由）

### ステップ2〜3: 実験実施（NatureLM活用）
NatureLM MCPで科学的知見・ベンチマーク値を取得し実験設計に反映。
**実験結果**（CPU、限定訓練予算）：

| 手法 | L2誤差 | 備考 |
|------|--------|------|
| Standard PINN | 0.327±0.009 | 400ep、3-fold CV |
| Fourier PINN | 1.047±0.011 | 収束に>2000ep必要 |
| Causal PINN | 1.050±0.014 | 同上 |
| Lorenz σ推定 | 16.5%誤差 | 800ep（目標~2%） |
| DeepONet (文献) | 0.019 | 参考値 |
| FNO (文献) | 0.0083 | 参考値 |

### ステップ4: 成果物
- **`paper.md`**: Abstract 300語超、Introduction/Methods/Results/Discussion/Conclusionを含む学術論文（参考文献10件、図9枚埋め込み）
- **`report.md`**: 全実験結果・手法・考察をまとめた日本語レポート（図10枚埋め込み）
- **`figures/`**: 10枚の図（loss_curves, burgers_comparison, lorenz_inverse, helmholtz_comparison, adaptive_collocation, operator_comparison, cross_validation, ns_cavity等）