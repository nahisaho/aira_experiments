I have enough literature. Now let me run the experiments.Now let me write and run the full experimental code:I notice a geometry issue - PennyLane's DHF uses Bohr, not Angstrom. Let me fix and rerun:Excellent results! Now let me create the two required output files:全ファイルが正常に作成されました。以下に結果をまとめます。

---

## 完了: VQEノイズ耐性研究

### ステップ1: 先行研究調査 (10件特定)

OpenAlex/Crossref経由で主要文献を特定:
- Cerezo et al. (2021) *Nature Reviews Physics* — VQAレビュー (2885引用)
- Bharti et al. (2022) *Rev. Mod. Phys.* — NISQアルゴリズム (1610引用)
- Cerezo et al. (2021) *Nat. Commun.* — バレンプラトー理論 (1031引用)
- Fedorov et al. (2022) *Materials Theory* — VQE ansatzレビュー
- LaRose et al. (2022) *Quantum* — Mitiqエラー軽減ライブラリ
- Gard et al. (2020) *npj Quantum Inf.* — H₂/LiH VQE対称性回路

### ステップ2-3: 実験結果 (PennyLane 0.45, Python 3.11)

| 系 | Ansatz | エネルギー (Ha) | 誤差 (mHa) | 化学精度 |
|---|---|---|---|---|
| H₂ (4q) | HEA (24p) | −1.137272 ± 0.000008 | 0.004 | ✓ |
| H₂ (4q) | UCCSD (1p) | −1.137276 ± 0.000000 | <0.001 | ✓ |
| LiH (6q) | HEA (36p) | −7.862016 ± 0.000011 | 1.047 | ✓ |
| LiH (6q) | UCCSD (6p) | −7.863063 ± 0.000000 | <0.001 | ✓ |

**エラー軽減** (p=0.01): 未軽減 18.46 mHa → ZNE 0.67 mHa (28×) → CDR 0.46 mHa (40×)

**⚠️ 自己批判**: 全結果はノイズレスシミュレーターに基づく。実機では±10×程度の性能劣化が予想される。

### 成果物
- **`paper.md`** — 英語学術論文 (Abstract 250語、References 10件、定量的結果表付き)
- **`report.md`** — 日本語実験レポート (図5枚埋め込み済み)
- **`figures/`** — 5つのグラフ (収束・エラー軽減・バレンプラトー・解離曲線・精度比較)