MCP tools are available. Running parallel literature searches now.Good results from ArXiv! Let me search for more specific papers.I have sufficient literature. Now let me implement the quantum network simulation.All experiments completed successfully. Now let me write the tests, then report.md and paper.md.All 28 tests pass. Now let me create the report.md and paper.md.All 28 tests pass, workspace is clean. Here is the summary:

---

## 実験完了: 量子インターネットQKDネットワーク設計シミュレーション

### 主要な科学的知見

1. **BB84有限鍵率 (5-fold CV)**: 距離50 km で **2.923×10⁻² ± 1.221×10⁻³ bits/pulse** (N=10¹⁰, CV=4.2%)、100 km で **2.951×10⁻³ ± 5.853×10⁻⁵ bits/pulse**
2. **量子リピータ最適点**: NVセンターメモリで500 km/4セグメントが最適（29.5 Hz、忠実度0.990、ノードあたり最低2量子ビット）
3. **DEJMPS蒸留収束**: F₀=0.70 から4〜6ラウンドでF=0.95+達成（16〜64ペア消費）；ゲート忠実度0.98では上限~0.92
4. **東京QKDネットワーク**: 最大忠実度ルーティングで全21ペア中、忠実度0.797〜0.940を達成
5. **チャネルロスが主要ボトルネック**: 標準SMF 50 km で光子伝送成功率~10%

### 最重要図

![BB84/E91 鍵率対距離](figures/fig1_key_rate_vs_distance.png)

### 成果物
- **ソースコード**: 4モジュール（~1,100行）
- **テスト**: 28件（全件パス）
- **図**: 9枚（`figures/`）
- **report.md**: 1,041語
- **paper.md**: 3,964語（18件参考文献、全DOI付き、60%以上が2020年以降）

### 制限事項
- Werner状態近似（精密密度行列追跡なし）
- 東京ネットワークトポロジーは2011年論文の近似
- マルチユーザ競合シナリオは未考慮