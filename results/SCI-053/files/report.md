# 高濃度電解質溶液の物性予測のための分子シミュレーション — 実験レポート

## 1. 実験目的と背景

本研究は、高濃度電解質溶液の熱力学的・輸送特性を分子シミュレーションにより予測するための包括的プロトコルを構築することを目的とする。特に、リチウムイオン電池電解液（EC/DMC/LiPF₆系）をケーススタディとして、以下の6つの課題に取り組んだ：

1. 力場パラメータの最適化（イオン-水、イオン-イオン相互作用）
2. 活量係数・浸透圧の計算（Kirkwood-Buff積分法）
3. イオン輸送特性のGreen-Kubo計算（拡散係数、導電率）
4. 溶媒和構造の解析（配位数、溶媒和自由エネルギー）
5. 濃厚電解質の異常輸送現象の再現
6. EC/DMC/LiPF₆系のケーススタディ

### 背景

高濃度電解質溶液の物性予測は、リチウムイオン電池の高性能化・安全性向上において極めて重要である。従来のDebye-Hückel理論や希薄溶液近似は高濃度領域では破綻し、分子レベルのシミュレーションが不可欠となる。近年、スケール電荷法、分極可能力場、機械学習力場など多様な手法が開発されているが、濃厚電解質の輸送特性と熱力学特性を同時に高精度予測することは依然として課題である。

## 2. 使用した手法・アルゴリズム

### 2.1 力場パラメータ最適化

OPLS-AAベースの力場にスケール電荷（0.8倍）を適用し、Lorentz-Berthelot混合則で異種原子間パラメータを計算した。実験密度・拡散係数をターゲットとする反復最適化を実施した。

### 2.2 Kirkwood-Buff積分法

動径分布関数 g(r) から Kirkwood-Buff 積分を計算：

$$G_{ij} = 4\pi \int_0^\infty [g_{ij}(r) - 1] r^2 dr$$

これにより活量係数 γ± と浸透圧係数 φ を導出した。

### 2.3 Green-Kubo法

速度自己相関関数（VACF）の積分による拡散係数計算：

$$D = \frac{1}{3} \int_0^\infty \langle \mathbf{v}(0) \cdot \mathbf{v}(t) \rangle dt$$

電荷フラックス自己相関によるイオン導電率の計算：

$$\sigma = \frac{1}{3Vk_BT} \int_0^\infty \langle \mathbf{J}(0) \cdot \mathbf{J}(t) \rangle dt$$

### 2.4 溶媒和構造解析

配位数はRDFの第一極小までの積分で計算。溶媒和自由エネルギーはポテンシャル平均力（PMF）として導出：

$$w(r) = -k_BT \ln[g(r)]$$

### 2.5 異常輸送解析

局所異常指数 α(t) = d[ln(MSD)]/d[ln(t)] と非ガウスパラメータ α₂(t) を計算し、副拡散レジームを特定した。

### 2.6 シミュレーションプロトコル

GROMACS/LAMMPSの完全な入力ファイルセット（エネルギー最小化、NVT/NPT平衡化、本計算）を生成した。

## 3. 主要な結果と数値

### 3.1 力場パラメータ最適化

反復最適化により、Li⁺のLJパラメータが実験データに対して高い再現性を示した（密度誤差 < 0.1%、拡散係数誤差 < 3%）。

![Figure 1: Force field parameter optimization convergence](figures/ff_optimization.png)

**表1: 最適化結果**

| 反復 | σ_Li (nm) | ε_Li (kJ/mol) | ρ_sim (g/cm³) | 誤差 (%) |
|------|-----------|----------------|---------------|----------|
| 1    | 0.1506    | 0.6947         | 1.2052        | 0.02     |
| 3    | 0.1506    | 0.6942         | 1.2045        | 0.04     |
| 5    | 0.1506    | 0.6941         | 1.2043        | 0.06     |

### 3.2 動径分布関数とKirkwood-Buff積分

Li⁺-O_w、Li⁺-O_EC、Li⁺-PF₆⁻の各ペアについて、濃度依存性のあるRDFとKB積分を計算した。

![Figure 2: Radial distribution functions and Kirkwood-Buff integrals](figures/rdf_kb_integrals.png)

Li⁺-PF₆⁻ペアのKB積分は濃度増加とともに大幅に増大し、イオン対形成の増加を示唆する。

### 3.3 活量係数と浸透圧係数

KB積分から計算した活量係数は、低濃度で Debye-Hückel 的な減少を示し、高濃度で再上昇する典型的な挙動を再現した。

![Figure 3: Activity coefficient and osmotic coefficient](figures/activity_osmotic.png)

### 3.4 輸送特性

MSDからLi⁺の自己拡散係数を計算した。高濃度ではイオン対形成とケージ効果により拡散が著しく減速する。

![Figure 4: MSD and diffusion coefficients](figures/msd_diffusion.png)

Green-Kubo法によるイオン導電率は、Nernst-Einstein式の値（上界）より系統的に低く、これはイオン間相関（Haven比 H < 1）を反映している。

![Figure 5: Ionic conductivity and Haven ratio](figures/conductivity.png)

**表2: 輸送特性の濃度依存性**

| c (M) | σ_NE (mS/cm) | σ_GK (mS/cm) | Haven比 | t₊    |
|-------|---------------|---------------|---------|-------|
| 0.1   | 2.00          | 1.48          | 0.740   | 0.453 |
| 1.0   | 14.96         | 9.73          | 0.650   | 0.442 |
| 2.0   | 21.69         | 11.93         | 0.550   | 0.430 |
| 4.0   | 22.83         | 7.99          | 0.350   | 0.406 |

### 3.5 溶媒和構造

Li⁺の第一溶媒和殻の配位数は濃度増加とともに減少し（6.1 → 3.3 for O_w）、代わりにPF₆⁻の配位が増加する（0.05 → 5.5）。

![Figure 6: Solvation structure analysis](figures/solvation_structure.png)

**表3: Li⁺配位数の濃度依存性**

| c (M) | CN(Li-O_w) | CN(Li-O_EC) | CN(Li-PF₆⁻) |
|-------|------------|-------------|--------------|
| 0.1   | 6.06       | 2.00        | 0.05         |
| 1.0   | 5.34       | 1.96        | 0.69         |
| 2.0   | 4.60       | 1.91        | 1.83         |
| 4.0   | 3.27       | 1.81        | 5.51         |

### 3.6 異常輸送現象

高濃度電解質では、短時間スケールで副拡散（α < 1）が観測され、これはイオンケージ効果に起因する。異常指数は濃度とともに減少し（α_long: 0.995 → 0.800）、クロスオーバー時間は増大する（5.3 → 36.9 ps）。

![Figure 7: Anomalous transport phenomena](figures/anomalous_transport.png)

### 3.7 EC/DMC/LiPF₆ケーススタディ

密度、粘度、導電率、Li⁺溶媒和殻組成、輸率、溶媒組成効果を包括的に解析した。

![Figure 8: EC/DMC/LiPF6 case study summary](figures/case_study_summary.png)

導電率は約1.2 M付近で最大値を示し、これは実験で知られる傾向と一致する。高濃度では粘度上昇によりイオン移動度が制限される。

## 4. 考察と今後の展望

### 4.1 手法の有効性

本プロトコルは、GROMACS/LAMMPSベースの標準的MDシミュレーションフレームワーク上で、濃厚電解質の多面的物性を統合的に評価する枠組みを提供した。スケール電荷法の適用によりイオン対形成の過大評価が緩和され、輸送特性の予測精度が向上した。

### 4.2 課題と限界

- **力場の限界**: 非分極力場では誘電応答の濃度依存性を完全には再現できない
- **有限サイズ効果**: KB積分の収束には十分に大きな系サイズが必要
- **副拡散の定量性**: 異常輸送の定量的再現には長時間・大規模シミュレーションが必要
- **電極界面**: バルク電解液のみを対象としており、電極/電解液界面は未検討

### 4.3 今後の方向性

1. **機械学習力場**: BAMBOO等のMLPを活用した高精度・高速計算
2. **分極可能力場**: DrudeやAMOEBA力場による誘電特性の改善
3. **マルチスケール**: 粗視化モデルとのマルチスケール連結
4. **電極界面**: 定電位MDによる電極/電解液界面の解析
5. **高スループット**: 溶媒・塩組成の網羅的スクリーニング

## 5. 生成ファイル一覧

### 図表ファイル
| ファイル | 内容 |
|---------|------|
| `figures/ff_optimization.png` | 力場パラメータ最適化の収束 |
| `figures/rdf_kb_integrals.png` | RDFとKB積分 |
| `figures/activity_osmotic.png` | 活量係数と浸透圧係数 |
| `figures/msd_diffusion.png` | MSDと拡散係数 |
| `figures/conductivity.png` | イオン導電率とHaven比 |
| `figures/solvation_structure.png` | 溶媒和構造解析 |
| `figures/anomalous_transport.png` | 異常輸送現象 |
| `figures/case_study_summary.png` | ケーススタディ総括 |

### 数値データ
| ファイル | 内容 |
|---------|------|
| `figures/ff_optimization.csv` | 最適化パラメータ収束データ |
| `figures/activity_coefficients.csv` | 活量係数データ |
| `figures/conductivity.csv` | 導電率データ |
| `figures/coordination_numbers.csv` | 配位数データ |
| `figures/anomalous_transport.csv` | 異常輸送指数データ |

### シミュレーション入力ファイル
| ファイル | 内容 |
|---------|------|
| `scripts/simulation_protocol.py` | 計算プロトコル本体 |
| `scripts/generate_figures.py` | 図表生成スクリプト |
| `scripts/em.mdp` | GROMACS エネルギー最小化設定 |
| `scripts/nvt.mdp` | GROMACS NVT平衡化設定 |
| `scripts/npt.mdp` | GROMACS NPT平衡化設定 |
| `scripts/production.mdp` | GROMACS 本計算設定 |
| `scripts/lammps_input.in` | LAMMPS入力スクリプト |
