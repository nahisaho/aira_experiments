# Quantum Machine Learning Benchmark — Experiment Report

**Date**: 2026-05-31  
**Research Theme**: Quantum Machine Learning Model Expressibility and Classical Comparison Framework  
**Notebook**: `qml_benchmark.ipynb`  
**Environment**: PennyLane 0.45.0, Python 3.11.2, `default.qubit` / `default.mixed` backends

---

## 1. 実験目的と背景

### 目的
本実験は、量子機械学習（QML）モデルの表現力（expressibility）・エンタングルメント能力・量子カーネル性能・バレンプラトー問題・ノイズ耐性を系統的に定量化し、古典的機械学習手法と比較するベンチマークフレームワークを構築することを目的とする。

### 背景
量子コンピューティングと機械学習の融合であるQMLは、指数的状態空間を活用した特徴写像により古典手法を超える可能性を秘めているが、以下の課題が実用化を妨げている：

1. **バレンプラトー（Barren Plateaus）**: 深い量子回路では勾配が指数的に消失し、学習が困難
2. **データエンコーディング依存性**: angle/amplitude/IQPエンコーディングによりカーネル構造が変化
3. **NISQハードウェアノイズ**: 脱分極ノイズが回路忠実度を劣化させる
4. **量子優位性の条件**: どのデータセット・タスクで量子が古典を超えるかが不明

---

## 2. 使用手法・アルゴリズムの概要

### 2.1 PQCアンザッツ比較（4種類）

| アンザッツ | 構成 | パラメータ数 |
|----------|------|------------|
| **Shallow** | 1層Ryゲート、エンタングルなし | 4 |
| **HWE-2L** | 2層(Ry+Rz)+線形CNOT | 16 |
| **SE-2L** | 2層Rot(Rz-Ry-Rz)+円形CNOT | 24 |
| **Deep-4L** | 4層Rot+全結合CZ | 48 |

### 2.2 表現力測定（KL divergence from Haar）

Sim et al. [2019]に基づき、PQCの出力状態ペアの忠実度分布とHaar乱数ユニタリの理論分布のKLダイバージェンスで表現力を定量化。

$$\text{Expr} = D_{KL}(\hat{P}_{\mathcal{U}}(F) \| P_{\text{Haar}}(F))$$

- $P_{\text{Haar}}(F) = (d-1)(1-F)^{d-2}$, $d = 2^n$
- 低い値ほど高表現力（Haar乱数に近い）

### 2.3 エンタングルメント能力（Meyer-Wallach測度）

$$Q = \frac{4}{n} \sum_{j=0}^{n-1} \left(1 - \text{Tr}(\rho_j^2)\right)$$

$Q \in [0,1]$（0:積状態、1:最大エンタングル）

### 2.4 量子カーネル法

3種のエンコーディング戦略：
- **Angle encoding**: $k(x_1,x_2) = \prod_k \cos^2((x_k^{(1)} - x_k^{(2)})/2)$
- **Amplitude encoding**: 正規化ベクトル内積の2乗
- **IQP encoding**: ZZ相互作用項を含む指数型カーネル

評価指標：AUROC（5分割交差検証）＋カーネルターゲットアライメント（KTA）

### 2.5 バレンプラトー解析

パラメータシフト法で勾配分散を測定：
$$\frac{\partial C}{\partial\theta_0} = \frac{C(\theta_0 + \pi/2) - C(\theta_0 - \pi/2)}{2}$$

グローバルコスト（$Z_0 \otimes Z_1$）とローカルコスト（$Z_0$のみ）を比較。

### 2.6 ノイズシミュレーション

IBMデバイスを模擬した脱分極ノイズ：

$$\mathcal{D}_p(\rho) = (1-p)\rho + \frac{p}{3}(X\rho X + Y\rho Y + Z\rho Z)$$

忠実度 $F = \langle\psi_\text{ideal}|\rho_\text{noisy}|\psi_\text{ideal}\rangle$ を評価。

---

## 3. ToolUniverse MCP 使用状況

### Semantic Scholar 検索（成功）
以下の論文を収集：
- Sim et al. 2019 (1066 citations)
- Havlíček et al. 2019 (2551 citations)
- Cerezo et al. 2020/2021 (134 citations)
- Pesah et al. 2021 (377 citations)
- Bowles et al. 2024 (161 citations)

### NatureLM MCP（接続失敗）
- **試行ツール名**: NatureLM MCP (`ask_naturelm`)
- **エラー内容**: ToolUniverseで0件一致（ツール未登録）
- **代替手段**: PennyLane シミュレーションによる定量的予測

### GALACTICA MCP（接続失敗）
- **試行ツール名**: GALACTICA MCP (`scientific_qa`, `predict_citations`)
- **エラー内容**: ToolUniverseで0件一致（ツール未登録）
- **代替手段**: Semantic Scholar文献検索 + 理論的検証

---

## 4. 主要な結果と数値

### 4.1 PQC表現力とエンタングルメント

| アンザッツ | KL発散（↓良） | Q（↑良） | パラメータ数 |
|----------|-------------|--------|----------|
| Shallow | **0.6231** | 0.0000±0.0000 | 4 |
| HWE-2L | 0.0953 | 0.7035±0.1556 | 16 |
| SE-2L | 0.0578 | **0.8509±0.0819** | 24 |
| Deep-4L | **0.0348** | 0.7521±0.1203 | 48 |

[cell:2], [cell:3]

**注目点**: SE-2LはDeep-4Lよりパラメータ数が少ないにもかかわらず、高いエンタングルメントを示す（Q=0.851 vs 0.752）。

![PQC Expressibility and Entanglement](figures/fig1_expressibility_entanglement.png)

### 4.2 量子カーネル vs 古典SVM（AUROC, 5分割CV）

| データセット | Classical RBF | Classical Linear | Q-Angle | Q-IQP |
|------------|---------------|-----------------|---------|-------|
| Linear | 0.972±0.056 | **0.997±0.006** | 0.972±0.056 | 0.966±0.069 |
| Moons | 0.947±0.040 | 0.959±0.025 | 0.938±0.046 | **0.972±0.025** |
| Circles | **0.997±0.006** | 0.409±0.135 | 0.994±0.008 | 0.991±0.013 |

[cell:5]

![Quantum vs Classical Kernel Comparison](figures/fig2_kernel_comparison.png)

### 4.3 データエンコーディング戦略の影響

| データセット | Angle AUROC | Amplitude AUROC | IQP AUROC | IQP KTA |
|------------|-------------|-----------------|-----------|---------|
| Linear-Sep | 0.942 | 0.729 | 0.973 | 0.2764 |
| Non-linear | 0.920 | 0.942 | **0.978** | 0.3166 |
| High-noise | 0.862 | 0.902 | 0.871 | 0.1927 |
| Random | 0.606 | 0.606 | 0.661 | 0.0779 |

[cell:7]

**発見**: IQPエンコーディングが最も高いKTAを示し、非線形データで最高性能（0.978 AUROC）。Amplitude encodingは線形データで低性能（0.729）。

### 4.4 量子優位性がある/ないデータセット

| データセット | RBF AUROC | Q-IQP AUROC | **量子優位性** |
|------------|-----------|-------------|-------------|
| Linear | 0.972 | 0.966 | −0.006 |
| **Quadratic** | 0.721 | **0.900** | **+0.179** ✓ |
| XOR | 0.835 | 0.705 | −0.130 ✗ |
| Checkerboard | 0.609 | 0.447 | −0.162 ✗ |
| High-D Interaction | 0.893 | 0.862 | −0.032 |

[cell:12]

**量子優位性が現れた条件**: 二次境界（放射状分離）を持つデータセットで、IQPカーネルの暗黙的ガウス積構造と一致する場合。

![Quantum Advantage by Dataset](figures/fig5_quantum_advantage.png)

### 4.5 バレンプラトー解析

| 量子ビット数 | Var(Global) | Var(Local) |
|-----------|-------------|------------|
| 2 | 1.279×10⁻¹ | 3.151×10⁻¹ |
| 4 | 3.225×10⁻² | 4.319×10⁻¹ |
| 6 | 6.954×10⁻³ | 4.279×10⁻¹ |
| 8 | 2.672×10⁻³ | 3.049×10⁻¹ |

[cell:8]

**指数減衰フィット**: Var(Global) ∝ exp(−**0.598** × n)，R² = **0.871**

Local costはほぼ一定（減衰率 α = 0.004, R² = 0.004）→ バレンプラトーを効果的に回避。

![Barren Plateau Analysis](figures/fig3_barren_plateau.png)

### 4.6 IBM Quantumノイズ下での実用性

| ノイズ率 p | 回路忠実度 | Std |
|----------|----------|-----|
| 0.000 (理想) | 1.0000 | 0.0000 |
| 0.001 | 0.9835 | 0.0010 |
| **0.010 (IBM典型)** | **0.8460** | **0.0072** |
| 0.050 | 0.4475 | 0.0190 |
| 0.100 | 0.2195 | 0.0173 |

[cell:10]

| 量子ビット数 | 忠実度（p=0.01） |
|-----------|----------------|
| 2 | 0.941±0.004 |
| 4 | 0.847±0.008 |
| 6 | 0.759±0.008 |
| 8 | 0.682±0.007 |

[cell:11]

**結論**: IBM現行ハードウェア（p≈0.001–0.01）では4量子ビット回路は実用的だが、8量子ビット超では忠実度が急速に劣化。

![Noise Analysis](figures/fig4_noise_analysis.png)

---

## 5. 総合ダッシュボード

![Comprehensive Overview](figures/fig0_overview.png)

---

## 6. 考察と今後の展望

### 6.1 主要な発見

1. **SE-2L が最適バランス**: 24パラメータで KL=0.058, Q=0.851 を達成。効率的なアンザッツ選択の指針。
2. **IQP優位性の条件**: 放射状/非線形構造を持つデータセットでのみ量子優位性が現れる（+0.179 AUROC）。XOR・チェッカーボードでは古典に劣る。
3. **バレンプラトー確認**: グローバルコストの勾配分散は exp(-0.598n) で指数的減少（R²=0.871）。ローカルコストへの切り替えが本質的解決策。
4. **ノイズ限界**: IBM典型ノイズ(p=0.01)で4量子ビット84.6%, 8量子ビット68.2%の忠実度。12量子ビット超で実用性が危うい。

### 6.2 自己批判的評価

- **合成データへの依存**: 実世界データでの汎化は未検証
- **カーネル近似**: 完全量子回路シミュレーションではなく解析的近似を使用
- **小サンプル問題**: n=80での5分割CVは高分散（最大std=0.20）
- **ノイズモデル簡略化**: コヒーレントエラーやクロストーク未考慮

### 6.3 今後の展望

1. **エラー緩和技術**: Zero-noise extrapolation, 確率的エラーキャンセル
2. **実機実装**: IBM Eagle (127 qubits) での実験的検証
3. **タスク適応型アンザッツ設計**: データ幾何学に基づく回路設計
4. **大規模ベンチマーク**: 10以上の量子ビット、実世界データセット適用
5. **量子優位性の理論的特徴付け**: どのデータ分布がIQPカーネルに有利かの数学的解析

---

## 7. 生成ファイル一覧

### 図表
| ファイル名 | 内容 |
|----------|------|
| `figures/fig0_overview.png` | 総合ダッシュボード（6パネル） |
| `figures/fig1_expressibility_entanglement.png` | PQC表現力とエンタングルメント |
| `figures/fig2_kernel_comparison.png` | 量子 vs 古典SVM比較 |
| `figures/fig3_barren_plateau.png` | バレンプラトー解析 |
| `figures/fig4_noise_analysis.png` | ノイズ下での忠実度 |
| `figures/fig5_quantum_advantage.png` | データセット別量子優位性 |

### データ
| ファイル名 | 内容 |
|----------|------|
| `data/raw/pip_freeze.txt` | Python依存パッケージ全一覧 |

### コード
| ファイル名 | 内容 |
|----------|------|
| `qml_benchmark.ipynb` | Jupyter実験ノートブック |
| `paper.md` | 学術論文形式文書 |
| `report.md` | 本実験レポート |

---

## 8. 再現性情報

| 項目 | 値 |
|-----|---|
| ランダムシード | 42 (`np.random.seed(42)`) |
| Pythonバージョン | 3.11.2 |
| PennyLane | 0.45.0 |
| NumPy | 2.4.6 |
| SciPy | 1.17.1 |
| scikit-learn | 1.8.0 |
| matplotlib | 3.10.9 |
| pandas | 3.0.3 |

---

*Report generated: 2026-05-31 | Framework: PennyLane QML Benchmark Suite*
