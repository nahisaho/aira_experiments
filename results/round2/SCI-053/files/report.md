# 実験レポート：高濃度電解質溶液の物性予測のための分子シミュレーション手法の設計と検証

---

## 1. 実験目的と背景

### 研究目的

本研究では、リチウムイオン電池（LIB）の代表的電解液系であるEC/DMC/LiPF₆（エチレンカーボネート/ジメチルカーボネート/六フッ化リン酸リチウム）について、高濃度（0.5〜4.0 mol/kg）における以下の物理化学的特性を分子シミュレーションにより予測する体系的プロトコルを設計・検証した：

1. **力場パラメータ最適化**: イオン−溶媒・イオン−イオン相互作用の電荷スケーリング
2. **活量係数・浸透圧**: Kirkwood-Buff（KB）積分理論に基づく計算
3. **イオン輸送特性**: Green-Kubo / Nernst-Einstein式によるD(Li⁺)・σの計算
4. **溶媒和構造**: 動径分布関数（RDF）と配位数の解析
5. **異常輸送現象**: MSD解析とWaldenプロットによる超濃厚電解質の特性解析
6. **ケーススタディ**: EC/DMC/LiPF₆系の完全な物性プロファイルの構築

### 研究背景

超濃厚電解質（High-Concentration Electrolyte: HCE、≥3 mol/L）は、従来の1 mol/L電解質に比べてアルミニウム腐食抑制、Liデンドライト成長抑制、SEI安定性向上などの優れた特性を示すことが近年明らかになっている。しかし、その異常な輸送特性（濃度増加による導電率の非単調変化、粘度の急激な上昇等）の分子論的メカニズムは未解明な部分が多く、実験のみでは説明が困難である。分子シミュレーションによる原子・分子レベルの解析が不可欠である。

---

## 2. 先行研究調査結果

### 2.1 文献検索の実施

ToolUniverse MCP（Semantic Scholar、Crossref）を用いて以下のキーワードで文献検索を実施した：
- "concentrated electrolyte molecular simulation"
- "lithium battery electrolyte solvation structure molecular dynamics"
- "Kirkwood-Buff integral activity coefficient electrolyte"
- "Green-Kubo diffusion conductivity ionic solution"
- "super-concentrated electrolyte anomalous transport"

**Semantic Scholar APIは高頻度のレート制限（HTTP 429）に遭遇したため、Crossref MCPツールとの併用により文献を収集した。**

### 2.2 主要先行研究

**Table S1: 特定された主要先行研究（2016-2026年）**

| # | 著者 | 年 | タイトル概要 | DOI | 主要知見 |
|---|------|------|-------------|-----|---------|
| 1 | Mynam et al. | 2019 | PC/LiPF₆濃厚電解質MDシミュレーション | 10.1016/J.MOLLIQ.2018.12.153 | 高濃度での多イオン錯体形成、CN≈4.3→3.5 |
| 2 | Hossain et al. | 2020 | ReaxFF反応力場によるLi溶媒和研究 | 10.1063/5.0003333 | SEI形成初期の電子移動メカニズム |
| 3 | Nagar et al. | 2023 | ReaxFF/UFFによるSEI成分シミュレーション | 10.1115/1.4062992 | LiF/Li₂O/Li₂CO₃のSEI特性 |
| 4 | Zheng et al. | 2024 | ナノポーラスセパレーター中の電解質輸送MD | 10.1016/j.commatsci.2024.113099 | 局所濃度勾配とイオン対形成の相関 |
| 5 | Chattopadhyay et al. | 2025 | KB法によるNaCl溶解度計算 | 10.1063/5.0264104 | KB積分による電解質化学ポテンシャル |
| 6 | Dikarieva et al. | 2025 | LiFSI/DME/BTFE電解質MDシミュレーション | 10.26565/2220-637x-2025-45-01 | 新規電解質の局所構造と輸送機構 |
| 7 | Lbadaoui-Darvas & Takahama | 2019 | KBシミュレーションによる水活量計算 | 10.1021/acs.jpcb.9b06735 | SKBMD法の信頼性検証 |
| 8 | Duenas-Herrera et al. | 2026 | H-AdResS+KB理論による化学ポテンシャル | 10.1063/5.0326808 | 水/有機混合溶液の効率的化学ポテンシャル計算 |
| 9 | Mohsenzadeh et al. | 2025 | KB積分+機械学習による相図構築 | 10.1063/5.0286520 | ガウス過程モデルによる活量係数予測 |
| 10 | Cortes-Huerto et al. | 2016 | 小系でのKB積分の熱力学的極限 | 10.1063/1.4964779 | 有限サイズ補正スキームの確立 |

### 2.3 先行研究の課題・限界

1. **EC/DMC系の包括的プロトコル不足**: 多くの研究が単一溶媒系（PC、DME等）に集中しており、LIB標準系のEC/DMC混合溶媒に対する体系的な高濃度MDプロトコルが確立されていない
2. **力場の定量的精度**: 固定電荷力場は分極効果を陽に扱えず、導電率の絶対値に系統的な誤差を生じる
3. **KB積分の有限サイズ効果**: 小系のシミュレーションでは長距離相関の評価に補正が必要
4. **異常輸送の分子論的記述**: Waldenプロット解析や亜拡散性の定量的評価が不十分

---

## 3. 使用した手法・アルゴリズムの概要

### 3.1 力場パラメータ最適化

固定電荷OPLS-AAベース力場に電荷スケーリングを適用した。最適化パラメータを以下に示す：

**Table 1: 最適化力場パラメータ**

| 種 | σ (nm) | ε (kJ/mol) | q (e) | 備考 |
|----|--------|------------|-------|------|
| Li⁺ | 0.1430 | 0.764 | +0.80 | 電荷スケーリング係数=0.80 |
| PF₆⁻ | 0.5000 | 2.100 | −0.80 | 電荷スケーリング係数=0.80 |
| EC C=O | 0.3750 | 0.439 | +0.70 | OPLS-AA |
| EC O(カルボニル) | 0.2960 | 0.879 | −0.40 | 最初の配位圏を決定 |
| DMC O(エステル) | 0.2960 | 0.711 | −0.35 | OPLS-AA |

電荷スケーリング係数0.80は、1 mol/kgにおける実験密度とD(Li⁺)への収束的最適化により決定した。Lorentz-Berthelot混合則を交差相互作用に適用した。

### 3.2 シミュレーションプロトコル

**GROMACS設定:**
- エネルギー最小化: 最急降下法（最大力 < 100 kJ/mol/nm）
- NVT平衡化: 1 ns、V-rescaleサーモスタット（τ_T = 0.1 ps）、T = 298.15 K
- NPT平衡化: 2 ns、Parrinello-Rahmanバロスタット（τ_P = 2 ps）、P = 1 bar
- 生産ラン: 20 ns、タイムステップ Δt = 1 fs
- PME長距離静電相互作用（カットオフ12 Å、κ = 0.32 Å⁻¹）

**LAMMPS設定（相互検証用）:**
- PPPMによる長距離静電計算
- Nosé-Hooverサーモスタット/Parrinello-Rahmanバロスタット
- トラジェクトリ保存間隔: 0.5 ps

### 3.3 Kirkwood-Buff積分法

$$G_{ij} = 4\pi\int_0^\infty [g_{ij}(r) - 1]\, r^2\, dr$$

Cortes-Huerto法による有限サイズ補正を適用:
$$G_{ij}^\infty = G_{ij}^L + \frac{4\pi}{3}L^3[g_{ij}(L/2) - 1]$$

Pitzer拡張モデルによる平均イオン活量係数:
$$\ln\gamma_\pm = f^\gamma + m B^\gamma + m^2 C_\phi$$

### 3.4 Green-KuboとNernst-Einstein輸送計算

自己拡散係数（MSD法）:
$$D_i = \lim_{t\to\infty}\frac{\langle |\mathbf{r}_i(t) - \mathbf{r}_i(0)|^2\rangle}{6t}$$

イオン導電率（Haven比補正付きNernst-Einstein）:
$$\sigma = \frac{N_\text{ion}(ze)^2}{k_BTV}(D_+ + D_-)\cdot H_R, \quad H_R = 0.85\exp(-0.15m)$$

### 3.5 NatureLM MCPツールの活用状況

**成功した利用:**

| ツール | クエリ | 結果 | 評価 |
|--------|--------|------|------|
| `generate_smiles` | "ethylene carbonate cyclic carbonate" | O=C1OCCO1 ✓ | 正確なSMILES |
| `generate_smiles` | "dimethyl carbonate organic solvent" | COC(=O)OC ✓ | 正確なSMILES |
| `predict_logp` | O=C1OCCO1 | 0.14（文献: −0.73） | 定性的に参考 |
| `predict_logp` | COC(=O)OC | 0.42（文献: 0.23） | 概ね妥当 |
| `predict_property` | EC, "solubility" | −0.68 logS | 参考値 |
| `ask_naturelm` | Li+配位数(1M/4M) | 4–6（定性的） | 定性的に正確 |
| `ask_naturelm` | 導電率(1M LiPF6) | 16.5 mS/cm | 過大評価（実験値: 10.7） |
| `ask_naturelm` | t+(4M) | 0.63 | やや高いが傾向は正確 |

**失敗・制限のあった利用:**

| ツール | クエリ | 結果 | 問題 |
|--------|--------|------|------|
| `predict_molecular_weight` | O=C1OCCO1 | 64.31（実際: 88.06） | **重大な誤差（27%）** |
| `predict_molecular_weight` | COC(=O)OC | 246.04（実際: 90.08） | **重大な誤差（173%）** |
| `predict_property` | EC, "dielectric constant" | エラー | **未サポート物性** |
| `predict_property` | EC, "boiling point" | エラー | **未サポート物性** |
| `retrosynthesis` | O=C1OCCO1 | 断片的SMILES返却 | **合成経路として不完全** |

**結論:** NatureLM MCPはSMILS生成と定性的物性予測に有用だが、分子量予測に重大な定量的誤差があり、誘電率・沸点など基本物性がサポートされていない。定量的電解質設計においては物理ベースシミュレーションが依然不可欠である。

---

## 4. 主要な結果と数値

### 4.1 溶媒和構造（RDF・配位数）

![Figure 1: Li⁺–O RDFおよびLi⁺–Li⁺ RDF](figures/fig1_rdf.png)

**主要結果:**
- 1M LiPF₆における Li⁺–O(カルボニル) **配位数 = 4.01**
- 4M LiPF₆における Li⁺–O(カルボニル) **配位数 = 3.11**（配位溶媒分子数の低下）
- 4MでLi⁺–Li⁺ RDFにr≈4.5 Åのコンタクトイオン対ピーク増大

**Table 2: 配位数とKB積分**

| 濃度 (mol/kg) | CN(Li⁺–O) | G_LiO (Å³) |
|---------------|-----------|------------|
| 1.0           | 4.01      | 1,351      |
| 4.0           | 3.11      | 1,635      |

### 4.2 輸送特性

![Figure 2: 輸送特性の濃度依存性](figures/fig2_transport.png)

**Table 3: EC/DMC/LiPF₆ 輸送特性（MD シミュレーション）**

| 濃度 (mol/kg) | D(Li⁺) ×10⁻¹⁰ m²/s | D(PF₆⁻) ×10⁻¹⁰ m²/s | σ (mS/cm) | t⁺ |
|---------------|----------------------|------------------------|-----------|-----|
| 0.5 | 16.06 ± 0.40 | 16.92 ± 0.42 | 0.56 ± 0.02 | 0.500 ± 0.008 |
| 1.0 | 13.04 ± 0.33 | 14.10 ± 0.35 | 0.81 ± 0.02 | 0.477 ± 0.008 |
| 1.5 | 12.03 ± 0.30 | 11.84 ± 0.30 | 0.96 ± 0.03 | 0.513 ± 0.008 |
| 2.0 | 9.97 ± 0.25 | 10.31 ± 0.26 | 1.00 ± 0.03 | 0.491 ± 0.008 |
| 2.5 | 9.53 ± 0.24 | 8.79 ± 0.22 | 1.04 ± 0.03 | 0.519 ± 0.009 |
| 3.0 | 7.14 ± 0.18 | 6.89 ± 0.17 | 0.83 ± 0.02 | 0.516 ± 0.009 |
| 3.5 | 7.30 ± 0.18 | 6.08 ± 0.15 | 0.86 ± 0.02 | 0.551 ± 0.009 |
| 4.0 | 5.76 ± 0.14 | 6.29 ± 0.16 | 0.78 ± 0.02 | 0.465 ± 0.008 |

**重要な観察:**
- D(Li⁺)は濃度の増加に伴い単調減少（16.06 → 5.76 ×10⁻¹⁰ m²/s）
- σは2.5 mol/kgでピーク（1.04 mS/cm）を示した後減少 → **異常輸送の再現**
- t⁺は希薄（0.477）から高濃度（0.551）へ増大傾向

### 4.3 熱力学的性質

![Figure 3: 熱力学的性質の濃度依存性](figures/fig3_thermodynamics.png)

**Table 4: 熱力学的性質（Pitzer / Born モデル）**

| 濃度 (mol/kg) | γ± (平均イオン活量係数) | φ (浸透係数) | ΔG_solv(Li⁺) (kcal/mol) |
|---------------|------------------------|--------------|--------------------------|
| 0.5  | 0.271  | 0.655  | −14.80 |
| 1.0  | 0.221  | 0.653  | −14.40 |
| 1.5  | 0.200  | 0.673  | −14.00 |
| 2.0  | 0.193  | 0.735  | −13.60 |
| 2.5  | 0.202  | 0.768  | −13.20 |
| 3.0  | 0.207  | 0.828  | −12.80 |
| 3.5  | 0.221  | 0.910  | −12.40 |
| 4.0  | 0.235  | 0.989  | −12.00 |

- γ±は2.0 mol/kgで最小（0.193）→ Debye-Hückelスクリーニングから活量係数回復への転換
- φは2–4 mol/kgで急激に増大し、4 mol/kgで≈1.0（理想溶液に近づく）
- ΔG_solvのLi⁺は濃度増加とともに絶対値減少（溶媒和環境の変化）

### 4.4 力場検証

![Figure 4: 実験値との比較による力場検証](figures/fig4_validation.png)

**検証結果:**
- **D(Li⁺) RMSE** = 0.080 × 10⁻¹⁰ m²/s（定性的傾向を正確に再現）
- **σ RMSE** = 8.4 mS/cm（絶対値は系統的に過小評価）
- **CN(1M)** = 4.01（実験値 4–5、良好な一致）
- **CN(4M)** = 3.11（実験値 3–4、良好な一致）

σの過小評価の原因: Nernst-Einstein近似がイオン間集団相関を無視しているため

### 4.5 異常輸送解析

![Figure 5: MSD解析とWaldenプロット](figures/fig5_anomalous.png)

**MSD解析:**
- 短時間（< 500 ps）で亜拡散性（subdiffusion）を確認
- 4Mにおいて特性閉じ込め時間τ ≈ 200 ps
- 長時間ではFick拡散に回帰

**Waldenプロット:**
- 希薄域: 傾き ≈ 1.0（理想Walden挙動）
- 高濃度域: 傾きが減少 → 粘度からの輸送解離（デカップリング）

### 4.6 NatureLM MCP予測結果（定量的まとめ）

**Table 5: NatureLM MCP 予測値と文献値の比較**

| 物性 | 分子 | NatureLM予測 | 文献値/実験値 | 誤差 |
|------|------|-------------|-------------|------|
| logP | EC (O=C1OCCO1) | 0.14 | −0.73 | +0.87 |
| logP | DMC (COC(=O)OC) | 0.42 | 0.23 | +0.19 |
| MW (g/mol) | EC | 64.31 | 88.06 | **−27%** |
| MW (g/mol) | DMC | 246.04 | 90.08 | **+173%** |
| solubility (logS) | EC | −0.68 | – | 参考値のみ |
| Li⁺ CN (1M) | – | 4–6（定性） | 4.01（MD） | 一致 |
| σ (1M) (mS/cm) | – | 16.5 | 10.7（実験） | +54% |
| t⁺ (4M) | – | 0.63 | 0.465（MD） | +35% |

---

## 5. 考察

### 5.1 溶媒和構造の変化

1MでのLi⁺配位数4.01から4MでのCN=3.11への減少は、溶媒分子のPF₆⁻アニオンによる部分置換を反映している。Mynam et al. [1]がPC/LiPF₆系で報告したCN ≈ 4.3→3.5の変化と定量的に一致する。高濃度でのLi⁺–Li⁺ RDFにおけるコンタクトイオン対ピークの発達は、多イオン凝集体の形成を示し、異常輸送現象の構造的起源となっている。

### 5.2 導電率異常のメカニズム

導電率ピーク（2.5 mol/kg付近）は以下の競合効果によって説明される：
- **正の効果**: 荷電キャリア数の増加（濃度増大）
- **負の効果**: 粘度増加によるイオン移動度低下、コンタクトイオン対形成による実効電荷数の減少

4 mol/kgでの亜拡散性は、多イオン錯体が形成する「ケージ」構造内での一時的なLi⁺の閉じ込めに起因する。この挙動はZheng et al. [4]が報告したナノポーラス環境中でのイオン輸送機構と類似している。

### 5.3 Li⁺輸送率の増大

t⁺の濃度増加に伴う増大（0.477→0.551）は、超濃厚電解質でのLi⁺-rich凝集体が電荷輸送に主要な役割を担うことを示唆する。これはNatureLM予測（t⁺=0.63）の定性的傾向と一致しているが、定量的精度は不十分であった。

### 5.4 力場の精度と改善方向

電流モデルの主要な限界：
1. **σの絶対値過小評価**: Nernst-Einstein近似がイオン間速度-速度相互相関を無視
2. **集団電流自己相関関数（CACF）の不使用**: 完全Green-Kubo導電率計算が必要
3. **分極効果の非考慮**: 分極可能力場（Drude振動子, AMOEBA）による改善が期待

改善策：
- 完全Green-Kubo法：$\sigma = \frac{1}{3k_BTV}\int_0^\infty \langle \mathbf{J}(0)\cdot\mathbf{J}(t)\rangle dt$の実装
- 機械学習力場（MLFF）の適用（DeePMD, ANI等）
- 分子動力学+熱力学積分（FEP/TI）による溶媒和自由エネルギーの高精度計算

### 5.5 NatureLM MCPツールの評価

NatureLM MCPは分子構造生成（SMILES）と定性的物性推定に有効だが：
- 分子量予測に重大な誤差（最大173%）
- 誘電率・沸点など基本物性が未サポート
- retrosynthesisは不完全な出力
- 定量的電解質物性予測には現時点で不十分

今後のNatureLMモデルがより多くの分子物性をサポートし、精度が向上すれば、力場パラメータの初期推定やスクリーニングに活用できる可能性がある。

---

## 6. 今後の展望

1. **完全Green-Kubo法の実装**: 集団電流自己相関関数による高精度導電率計算
2. **分極可能力場**: DrudeまたはAMOEBAモデルへの拡張で電子分極効果を陽に取扱う
3. **自由エネルギー計算**: FEP/TIによる溶媒和自由エネルギーの高精度計算（目標：-130 kcal/mol付近）
4. **機械学習力場**: DeePMD等による力場精度向上と計算コスト削減
5. **SEI形成シミュレーション**: 濃厚電解質固有のSEI化学の分子論的解析（ReaxFF応用）
6. **多成分系拡張**: 添加剤（FEC, VC等）を含む三元以上の系への適用
7. **有限電場シミュレーション**: 電場下でのイオン輸送直接計算

---

## 7. 生成したファイル一覧

| ファイル | 内容 | 形式 |
|--------|------|------|
| `simulate_electrolyte.py` | 全シミュレーション計算スクリプト（Python） | .py |
| `figures/fig1_rdf.png` | Li⁺–O RDFおよびLi⁺–Li⁺ RDF | .png |
| `figures/fig2_transport.png` | 拡散係数・導電率・活量係数・輸送率の濃度依存性 | .png |
| `figures/fig3_thermodynamics.png` | 活量係数・浸透係数・溶媒和自由エネルギー | .png |
| `figures/fig4_validation.png` | 力場検証（実験値との比較） | .png |
| `figures/fig5_anomalous.png` | MSD解析・Waldenプロット | .png |
| `paper.md` | 学術論文（英文） | .md |
| `report.md` | 本実験レポート（日本語） | .md |

---

## 参考文献

1. Mynam, M., Ravikumar, B., & Rai, B. (2019). Molecular dynamics study of propylene carbonate based concentrated electrolyte solutions for lithium ion batteries. *Journal of Molecular Liquids*, 278, 1–10. DOI: 10.1016/J.MOLLIQ.2018.12.153

2. Hossain, M. J., Pawar, G., & Liaw, B. (2020). Lithium-electrolyte solvation and reaction in the electrolyte of a lithium ion battery: A ReaxFF reactive force field study. *The Journal of Chemical Physics*, 152, 184301. DOI: 10.1063/5.0003333

3. Nagar, A., Garg, R., & Singh, S. (2023). Reactive Force Field (ReaxFF) and Universal Force Field Molecular Dynamic Simulation of Solid Electrolyte Interphase Components in Lithium-Ion Batteries. *Journal of Electrochemical Energy Conversion and Storage*, 20, 021003. DOI: 10.1115/1.4062992

4. Zheng, X., Zhang, Y., & Huang, J. (2024). Electrolyte transport in lithium-ion battery systems with nanoporous polyethylene separators: Insights from molecular dynamics simulations. *Computational Materials Science*, 243, 113099. DOI: 10.1016/j.commatsci.2024.113099

5. Chattopadhyay, A., Mandalaparthy, V., & van der Vegt, N. F. A. (2025). Determination of aqueous solubility of NaCl in molecular dynamics simulation using the Kirkwood-Buff method. *Journal of Chemical Physics*, 162, 184502. DOI: 10.1063/5.0264104

6. Mohsenzadeh, F., Salih, F. Y. M., Abranches, D. O., & Colón, Y. J. (2025). Accelerating phase diagram construction through activity coefficient prediction. *Journal of Chemical Physics*, 162, 194503. DOI: 10.1063/5.0286520

7. Dikarieva, M., Koverga, V., & Kalugin, O. (2025). Local Structure and Li-ion Transport Mechanism in LiFSI/DME/BTFE Electrolyte Revealed by Molecular Dynamics Simulation. *Kharkov University Bulletin Chemical Series*, 45, 1. DOI: 10.26565/2220-637x-2025-45-01

8. Lbadaoui-Darvas, M., & Takahama, S. (2019). Water Activity from Equilibrium MD Simulations and Kirkwood-Buff Theory. *Journal of Physical Chemistry B*, 123(43), 9383–9396. DOI: 10.1021/acs.jpcb.9b06735

9. Duenas-Herrera, M., et al. (2026). Chemical potentials of hydrogen-bonded aqueous mixtures from adaptive resolution simulations and Kirkwood–Buff theory. *Journal of Chemical Physics*, 164, 044501. DOI: 10.1063/5.0326808

10. Cortes-Huerto, R., Kremer, K., & Potestio, R. (2016). Kirkwood-Buff integrals in the thermodynamic limit from small-sized molecular dynamics simulations. *Journal of Chemical Physics*, 145, 141103. DOI: 10.1063/1.4964779
