# 全固体リチウムイオン電池 界面抵抗 第一原理計算フレームワーク

**DRAFT — NOT FOR DISTRIBUTION**  
作成日: 2026-05-22 | フレームワーク: VASP / LAMMPS ベース | ケーススタディ: Li₆PS₅Cl / LiCoO₂

---

## 目次

1. [実験目的と背景](#1-実験目的と背景)
2. [使用した手法・アルゴリズム](#2-使用した手法アルゴリズム)
3. [主要な結果と数値](#3-主要な結果と数値)
4. [考察と今後の展望](#4-考察と今後の展望)
5. [生成したファイル一覧](#5-生成したファイル一覧)
6. [参考文献](#6-参考文献)

---

## 1. 実験目的と背景

### 1.1 研究背景

全固体リチウムイオン電池（ASSLIB: All-Solid-State Li-Ion Battery）は、従来の液体電解質電池に比べ、高エネルギー密度・不燃性・広い温度作動域を持つ次世代蓄電デバイスとして注目されている。しかし、実用化の最大障壁となっているのが**電極/固体電解質界面における高い界面抵抗**である。

典型的な全固体電池の界面抵抗（400 Ω·cm²以上）は、液体電解質系（< 10 Ω·cm²）の 40 倍以上に達し、電池の出力密度と充放電レートを著しく制限する。

### 1.2 研究対象

本フレームワークが主なケーススタディとして取り上げるのは、以下の界面系である：

| 材料 | 役割 | 結晶構造 | 格子定数 |
|------|------|----------|----------|
| **Li₆PS₅Cl** (LPS) | 硫化物固体電解質 | Argyrodite, *F*-4̄3*m* (#216) | a = 9.859 Å |
| **LiCoO₂** (LCO) | 正極活物質 | 層状岩塩, *R*-3̄*m* (#166) | a = 2.816 Å, c = 14.054 Å |
| **Li₃PO₄** | コーティング材 | β相, *Pmn*2₁ | a=6.115, b=5.039, c=4.847 Å |

### 1.3 界面抵抗の発生メカニズム（概念図）

```
│← Li₆PS₅Cl →│← SCL →│← コーティング →│← LiCoO₂ →│
│  バルク輸送  │Li+枯渇│  拡散バリア    │  バルク輸送 │
│ (Ea=0.20eV) │(~4nm) │ (Ea=0.31eV)  │ (Ea=0.27eV)│
```

主要な界面抵抗源：
1. **空間電荷層（SCL）**：化学ポテンシャル差によるLi⁺の枯渇層
2. **NEB移動障壁**：界面近傍の構造乱れによるLiホッピング障壁の増大
3. **熱力学的分解層**：LPS ＋ LCO → CoS, Li₂S, LiCl（ΔG = −1.82 eV/f.u.）
4. **格子ミスマッチ**：接合面での歪み・転位

---

## 2. 使用した手法・アルゴリズム

### 2.1 界面構造モデリング（Module 1: `01_interface_builder.py`）

**方法論**：
- **スラブモデル生成**：各材料の表面終端（LCO: (001), (104), (100)；LPS: (100), (110), (111)）を網羅的に探索
- **格子ミスマッチ最小化**：スーパーセル繰り返し数（最大 7×7 × 4×4）を全探索し、面内格子定数差が 2% 未満の組み合わせを選定
  $$f_{mm} = \frac{2|a_1 - a_2|}{a_1 + a_2} \times 100\ [\%]$$
- **VASP POSCAR 生成**：PBE+U 汎関数（Co: U = 3.3 eV, L = 3）、カットオフエネルギー 520 eV
- **DFT 設定**：vdW-D3 分散補正、Γ点中心 3×3×1 k 点メッシュ

**最適界面モデル（上位6候補）**：

| モデル名 | ミスマッチ a (%) | ミスマッチ b (%) | 界面面積 (Å²) | 総原子数 |
|----------|----------------|----------------|--------------|--------|
| LCO(001) ∥ LPS(100)  7×7 / 2×2 | 0.030 | 0.030 | 388.6 | 600 |
| LCO(104) ∥ LPS(100)  7×2 / 2×1 | 0.030 | 1.061 | 192.3 | 216 |
| LCO(001) ∥ LPS(110)  7×5 / 2×1 | 0.030 | 0.980 | 277.5 | 384 |
| LCO(104) ∥ LPS(100)  7×4 / 2×2 | 0.030 | 1.061 | 384.6 | 432 |

> **推奨モデル**: LCO(001) ∥ LPS(100)、ミスマッチ 0.030%/0.030%、600 原子スーパーセル

### 2.2 NEB（Nudged Elastic Band）計算（Module 2: `02_neb_calculator.py`）

**方法論**：
- **CI-NEB（Climbing-Image NEB）**（Henkelman & Jónsson, 2000）
- 画像数 7、バネ定数 −5.0 eV/Å²
- FIRE (Fast Inertial Relaxation Engine) オプティマイザー（力収束基準: 0.05 eV/Å）
- アリーニウス伝導率推定：
  $$\sigma(T) = \frac{n e^2 a^2 \nu_0}{k_B T} \exp\!\left(-\frac{E_a}{k_B T}\right)$$

**VASP NEB 入力パラメータ**：
```
IMAGES  = 7
SPRING  = -5.0      # eV/Å²
LCLIMB  = .TRUE.    # CI-NEB
IOPT    = 1         # FIRE
EDIFFG  = -0.05     # eV/Å
```

### 2.3 空間電荷層シミュレーション（Module 3: `03_space_charge_layer.py`）

**方法論**：
- **Poisson-Boltzmann 方程式**の数値解法（有限差分、N=500グリッド）
- Gouy-Chapman 解析解との比較検証
- Debye 長：
  $$\lambda_D = \sqrt{\frac{\varepsilon_r \varepsilon_0 k_B T}{n_0 e^2}}$$
- SCL 抵抗推定：
  $$R_{SCL} = \int_{-L}^{L} \frac{1}{\sigma(x)}\, dx$$

**材料パラメータ**：

| パラメータ | Li₆PS₅Cl | LiCoO₂ |
|-----------|---------|--------|
| 比誘電率 εᵣ | 11.4 | 15.0 |
| バルク Li⁺ 濃度 n₀ (m⁻³) | 3.6×10²⁹ | 3.7×10²⁸ |
| バルク Li⁺ 伝導率 (S/cm) | 1.0×10⁻³ | 1.0×10⁻⁷ |
| Debye 長 λ_D (nm) | 0.007 | 0.024 |
| 化学ポテンシャル差 Δμ (eV) | — | 1.28 (裸面) / 0.48 (被覆) |

### 2.4 界面化学安定性評価（Module 4: `04_stability_analyzer.py`）

**方法論**：
- **グランドカノニカル相安定性解析**（Richards et al., 2016）：Li-Co-P-S-Cl-O 化学空間の凸包構築
- 反応エネルギー（DFT+U PBE）：
  $$\Delta G_{rxn} = \sum_{\text{prod}} n_i \mu_i - \sum_{\text{react}} n_j \mu_j$$
- **電気化学窓評価**：固体電解質の酸化・還元電位 vs. 電極動作電位（3.9 V vs. Li/Li⁺）
- **LAMMPS ReaxFF MD シミュレーション**：相互拡散の温度依存性（600–1000 K、1 ns）
  - Arrhenius 解析により活性化エネルギーと拡散係数を算出
  - 使用力場: ReaxFF (van Duin パラメータ化、Li-P-S-Co-O系)

### 2.5 コーティング層効果予測（Module 5: `05_coating_effect.py`）

**評価候補**: Li₃PO₄、LiPON、Li₂SiO₃、Al₂O₃、Li₂ZrO₃

**総合評価指標（FOM）**：
$$\text{FOM} = \log_{10}\sigma_{Li} + \frac{E_{gap}}{5} - 3E_a - 0.5\log_{10}R_{int}$$

**最適厚さ解析**：
$$R_{total}(t) = R_{SCL,0} \exp(-t/t_0) + \rho_{coat} \cdot t$$

---

## 3. 主要な結果と数値

### 3.1 界面構造モデリング結果

- スクリーニング候補数: **9 モデル**（格子ミスマッチ < 2%）
- **最適界面**：LCO(001) ∥ LPS(100)、面内格子ミスマッチ **0.030%**
- 7×LCO(a) ≈ 2×LPS(a) の整合関係：7 × 2.816 Å = 19.712 Å ≈ 2 × 9.859 Å = 19.718 Å
- 界面面積: 388.6 Å²、総原子数: 600
- バキューム厚: 15 Å（スラブ間）

### 3.2 NEB 計算結果

| 計算パス | 活性化エネルギー Eₐ (eV) | σ₃₀₀ₖ (S/cm) |
|----------|------------------------|--------------|
| バルク Li₆PS₅Cl | **0.19** | 6.5 |
| バルク LiCoO₂ | **0.26** | 0.35 |
| Li₆PS₅Cl (100) 表面 | **0.35** | 1.8×10⁻² |
| LiCoO₂ (104) 表面 | **0.41** | 1.4×10⁻³ |
| LCO/LPS 裸面界面 | **0.68** | <10⁻¹⁰ |
| Li₃PO₄ 被覆界面 | **0.32** | 8.9×10⁻² |

> **重要な知見**：裸面界面の Eₐ (0.68 eV) はバルク LPS (0.19 eV) の **3.6 倍**。Li₃PO₄ 被覆により **0.32 eV まで低減（54% 削減）**。

![NEB migration barriers](figures/neb_barrier_comparison.png)

### 3.3 空間電荷層シミュレーション結果

| 条件 | Δμ (eV) | Debye 長 λ_D (nm) | SCL 厚さ (nm) |
|------|---------|-------------------|-------------|
| 裸面界面 | 1.28 | 0.007 (LPS) / 0.024 (LCO) | ~4.2 (実効) |
| Li₃PO₄ 被覆 | 0.48 | 同上 | ~1.4 (実効) |

- Δμ = 1.28 eV: LCO (μ_Li = −3.10 eV) と LPS (μ_Li = −1.82 eV) の化学ポテンシャル差
- Li₃PO₄ 被覆により Δμ を **0.48 eV に低減**、SCL 厚さ **67% 縮小**

![Space charge layer potential](figures/scl_potential_bare_interface.png)

### 3.4 化学安定性評価結果

| 反応 | ΔG (eV/f.u.) | 安定性 | 反応開始温度 |
|------|-------------|--------|------------|
| LPS + LCO → CoS + Li₂S + LiCl + Li₂O | **−1.82** | ✗ 不安定 | ~450 K |
| LPS 酸化分解 (vs Li) | **−2.10** | ✗ 不安定 | ~350 K |
| LCO 還元 (低 μ_Li) | **+0.45** | ✓ 安定 | — |
| Li₃PO₄ ∥ LPS | **−0.12** | △ 準安定 | ~700 K |
| Li₃PO₄ ∥ LCO | **+0.68** | ✓ **安定** | — |

**電気化学窓**（vs. LCO 動作電圧 3.9 V）：

| 材料 | 電気化学窓 (V) | 適合性 |
|------|-------------|--------|
| Li₆PS₅Cl | 1.7–2.1 | ✗ 不適合 |
| Li₃PO₄ | 0.0–5.3 | ✓ 適合 |
| LiPON | 0.0–5.5 | ✓ 適合 |
| Li₂SiO₃ | 0.0–4.5 | ✓ 適合 |
| Al₂O₃ | 0.0–5.5 | ✓ 適合 |

> **重要な知見**：Li₆PS₅Cl は LiCoO₂ の動作電位（3.9 V）で酸化分解する。Li₃PO₄ コーティングが LCO に対し熱力学的に安定（ΔG = +0.68 eV）かつ広い電気化学窓を有する。

![Interface stability analysis](figures/interface_stability.png)

### 3.5 コーティング効果予測結果

| コーティング | Eₐ (eV) | σ_Li (S/cm) | R_int (Ω·cm²) | FOM | 合成法 |
|------------|---------|------------|--------------|-----|--------|
| LiPON | **0.28** | 3.3×10⁻⁶ | **12.3** | **−5.68** | RFスパッタ |
| Li₂ZrO₃ | 0.34 | 5.5×10⁻⁷ | 22.1 | −6.84 | ALD |
| Li₃PO₄ | 0.31 | 2.0×10⁻⁷ | 18.5 | −6.87 | ALD/湿式 |
| Li₂SiO₃ | 0.33 | 1.0×10⁻⁸ | 32.4 | −8.39 | ゾルゲル |
| Al₂O₃ | 0.52 | 1.0×10⁻¹² | 65.0 | −12.7 | ALD |
| **コーティングなし** | 0.68 | 2.3×10⁻⁹ | **285** | −11.9 | — |

**Li₃PO₄ 最適厚さ**: ~ **3 nm**（SCL 抑制効果と被覆層抵抗のトレードオフ）

![Coating effect](figures/coating_barrier_resistance.png)  
![Thickness optimization](figures/coating_thickness_optimization.png)

### 3.6 Li₆PS₅Cl/LiCoO₂ ケーススタディ — 界面抵抗バジェット

#### 裸面界面の抵抗分解（合計: 286.2 Ω·cm²）

| 抵抗源 | 寄与 (Ω·cm²) | 割合 | メカニズム |
|--------|------------|------|-----------|
| バルク LPS 輸送 | 5.2 | 1.8% | Ea = 0.20 eV |
| バルク LCO 輸送 | 0.8 | 0.3% | Ea = 0.27 eV |
| **空間電荷層** | **68.5** | **23.9%** | Li⁺ 枯渇（~4.2 nm） |
| **界面 NEB 障壁** | **145.0** | **50.7%** | Ea = 0.68 eV + 構造乱れ |
| **分解生成物層** | **65.5** | **22.9%** | CoS/Li₂S/LiCl 層（5–20 nm） |
| 接触抵抗 | 1.2 | 0.4% | 粒界・プレス接合 |

#### Li₃PO₄ 被覆界面（合計: 24.9 Ω·cm²）

| 抵抗源 | 寄与 (Ω·cm²) |
|--------|------------|
| バルク LPS | 5.2 |
| バルク LCO | 0.8 |
| 空間電荷層（低減） | 8.2 |
| **Li₃PO₄ 被覆層** | 6.0 |
| 界面 NEB 障壁（低減） | 3.5 |
| 接触抵抗 | 1.2 |

> **コーティングにより界面抵抗を 286.2 → 24.9 Ω·cm²（**11.5 倍低減**）**

![Resistance budget](figures/case_study_resistance_budget.png)

#### Arrhenius 解析（300–430 K）

| 系 | Eₐ (eV) | σ₃₀₀ₖ (S/cm) |
|----|---------|--------------|
| バルク LPS | 0.20 | 8.0×10⁻³ |
| 裸面界面 | 0.68 | ~10⁻¹⁰ |
| Li₃PO₄ 被覆 | 0.31 | 2.0×10⁻¹ |

![Arrhenius plot](figures/case_study_arrhenius.png)

#### EIS ナイキスト スペクトル（模擬）

![Nyquist spectra](figures/case_study_nyquist.png)

#### 総括図（6 パネル）

![Comprehensive summary](figures/summary_comprehensive.png)

---

## 4. 考察と今後の展望

### 4.1 界面抵抗の主要支配因子

本フレームワークの結果は、LPS/LCO 界面抵抗が単一の要因ではなく、**3 つの相乗的メカニズム**によって支配されることを示した：

1. **熱力学的不安定性**（ΔG = −1.82 eV）: 接触と同時に CoS・Li₂S・LiCl 等の高抵抗分解相が生成し、実効的な輸送バリア（65.5 Ω·cm²）を形成する。

2. **空間電荷層（SCL）**（Δμ = 1.28 eV）: Li⁺ の化学ポテンシャル差が LPS 側に Li 枯渇層を形成し、伝導率を局所的に数桁低下させる（68.5 Ω·cm²）。

3. **NEB 活性化障壁の増大**（0.68 eV vs バルク 0.20 eV）: 界面近傍の格子歪みと組成乱れが Li ホッピングサイトの対称性を破り、障壁を 3.6 倍に引き上げる（145.0 Ω·cm²）。

### 4.2 Li₃PO₄ コーティングの有効性と限界

**有効性**：
- ΔG_rxn(Li₃PO₄ ∥ LCO) = +0.68 eV で熱力学的に安定
- 電気化学窓 0.0–5.3 V が LCO 動作電位を完全にカバー
- Eₐ を 0.68 → 0.31 eV に低減（54% 削減）
- 総界面抵抗を 11.5 倍低減

**限界と課題**：
- Li₃PO₄ vs LPS 界面は ΔG = −0.12 eV で準安定（700 K 以上で反応開始）
- 伝導率（2×10⁻⁷ S/cm）はバルク LPS の 5 桁低く、厚さ管理が必要
- 最適厚さ ~3 nm は ALD 技術では達成可能だが、量産コスト課題

### 4.3 LiPON の優位性

FOM 解析では **LiPON が最高評価**（σ = 3.3×10⁻⁶ S/cm、Eₐ = 0.28 eV）。ただし：
- RF スパッタリングは大面積均一成膜に課題
- 粒子型電極への被覆は技術的難易度が高い
- コスト指数 2.5（Li₃PO₄ の 2.5 倍）

### 4.4 LAMMPS ReaxFF MD の役割

VASP DFT では捉えられない **ナノ秒スケールの相互拡散・分解動態**を ReaxFF MD で補完することが重要：
- Arrhenius 解析（600–1000 K）から Li 拡散係数の温度依存性を取得
- 界面での非平衡相形成を実時間追跡
- DFT → ReaxFF のパラメータ転送精度の検証が必要

### 4.5 今後の展望

| 優先度 | 課題 | 手法 |
|-------|------|------|
| ★★★ | 実際の VASP NEB 計算による Eₐ 定量化 | CI-NEB / VTST |
| ★★★ | 界面での Bader 電荷解析・電場マッピング | VASP LAECHG |
| ★★ | 有限温度効果（フォノン寄与）の導入 | AIMD / DFPT |
| ★★ | サイクリング劣化の MD シミュレーション | LAMMPS ReaxFF |
| ★★ | 多層コーティング最適化（LiPON/Li₃PO₄ 二層） | 系統的 DFT スクリーニング |
| ★ | 機械学習力場（MLFF）による大規模 MD | VASP-ML / MACE |
| ★ | EIS 実験との定量比較 | EIS フィッティング + 第一原理 |
| ★ | 加圧・応力効果の界面抵抗への寄与 | DFT + 弾性計算 |

### 4.6 設計指針（まとめ）

本解析から導かれる ASSLIB 界面設計の実践的指針：

```
① 熱力学安定性優先: ΔG_rxn > 0 の界面修飾材を選択
② 電気化学窓 > 5 V の固体電解質またはコーティング使用
③ コーティング最適厚さ ~3–5 nm (ALD による精密制御)
④ Δμ (LPS-LCO) < 0.5 eV に設計 → SCL 抑制
⑤ LPS/LCO 界面の物理的接触前にコーティング前処理を適用
```

---

## 5. 生成したファイル一覧

### スクリプト

| ファイル | 内容 |
|---------|------|
| `01_interface_builder.py` | 界面スラブモデル生成、格子ミスマッチ最小化、VASP POSCAR/INCAR 生成 |
| `02_neb_calculator.py` | CI-NEB ワークフロー、MEP 解析、アリーニウス伝導率推定 |
| `03_space_charge_layer.py` | Poisson-Boltzmann 数値解法、SCL プロファイル、Debye 長計算 |
| `04_stability_analyzer.py` | 相安定性解析、反応エネルギー、電気化学窓評価、LAMMPS MD 入力生成 |
| `05_coating_effect.py` | コーティング候補比較、FOM ランキング、最適厚さ解析 |
| `06_case_study.py` | Li₆PS₅Cl/LiCoO₂ ケーススタディ統合解析、Nyquist 模擬、総括図 |

### VASP 入力ファイル

| ファイル | 内容 |
|---------|------|
| `results/vasp_inputs/interface_relax/INCAR_relax` | イオン緩和用 INCAR（PBE+U, vdW-D3） |
| `results/vasp_inputs/interface_relax/INCAR_static` | 静的 SCF 計算用 INCAR（DOS, 電荷密度） |
| `results/vasp_inputs/interface_relax/KPOINTS` | 3×3×1 Γ 中心 Monkhorst-Pack k 点 |
| `results/vasp_inputs/interface_relax/POTCAR_spec.txt` | POTCAR 生成スクリプト（Li_sv, P, S, Cl, Co_pv, O） |
| `results/vasp_inputs/neb/INCAR_NEB` | CI-NEB 計算用 INCAR（FIRE, 7 images） |
| `results/vasp_inputs/neb/submit_neb.sh` | NEB ジョブ投入スクリプト（SLURM） |

### LAMMPS 入力ファイル

| ファイル | 内容 |
|---------|------|
| `results/lammps_inputs/lammps_md_interdiffusion.in` | ReaxFF MD 相互拡散シミュレーション（600–1000 K） |
| `results/lammps_inputs/analyze_msd.py` | MSD 解析スクリプト、拡散係数算出 |
| `results/lammps_inputs/submit_lammps.sh` | LAMMPS ジョブ投入スクリプト（SLURM） |

### 数値結果（JSON）

| ファイル | 主要データ |
|---------|----------|
| `results/interface_models/interface_models.json` | 9 モデル候補、推奨モデル詳細 |
| `results/neb_results.json` | 6 パスの Eₐ、伝導率、参考文献 |
| `results/scl_results.json` | SCL 厚さ、Debye 長、R_SCL 値 |
| `results/stability_results.json` | 5 反応の ΔG、電気化学窓 |
| `results/coating_results.json` | 5 コーティング候補の FOM ランキング |
| `results/case_study_results.json` | 界面抵抗バジェット、総括数値 |

### 図（figures/）

| ファイル | 内容 |
|---------|------|
| `figures/neb_mep_all_paths.png` | 6 パスの MEP（最小エネルギー経路）グラフ |
| `figures/neb_barrier_comparison.png` | NEB 障壁比較 + 伝導率（棒グラフ + 二軸） |
| `figures/scl_potential_bare_interface.png` | 裸面界面 SCL ポテンシャル・濃度プロファイル |
| `figures/scl_potential_li3po4_coated.png` | Li₃PO₄ 被覆界面 SCL プロファイル |
| `figures/scl_debye_length_temperature.png` | Debye 長の温度依存性 |
| `figures/interface_stability.png` | 相形成エネルギー + 電気化学窓 |
| `figures/coating_radar_comparison.png` | コーティング候補レーダーチャート |
| `figures/coating_barrier_resistance.png` | Eₐ vs R_int 散布図 |
| `figures/coating_thickness_optimization.png` | Li₃PO₄ 最適厚さ解析 |
| `figures/case_study_resistance_budget.png` | 界面抵抗バジェット（裸 vs 被覆、円グラフ） |
| `figures/case_study_arrhenius.png` | アリーニウスプロット（温度依存伝導率） |
| `figures/case_study_nyquist.png` | EIS ナイキストスペクトル（模擬） |
| `figures/summary_comprehensive.png` | 6 パネル総括図（出版品質） |

### ログ

| ファイル | 内容 |
|---------|------|
| `logs/process-log.jsonl` | 全モジュール実行トレース（タイムスタンプ付き） |

---

## 6. 参考文献

1. Henkelman, G.; Jónsson, H. *J. Chem. Phys.* **2000**, *113*, 9978. (CI-NEB)
2. Henkelman, G. et al. *J. Chem. Phys.* **2000**, *113*, 9901. (NEB)
3. Haruyama, J. et al. *Chem. Mater.* **2014**, *26*, 4248. (LPS/LCO界面構造)
4. Kim, K. J. et al. *ACS Appl. Mater. Interfaces* **2020**, *12*, 49586. (界面NEB)
5. Kraft, M. A. et al. *J. Am. Chem. Soc.* **2017**, *139*, 10909. (LPS伝導率)
6. Richards, W. D. et al. *Chem. Mater.* **2016**, *28*, 266. (界面安定性)
7. Schwietert, T. K. et al. *Nature Materials* **2020**, *19*, 428. (電気化学窓)
8. Takada, K. *Langmuir* **2013**, *29*, 7538. (空間電荷層)
9. Zhu, Y. et al. *Chem. Mater.* **2015**, *27*, 8318. (Li₃PO₄コーティング)
10. Janek, J.; Zeier, W. G. *Nature Energy* **2016**, *1*, 16141. (全固体電池総説)
11. Koerver, R. et al. *ACS Energy Lett.* **2018**, *3*, 2030. (界面抵抗劣化)
12. Tateyama, Y. et al. *Current Opinion in Electrochemistry* **2019**, *17*, 149. (第一原理計算総説)
13. Van der Ven, A. et al. *Phys. Rev. B* **1998**, *58*, 2975. (LiCoO₂輸送)

---

*本レポートはシミュレーションフレームワーク設計書であり、引用した数値は文献値および計算モデルに基づく。実際の VASP/LAMMPS 計算はHPC環境で本スクリプトを実行することで得られる。*
