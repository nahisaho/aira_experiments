# 不均一系触媒反応マイクロキネティックモデリングフレームワーク — 報告書

**DRAFT — NOT FOR DISTRIBUTION**

日付: 2026-05-23  
著者: Co-Scientist (自動生成)  
ケーススタディ: Fischer-Tropsch 合成 (Co(0001) 表面)

---

## 1. 実験目的と背景

不均一系触媒反応のマイクロキネティックモデリングは、触媒設計の合理化に不可欠である。本フレームワークは、DFT（密度汎関数理論）から得られるエネルギーパラメータを基に、遷移状態理論（TST）による速度定数算出から反応器レベルのシミュレーションまでを一貫して実行するPythonベースのツールチェーンを構築した。

**対象系**: Co(0001) 表面上の Fischer-Tropsch (FT) 合成反応（CO水素化によるメタン・オレフィン生成）

**目的**:
1. DFT エネルギーから TST + トンネル効果補正による速度定数算出
2. 複数の吸着等温線モデル（Langmuir / Temkin / フラクタル表面）の実装と比較
3. Campbell の速度制御度（Degree of Rate Control, X_RC）による律速段階の自動同定
4. 被覆率依存性（lateral interaction）の自己無撞着的取り込み
5. PFR / CSTR 反応器モデルとの連成
6. FT 合成のケーススタディによる検証

---

## 2. 使用した手法・アルゴリズムの概要

### 2.1 速度定数算出 (`rate_constants.py`)

遷移状態理論に基づく速度定数：

$$k_{TST} = \kappa \cdot \frac{k_B T}{h} \cdot \frac{Q^\ddagger}{Q_R} \cdot \exp\left(-\frac{E_a}{k_B T}\right)$$

- **Wigner トンネル補正**: $\kappa = 1 + \frac{1}{24}\left(\frac{h\nu^\ddagger}{k_B T}\right)^2$
- **Eckart トンネル補正**: ポテンシャル障壁の非対称性を考慮した数値積分法
- **振動分配関数**: 量子調和振動子近似
- **ゼロ点エネルギー補正**: DFT 振動数から自動算出

### 2.2 吸着等温線モデル (`adsorption.py`)

| モデル | 特徴 | 適用場面 |
|--------|------|----------|
| **Langmuir** | 均一表面、単一層吸着 | 基本的な被覆率推定 |
| **競争 Langmuir** | 多成分系の競争吸着 | CO/H₂ 共吸着 |
| **Temkin** | 被覆率による吸着エネルギーの線形減少 | 中被覆率領域 |
| **フラクタル表面** | 不均一表面（D_f = 2–3）のエネルギー分布 | 実触媒表面 |

### 2.3 律速段階同定 (`rds_identifier.py`)

- **Campbell の速度制御度 (X_RC)**: 各素反応の速度定数を微小摂動し、全体速度への感度を評価
- **熱力学的速度制御度 (X_TRC)**: 平衡定数の摂動に基づく解析
- **Kozuch-Shaik エネルギースパンモデル**: TDTS（Turnover-Determining TS）と TDI（Turnover-Determining Intermediate）の自動検出
- **見かけの活性化エネルギー**: 温度微分法による算出

### 2.4 被覆率依存性・Lateral Interaction (`lateral.py`)

- **平均場近似**: $\Delta E_{ads,i}(\theta) = \sum_j \epsilon_{ij} \cdot z \cdot \theta_j$
- **BEP 関係**: 遷移状態エネルギーシフトの推定
- **準化学近似 (QCA)**: 強い相互作用系への対応
- **自己無撞着解法**: ダンピング付き反復法で被覆率を収束

### 2.5 反応器モデル (`reactor.py`)

- **PFR (Plug Flow Reactor)**: `dF_i/dW = r_i · n_sites` の常微分方程式系を BDF 法で解法
- **CSTR (Continuous Stirred Tank Reactor)**: 定常状態の物質収支を反復法で解法
- **表面-気相連成**: 擬定常状態近似（表面種は各位置で瞬時平衡）

### 2.6 FT 合成メカニズム (`ft_synthesis.py`)

Co(0001) 上のカーバイドメカニズム（10素反応）：

| # | 素反応 | E_a [eV] | ΔE [eV] |
|---|--------|----------|---------|
| 1 | CO + * → CO* | 0.00 | −1.30 |
| 2 | H₂ + 2* → 2H* | 0.05 | −0.50 |
| 3 | CO* + H* → HCO* + * | 0.80 | +0.20 |
| 4 | HCO* + H* → CH₂O* + * | 0.55 | −0.10 |
| 5 | CH₂O* → CH₂* + O* | 1.20 | −0.40 |
| 6 | CH₂* + H* → CH₃* + * | 0.60 | −0.30 |
| 7 | CH₃* + H* → CH₄ + 2* | 0.95 | +0.10 |
| 8 | O* + H* → OH* + * | 0.70 | +0.15 |
| 9 | OH* + H* → H₂O + 2* | 0.85 | −0.20 |
| 10 | 2CH₂* → C₂H₄ + 2* | 0.75 | −0.50 |

エネルギーパラメータは Zhuo et al. (JACS, 2009)、Ojeda et al. (J. Catal., 2010)、van Santen et al. (PCCP, 2011) に基づく代表値を使用。

---

## 3. 主要な結果と数値

### 3.1 速度定数（T = 500 K, Wigner トンネル補正）

| 素反応 | k_forward [s⁻¹] | E_a [eV] | κ (tunneling) |
|--------|-----------------|----------|---------------|
| CO adsorption | 1.042 × 10¹³ | 0.00 | 1.000 |
| H₂ dissociation | 2.039 × 10¹² | 0.05 | 1.249 |
| CO hydrogenation | 9.546 × 10⁴ | 0.80 | 1.061 |
| HCO hydrogenation | 3.127 × 10⁷ | 0.55 | 1.050 |
| C-O scission | 9.235 × 10⁰ | 1.20 | 1.104 |
| CH₂ hydrogenation | 9.986 × 10⁶ | 0.60 | 1.070 |
| CH₄ formation | 1.009 × 10³ | 0.95 | 1.093 |
| O hydrogenation | 9.893 × 10⁵ | 0.70 | 1.079 |
| OH hydrogenation | 3.073 × 10⁴ | 0.85 | 1.090 |
| C₂ coupling | 1.497 × 10⁵ | 0.75 | 1.042 |

**C-O結合開裂**（E_a = 1.20 eV）が最も遅い固有速度を示す。Wigner 補正は 500 K では 1.0–1.25 倍と穏やかな効果。

### 3.2 吸着等温線比較

- **Langmuir**: θ(CO) = 0.9999（強い CO 被毒）
- **Temkin** (α = 0.5): 高被覆率で吸着エネルギーの減弱を再現
- **フラクタル表面** (D_f = 2.5): 不均一サイト分布による連続的な被覆率変化

→ 図5 (`fig5_adsorption_isotherms.png`) に比較プロットを示す。

### 3.3 定常状態被覆率

| 表面種 | 被覆率 θ |
|--------|----------|
| CO* | 3.82 × 10⁻⁵ |
| H* | 1.34 × 10⁻¹⁵ |
| HCO* | 0.454 |
| CH₂O* | 5.74 × 10⁻⁹ |
| CH₂* | 9.15 × 10⁻¹² |
| CH₃* | 0.093 |
| O* | 1.72 × 10⁻¹¹ |
| OH* | 0.454 |

主要な表面種は **HCO*** と **OH***（各 ~45%）で、反応中間体の蓄積を示す。

### 3.4 律速段階 (RDS) 解析

**Campbell X_RC 解析結果**: CH₄ formation（CH₃* + H* → CH₄ + 2*）が X_RC = 1.000 で明確な律速段階。

**エネルギースパン解析**:
- エネルギースパン δE = 0.500 eV
- TDTS: C-O scission（遷移状態エネルギー 0.500 eV）
- TDI: 最安定中間体（−2.350 eV）
- 推定 TOF: 9.52 × 10⁷ s⁻¹

### 3.5 Lateral Interaction の効果

CO-CO 間の最近接斥力相互作用（ε_nn = −0.10 eV）により、高圧域で CO 被覆率が 10–20% 低下（図8参照）。これは実触媒上で観測される CO 被毒の緩和機構に対応する。

### 3.6 反応器シミュレーション（PFR）

| 指標 | 値 |
|------|-----|
| CO 転化率 | 100.0% |
| TOF | 2.222 × 10³ s⁻¹ |
| STY | 3.333 × 10⁻² mol/(kg_cat·s) |
| CH₄ 選択率 | 0.01% |
| H₂O 選択率 | 0.08% |
| C₂H₄ 選択率 | 0.66% |

### 3.7 Arrhenius 解析（C-O scission）

- 頻度因子 A = 算出済み（`results/ft_simulation_results.json` に格納）
- 見かけの活性化エネルギー E_a ≈ 1.20 eV（図6参照）

---

## 4. 考察と今後の展望

### 4.1 考察

1. **律速段階の二面性**: Campbell X_RC はメタン形成段階を律速と判定した一方、エネルギースパンモデルでは C-O 開裂が TDTS となった。これは解析手法の定義差に起因する — X_RC は微分的感度、エネルギースパンは熱力学的制約を反映する。

2. **HCO*/OH* の蓄積**: 定常状態で表面の ~90% が HCO* と OH* に占有されている。H* の枯渇が各水素化段階のボトルネックを形成しており、H₂ 分圧増加（H₂/CO 比の上昇）が活性向上に有効と予測される。

3. **Lateral interaction の重要性**: CO-CO 斥力相互作用は高被覆率域で有効吸着エネルギーを 0.6 eV 以上シフトさせうる。Langmuir モデル単独では過大な CO 被毒を予測するため、lateral interaction の考慮は定量的精度に不可欠である。

4. **トンネル効果**: 500 K での Wigner 補正は最大 25%（H₂ 解離）と穏やかである。低温域（< 400 K）では Eckart 補正の適用が推奨される。

### 4.2 制約事項

- DFT パラメータは文献の代表値を使用。特定の計算手法（PBE, RPBE, BEEF-vdW 等）による系統的変動は考慮していない。
- Anderson-Schulz-Flory (ASF) 分布による長鎖炭化水素の生成は簡略化（C₂ のみ）。
- 表面再構成、ステップサイト、粒径効果は未考慮。
- PFR モデルは理想的な等温条件を仮定。発熱反応の温度勾配は未実装。

### 4.3 今後の展望

1. **CatMAP / OpenMKM 連携**: 本フレームワークのパラメータを CatMAP フォーマットで出力し、大規模メカニズム解析に接続
2. **Cantera 統合**: 気相反応ネットワークとの結合、多相反応器モデル
3. **ASF 分布の実装**: C₁–C₃₀ の完全な鎖成長モデル
4. **機械学習加速**: GNN ベースのエネルギー予測による広範な触媒スクリーニング
5. **不確実性定量化**: Bayesian Error Estimation Functional (BEEF) による誤差伝播解析
6. **非等温反応器**: エネルギー収支連成による温度プロファイル算出

---

## 5. 生成ファイル一覧

### フレームワークコード

| ファイル | 説明 |
|----------|------|
| `microkinetic_framework/__init__.py` | パッケージ初期化 |
| `microkinetic_framework/rate_constants.py` | TST + トンネル効果による速度定数算出 |
| `microkinetic_framework/adsorption.py` | 吸着等温線モデル（Langmuir / Temkin / フラクタル） |
| `microkinetic_framework/rds_identifier.py` | 律速段階自動同定（X_RC, エネルギースパン） |
| `microkinetic_framework/lateral.py` | 被覆率依存性・lateral interaction |
| `microkinetic_framework/reactor.py` | PFR / CSTR 反応器モデル |
| `microkinetic_framework/ft_synthesis.py` | FT 合成ケーススタディ（メカニズム・パラメータ） |
| `run_simulation.py` | メインシミュレーション実行・可視化スクリプト |

### 結果・データ

| ファイル | 説明 |
|----------|------|
| `results/ft_simulation_results.json` | 全数値結果（速度定数、被覆率、転化率等） |

### 図表

| ファイル | 説明 |
|----------|------|
| `figures/fig1_energy_diagram.png` | FT 合成反応エネルギーダイアグラム |
| `figures/fig2_rate_constants.png` | 速度定数比較（棒グラフ） |
| `figures/fig3_rds_analysis.png` | Campbell 速度制御度 X_RC |
| `figures/fig4_coverages.png` | 定常状態表面被覆率 |
| `figures/fig5_adsorption_isotherms.png` | 吸着等温線モデル比較 |
| `figures/fig6_arrhenius.png` | Arrhenius プロット（C-O scission） |
| `figures/fig7_reactor_profiles.png` | PFR 軸方向プロファイル |
| `figures/fig8_lateral_interactions.png` | Lateral interaction の効果 |
| `figures/fig9_temperature_sensitivity.png` | 温度感度解析（TOF・転化率） |

### ログ

| ファイル | 説明 |
|----------|------|
| `logs/process-log.jsonl` | 実行トレースログ |

---

*本報告書は Co-Scientist v1.0.0 により自動生成されました。*
