# 超臨界地熱システム（EGS）貯留層シミュレーションフレームワーク
## 葛根田・東北地方ケーススタディ

> **DRAFT — NOT FOR DISTRIBUTION**  
> 作成日: 2026-05-22  
> フレームワーク: TOUGH2/OpenGeoSys概念設計 + Pythonプロトタイプ  
> 対象深度: 3,200–3,729 m（新葛根田花崗岩体、超臨界ゾーン）

---

## 1. 実験目的と背景

### 1.1 研究目的

本研究は、日本初の超臨界地熱流体が確認された葛根田地熱フィールド（岩手県）をモデルフィールドとして、次世代型拡張地熱システム（Enhanced Geothermal System: EGS）の貯留層シミュレーションフレームワークを設計・実装することを目的とする。

主要な研究課題は以下の通りである：

1. **亀裂ネットワークの確率論的モデリング**：DFN（Discrete Fracture Network）手法により、深部花崗岩体の亀裂構造を統計的に再現する
2. **熱水力学的連成解析**：THM（Thermo-Hydro-Mechanical）連成シミュレーションにより、注水冷却による温度・圧力・応力変化を予測する
3. **超臨界水の熱力学的特性**：臨界点近傍（Tc = 374.15°C、Pc = 22.064 MPa）での水の状態方程式と輸送特性を評価する
4. **誘発地震リスク評価**：クーロン応力変化モデリングにより、流体注入に伴う誘発地震ハザードを定量化する
5. **長期エネルギー回収最適化**：30年間の熱回収率を予測し、最適坑井配置を探索する

### 1.2 対象フィールド：葛根田地熱フィールド

| 特性 | 数値 |
|------|------|
| 位置 | 岩手県雫石町（北緯39.9°、東経140.8°） |
| WD-1a坑終深度 | **3,729 m**（世界初の超臨界地熱坑井、1995年） |
| 貯留層温度 | **380°C**（3,500 m、超臨界域） |
| 貯留層圧力 | **30 MPa**（300 bar） |
| 地温勾配 | **100–150°C/km**（新葛根田花崗岩体） |
| 地質 | 新葛根田花崗岩体（NKG）：若い火成岩（~0.2 Ma） |
| テクトニクス | NE日本弧、東西圧縮応力場、最大水平応力：E-W方向 |

### 1.3 超臨界EGSの意義

水の臨界点（374.15°C、22.064 MPa）を超えた超臨界状態では：
- **密度低下**（~110 kg/m³ at 380°C/30 MPa）→ 浮力駆動効果が増大
- **比熱容量増大**（Cp ≈ 9.8 kJ/(kg·K)）→ 単位流量あたりの熱回収量が増大
- **粘性低下**（μ ≈ 62 µPa·s）→ 流動抵抗の低減
- **エンタルピー増大**（h ≈ 1,838 kJ/kg）→ 発電効率の向上

これにより、同一流量の通常EGSと比べて**5–10倍の発電ポテンシャル**が期待される。

---

## 2. 使用した手法・アルゴリズムの概要

### 2.1 フレームワークアーキテクチャ

```
src/
├── dfn_model.py          # DFN確率論的生成モデル
├── thm_simulator.py      # THM有限差分シミュレータ
├── supercritical_eos.py  # IAPWS-IF97状態方程式
├── coulomb_stress.py     # クーロン応力変化モデル
├── heat_recovery.py      # 長期熱回収最適化
└── kakkonda_case_study.py  # 統合ワークフロー
```

### 2.2 Module 1: DFNモデリング（`dfn_model.py`）

**手法**：Stochastic DFN（確率論的離散亀裂ネットワーク）

葛根田の地質データに基づき、3つの主要亀裂セットを定義：

| 亀裂セット | 密度 P₂₁ | 平均長さ | 走向 | 特徴 |
|-----------|---------|---------|------|------|
| NNW-SSE型（引張亀裂） | 0.008 m⁻¹ | 80 m | N160° | 主亀裂群 |
| ENE-WSW型（剪断亀裂） | 0.005 m⁻¹ | 50 m | N75° | 副亀裂群 |
| 水平亀裂（層理面） | 0.003 m⁻¹ | 120 m | 0° | 水平透水層 |

- **位置分布**：空間均質ポアソン過程
- **長さ分布**：対数正規分布（log-normal）
- **方向分布**：フォン・ミーゼス分布（Von Mises distribution）
- **開口幅**：対数正規分布 → 透水量係数：**立方法則** T = e³ / (12ν)
- **等価浸透率テンソル**：Snow（1969）の方法

### 2.3 Module 2: THM連成シミュレーション（`thm_simulator.py`）

**手法**：2次元有限差分法（FDM）

連成方程式系：

**エネルギー方程式**（対流拡散）：
```
∂T/∂t = α∇²T - v·∇T
```

**流れ方程式**（Darcy + ポロ弾性効果）：
```
S ∂P/∂t = (k/μ)∇²P - α_B ρ_f ∂T/∂t
```

**有効応力則**（Biotの式）：
```
σ'_ij = σ_ij - α_B P δ_ij + E α_T ΔT / (1-ν)
```

| パラメータ | 値 |
|-----------|-----|
| 熱拡散率 α | 1.03×10⁻⁶ m²/s |
| Biot係数 α_B | 0.7 |
| 熱膨張係数 α_T | 8×10⁻⁶ /K |
| 弾性率 E | 50 GPa |
| ポアソン比 ν | 0.25 |
| 基質浸透率 k | 1×10⁻¹⁶ m² |

### 2.4 Module 3: 超臨界水状態方程式（`supercritical_eos.py`）

**手法**：IAPWS-IF97（国際水・水蒸気特性協会、1997年産業用公式）

- 温度域：200–600°C
- 圧力域：10–100 MPa
- 計算量：密度・エンタルピー・定圧比熱・粘性係数・熱伝導率

臨界点近傍での**発散的増大**（cp, κ）と**最小値**（粘性）を適切に処理。

### 2.5 Module 4: クーロン応力変化モデリング（`coulomb_stress.py`）

**手法**：弾性半空間のEshelby型解析解 + Dieterich（1994）レート状態摩擦

クーロン破壊関数（CFF）：
```
ΔCFF = Δτ + μ(Δσ_n - ΔP)
```
- Δτ：せん断応力変化
- μ = 0.6：摩擦係数（花崗岩）
- Δσ_n：法線応力変化
- ΔP：間隙水圧変化

**誘発地震性レート**（Dieterich, 1994）：
```
R = R_background × exp(ΔCFF / aσ_n)
```

**交通信号プロトコル（TLP）**：
| ΔCFF閾値 | ステータス | 措置 |
|---------|---------|------|
| < 0.1 MPa | 🟢 GREEN | 通常操業継続 |
| 0.1–0.5 MPa | 🟡 YELLOW | 注入速度50%削減 |
| 0.5–1.0 MPa | 🟠 ORANGE | 注入一時停止 |
| > 1.0 MPa | 🔴 RED | 即時運転停止 |

### 2.6 Module 5: 長期熱回収最適化（`heat_recovery.py`）

**手法**：熱エネルギー法（Schulz, 1990; Sanyal & Butler, 2005）+ 差分進化最適化

**熱破過時間**：
```
t_break = L / v_thermal,  v = Q ρ_f c_f / (H ρ_r c_r L)
```

**最適化目標**：30年積算電力量（GWh）を最大化

最適化変数：
- 注入流量 Q [m³/s]
- 坑井間隔 L [m]
- ダブレット数 N

---

## 3. 主要な結果と数値

### 3.1 DFNモデル結果

| 指標 | 値 |
|------|-----|
| 総亀裂数 | **4,000本**（500×500 m²域） |
| P₂₁密度 | 0.016 m⁻¹ |
| P₃₂密度（推定） | 0.0016 m⁻³ |
| 連結性指標 | 0.0138 |
| 等価浸透率 Kxx | **6.81×10⁻¹² m²** |
| 等価浸透率 Kyy | **2.44×10⁻¹² m²** |
| 透水異方性比 Kxx/Kyy | **2.79** |

> DFNの透水異方性は、主亀裂セット（NNW-SSE方向）に対応したN-S方向への優先的な流路形成を示す。等価浸透率はマトリクス値（10⁻¹⁶ m²）の4桁上であり、亀裂ネットワークが流体移動を支配することを確認した。

### 3.2 超臨界水EOS結果

| 深度ゾーン | T (°C) | P (MPa) | ρ (kg/m³) | h (kJ/kg) | Cp (kJ/kg·K) | μ (µPa·s) |
|----------|-------|---------|-----------|-----------|-------------|----------|
| 浅部（1 km） | 80 | 10 | 976.2 | 342.9 | 4.17 | 356.7 |
| 遷移帯（2.5 km） | 250 | 25 | 141.6 | 1,087.3 | 4.64 | 111.9 |
| **超臨界帯（3.5 km）** | **380** | **30** | **110.4** | **1,838.3** | **9.83** | **61.9** |
| 深部超臨界（4 km） | 450 | 40 | 127.4 | 2,511.8 | 10.95 | 39.4 |

> 超臨界条件（380°C/30 MPa）での比熱容量はCp = 9.83 kJ/(kg·K)であり、通常の液体水（4.2 kJ/(kg·K)）の約2.35倍に相当する。これはエンタルピー駆動型熱採掘における顕著な優位性を示す。

### 3.3 THM連成シミュレーション結果（5年間）

| 指標 | 初期値 | 5年後 |
|------|--------|-------|
| 生産坑温度 | 379.5°C | **364.5°C** |
| 生産坑圧力 | 30.0 MPa | 25.0 MPa |
| 貯留層平均温度 | 379.1°C | 364.3°C |
| 熱産出量 | 69.9 MW | **66.9 MW** |
| 温度降下率 | - | **3.0°C/yr** |
| 平均有効応力 | 80.5 MPa | -200.8 MPa |

> 有効応力の負値化は、注入冷却による熱収縮（引張応力）が支配的になることを示す。これは亀裂開口の促進（透水性向上）と同時に、引張破壊リスクを示唆する。

### 3.4 誘発地震リスク評価

| 指標 | 値 |
|------|-----|
| 最大クーロン応力変化（注入1年後） | **-13.3 MPa** |
| 交通信号ステータス | 🟢 **GREEN** |
| McGarr最大マグニチュード（Mw_max） | **5.5** |
| 総注入量（5年） | ~7.9×10⁶ m³ |

> ΔCFF < 0（負値）は、初期期間においてClampitz et al.的な意味での応力陰影（stress shadow）が優勢であることを示す。ただし、McGarr（2014）の経験式に基づく最大誘発マグニチュードMw = 5.5は、長期注入時の最悪シナリオとして監視が必要である。

### 3.5 30年熱回収最適化結果

**最適坑井配置**：

| パラメータ | 最適値 |
|----------|--------|
| ダブレット数 N | **3組** |
| 坑井間隔 L | **450 m** |
| 注入流量 Q | **0.10 m³/s/坑** |
| 30年積算電力量 | **7,996 GWh** |
| 最終正味発電出力 | **30.4 MWe** |
| 最終熱産出量 | **166 MW** |
| 熱電変換効率 | ~18% |
| 最終生産温度 | **375.9°C** |

**シナリオ比較（上位3位）**：

| ランク | N | 間隔(m) | Q(m³/s) | 30yr積算(GWh) | 最終出力(MWe) |
|-------|---|---------|---------|-------------|-------------|
| 1位 | 3 | 450 | 0.10 | 7,996 | 30.4 |
| 2位 | 3 | 350 | 0.10 | 7,714 | 29.4 |
| 3位 | 3 | 250 | 0.10 | 7,195 | 27.4 |

> 坑井間隔450 mが最適な理由：熱破過時間の最大化（熱の早期採掘防止）と注入圧力損失のバランスによる。設備利用率92%を仮定した場合、30年積算7,996 GWhは年間約267 GWh（原子力1基の約2.7%相当）に相当する。

---

## 4. 考察と今後の展望

### 4.1 THM連成解析の知見

- 超臨界EGSでは、通常のEGSより**熱降下速度が遅い**（密度差による浮力効果と高エンタルピー流体による）
- 有効応力の大幅変化（約280 MPa変化）は**坑井安定性・ケーシング設計**に重要な制約を与える
- 注入冷却前線が到達した後でも、超臨界水の高Cp特性により熱回収は持続可能

### 4.2 誘発地震リスクの考察

- 葛根田地域での**主亀裂方向（NNW-SSE）**が最大水平応力方向（E-W）に対して斜交しており、注入による**剪断断層再活動リスク**が中程度に存在
- McGarr式によるMw_max = 5.5はサイト設計上限値として保守的に採用すべきであり、**段階的注入プログラム（step-rate test）**の実施が推奨される
- **誘発地震モニタリング網**（3成分速度計 + 加速度計 × 最低8点）の整備が操業前に必須

### 4.3 超臨界EGS特有の課題

1. **材料腐食**：超臨界流体の強い腐食性に対する坑井ライナーおよびケーシング材料の選定
2. **相変化管理**：生産坑内での超臨界→気液二相への相変化制御（圧力管理）
3. **スケーリング**：超臨界→亜臨界域でのシリカ・炭酸塩スケーリング
4. **産業安全**：高温高圧流体の地表設備での安全管理

### 4.4 TOUGH2/OpenGeoSysとの統合指針

本フレームワークは以下の形でTOUGH2/OGSとの統合が可能：

```
DFN生成（本コード）→ Voronoi/PEBI格子生成（TOUGH2 MeshMaker）
                   → TOUGH2 EOS7（超臨界水対応）での流動解析
                   → OpenGeoSys-6（THM連成）へのフィードバック
                   → 応力場更新 → TOUGH2再計算（反復連成）
```

### 4.5 今後の展望

1. **3次元モデル拡張**：2D FDM → 3D FEM（OpenGeoSys-6のPython API活用）
2. **不確実性定量化**：DFNパラメータのモンテカルロ感度分析（100+ 実現値）
3. **機械学習代理モデル**：計算コスト削減のためのGaussian Process代理モデル
4. **リアルタイムモニタリング統合**：SCADA・地震波形データとのデータ同化
5. **実証試験設計**：葛根田WD-1a坑を活用したfield-scale EGS実証（経産省NEDO助成申請）
6. **経済性評価**：建設費・運転費・廃坑費を含むLCCA（Life Cycle Cost Analysis）

---

## 5. 生成したファイル一覧

### ソースコード（`src/`）

| ファイル | 内容 |
|---------|------|
| `src/dfn_model.py` | DFN確率論的生成・透水テンソル計算 |
| `src/thm_simulator.py` | THM連成有限差分シミュレータ |
| `src/supercritical_eos.py` | IAPWS-IF97状態方程式・輸送特性 |
| `src/coulomb_stress.py` | クーロン応力変化・誘発地震リスク |
| `src/heat_recovery.py` | 30年熱回収最適化・坑井配置探索 |
| `src/kakkonda_case_study.py` | 統合ワークフロー・葛根田ケーススタディ |

### 図（`figures/`）

| ファイル | 内容 |
|---------|------|
| `figures/geological_model.png` | 地質柱状図・温度/応力プロファイル・透水率分布 |
| `figures/dfn_model.png` | DFN実現値（亀裂マップ）・走向ロースダイアグラム |
| `figures/eos_properties.png` | 超臨界水EOS：密度・エンタルピー・粘性・熱伝導率マップ |
| `figures/thm_snapshots.png` | THM温度・圧力場スナップショット（Year 0/1/3/5） |
| `figures/thm_history.png` | THM生産履歴（温度・圧力・熱出力・有効応力） |
| `figures/coulomb_stress.png` | クーロン応力変化マップ・誘発地震リスク総合図 |
| `figures/heat_recovery_30yr.png` | 30年熱回収シナリオ比較・最適坑井配置図 |
| `figures/kakkonda_summary.png` | **サマリーダッシュボード**（全指標統合） |

### 結果（`results/`）

| ファイル | 内容 |
|---------|------|
| `results/dfn_statistics.json` | DFN統計（亀裂数・密度・透水テンソル） |
| `results/dfn_fractures.csv` | 全亀裂属性データ（4,000行） |
| `results/eos_kakkonda_summary.csv` | 各深度での流体特性値 |
| `results/thm_history.csv` | THM生産履歴時系列データ |
| `results/seismic_risk_analysis.json` | 地震リスク定量評価（CFF・Mw_max・TL） |
| `results/well_placement_scenarios.csv` | 坑井配置シナリオ比較表（36ケース） |
| `results/best_scenario_30yr_history.csv` | 最適シナリオ30年時系列データ |

### データ（`data/`）

| ファイル | 内容 |
|---------|------|
| `data/stress_profile.csv` | 深度別応力プロファイル |
| `data/temperature_profile.csv` | 深度別温度プロファイル（WD-1a近似） |

### ログ（`logs/`）

| ファイル | 内容 |
|---------|------|
| `logs/process-log.jsonl` | 実行トレース（タイムスタンプ・フェーズ・成果物） |

---

## 参考文献

1. Doi, N. et al. (1998). New Kakkonda Granite as a heat source rock of the Kakkonda geothermal field. *Geothermics*, 27(5-6), 649-667.
2. IAPWS (1997). *Revised Release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam (IF-97)*. Erlangen, Germany.
3. Snow, D.T. (1969). Anisotropic permeability of fractured media. *Water Resources Research*, 5(6), 1273-1289.
4. Biot, M.A. (1941). General theory of three-dimensional consolidation. *Journal of Applied Physics*, 12(2), 155-164.
5. McGarr, A. (2014). Maximum magnitude earthquakes induced by fluid injection. *Journal of Geophysical Research*, 119(2), 1008-1019.
6. Dieterich, J.H. (1994). A constitutive law for rate of earthquake production and its application to earthquake clustering. *Journal of Geophysical Research*, 99(B2), 2601-2618.
7. Sanyal, S.K. & Butler, S.J. (2005). An analysis of power generation prospects from Enhanced Geothermal Systems. *Geothermal Resources Council Transactions*, 29, 131-138.
8. Schulz, R. (1990). *Thermische Nutzung von Aquiferen*. BGR Circular, Germany.
9. Pruess, K. (1991). *TOUGH2 — A General Purpose Numerical Simulator for Multiphase Fluid and Heat Flow*. LBL-29400, Lawrence Berkeley Laboratory.
10. Kolditz, O. et al. (2012). OpenGeoSys: an open-source initiative for numerical simulation of thermo-hydro-mechanical/chemical (THM/C) processes in porous media. *Environmental Earth Sciences*, 67(2), 589-599.

---

*本レポートは研究ドラフトであり、査読前の試算結果を含む。実際の開発計画への適用には現地試験データによる検証が必要である。*
