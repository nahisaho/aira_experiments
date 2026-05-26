# マイクロキネティックモデリングフレームワーク：不均一系触媒反応への適用
## Fischer-Tropsch合成のケーススタディ

---

## 1. 実験目的と背景

不均一系触媒反応のマイクロキネティックモデリングは、DFT計算から得られるエネルギー論に基づき、触媒反応の速度論を定量的に予測するための重要な手法である。本研究では、以下の機能を統合したPythonベースのマイクロキネティックモデリングフレームワークを開発した：

1. **DFTエネルギーからの速度定数算出**（遷移状態理論＋Wignerトンネル効果補正）
2. **吸着等温線モデル**（Langmuir/Temkin/フラクタル表面）
3. **反応速度支配段階の自動同定**（Degree of Rate Control, DRC解析）
4. **被覆率依存性**（lateral interaction を考慮した平均場モデル）
5. **反応器モデルとの連成**（PFR/CSTR）
6. **Fischer-Tropsch合成**のケーススタディ（Co(0001)上）

本フレームワークは CatMAP/Cantera/OpenMKM の設計思想を参考にしつつ、モジュラーな Python 実装により、各要素を独立に拡張可能な構成とした。

---

## 2. 使用した手法・アルゴリズムの概要

### 2.1 速度定数の算出

遷移状態理論（TST）に基づくEyring式：

$$k = \frac{k_B T}{h} \exp\left(-\frac{E_a}{k_B T}\right)$$

Wignerトンネル効果補正：

$$\kappa = 1 + \frac{1}{24}\left(\frac{h\nu^{\ddagger}}{k_B T}\right)^2$$

### 2.2 吸着等温線モデル

- **Langmuir**: $\theta = KP/(1+KP)$
- **Temkin**: $K(\theta) = K_0 \exp(-\alpha\theta)$
- **フラクタル表面**: $\theta = (KP)^{d/D}$

### 2.3 DRC解析

Campbell の Degree of Rate Control：

$$X_{RC,i} = \frac{k_i}{r}\left(\frac{\partial r}{\partial k_i}\right)_{K_{eq,i}, k_{j\neq i}}$$

### 2.4 Lateral Interaction

平均場近似による被覆率依存性：

$$\Delta E_i = z \sum_j \epsilon_{ij} \theta_j$$

BEP関係を通じた活性化エネルギー修正：

$$E_a(\theta) = E_{a,0} + \alpha_{BEP} \cdot \Delta E$$

### 2.5 反応器モデル

- **PFR**: $dF_i/dW = r_i \cdot \rho_{sites}$
- **CSTR**: $0 = F_{i,0} - F_i + r_i \cdot \rho_{sites} \cdot W$

---

## 3. 主要な結果と数値

### 3.1 Arrhenius プロットとトンネル効果

![Arrhenius plots for FT elementary steps](figures/arrhenius_plots.png)

**図1**: (a) FT素反応のArrhenius プロット。実線はWignerトンネル補正付き、破線は古典TST。CO解離が最も高い活性化障壁（1.60 eV）を持つ。(b) Wignerトンネル補正係数。低温域でC水素化（ν‡=1100 cm⁻¹）の補正が最大となる。

### 3.2 吸着等温線の比較

![Adsorption isotherm comparison](figures/adsorption_isotherms.png)

**図2**: (a) Langmuir等温線：K値増加に伴い飽和被覆率への到達が速くなる。(b) Temkin等温線：α増加により高被覆率域での吸着が抑制される。(c) フラクタル表面等温線：フラクタル次元Dの増加に伴い中間被覆率が増加する。

### 3.3 Lateral Interaction

![Lateral interaction analysis](figures/lateral_interactions.png)

**図3**: (a) CO結合エネルギーの被覆率依存性（3つのモデル比較）。piecewiseモデルでは高被覆率域で急激な反発増大を表現。(b) Co(0001)上のFT吸着種間の lateral interaction matrix。O*-O* 間の反発（0.15 eV）が最も強い。

### 3.4 表面被覆率の解析

![Surface coverage analysis](figures/surface_coverages.png)

**図4**: T=500K, P=20 bar, H₂/CO=2 における (a) 被覆率の時間発展と (b) 定常状態被覆率。CO*が支配的（θ_CO=0.613）、H*=0.160、空きサイト *=0.227。

**定常状態被覆率（T=500K）**:

| 吸着種 | 被覆率 θ |
|--------|----------|
| CO*    | 0.6133   |
| H*     | 0.1595   |
| C*     | ~10⁻⁹   |
| O*     | ~10⁻⁸   |
| *      | 0.2272   |

### 3.5 温度依存性

![Temperature study results](figures/temperature_study.png)

**図5**: (a) ターンオーバー頻度（TOF）の温度依存性：見かけの活性化エネルギー約100 kJ/mol。(b) PFRとCSTRでのCO転化率比較：500Kで PFR=36.6%, CSTR=39.0%。(c) 主要吸着種の温度依存性：高温でCO*が減少しH*が増加。(d) CH₄選択率が温度とともに変化。

### 3.6 速度支配段階の同定

![Degree of Rate Control analysis](figures/degree_of_rate_control.png)

**図6**: DRC解析結果。CO解離反応の X_RC = 1.000 であり、明確な律速段階として同定された。これはCo触媒上FTの既知の実験結果と一致する。

**DRC値**:

| 素反応 | X_RC |
|--------|------|
| CO解離 | 1.000 |
| OH+H   | 0.000 |
| その他 | ~0.000 |

### 3.7 PFR濃度プロファイル

![PFR concentration profiles](figures/pfr_profiles.png)

**図7**: (a) PFR内のモル流量プロファイル。COとH₂が触媒重量に沿って消費され、CH₄とH₂Oが生成される。(b) CO転化率プロファイル：触媒出口で36.6%に到達。

### 3.8 Lateral Interaction の効果

![Lateral interaction effect on rates](figures/lateral_interaction_effect.png)

**図8**: (a) lateral interaction の有無による反応速度の比較。lateral interaction により高被覆率域で速度が低下する。(b) 見かけの活性化エネルギーの変化：lateral interaction により Ea_app が増加する。

### 3.9 エネルギーダイアグラム

![Potential energy surface](figures/energy_diagram.png)

**図9**: Co(0001)上FTのポテンシャルエネルギー面。CO解離の遷移状態（TS₁, +0.25 eV）が最高エネルギー点であり、律速段階であることを示す。

---

## 4. 考察と今後の展望

### 4.1 主要な知見

1. **CO解離がFTの律速段階**: DRC解析により X_RC(CO diss.)=1.0 と確認。これは Motagamwala & Dumesic (2021) や Xie et al. (2022) の報告と一致する。

2. **表面はCO*が支配的**: θ_CO=0.61 は Co(0001) 上の実験値（0.4-0.7）の範囲内。

3. **Lateral interaction の影響**: 見かけの活性化エネルギーが増大し、高温域での反応速度予測に重要な補正となる。

4. **反応器比較**: CSTR がPFR よりやや高い転化率を示す。これは CO解離が律速でありCO分圧低下が有利に働くためである。

### 4.2 フレームワークの特長

- **モジュラー設計**: 各モジュール（速度定数、吸着、lateral interaction、反応器）が独立に交換可能
- **解析的PSS**: 極めて剛性の高い系に対して安定な定常状態解を提供
- **DRC自動解析**: 任意の反応ネットワークに対して律速段階を自動同定

### 4.3 今後の展望

1. **kMC連成**: 平均場近似を超えるkMCシミュレーションとの統合
2. **機械学習加速**: Neural Network Potential による lateral interaction の高速予測
3. **多サイトモデル**: terrace/step/kink サイトの区別
4. **長鎖生成物**: Anderson-Schulz-Flory分布の組み込み
5. **実験データとの直接比較**: TPD/TPR スペクトルのシミュレーション

---

## 5. 生成したファイル一覧

### ソースコード
| ファイル | 説明 |
|----------|------|
| `src/__init__.py` | パッケージ初期化 |
| `src/rate_constants.py` | TST速度定数 + トンネル効果補正 |
| `src/adsorption.py` | Langmuir/Temkin/フラクタル等温線 |
| `src/rate_control.py` | DRC解析・感度解析 |
| `src/lateral_interactions.py` | 平均場lateral interactionモデル |
| `src/reactor_models.py` | PFR/CSTR反応器モデル + MKMソルバー |
| `src/fischer_tropsch.py` | FT合成ケーススタディ |
| `run_simulation.py` | メインシミュレーション実行スクリプト |

### 生成図
| ファイル | 内容 |
|----------|------|
| `figures/arrhenius_plots.png` | Arrheniusプロット + トンネル補正 |
| `figures/adsorption_isotherms.png` | 吸着等温線比較 |
| `figures/lateral_interactions.png` | Lateral interaction解析 |
| `figures/surface_coverages.png` | 表面被覆率 |
| `figures/temperature_study.png` | 温度依存性研究 |
| `figures/degree_of_rate_control.png` | DRC解析 |
| `figures/pfr_profiles.png` | PFR濃度プロファイル |
| `figures/lateral_interaction_effect.png` | Lateral interaction効果 |
| `figures/energy_diagram.png` | ポテンシャルエネルギー面 |

### ドキュメント
| ファイル | 内容 |
|----------|------|
| `report.md` | 実験レポート（本ファイル） |
| `paper.md` | 学術論文形式の文書 |
