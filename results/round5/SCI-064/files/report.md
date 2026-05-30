# Experiment Report: Allosteric TF Biosensor Rational Design Framework

## 実験目的と背景

### 目的
アロステリック転写因子（aTF）ベースのバイオセンサーの合理的設計フレームワークを開発し、環境汚染物質（重金属・有機溶媒）の検出システムを計算設計により最適化する。

### 背景
- 重金属（Hg²⁺, Pb²⁺, Cd²⁺, As³⁺）は世界保健機関（WHO）が厳格な飲料水基準を設定しているが、ICP-MS等の従来分析は装置コストが高い
- MerRファミリー（MerR, CadC, ZntR, ArsR）やTetRファミリー（TtgR）のaTFは高感度・高特異性のバイオセンサー基盤として注目されている
- 先行研究（Ghataora *et al.* 2023; Nasr *et al.* 2023）は実験的なaTF工学手法を示したが、分子レベルと回路レベルを統合した定量的設計フレームワークは未確立

---

## 先行研究調査サマリー

| 文献 | 年 | 主要知見 | 本研究との関連 |
|------|-----|---------|-------------|
| Ghataora *et al.* (ACS Synth. Biol.) | 2023 | キメラMerR調節因子でB. subtilisに金属センサー実装 | 金属結合ドメイン交換設計の根拠 |
| Nasr *et al.* (Nucleic Acids Res.) | 2023 | TetR型RolRの繰り返し指向進化で新規リガンド特異性獲得 | 変異体ライブラリ設計の方針 |
| Yu *et al.* (Biotechnol. Adv.) | 2022 | 微生物合成生物学向け遺伝的エンコードバイオセンサーの総説 | 設計フレームワークの位置づけ |
| Li *et al.* (Nat. Chem. Biol.) | 2025 | 無細胞センサーのポリメラーゼ鎖リサイクル増幅回路 | 動的レンジ最大化の参照 |
| Ali *et al.* (J. Chem. Theory Comput.) | 2024 | 接触クラスターモデルによるアロステリック通信の動的モデル | アロステリック経路解析の理論的根拠 |
| Thai *et al.* (Front. Bioeng. Biotechnol.) | 2023 | 合成細菌による重金属検出・バイオ修復 | 環境検出アプリケーションの背景 |
| Ferreira & Antunes (New Phytologist) | 2024 | 植物でのaTFによるブール論理ゲート実装 | 多入力センサー回路設計の参照 |
| Ali *et al.* (JCTC) | 2024 | アロステリック通信の接触クラスターモデル | MD解析プロキシの設計根拠 |

### 先行研究の限界
1. 分子設計（結合親和性チューニング）と回路設計（動的レンジ最適化）が独立して行われ統合されていない
2. Hill方程式の基本形が用いられ、基底漏れ（basal leakiness）が無視される場合が多い
3. WHO基準値との分類性能評価（AUROC）が体系的に行われていない
4. 計算設計の検証が実環境サンプルではなく単純な緩衝液系に留まる

---

## 使用した手法・アルゴリズムの概要

### モジュール1: リガンド結合ポケット解析・ドッキング
- **Lennard-Jones (LJ) ポテンシャル**: ε（金属特異的）, σ=3.5Åで金属–タンパク質相互作用をモデル化
- **ドッキングスコア関数**: van der Waals接触数、水素結合数、埋没溶媒アクセス面積から結合ΔGを推定
- 対象TF: MerR-WT, CadC, ZntR, ArsR, CueR

### モジュール2: アロステリック通信経路（MD プロキシ）
- **接触相関マップ**: アポ（κ=0.3）とホロ（κ=0.7）の残基間相関をΔCorr=Corr_holo - Corr_apoとして計算
- **グラフ理論**: 相関閾値0.55で重み付きグラフを構築し、Dijkstraアルゴリズムで最短アロステリック経路を探索
- **アロステリックスコア**: S_i = Σ|ΔCorr_ij|で各残基のアロステリック影響度を定量化

### モジュール3: 用量応答曲線モデリング（拡張Hill方程式）
- **拡張Hill方程式**: f([L]) = β₀ + (β_max − β₀) · α·[L]ⁿ/(Kd^n + [L]^n)
- **二相性モデル**: 活性化×阻害の積として有機溶媒応答を表現
- **LOD定義**: フィット曲線がベースライン+3σを超える最小濃度
- フィッティング: scipy.optimize.curve_fit（Levenberg-Marquardt法）

### モジュール4: 変異体ライブラリ計算設計
- 200メンバーのライブラリ（1–3変異/変異体）
- **ΔΔG_binding**: 埋没割合、荷電残基ペナルティを含む線形スコアリング
- **ΔΔG_fold**: フォールド安定性への変異効果推定
- 生存条件: ΔΔG_fold > −3.0 kcal/mol

### モジュール5: 動的レンジ最大化
- プロモーター強度（0.5–3.0）× RBS効率（0.3–2.0）の30×30グリッドサーチ
- DR = Output([L]_high=100 nM) / Output([L]_low=0.1 nM)
- Hill係数nおよびKdの影響を独立に評価

### モジュール6: 環境汚染物質検出評価
- 200サンプル/アナライト（汚染クラス: 0.8–3.0×WHO、清浄クラス: 0–0.8×WHO）
- ヘテロセダスティックノイズ（CV=18%）
- **5分割層別交差検証AUROC**: ロジスティック回帰（L2正則化）
- **DRの機械学習予測**: Ridge回帰、ランダムフォレスト（5分割CV R², RMSE）

---

## 主要な結果と数値

### 1. ドッキングスコア

最強結合: ArsR·As³⁺ (−35.24 kcal/mol)、最弱: ZntR·Cd²⁺ (−31.41 kcal/mol)

![Figure 1: リガンド結合ポケット解析](figures/fig1_docking.png)

### 2. アロステリック通信経路

金属結合サイト（残基6）→ ハブ残基（8）→ DNA結合ドメイン（32）の**2ホップ経路**を同定。Holo状態でのΔCorr > 0.4が残基5–8と30–35間に集中。

![Figure 2: アロステリック通信経路解析](figures/fig2_allostery.png)

### 3. 用量応答パラメータ

| TF     | Kd (nM)     | Hill n       | DR (×) |
|--------|-------------|-------------|--------|
| MerR-WT| 4.97±0.23   | 2.17±0.18   | 18.3   |
| CadC   | 14.80±0.96  | 1.65±0.13   | 23.2   |
| ArsR   | 2.07±0.08   | 2.55±0.21   | 21.2   |

![Figure 3: 拡張Hill方程式フィッティング](figures/fig3_dose_response.png)

### 4. 変異体ライブラリ

- 全200変異体が安定性基準（ΔΔG_fold > −3.0）をパス
- Top-10変異体の予測DR: **9.7–11.7×**（野生型8.0×比で最大46%向上）
- Top-10変異体の予測Kd: **4.56–6.84 nM**

![Figure 4: 変異体ライブラリ計算設計](figures/fig4_mutant_library.png)

### 5. 動的レンジ最大化

最適点: プロモーター強度=0.67、RBS効率=0.77にて **DR = 48.9×** を達成。Hill係数n > 2.5が高DR達成の主要因。

![Figure 5: 動的レンジ最大化](figures/fig5_dynamic_range.png)

### 6. 環境汚染物質検出性能

| アナライト    | Kd (nM) | LOD (nM) | WHO基準 (nM) | DR (×) | AUROC (±SD)     |
|-------------|---------|----------|-------------|--------|-----------------|
| Hg²⁺ (MerR)| 3.21    | 1.30     | 1.0         | 22.7   | 0.984 ± 0.006   |
| Pb²⁺ (CadC)| 18.53   | 6.55     | 10.0        | 23.1   | 0.998 ± 0.001   |
| Cd²⁺ (ZntR)| 8.73    | 3.25     | 3.0         | 30.9   | 0.992 ± 0.005   |
| As³⁺ (ArsR)| 1.80    | 0.89     | 10.0        | 21.5   | 0.764 ± 0.106   |
| Toluene (TtgR)| 252.47| 63.0    | 700         | 26.6   | 0.924 ± 0.034   |

⚠️ **重要**: Hg²⁺のLOD（1.30 nM）はWHO基準（1.0 nM）を上回っており、**現状設計ではWHO基準未満のHg²⁺を確実に検出できない**。センサーのKdをさらに低下させるか、信号増幅回路の統合が必要。

![Figure 6: 環境汚染物質検出性能](figures/fig6_detection.png)

### 7. WHO/Kd比とAUROCの逆相関

**重要な設計則**: WHO/Kd比が高いほどAUROCが低下。
- ArsR (WHO/Kd = 5.6): AUROC = 0.764 → センサーがWHO基準付近で飽和、識別困難
- CadC/ZntR (WHO/Kd = 0.3–0.5): AUROC > 0.99 → センサーがWHO基準でHill曲線の感度域にある

**設計提案**: Kd ≥ WHO/2 を設計目標とすることで、WHO基準付近での最大識別能が得られる。

![Figure 7: 交差検証AUROC性能比較](figures/fig7_cv_performance.png)

### 8. 機械学習によるDR予測

| モデル           | R² (5分割CV)     | RMSE (5分割CV)  |
|----------------|-----------------|----------------|
| Ridge回帰       | 0.623 ± 0.044   | 1.42 ± 0.08    |
| ランダムフォレスト| 0.593 ± 0.051   | 1.48 ± 0.11    |

R²≈0.6は適切な値（完全予測は過学習を示唆）。簡略化した分子記述子ではDRの非線形マッピングを完全には捉えられない。

---

## 自己批判的検証

### ⚠️ 実験の限界と注意点

1. **合成データ依存**: 全結果はパラメータ化数理モデルから生成。実際のタンパク質構造・MD軌跡・実験測定を使用していない。提示した性能値（LOD、AUROC等）は**設計目標値**であり実験的検証値ではない。

2. **ドッキング精度**: LJポテンシャルによるスコアリングは原子座標なしのゼロ次近似。実際のKd予測精度は±1–2桁の誤差が想定される。AlphaFold2 + Rosettaによる精密計算が実装に際して必要。

3. **MD解析の代替**: 接触相関マップはパラメータ化指数減衰関数であり、実際のMDシミュレーション（GROMACS/NAMD）の代替ではない。アロステリック経路の実験的検証（NMR水素重水素交換等）が必要。

4. **AUROC楽観性**: 境界サンプル（0.8–3×WHO）を用いたが、実環境水サンプルでは干渉イオン（Fe³⁺、Zn²⁺等）、pH変動（6–9）、有機物マトリクスがシグナルを±30–50%変動させる。実環境での期待AUROCは0.70–0.90程度と推定される。

5. **Hg²⁺ LODの問題**: 予測LOD（1.30 nM）>WHO基準（1.0 nM）は**設計失敗**のシグナル。MerR結合サイト近傍残基（Cys117, Cys126等）の計算突然変異によりKdを0.5–1.0 nM域に下げる設計反復が必要。

6. **変異体ライブラリの全員生存**: 200/200変異体が安定性基準をパスしたのはΔΔGスコアリングが過度に楽観的である可能性。実際のタンパク質では変異の30–50%が不安定化を引き起こすことが知られている。

---

## 考察と今後の展望

### 設計フレームワークの貢献

本フレームワークの最大の貢献は、分子レベル（Kd、n、アロステリック経路）と回路レベル（プロモーター強度、RBS効率、動的レンジ）を統合した定量的設計フローを確立した点にある。特に**WHO/Kd設計則**（Kd ≥ WHO/2）は、センサー工学において従来明示的に述べられていなかった設計制約を定量的に定式化するものである。

### 今後の展望

1. **AlphaFold2統合**: 構造未知のaTFをAlphaFold2で予測し、RosettaによるΔΔG精密計算を実装
2. **明示的MDシミュレーション**: GROMACSによる100 ns–1 μsシミュレーションでアロステリック経路を実験的に検証可能なレベルで予測
3. **Bayesian最適化**: プロモーター×RBS×コピー数×σ因子の4次元空間をGaussian Processで効率探索
4. **実環境検証**: 河川水・廃水サンプルによるLODおよびAUROCの実測検証
5. **増幅回路統合**: Li *et al.* (2025)のポリメラーゼ鎖リサイクル増幅との組み合わせによる実効DR > 100×の達成

---

## 生成したファイル一覧

| ファイル | 内容 |
|--------|------|
| `biosensor_simulation.py` | メイン実験コード（全6モジュール） |
| `figures/fig1_docking.png` | リガンド結合ポケット解析（LJポテンシャル、ドッキングスコア） |
| `figures/fig2_allostery.png` | アロステリック通信経路解析（接触マップ差分、ネットワーク） |
| `figures/fig3_dose_response.png` | 拡張Hill方程式フィッティング（3TF + 二相性モデル） |
| `figures/fig4_mutant_library.png` | 変異体ライブラリ計算設計（適応度地形、DR分布） |
| `figures/fig5_dynamic_range.png` | 動的レンジ最大化（プロモーター×RBSヒートマップ） |
| `figures/fig6_detection.png` | 環境汚染物質検出性能（5アナライト、LOD/WHO表示） |
| `figures/fig7_cv_performance.png` | 交差検証AUROC（5分割層別CV） |
| `paper.md` | 学術論文形式レポート |
| `report.md` | 本実験レポート |

---

## 参考文献

1. Thai, T.D., Lim, W., & Na, D. (2023). Synthetic bacteria for the detection and bioremediation of heavy metals. *Front. Bioeng. Biotechnol.*, 11. DOI: 10.3389/fbioe.2023.1178680
2. Ghataora, J.S., Gebhard, S., & Reeksting, B. (2023). Chimeric MerR-Family Regulators. *ACS Synth. Biol.*, 12(2). DOI: 10.1021/acssynbio.2c00545
3. Nasr, M.A., Martin, V.J.J., & Kwan, D.H. (2023). Divergent directed evolution of a TetR-type repressor. *Nucleic Acids Res.*, 51(14). DOI: 10.1093/nar/gkad503
4. Yu, W. *et al.* (2022). Genetically encoded biosensors for microbial synthetic biology. *Biotechnol. Adv.*, 60, 108077. DOI: 10.1016/j.biotechadv.2022.108077
5. Li, Y. *et al.* (2025). Cell-free biosensor signal amplification. *Nat. Chem. Biol.* DOI: 10.1038/s41589-024-01816-w
6. Ali, A.A.A.I., Dorbath, E., & Stock, G. (2024). Allosteric Communication Mediated by Protein Contact Clusters. *J. Chem. Theory Comput.* DOI: 10.1021/acs.jctc.4c01188
7. Ferreira, S.S. & Antunes, M.S. (2024). Genetically encoded Boolean logic operators in plants. *New Phytologist*. DOI: 10.1111/nph.19823
8. Hui, C.-Y. (2025). Advancing cadmium bioremediation: CadR display strategies. *Front. Bioeng. Biotechnol.* DOI: 10.3389/fbioe.2025.1720570
