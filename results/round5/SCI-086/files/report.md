# 実験レポート：患者個別心臓デジタルツイン（CardioTwin）フレームワーク

---

## 1. 実験目的と背景

### 目的
本実験は、心臓MRIから始まり心房細動アブレーション効果予測に至るまでの、患者固有の心臓デジタルツイン（CDT）を構築するための統合的なコンピューテーショナルフレームワーク「CardioTwin」を設計・実証することを目的とする。

### 背景
- 心臓不整脈（心室細動、心房細動等）は世界で年間300万人以上が死亡する主要疾患
- 現在の臨床リスク層別化（LVEF、QRS幅）は不十分
- 計算科学的心臓デジタルツインは個別化医療の実現に向けた革新的アプローチ
- OpenCARP（電気生理学）とFEBio（力学）を基盤とした統合フレームワークが必要

---

## 2. 先行研究調査結果

ToolUniverse MCP（Semantic Scholar, OpenAlex, Crossref）を用いた文献調査により、以下の主要論文を特定した。

### 主要先行研究

| # | タイトル | 著者 | 年 | DOI | 主要知見 |
|---|---------|------|----|-----|---------|
| 1 | A Framework for digital twins of cardiac EP from 12-leads ECGs | Camps et al. | 2021 | 10.1016/j.media.2021.102080 | 12誘導ECGのみからCDTを非侵襲的に生成する初の統合フレームワーク |
| 2 | A comprehensive model of whole human heart electromechanics | Piersanti et al. | 2023 | 10.1016/j.cma.2023.115983 | 心房収縮を含む4腔全心臓電気力学モデルの初の実用実装 |
| 3 | Personalized ablation vs. conventional ablation in AF | Azzolin et al. | 2022 | 10.1093/europace/euac116 | 29例患者特異的デジタル双子で個別化HDF戦略が標準PVIより優れることを実証 |
| 4 | Computational modeling of cardiac arrhythmogenesis | Trayanova et al. | 2023 | 10.1152/physrev.00017.2023 | 心臓計算モデルから臨床応用への包括的レビュー（引用数98件） |
| 5 | Solving the Inverse Problem of ECG for CDTs: A Survey | Li et al. | 2024 | 10.1109/rbme.2024.3486439 | 決定論的・確率論的ECG逆問題手法の体系的レビュー |
| 6 | PINNs for Cardiac Activation Mapping | Sahli Costabal et al. | 2020 | 10.3389/fphy.2020.00042 | 物理インフォームドNNで心房内活性化マップを精度よく推定（引用数400件） |
| 7 | Precision medicine in human heart modeling | Peirlinck et al. | 2021 | 10.1007/s10237-021-01421-z | 心臓精密医療のための人口ライブラリ型パーソナライゼーション戦略 |
| 8 | ML for AF driver location prediction | Luongo et al. | 2021 | 10.1016/j.cvdhj.2021.03.002 | 機械学習で12誘導ECGからAFドライバー位置を予測（感度73.9%、特異度82.6%） |
| 9 | Whole-heart EM-driven CFD simulation | Zingaro et al. | 2024 | 10.1016/j.jcp.2024.112885 | 電気力学-流体連成シミュレーションで左脚ブロック病態の血行動態変化を再現 |
| 10 | Cyclical FIM for arrhythmia in openCARP | Barrios Espinosa et al. | 2025 | 10.1007/s00366-024-02094-9 | DREAM法でモノドメインより87倍高速な不整脈シミュレーションを実現 |

### 先行研究の課題・限界
1. **統合性の欠如**: 多くの研究は電気生理学 OR 力学のどちらかに特化し、両者の連成が不完全
2. **臨床コホートの小ささ**: Azzolin et al. は n=29 と限られており、外的妥当性が未確認
3. **計算コスト**: フル3Dモノドメインシミュレーション（500,000要素、dt=0.02ms）は48コアで4〜8時間
4. **パラメータ同定可能性**: ECG逆問題は本質的に不適切問題（ill-posed）であり、一意的な解が存在しない
5. **実世界への転換**: 合成データ/シミュレーション由来のモデルは臨床データに適用した際の性能劣化が大きい

---

## 3. 実験計画

### フレームワーク設計（6モジュール）

```
Module 1: 心臓MRIセグメンテーション + メッシュ生成
    ↓
Module 2: 電気生理学シミュレーション（openCARP）
    - Aliev-Panfilov モデル（高速、2変数）
    - ten Tusscher-Panfilov TP06（詳細イオンチャネル）
    ↓
Module 3: 電気-力学連成モデル（FEBio）
    - 能動応力定式化
    - Holzapfel-Ogden 受動弾性
    ↓
Module 4: 患者固有パラメータ逆問題推定（MCMC）
    ↓
Module 5: 不整脈リスク評価
    - APD restitution 解析
    - 脆弱性ウィンドウ評価
    ↓
Module 6: 心房細動アブレーション予測
    - PersonAL（個別化アブレーション）戦略
```

### 評価計画
- 合成患者コホート n=200 による5分割交差検証
- 評価指標: AUROC, F1, nRMSE
- 批判的自己評価: 合成データ前提への依存度を明示

---

## 4. 実験結果

### 4.1 フレームワーク概要図

![Figure 1: CardioTwin Framework](figures/fig1_framework_overview.png)

CardioTwinフレームワークの全体構造。6つのモジュールが心臓MRI画像入力から臨床アウトカム予測まで接続される。

---

### 4.2 活動電位モデル

![Figure 2: Action Potential Models](figures/fig2_action_potentials.png)

**左パネル（Aliev-Panfilov）**: 2変数モデルによる正規化活動電位。除分極・再分極ダイナミクスを効率的に再現。

**右パネル（ten Tusscher）**: 心内膜、中層筋、心外膜の膜電位波形。APD₉₀の透壁不均一性（内膜>外膜）が臨床的に観察される値と一致。

| 細胞層 | APD₉₀ (ms) | 文献値 (ms) |
|-------|-----------|-----------|
| 心内膜 | 275 | 260–290 |
| 中層筋 | 285 | 275–320 |
| 心外膜 | 265 | 240–270 |

---

### 4.3 2D電気伝播シミュレーション

![Figure 3: 2D Propagation](figures/fig3_propagation.png)

60×60グリッド（dx=0.025 cm）でのAliev-Panfilov伝播シミュレーション。コーナーから開始した刺激が組織全体に伝播する様子を4時点で示す。

- **見掛け伝導速度**: ~0.47 m/s（等方性）
- **文献値**: 横断方向 0.3–0.5 m/s ✓
- 渦状再入波（reentry）のシミュレーションは次フェーズで実施予定

---

### 4.4 電気-力学連成モデル

![Figure 4: Electro-Mechanical Coupling](figures/fig4_electromechanical.png)

3種の疾患表現型に対する能動応力波形（左）、左室圧-容積ループ（中央）、および透壁線維応力分布（右）。

**LV圧-容積ループ結果**:

| 表現型 | EDV (mL) | ESV (mL) | EF (%) | 最大能動応力 (kPa) |
|-------|---------|---------|--------|-----------------|
| 正常 | 120 | 45 | 62.5 | 50.0 |
| 拡張型心筋症 | 160 | 104 | 35.0 | 35.0 |
| 肥大型心筋症 | 100 | 30 | 70.0 | 58.0 |

心内膜（endocardium）での線維応力が心外膜より約70%高く、ラプラス則と一致する透壁応力分布が再現された。

---

### 4.5 逆問題：パラメータ推定

![Figure 5: Inverse Problem MCMC](figures/fig5_inverse_problem.png)

MCMC（Metropolis-Hastings法）によるベイズパラメータ推定の事後分布。バーンイン1,000サンプル後の4,000サンプルを使用。

**推定結果**:

| パラメータ | 真値 | 推定平均 | 推定SD | 相対誤差 (%) |
|---------|-----|-------|------|-----------|
| σ_t (cm/ms) | 0.0800 | 0.0793 | 0.0041 | 0.9 |
| g_Na (nS/pF) | 14.80 | 14.67 | 0.78 | 0.9 |
| g_CaL (μS/μF) | 3.98×10⁻⁵ | 3.91×10⁻⁵ | 2.1×10⁻⁶ | 1.8 |
| g_Kr (nS/pF) | 0.153 | 0.149 | 0.011 | 2.6 |

**nRMSE（5分割CV）: 0.044 ± 0.010**

> ⚠️ **重要注記**: これらの誤差は楽観的であり、推定に使用したフォワードモデルがデータ生成モデルと同一であるため、完全な同定可能性が保証されている。実臨床ECGデータでは、モデルミスマッチによりRMSEが5〜15倍増加することが予想される。

---

### 4.6 不整脈リスク評価 + 心房細動アブレーション

![Figure 6: Arrhythmia Risk and AF Ablation](figures/fig6_arrhythmia_ablation.png)

**APD restitutionカーブ**: 心不全では傾きdAPD/dDI > 1の臨界DIが通常の<10 msから35 msまで拡大し、alternansリスクが増大。

**アブレーション成功率（線維化度別）**:

| 線維化 (%) | PVI成功率 (%) | PersonAL成功率 (%) | 改善幅 (%) |
|---------|------------|-----------------|----------|
| 0–10 | 68.2 ± 8.1 | 76.4 ± 7.3 | +8.2 |
| 10–20 | 58.6 ± 9.4 | 70.1 ± 8.8 | +11.5 |
| 20–30 | 44.3 ± 11.2 | 61.8 ± 10.5 | +17.5 |
| 30–40 | 31.7 ± 12.8 | 54.2 ± 11.9 | +22.5 |

線維化が高いほど個別化戦略（PersonAL）の優位性が顕著（最大+22.5ポイント改善）。

---

### 4.7 交差検証パフォーマンス

![Figure 7: Cross-Validation Performance](figures/fig7_cv_performance.png)

**5分割交差検証結果（n=200 仮想患者）**:

| タスク | 評価指標 | Mean ± SD |
|------|---------|-----------|
| 不整脈リスク予測 | AUROC | **0.902 ± 0.026** |
| リスク分類 | F1-score | **0.830 ± 0.031** |
| パラメータ推定 | nRMSE | **0.044 ± 0.010** |
| AFアブレーション予測 | AUROC | **0.811 ± 0.020** |

**学習曲線の解釈**: 臨床的に許容可能なAUROC 0.80に到達するには約80〜100名の学習データが必要。

---

## 5. 考察と批判的検証

### 5.1 合成データへの依存

本実験の最大の限界は、すべての定量的性能指標が**合成的に生成された仮想患者**から得られていることである。アウトカムラベル（不整脈有無）は同じ数理モデル（ロジスティック回数）から生成されており、実質的に学習-評価間での*循環依存*が存在する。

実臨床データへの適用を想定した場合の性能劣化の要因：
1. 生物学的多様性（イオンチャネル密度の個人差、繊維化の空間的不均一性）
2. 画像セグメンテーション誤差（MRI解像度限界、呼吸・心拍アーチファクト）
3. ECG逆問題の非一意性（異なるパラメータ組合せが同一ECGを生成）
4. 時間的安定性の欠如（心室リモデリング、自律神経変調）

**推定実世界 AUROC: 0.70–0.78**（Azzolin 2022, Luongo 2021 の臨床検証値に基づく）

### 5.2 モデル仮定への依存

| モデル仮定 | 実験での設定 | 実世界での影響 |
|----------|-----------|-------------|
| 等方性拡散 (AP) | D = 0.001 cm²/ms | 実際の線維走向依存異方性比は10:1 |
| 固定パラメータ (TP06) | 一様な心室壁 | LGEで検出される線維化では局所的GK/GNa変化 |
| 2D伝播 | 60×60グリッド | 3D解剖学的構造・Purkinje系が必要 |
| 単純フォワードモデル | 代数的ECG特徴量 | 完全ECG合成には前胸部誘導の空間統合が必要 |

### 5.3 自己批判的評価のまとめ

> ⚠️ **本フレームワークの現状は、臨床使用より研究・方法論開発目的**に適している。AUROC 0.90 という数値は合成データ上での評価に過ぎず、実患者への直接適用を保証するものではない。臨床有用性を確立するには、前向き多施設コホート試験による外的妥当性検証が不可欠である。

---

## 6. 今後の展望

1. **臨床コホート検証**: 心房細動アブレーション患者100名以上による前向きコホート研究（12ヶ月フォローアップ）
2. **代理モデル加速**: Neural Operator（FNO, DeepONet）によるEPシミュレーションの秒単位加速
3. **マルチモーダルデータ統合**: 心臓MRI + 12誘導ECG + 体表面電位マッピング + LGE線維化
4. **不確実性定量化**: パラメータ不確実性から臨床アウトカム不確実性までのベイズ伝播
5. **規制対応**: FDA/CEマーク取得に向けたASME V&V 40に準拠したデジタルツイン検証戦略

---

## 7. 生成ファイル一覧

| ファイル名 | 説明 | サイズ |
|----------|-----|------|
| `figures/fig1_framework_overview.png` | CardioTwinフレームワーク概要図 | 113 KB |
| `figures/fig2_action_potentials.png` | Aliev-PanfilovとTP06活動電位モデル比較 | 108 KB |
| `figures/fig3_propagation.png` | 2D電気伝播シミュレーション（4時点） | 56 KB |
| `figures/fig4_electromechanical.png` | 電気-力学連成：能動応力・PVループ・透壁応力 | 169 KB |
| `figures/fig5_inverse_problem.png` | MCMCベイズパラメータ推定の事後分布 | 145 KB |
| `figures/fig6_arrhythmia_ablation.png` | APD restitution + AFアブレーション成功率 | 201 KB |
| `figures/fig7_cv_performance.png` | 5分割CV性能評価 + 学習曲線 | 134 KB |
| `paper.md` | 学術論文形式の詳細報告書 | 29 KB |
| `report.md` | 本実験レポート | — |

---

## 参考文献

1. Camps, J. et al. (2021). Medical Image Analysis, 71, 102080. https://doi.org/10.1016/j.media.2021.102080
2. Piersanti, R. et al. (2023). Comput. Methods Appl. Mech. Eng., 410, 115983. https://doi.org/10.1016/j.cma.2023.115983
3. Azzolin, L. et al. (2022). EP Europace, euac116. https://doi.org/10.1093/europace/euac116
4. Trayanova, N.A. et al. (2023). Physiological Reviews. https://doi.org/10.1152/physrev.00017.2023
5. Li, L. et al. (2024). IEEE Rev. Biomed. Eng. https://doi.org/10.1109/rbme.2024.3486439
6. Sahli Costabal, F. et al. (2020). Frontiers in Physics, 8, 42. https://doi.org/10.3389/fphy.2020.00042
7. Luongo, G. et al. (2021). Cardiovasc. Digital Health J. https://doi.org/10.1016/j.cvdhj.2021.03.002
8. Zingaro, A. et al. (2024). J. Comput. Phys., 504, 112885. https://doi.org/10.1016/j.jcp.2024.112885
9. Barrios Espinosa, C. et al. (2025). Eng. Comput. https://doi.org/10.1007/s00366-024-02094-9
10. Peirlinck, M. et al. (2021). Biomech. Model. Mechanobiol. https://doi.org/10.1007/s10237-021-01421-z
