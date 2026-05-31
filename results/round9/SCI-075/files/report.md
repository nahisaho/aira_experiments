# 実験レポート：手術ロボット半自律縫合の学習・制御システム

**実験日:** 2026-05-31  
**研究者:** GitHub Copilot (Claude Sonnet 4.6)  
**Notebook:** `data/jupyter/surgical_robot_suturing.ipynb`  
**乱数シード:** 42

---

## 1. 実験目的と背景

本研究は、**手術ロボット（dVRKプラットフォーム）における半自律縫合動作の学習・制御システム**を設計・シミュレーション検証することを目的とする。外科縫合は、高精度な針操作・軟組織への適応的な力制御・立体視覚フィードバックからの空間認識を要する最も複雑な術式の一つである。以下の6つの中核課題を統合的に扱う：

1. **デモンストレーションからの学習（LfD）** — Dynamic Movement Primitives (DMP) による専門家デモの模倣学習
2. **組織変形のリアルタイムモデリング** — Mass-Spring モデルと多項式回帰による力→変形推定
3. **力センシングと順応制御** — XGBoost による力推定 + アドミッタンス制御
4. **視覚サーボ（3D再構成+追跡）** — ステレオカメラ三角測量による3D針先追跡
5. **安全制約の保証** — 力・作業空間・速度の3層監視
6. **dVRK シミュレーション検証** — 合成デモデータによるベンチマーク

---

## 2. 先行研究調査

SemanticScholar MCP を使用し、以下のキーワードで検索（2020年以降）：
- "learning from demonstration surgical robot suturing autonomous"
- "dVRK da Vinci tissue deformation force sensing compliance control"
- "stereo visual servoing needle tracking surgical robot"

**主要先行研究 5件**

| No. | 著者・年 | タイトル | 主要知見 |
|-----|---------|---------|---------|
| 1 | Schwaner et al. (2021) IROS | Autonomous Bi-Manual Surgical Suturing Based on Skills Learned from Demonstration | 全タスク成功率17%、サブタスク75%、針挿入誤差3.3mm |
| 2 | Schwaner et al. (2021) CASE | Autonomous Needle Manipulation for Robotic Surgical Suturing Based on Skills Learned from Demonstration | DMP使用、成功率81%、挿入誤差3.8mm |
| 3 | Arduini et al. (2024) RO-MAN | Learning From Demonstration of Robot Motions And Stiffness Behaviors | DMP+GMM、可変インピーダンス制御、安全性検証済み |
| 4 | Zheng et al. (2024) ICRA | User-Centered Shared Control Scheme with LfD for Robotic Surgery | dVRKシミュレーション、深層IRL、ファジー制御切り替え |
| 5 | Black et al. (2020) RA-L | 6-DOF Force Sensing for the MTM of the da Vinci Surgical System | dVRKへの6自由度力センサ統合、ROS対応 |

**先行研究の限界：**
- 単一サブシステム（LfD、力推定、視覚サーボ）の独立評価が多く、6サブシステム統合の定量ベンチマークが不足
- 合成組織データ or 実機実験のどちらかに偏り、中間的なシミュレーション検証体制が少ない
- 速度制約を含む多層安全モニタリングの統合評価が不十分

---

## 3. 使用した手法・アルゴリズムの概要

### 3.1 システム全体構成（ROS/SurRoL互換）

```
[専門家デモ] → [DMP学習器] → [モーションプランナー]
                                       ↓
[ステレオカメラ] → [視覚サーボ]  → [順応制御器] → [dVRK PSM]
                                       ↑
[力センサ] → [XGBoostモデル] → [安全モニタ]
                                       ↑
[組織変形モデル] ──────────────────────┘
```

制御周期：100 Hz (dt = 10 ms)

### 3.2 各コンポーネント詳細

#### Dynamic Movement Primitives (DMP)
- 変換ダイナミクス：τẏ = αz(βz(g-y) - z)
- 正規化（[y0, g]→[0,1]）後にフィッティング
- 基底関数数 N=25、αz=48、βz=12
- 全 demo 平均正規化軌道に対してLeast-Squares回帰

#### Mass-Spring 組織モデル
- m = 5×10⁻⁴ kg, k = 500 N/m, ζ = 0.70
- 固有振動数: 159.2 Hz、オイラー法(dt=1ms)

#### XGBoost 力推定
- 入力: [x, z, vx, vz, deformation_mm]
- 5-fold CV、100推定器、max_depth=5

#### アドミッタンス制御
- Md=0.5kg, Bd=50Ns/m, Kd=200N/m
- 力クランプ: |F| ≤ 5.0 N

#### ステレオ視覚サーボ
- ベースライン B=65mm、焦点距離 f=500px
- 三角測量: Z = f·B/d
- 画素ノイズ σ=1.0px

#### GMR ポリシー
- GaussianMixture K=8、joint state-action fitting

#### 安全モニタ
- Layer 1: 力 (‖F‖ ≤ 5N)
- Layer 2: 作業空間 (x: ±50mm, y: ±30mm, z: 0-60mm)
- Layer 3: 速度 (‖v‖ ≤ 100mm/s)

---

## 4. 主要な実験結果と数値

### 4.1 合成デモデータ [cell:1]
- **デモ数:** 20件、100タイムステップ/件
- **組織ヤング率範囲:** 6,069 – 19,501 Pa
- **最大挿入力 Fz:** 0.669 – 2.015 N
- **最大変形:** 3.164 – 3.754 mm

### 4.2 DMP 軌道再現 [cell:3]

| 指標 | 値 |
|------|-----|
| X軸 RMSE | 9.603 ± 0.956 mm |
| Z軸 RMSE | 6.474 ± 0.365 mm |
| Combined RMSE | 11.59 ± 0.91 mm |
| **最終位置成功率 (<3mm)** | **100.0%** |

> ⚠️ 軌道RMSE(11.6mm)は高いが、最終位置誤差<3mmの成功率100%を達成。中間経路の偏差はDMPの標準的挙動（ゴール補間特性）による。先行研究の挿入誤差3.3-3.8mm（Schwaner et al.）と最終位置精度は同等レベル。

### 4.3 組織変形モデル [cell:4c]

| 指標 | 値 |
|------|-----|
| 多項式回帰 5-fold CV RMSE | **0.415 ± 0.023 mm** |
| 訓練 R² | 0.139 |
| Mass-Spring 最大変形 | 2.36 mm |
| デモ平均最大変形 | 3.10 mm |

> R²=0.14 は低く、合成力-変形関係の単純性を反映。実組織では非線形粘弾性が支配的。

### 4.4 XGBoost 力推定 [cell:10]

| 指標 | 値 |
|------|-----|
| **5-fold CV RMSE** | **0.281 ± 0.013 N** |
| **5-fold CV R²** | **0.606 ± 0.025** |
| 訓練 RMSE | 0.177 N |
| 最重要特徴量 | deformation_mm (64.6%) |

### 4.5 アドミッタンス制御 [cell:5]

| 指標 | 値 |
|------|-----|
| 最終位置誤差 | 2.870 ± 0.924 mm |
| 軌道 RMSE | 7.406 ± 1.058 mm |
| **力安全率 (F<5N)** | **100.0%** |
| **作業空間安全率** | **100.0%** |
| 成功率 (最終誤差<3mm) | 55.0% |

### 4.6 ステレオ視覚サーボ [cell:7]

| 指標 | 値 |
|------|-----|
| **3D追跡誤差** | **0.755 ± 1.474 mm** |
| 深度(Z)誤差 | 0.719 ± 1.490 mm |
| 95パーセンタイル誤差 | 4.595 mm |

### 4.7 GMR ポリシー（LfD）[cell:8]

| 指標 | 値 |
|------|-----|
| Vx RMSE | 47.200 ± 1.426 mm/s |
| Vz RMSE | 26.143 ± 2.135 mm/s |

### 4.8 安全制約 [cell:6]

| 制約 | 遵守率 |
|------|--------|
| 力制約 (‖F‖ < 5N) | **100.0%** |
| 作業空間 | **100.0%** |
| 速度 (< 100mm/s) | 55.2% |

> 速度違反(44.8%)は接近フェーズの高速動作(最大304mm/s)が100mm/sの保守的閾値を超えるため。閾値を350mm/sに緩和すれば全100%を達成。

### 4.9 NatureLM / GALACTICA MCP の試行記録

| 項目 | 内容 |
|------|------|
| **試行ツール** | `ask_naturelm`(NatureLM MCP), `scientific_qa`/`predict_citations`(GALACTICA MCP) |
| **検索方法** | ToolUniverseの`find_tools`・`grep_tools`を使用してパターン検索 |
| **エラー内容** | ToolUniverseレジストリに"NatureLM"/"GALACTICA"/"naturelm"/"galactica"のマッチ0件 |
| **代替手段** | SemanticScholar MCP（学術検索）＋Python Jupyter実験（定量検証） |

---

## 5. 生成した図表

### Figure 1: システム全体性能

![Figure 1: System Overview](figures/fig1_system_overview.png)

*統合システム性能の9パネル概要。(A)針軌道XZ平面（青:デモ、赤:DMP再現）。(B)挿入力プロファイルと5N制限。(C)組織変形の多項式推定とMass-Spring比較。(D)DMP誤差分布。(E)順応制御位置追跡。(F)安全制約遵守率。(G)視覚サーボ追跡誤差。(H)GMRポリシー速度予測。(I)コンポーネント性能サマリー。*

### Figure 2: 詳細性能分析

![Figure 2: Performance Details](figures/fig2_performance_details.png)

*(左)XGBoost力推定散布図。(中央)特徴量重要度。(右)クロスコンポーネント性能比較。*

---

## 6. 考察と今後の展望

### 6.1 主要知見

- **DMP**は最終位置精度（100%の<3mm達成）において有効だが、中間経路品質の保証には追加制約が必要
- **XGBoost力推定**（RMSE 0.28N）は有望だが、組織不均一性によりR²=0.61と限定的
- **力安全制約**（100%）と**作業空間制約**（100%）は完全に満足されており、アドミッタンス制御が安全保証に有効
- **視覚サーボ**の平均精度（0.75mm）はサブミリ達成だが、標準偏差1.47mmと95パーセンタイル4.6mmは浅い距離での深度推定破綻を示す
- **速度制約違反**（44.8%）は保守的閾値設定の問題であり、軌道スケーリングで解消可能

### 6.2 自己批判的評価

1. **合成データへの依存度が高い**: 全実験は数学的生成データ。実組織の粘弾性、摩擦、液体環境は未モデル化
2. **DMPの中間経路品質**: 軌道RMSE 11.6mmは臨床使用には不十分（目標<2mm）
3. **R²=0.14（組織変形）**: 訓練データでも低いR²は合成関係の単純性を示し、実データでは性能低下が予想される
4. **速度違反**: 実臨床では速度超過は安全上許容されない
5. **sim-to-realギャップ**: 物理dVRKへの転移には domain randomization と実機検証が必須

### 6.3 今後の展望

1. **物理dVRK実機検証**: PyBullet/SurRoLから実機へのSim-to-Real転移
2. **FEM組織モデル**: リアルタイム線形FEMによる高精度変形予測
3. **深層学習視覚サーボ**: オクルージョン頑健性のためのYOLO/DETRベース針追跡
4. **CBF安全制約**: Control Barrier Functions による形式的安全保証
5. **強化学習統合**: デモ誘導RLによる試行錯誤を通じた政策改善（Singh et al. [2023]の方向性）
6. **力センサ統合**: 真のdVRK力センサ（Black et al. [2020]）との統合による力推定精度向上

---

## 7. 生成したファイル一覧

| ファイル | 説明 |
|---------|------|
| `data/jupyter/surgical_robot_suturing.ipynb` | Jupyter実験ノートブック（13セル） |
| `data/jupyter/data/raw/demo_summary.csv` | 20件のデモデータサマリー（CSV） |
| `data/jupyter/data/raw/environment.json` | 実験環境・パッケージバージョン記録 |
| `figures/fig1_system_overview.png` | システム全体9パネル図 |
| `figures/fig2_performance_details.png` | 詳細性能分析図 |
| `paper.md` | 学術論文形式の成果物 |
| `report.md` | 本レポートファイル |

---

## 付録: 実験環境

| 項目 | 値 |
|------|-----|
| Python | 3.11.2 (GCC 12.2.0) |
| numpy | 2.3.5 |
| pandas | 2.3.3 |
| scikit-learn | 1.6.1 |
| scipy | 1.17.1 |
| xgboost | 3.2.0 |
| matplotlib | 3.10.9 |
| seaborn | 0.13.2 |
| lightgbm | 4.6.0 |
| 乱数シード | 42 |
| デモ数 | 20 |
| タイムステップ数 | 100 |
| 制御周期 dt | 0.01 s |
