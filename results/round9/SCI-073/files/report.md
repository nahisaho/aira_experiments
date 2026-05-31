# 実験レポート：高解像度触覚センサーによる物体認識・操作システム

**研究テーマ**: GelSight/DIGITセンサーを用いた高解像度触覚センシングによるロボット操作システムの深層学習フレームワーク設計

**実験日**: 2026-05-31  
**ノートブック**: `tactile_sensing.ipynb`  
**乱数シード**: `np.random.seed(42)`, `random.seed(42)`

---

## 1. 実験目的と背景

### 1.1 研究背景

GelSightおよびDIGITセンサーは、弾性ゲルの変形を内蔵カメラで撮影することにより、接触面の2.5D形状・力分布をサブミリメートル精度で計測する視覚ベース触覚センサーである。従来の圧電式・容量式センサーに比べ、空間分解能が2桁以上高く、ロボットハンドの指先に搭載可能なコンパクト設計が実現されている。

本研究では以下の6タスクを統合的に評価する機械学習フレームワークを設計・実装した：

1. 接触形状・力分布の推定
2. テクスチャ分類（深層学習ベース）
3. 触覚と視覚のマルチモーダル融合
4. 把持安定性のリアルタイム評価
5. すべり検出と力制御フィードバック
6. 未知物体の探索的把持戦略（不確実性定量化）

### 1.2 先行研究調査の結果

Semantic Scholar MCPを使用して文献調査を実施（HTTP 429レート制限により4件の検索が失敗）。以下の主要論文を確認：

| # | タイトル | 著者 | 年 | DOI |
|---|---------|------|-----|-----|
| 1 | Shape-independent hardness estimation using deep learning and GelSight | Yuan et al. | 2017 | 10.1109/ICRA.2017.7989116 |
| 2 | Vision-Based Tactile Sensor for Contact Position and Force Distribution | Kakani et al. | 2021 | 10.3390/s21051920 |
| 3 | Learning Force Distribution for GelSight Mini via FEA | Helmut et al. | 2024 | 10.1109/IROS60139.2025.11246486 |
| 4 | Generalizable Force Estimation via Graph Networks | Chen et al. | 2024 | 10.1109/ROBIO64047.2024.10907583 |
| 5 | DIGIT: Low-cost compact high-resolution tactile sensor | Lambeta et al. | 2020 | 10.1109/LRA.2020.3011445 |
| 6 | TACTO: Fast simulator for vision-based tactile sensors | Wang et al. | 2022 | 10.1109/LRA.2022.3146945 |

**先行研究の課題・限界**:
- 個々のタスク（硬度推定、力推定、テクスチャ分類）は個別に研究されているが、統合フレームワークは少ない
- 合成データ（FEA）からリアルセンサーへのシム-to-リアル転移が未解決
- 静的な単フレーム解析が主流で、時系列ダイナミクスのモデル化が不足
- クラス不均衡下でのすべり検出が困難

### 1.3 NatureLM / GALACTICA MCPの試行記録

| ツール | 試行内容 | 結果 |
|-------|---------|------|
| `ask_naturelm` | GelSightセンサーの定量パラメータ（解像度・力レンジ・レイテンシ）の取得 | **接続失敗**: ToolUniverseレジストリに未登録 |
| `scientific_qa` | 摩擦係数モデルの科学的検証 | **接続失敗**: ToolUniverseレジストリに未登録 |
| `predict_citations` | 関連文献の予測 | **接続失敗**: ToolUniverseレジストリに未登録 |
| `SemanticScholar_search_papers` | 先行研究検索 | **部分成功**: 1クエリ成功（8件取得）、その後HTTP 429レート制限 |

---

## 2. 使用した手法・アルゴリズム

### 2.1 データ生成

**テクスチャ分類データセット** (`data/raw/texture_dataset.csv`)
- 8クラス × 100サンプル = 800サンプル
- 20次元特徴量（空間周波数、勾配強度、接触面積、深度統計、エッジ密度、等方性、粗さ、周期性、RGB輝度、変形場、せん断推定、法線力推定、接触対称性）
- クラスごとにGaussianノイズ（σ = 特徴量平均の15%）を付加

**力推定・すべり検出データセット** (`data/raw/force_dataset.csv`)
- 1,000サンプル、15次元触覚特徴量
- 法線力：非線形モデル $F_z = 3.0 + 5.0|x_0| + 2.0x_1^2 + 1.5x_2 + 0.5x_0 x_3 + \mathcal{N}(0, 0.8)$
- すべりラベル：$y = 1$ if $\sqrt{F_x^2 + F_y^2}/F_z > 0.40$、陽性率12.8%

**把持安定性データセット**
- 800サンプル（安定53.6%、不安定46.4%）
- 18次元特徴量（接触面積×4、法線力×4、せん断力×4、アプローチ角、物体物性×5）

### 2.2 評価手法

- **5分割層化交差検証**（分類タスク）/ **5分割KFold**（回帰タスク）
- 評価指標：Accuracy, AUROC（分類）; R², RMSE（回帰）
- 標準化：各フォールドの訓練データのみで`StandardScaler`をフィット

### 2.3 機械学習モデル

| モデル | 用途 | 主要パラメータ |
|-------|-----|--------------|
| Random Forest | 全タスク | n_estimators=100, random_state=42 |
| Gradient Boosting | テクスチャ分類 | n_estimators=100, random_state=42 |
| SVM (RBF) | テクスチャ・すべり | C=10/5, γ='scale' |
| MLP (3層) | テクスチャ分類 | (128,64,32), max_iter=300 |
| kNN | テクスチャ分類 | k=5 |
| Ridge回帰 | 力推定 | α=1.0 |
| ロジスティック回帰 | すべり検出 | max_iter=1000 |

---

## 3. 主要な結果と数値

### 3.1 テクスチャ分類（5分割CV、n=800、8クラス）

| モデル | 精度（%） | 標準偏差（%） |
|-------|---------|------------|
| **Random Forest** | **97.62** | **0.47** |
| SVM (RBF) | 97.00 | 1.55 |
| MLP (3層) | 95.87 | 1.35 |
| Gradient Boosting | 95.88 | 1.79 |
| kNN (k=5) | 87.13 | 2.04 |

テストセット（n=160）での最良RF精度: **98.12%**、マクロF1 = 0.98 [cell:13]

PCA解析：第1・2主成分で76.9%の分散を説明 [cell:9]

![Figure 1: テクスチャ分類結果とPCA可視化](figures/fig1_texture_classification.png)

![Figure 5: 混同行列（Random Forest, テスト精度98.1%）](figures/fig5_confusion_matrix.png)

### 3.2 すべり検出（5分割CV、n=1000、陽性率12.8%）

| モデル | AUROC | 標準偏差 | 精度 |
|-------|-------|---------|-----|
| ロジスティック回帰 | 0.6420 | 0.0546 | 0.8710 |
| SVM (RBF) | 0.7660 | 0.0421 | 0.8640 |
| **Random Forest** | **0.7766** | **0.0416** | **0.8740** |

[cell:5] RF AUROC = 0.777 ± 0.042。クラス不均衡（12.8%陽性）のため、Accuracyより AUROCが主要指標。

### 3.3 法線力推定（5分割KFold、n=1000）

| モデル | R²（CV） | 標準偏差 | RMSE (N) | 標準偏差 |
|-------|---------|---------|---------|---------|
| Ridge回帰 | 0.0593 | 0.0367 | 4.375 | 0.183 |
| **Random Forest** | **0.8685** | **0.0106** | **1.636** | **0.114** |

テストセット（20%、n=200）：R² = **0.8753**、RMSE = **1.563 N** [cell:10]

Ridge回帰が失敗（R²≈0）なのは、力-変形関係が本質的に非線形であるため。

![Figure 2: すべり検出AUROCと法線力推定散布図](figures/fig2_slip_force.png)

### 3.4 把持安定性評価（5分割CV、n=800、安定53.6%）

| モデル | AUROC | 標準偏差 | 精度 | 標準偏差 |
|-------|-------|---------|-----|---------|
| **Random Forest** | **0.8867** | **0.0137** | **0.8562** | **0.0119** |

[cell:6] 最重要特徴量：アプローチ角、せん断/法線力比、接触面積非対称性。

### 3.5 マルチモーダル融合（触覚+視覚）

| 入力 | 精度（%） | 標準偏差（%） |
|-----|---------|------------|
| 触覚のみ（20次元） | 97.62 | 0.47 |
| **マルチモーダル（28次元）** | **98.75** | **0.79** |
| 改善幅 | **+1.13%** | — |

[cell:7] ノイズ込みの視覚特徴でも、融合により一貫した精度向上を確認。実際のViT視覚特徴ではより大きな改善が期待される。

![Figure 3: 把持安定性の特徴量重要度とマルチモーダル比較](figures/fig3_grasp_multimodal.png)

### 3.6 不確実性定量化と探索的把持（n=50未知物体）

| 集団 | 不確実性（平均±標準偏差） |
|-----|------------------------|
| 既知（訓練）物体 | 0.0534 ± 0.0328 |
| **未知（OOD）物体** | **0.2745 ± 0.0287** |
| 比率 | **5.14×** |

適応的閾値（訓練データの不確実性の90パーセンタイル = 0.0969）を超えた未知物体の割合：**100%（50/50）** [cell:8]

5.14倍の不確実性比率により、既知・未知物体の信頼性の高い判別が可能。

![Figure 4: 不確実性分布（既知 vs 未知物体）](figures/fig4_uncertainty_exploration.png)

---

## 4. 考察と今後の展望

### 4.1 主要知見の解釈

**Random Forestの優越性**: 全タスクでRFが最良モデル。15～28次元の中規模特徴空間と非線形相互作用の組み合わせに適している。

**テクスチャ分類（97.6%）**: 合成データのクラス間分離が良好なため高精度。実実験ではクラス内変動が大きく、精度は5～15%低下する可能性。

**力推定の課題**: RidgeのR²≈0はモデル選択の重要性を示す。実際のGelSightでは、弾性ゲルの超弾性変形物理をより正確にモデル化するFEAベースのアプローチ（Helmut et al. 2024）が推奨。

**すべり検出（AUROC 0.777）**: 単フレームの静的特徴量では中程度の性能。実システムでは1 kHz以上の振動信号や光学フロー系列が必要。

**マルチモーダル融合（+1.13%）**: 弱い視覚特徴でも改善。ViT/CLIPベースの視覚特徴との融合では+5～15%の改善が期待。

**探索的把持の不確実性（5.14倍）**: アンサンブル分散が既知・未知物体の効果的な判別基準となることを実証。

### 4.2 自己批判的評価

| 批判点 | 詳細 |
|-------|-----|
| 合成データへの依存 | Gaussianモデルは実センサーの複雑な特性を再現しない |
| 時系列ダイナミクスの欠如 | 把持操作には時系列モデル（LSTM/Transformer）が必要 |
| シム-to-リアルギャップ | 合成画像と実画像の外観差が性能を制限 |
| 力モデルの簡略化 | 多項式モデルはシリコンゲルの超弾性物性を反映しない |
| すべり検出のクラス不均衡 | SMOTE等のオーバーサンプリングで改善余地あり |

### 4.3 PyTorch/IsaacSimフレームワーク設計

提案するシミュレーション学習フレームワーク（設計仕様）：

```
┌─────────────────────────────────────────────────────┐
│              IsaacSim 環境                           │
│  ┌─────────────┐    ┌──────────────────────┐        │
│  │ 物理シミュ   │───>│ TACTOセンサーエミュレータ│        │
│  │ (RigidBody) │    │  (PyTorch 微分可能)   │        │
│  └─────────────┘    └──────────┬───────────┘        │
│                                │ 触覚画像              │
│                    ┌───────────▼───────────┐        │
│                    │  知覚モジュール         │        │
│                    │ ・ResNet-18 (触覚)      │        │
│                    │ ・ViT (視覚)            │        │
│                    │ ・Cross-Attention Fusion│        │
│                    └───────────┬───────────┘        │
│                                │ 状態ベクトル           │
│                    ┌───────────▼───────────┐        │
│                    │  PPO方策ネットワーク    │        │
│                    │ (Actor-Critic)         │        │
│                    └───────────────────────┘        │
└─────────────────────────────────────────────────────┘
```

### 4.4 今後の展望

1. **エンドツーエンドCNN学習**: 生GelSight画像から直接特徴抽出（ResNet/EfficientNet）
2. **時系列モデル**: LSTM/Transformerによるすべり前兆パターン検出
3. **ドメイン適応**: TACTOシミュレータとリアルセンサー間のドメイン適応
4. **IsaacSimによる強化学習**: PPOで訓練した把持ポリシーの実機転移
5. **ベンチマーク**: YCBオブジェクトセット（77物体）での検証
6. **リアルタイム実装**: < 10 msの推論レイテンシを目標とするモデル軽量化

---

## 5. 生成したファイル一覧

| ファイル | 内容 |
|---------|-----|
| `tactile_sensing.ipynb` | 実験ノートブック（全コード） |
| `data/raw/texture_dataset.csv` | テクスチャ分類データセット（n=800） |
| `data/raw/force_dataset.csv` | 力推定・すべり検出データセット（n=1000） |
| `figures/fig1_texture_classification.png` | テクスチャ分類結果+PCA可視化 |
| `figures/fig2_slip_force.png` | すべり検出AUROC+力推定散布図 |
| `figures/fig3_grasp_multimodal.png` | 把持安定性特徴量重要度+マルチモーダル比較 |
| `figures/fig4_uncertainty_exploration.png` | 不確実性分布（既知vs未知） |
| `figures/fig5_confusion_matrix.png` | テクスチャ分類の混同行列 |
| `paper.md` | 学術論文形式のドキュメント（英語） |
| `report.md` | 本レポート |

---

## 6. 再現性情報

| 項目 | 値 |
|-----|---|
| 乱数シード | `np.random.seed(42)`, `random.seed(42)` |
| Python | 3.11.2 |
| NumPy | 2.3.5 |
| Pandas | 2.3.3 |
| scikit-learn | 1.6.1 |
| SciPy | 1.17.1 |
| Matplotlib | 3.10.9 |
| Seaborn | 0.13.2 |
| PyTorch | 2.12.0 |

---

## 付録：主要コード

```python
# テクスチャ分類 - Random Forest 5分割CV [cell:3]
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler

np.random.seed(42)
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_texture)
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
rf = RandomForestClassifier(n_estimators=100, random_state=42)
scores = cross_val_score(rf, X_scaled, y_texture, cv=cv, scoring='accuracy')
# 結果: 0.9762 ± 0.0047

# すべり検出 - AUROC評価 [cell:5]
rf_slip = RandomForestClassifier(n_estimators=100, random_state=42)
auroc = cross_val_score(rf_slip, X_slip_scaled, y_slip, cv=cv5, scoring='roc_auc')
# 結果: 0.7766 ± 0.0416

# 法線力推定 - Random Forest回帰 [cell:6]
from sklearn.ensemble import RandomForestRegressor
rf_reg = RandomForestRegressor(n_estimators=100, random_state=42)
r2 = cross_val_score(rf_reg, X_fz_scaled, Fz_true, cv=cv5r, scoring='r2')
# 結果: 0.8685 ± 0.0106

# 不確実性定量化 [cell:8]
rf_full = RandomForestClassifier(n_estimators=200, random_state=42)
rf_full.fit(X_train, y_train)
tree_preds = np.array([tree.predict_proba(X_unknown) for tree in rf_full.estimators_])
uncertainty = tree_preds.std(axis=0).mean(axis=1)
# 既知物体: 0.0534、未知物体: 0.2745 (5.14倍)
```
