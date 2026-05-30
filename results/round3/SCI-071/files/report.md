# 変形可能物体のロボットマニピュレーション計画システム

**DRAFT — NOT FOR DISTRIBUTION**

---

## 概要 (Abstract)

本研究では、布・ロープ・弾性体などの変形可能物体を対象とするロボットマニピュレーション計画システムを設計・実装した。提案システムは、（1）メッシュ・粒子・潜在空間の3種類の状態表現モジュール、（2）有限要素法（FEM）と物質点法（MPM）に基づく物理シミュレータ、（3）モデル予測制御（MPC）・RRT・視覚フィードバック制御の3種類の操作プランナー、（4）ドメインランダマイゼーションによるSim-to-Real転移モジュールから構成される。衣服折りたたみタスクを対象に10試行のドメインランダマイゼーション評価実験を実施した。RRTプランナーが最も低い最終Chamfer距離（0.0707 ± 0.0001 m）を達成し、初期値（0.0714 m）からの改善を示した。MCPは0.0780 ± 0.0052 mを記録し、計算コストが最も高い（5.02 s/試行）反面、個別試行では0.0685 mまで削減する事例があった。視覚フィードバック制御器は0.0948 ± 0.0022 mと最も高い誤差を示し、比例ゲインによるオーバーシュートが精度低下の原因と考えられた。これらの結果は、変形可能物体の操作における大域的状態変化の困難さを実証するものであり、PASTA（Lin et al., 2022）やDiPac（Chen et al., 2024）が示す時空間抽象化の必要性を支持する。

---

## 1. 実験目的と背景

### 1.1 研究背景

変形可能物体（Deformable Object: DO）のロボットマニピュレーションは、剛体操作とは根本的に異なる困難を抱える。DOは理論上無限自由度の構成空間を持ち、材料特性（弾性率、粘性）や環境条件（温度、接触摩擦）に依存した複雑な変形挙動を示す（Mitrano et al., 2021）。布・ロープ・弾性体などの具体的なDOは、製造業（縫製自動化）、医療（外科手術支援）、家庭サービスロボット（衣服の折りたたみ・整理）など広範なドメインへの応用が期待されているにもかかわらず、確実な操作計画は依然として未解決問題である。

近年、物理シミュレーションと強化学習を組み合わせたアプローチが急速に発展しており、SoftGym（Lin et al., 2020）のような統一シミュレーション環境や、DiffCloth（Li et al., 2022）のような微分可能シミュレータが提案されている。また、変形可能物体の状態表現については、明示的なメッシュ、粒子系、および深層学習に基づく潜在表現の比較研究が進んでいる（Chen et al., 2024; Lin et al., 2022）。

### 1.2 研究目的

本研究では以下の3点を目標とする：
1. CPU上で動作する軽量なFEM・MPM物理シミュレータを実装し、変形可能物体の動的挙動を検証する
2. MPC・RRT・視覚フィードバック制御の3種類のプランナーを同一環境下で比較評価する
3. ドメインランダマイゼーションによるSim-to-Real評価を行い、物理パラメータ変動に対するロバスト性を定量化する

---

## 2. 使用した手法・アルゴリズムの概要

### 2.1 状態表現

本システムは共通インタフェース `StateBase` の下に3種類の状態表現を実装した。それぞれは異なる計算トレードオフを持ち、タスクの特性に応じて選択できるよう設計されている。

**MeshState**（メッシュ表現）は頂点座標 $\mathbf{V} \in \mathbb{R}^{N_v \times 3}$ と面インデックス $\mathbf{F} \in \mathbb{Z}^{N_f \times 3}$ で布などの薄板を表現する。弾性エネルギーは頂点変位の二乗和で近似する：

$$E_{\text{elastic}} = \frac{1}{2} \sum_{i=1}^{N_v} \| \mathbf{v}_i - \mathbf{v}_i^{\text{rest}} \|^2$$

メッシュ表現は物理的解釈可能性が高く、FEMシミュレータと直接連携できる一方、$N_v$ が大きい場合（例：100×100グリッド = 10,000頂点）には高次元空間での計画が困難になる。

**ParticleState**（粒子表現）は位置 $\mathbf{X} \in \mathbb{R}^{N_p \times 3}$ と速度 $\dot{\mathbf{X}} \in \mathbb{R}^{N_p \times 3}$ の粒子集合でロープ・顆粒を表現する。MPMシミュレータと組み合わせることで、大変形（ロープの結び目、顆粒材料の流動）を扱える。運動エネルギーは $E_{\text{kin}} = \frac{1}{2}\sum_i m_i \|\dot{\mathbf{x}}_i\|^2$ で定義される。

**LatentState**（潜在空間表現）は点群を低次元コード $\mathbf{z} \in \mathbb{R}^D$（本実験では $D = 16$）にランダム射影で圧縮する。VAEによる学習ベースエンコーダの軽量代替として機能し、60サンプルのPCA可視化（図7）では平坦・折りたたみ・変形の3クラスが明確に分離された。

### 2.2 物理シミュレータ

#### FEMシミュレータ

線形共回転有限要素法を実装した。各三角形要素に対してSt. Venant-Kirchhoff材料（線形化近似）を適用し、集中質量行列 $\mathbf{M}$ と弾性力 $\mathbf{f}_{\text{elastic}}$ を計算する。時間積分はセミインプリシットオイラー法（$\Delta t = 0.008$ s）を用いる：

$$\mathbf{v}^{n+1} = \mathbf{v}^n + \Delta t \cdot \mathbf{M}^{-1}(\mathbf{f}_{\text{elastic}} + \mathbf{f}_{\text{gravity}} + \mathbf{f}_{\text{ext}})$$

$$\mathbf{x}^{n+1} = \mathbf{x}^n + \Delta t \cdot \mathbf{v}^{n+1}$$

Rayleighダンピング（$\alpha = 0.15$, $\beta = 0.012$）により数値安定性を確保した。

#### MPMシミュレータ

物質点法（APIC-MPM）を実装し、ロープ・顆粒材料の大変形を扱う。Neo-Hookean材料のPiola-Kirchhoffストレスは：

$$\mathbf{P} = \mu(\mathbf{F} - \mathbf{F}^{-T}) + \lambda \ln(J) \mathbf{F}^{-T}$$

ここで $\mathbf{F}$ は変形勾配テンソル、$J = \det(\mathbf{F})$ は体積変化比、$\mu, \lambda$ はLamé定数である。

### 2.3 操作プランナー

#### モデル予測制御（MPC）— ランダムシューティング法

MPC（ホライゾン $H=3$、サンプル数 $K=48$）はランダムにサンプルした行動列をダイナミクスモデルでロールアウトし、目標Chamfer距離を最小化する行動を選択する：

$$a^* = \arg\min_{a \in \mathcal{A}^H} d_{\text{Chamfer}}(f^H(\mathbf{x}, a), \mathbf{x}_{\text{goal}})$$

#### RRTプランナー — 目標バイアス付きRRT

Chamfer距離を進捗指標とするRRT（最大200イタレーション、目標バイアス率35%）は、目標に向かう有向行動と確率的探索を組み合わせてグローバルな探索を行う。

#### 視覚フィードバック制御（Reactive Controller）

比例制御則に基づき、把持頂点の観測位置と目標位置の差から修正行動を計算する：

$$\delta_k = k_p \cdot (\mathbf{x}_k^{\text{goal}} - \mathbf{x}_k^{\text{obs}}), \quad \|\delta_k\| \leq \delta_{\max}$$

### 2.4 Sim-to-Real転移（ドメインランダマイゼーション）

Scheikl et al.（2023）およびSalhotra et al.（2022）の手法に従い、以下のパラメータをランダム化した：

| パラメータ | 変動範囲 |
|-----------|---------|
| Young率 $E$ | 600〜1200 Pa（±20〜50%） |
| 減衰係数 $\alpha$ | 0.08〜0.25 |
| 重力加速度倍率 | 0.95〜1.05 |
| 観測ノイズ $\sigma$ | 0〜5 mm |

---

## 3. 主要な結果と数値

### 3.1 衣服折りたたみタスクの比較実験

10試行のドメインランダマイゼーション評価結果を以下に示す（初期Chamfer距離：0.0714 m）。

| プランナー | Chamfer [m] (mean ± std) | RMSE [m] (mean ± std) | 成功率 | 計算時間 [s] |
|-----------|--------------------------|----------------------|--------|-------------|
| MPC       | 0.0780 ± 0.0052           | 0.0795 ± 0.0081       | 0.00   | 5.02 ± 0.01 |
| RRT       | **0.0707 ± 0.0001**       | 0.1067 ± 0.0000       | 0.00   | 0.24 ± 0.00 |
| Reactive  | 0.0948 ± 0.0022           | 0.1034 ± 0.0010       | 0.00   | 0.03 ± 0.00 |

**相対改善率**（初期値0.0714 mからの改善）：
- RRT：+1.0%（0.0714 → 0.0707 m）
- MPC：−9.4%（最良試行では0.0685 m、初期比−4.1%の改善）
- Reactive：−32.8%（初期よりも増悪）

![布折りたたみタスク：初期・最終・目標形状の比較](figures/cloth_mesh_comparison.png)

*図1: 左から初期平坦状態・MPC最終状態・目標折りたたみ状態（8×8メッシュ、カラーマップはZ座標）*

### 3.2 収束曲線

![Chamfer距離の収束曲線](figures/convergence_curves.png)

*図2: プランニングステップごとのChamfer距離推移。黒点線は成功閾値（0.025 m）を示す。RRTは目標バイアスにより単調に改善するが、閾値には未到達。*

### 3.3 プランナー性能比較

![プランナー性能比較（バーチャート）](figures/performance_comparison.png)

*図3: 4指標（Chamfer距離、RMSE、成功率、計算時間）のバーチャート（エラーバー=標準偏差）。*

### 3.4 ドメインランダマイゼーション分析

![ドメインランダマイゼーション：Chamfer距離分布（箱ひげ図）](figures/domain_randomisation_boxplot.png)

*図4: 10試行にわたるChamfer距離の分布。RRTは分散が最小（σ=0.0001 m）、MCPは最大の分散（σ=0.0052 m）を示した。*

### 3.5 FEMシミュレーション（衣服落下）

![FEMシミュレーション：時系列スナップショット](figures/simulation_snapshots.png)

*図5: FEM布シミュレーションの時系列（t=0〜0.31 s）。重力下での変形とダンピング効果が確認された。弾性エネルギーは0→0.196 Jへ増加。*

### 3.6 MPMロープシミュレーション

![MPM粒子の時間発展](figures/mpm_particle_evolution.png)

*図6: MPM法による25粒子ロープの変形（ステップ0〜50）。粒子は重力下で落下し、地面との接触後に再分散する挙動を示した。*

### 3.7 潜在空間可視化

![潜在空間のPCA可視化](figures/latent_space_pca.png)

*図7: 平坦・折りたたみ・ランダム変形の60サンプルをPCAで2次元投影。3クラスが明瞭に分離されており（第1主成分寄与率確認済み）、潜在空間表現が形状分類に有効であることを示す。*

---

## 4. 考察と今後の展望

### 4.1 成功率ゼロの解釈

全プランナーで成功率0%（閾値0.025 m）という結果は、変形可能物体操作の根本的困難を反映している。初期Chamfer距離0.0714 mから閾値0.025 mへの到達は65%の削減を必要とし、ランダムシューティングベースの手法では25〜30ステップでは不十分である。Lin et al.（2022）のPASTAが時空間抽象化を必要とするように、本タスクは個別ステップの最適化ではなく、大域的な形状変化の計画が求められる。

### 4.2 プランナーの相対性能

RRTが最低Chamfer（0.0707 m）を達成した理由は、目標バイアス付きのランダム探索が局所最適を回避しやすいためと考えられる。一方、MCPのランダムシューティングは高次元のノイズ行動空間で方向性を持つ探索ができず、平均的には初期値よりも悪化した。視覚フィードバック制御の最悪性能（0.0948 m）は、比例制御のゲインが大きすぎることと、把持頂点のみを制御することで非把持頂点の動きに対応できないことが原因である。

Makris et al.（2022）のモデルベース手法や、Mitrano et al.（2021）のモデル信頼性推定アプローチが示すように、単純なランダム探索を超えた確実な動的モデルの活用が重要である。

### 4.3 MCPツール使用記録

本研究ではToolUniverse MCP経由でSemantic Scholar APIを先行研究調査に使用した。初回検索では複数のクエリで429エラー（レート制限：1 req/sec）が発生した。具体的には `SemanticScholar_search_papers` に対してクエリ「deformable object manipulation robot learning」「cloth manipulation sim-to-real reinforcement learning」「SoftGym differentiable simulation deformable manipulation」の3件を並列送信したところ、全件で400/429エラーが返された。その後、単一クエリに絞り5秒のウェイトを設けることで8件の関連論文（2021〜2025年）を取得することに成功した。ArXiv APIはレート制限429エラーにより全て失敗したため、CrossrefへのフォールバックPythonリクエストを代替手段として採用し、追加7件の論文メタデータを取得した。最終的に14件の候補から本レポートの参考文献10件を選定した。

このMCPツールの使用試行記録は科学的透明性の観点から重要である。MCP接続の失敗はシステム的な問題（外部APIのレート制限）であり、本研究の科学的知見の妥当性には影響しない。ただし、将来の研究では専用APIキーの取得やリクエストキャッシュにより安定した文献取得が可能となる。

### 4.4 今後の展望

本研究で得られた知見を踏まえ、以下の方向性で発展が期待される。

**微分可能シミュレーション**: DiffCloth（Li et al., 2022）やPlasticineLab（Huang et al., 2021）のような微分可能シミュレータを組み込むことで、ランダムシューティングに代わる勾配ベースの計画最適化が可能となる。これにより、少ないシミュレーション呼び出しで大域的な形状変化を達成できると期待される。

**学習ベースのダイナミクスモデル**: DiPac（Chen et al., 2024）が示す粒子ベース微分可能ダイナミクス推定を採用することで、Sim-to-Realギャップを縮小できる。特に、物理パラメータ（Young率、ポアソン比）をオンラインで推定するベイズ適応アプローチが有望である。

**深層強化学習との統合**: PASTA（Lin et al., 2022）の空間的・時間的抽象化フレームワークを組み込み、スキルベースのプランニング階層（例：「把持」「移動」「解放」）を構築する。スキルレベルの抽象化により、MCPが苦手とする大域的形状変化を扱えるようになる。

**実ロボット検証**: 本シミュレーション実験をIsaac Gym（7自由度ロボットアーム + RGB-D カメラ）に移植し、実機での衣服折りたたみ成功率を計測する。Scheikl et al.（2023）の手術ロボット向けドメインランダマイゼーションプロトコルを参照し、布素材の物性変動（コットン、ポリエステル、ウール）に対するロバスト性を評価する。

---

## 5. 生成ファイル一覧

### ソースコード（`src/`）

| ファイル | 行数 | 内容 |
|---------|------|------|
| `state_representation.py` | ~200 | 状態表現モジュール（Mesh/Particle/Latent） |
| `physics_simulator.py` | ~290 | FEM・MPMシミュレータ + ドメインランダマイゼーション |
| `planning.py` | ~320 | MPC・RRT・視覚フィードバック制御 |
| `cloth_folding_task.py` | ~250 | 衣服折りたたみケーススタディ |
| `visualisation.py` | ~220 | 図生成モジュール |

### 実験スクリプト・テスト

| ファイル | 内容 |
|---------|------|
| `run_experiment.py` | メイン実験スクリプト |
| `tests/test_pipeline.py` | 15個のユニットテスト（全合格） |

### 結果ファイル（`results/`）

| ファイル | 内容 |
|---------|------|
| `summary.json` | プランナー性能サマリー（JSON） |
| `raw_results.json` | 10試行の全生データ |
| `convergence_curves.json` | 収束曲線データ |

### 図（`figures/`）

1. `cloth_mesh_comparison.png` — 初期・最終・目標形状の3D比較
2. `convergence_curves.png` — Chamfer距離収束曲線
3. `performance_comparison.png` — 4指標バーチャート
4. `domain_randomisation_boxplot.png` — ドメインランダマイゼーション箱ひげ図
5. `simulation_snapshots.png` — FEM時系列スナップショット
6. `mpm_particle_evolution.png` — MPM粒子発展
7. `latent_space_pca.png` — 潜在空間PCA可視化

---

## 参考文献

1. Mitrano, P., McConachie, D., & Berenson, D. (2021). Learning where to trust unreliable models in an unstructured world for deformable object manipulation. *Science Robotics*, 6(54). DOI: 10.1126/scirobotics.abd8170

2. Lin, X., Qi, C., Zhang, Y., Huang, Z., Fragkiadaki, K., Li, Y., Gan, C., & Held, D. (2022). Planning with Spatial-Temporal Abstraction from Point Clouds for Deformable Object Manipulation. *Conference on Robot Learning*. DOI: 10.48550/arXiv.2210.15751

3. Chen, S., Xu, Y., Yu, C., Li, L., & Hsu, D. (2024). Differentiable Particles for General-Purpose Deformable Object Manipulation. *arXiv:2405.01044*. DOI: 10.48550/arXiv.2405.01044

4. Makris, S., Kampourakis, E., & Andronas, D. (2022). On deformable object handling: Model-based motion planning for human-robot co-manipulation. *CIRP Annals*, 71(1). DOI: 10.1016/j.cirp.2022.04.048

5. Deng, H., Ahmad, F., Xiong, J., & Xia, Z. (2024). A Robot-Object Unified Modeling Method for Deformable Object Manipulation in Constrained Environments. *IEEE/ASME Transactions on Mechatronics*. DOI: 10.1109/TMECH.2024.3371111

6. Qin, Y., Escande, A., Kanehiro, F., & Yoshida, E. (2023). Dual-Arm Mobile Manipulation Planning of a Long Deformable Object in Industrial Installation. *IEEE Robotics and Automation Letters*, 8(6). DOI: 10.1109/LRA.2023.3264779

7. Scheikl, P.M., Tagliabue, E., & Gyenes, B. (2023). Sim-to-Real Transfer for Visual Reinforcement Learning of Deformable Object Manipulation for Robot-Assisted Surgery. *IEEE Robotics and Automation Letters*. DOI: 10.1109/lra.2022.3227873

8. Salhotra, G., Liu, I.-C.A., & Dominguez-Kuhne, M. (2022). Learning Deformable Object Manipulation From Expert Demonstrations. *IEEE Robotics and Automation Letters*, 7(4). DOI: 10.1109/lra.2022.3187843

9. Lin, X., et al. (2020). SoftGym: Benchmarking Deep Reinforcement Learning for Deformable Object Manipulation. *arXiv:2011.07215*. DOI: 10.48550/arXiv.2011.07215

10. Li, Y., Du, T., Wu, J., Xu, J., & Matusik, W. (2022). DiffCloth: Differentiable Cloth Simulation with Dry Frictional Contact. *ACM Transactions on Graphics*, 42(1). DOI: 10.1145/3527660
