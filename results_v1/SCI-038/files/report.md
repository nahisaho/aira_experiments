# Active Debris Removal (ADR) ミッション最適軌道設計システム報告書

**タイムスタンプ:** 2026-05-22T14:34:31+00:00

## 1. 実験目的と背景

本研究の目的は、低軌道（LEO）に存在する宇宙デブリを対象として、**候補デブリ生成 → 優先度評価 → 低推力軌道遷移評価 → 近傍運動解析 → タンブリング評価 → 捕獲機構評価 → ミッション系列最適化**を一体化した Python ベースの ADR ミッション設計系を構築することである。実運用の初期設計段階では、高忠実度ツールに先立って、物理整合的な簡略モデルで多数案を高速比較できることが重要である。本システムはそのためのトレードスタディ基盤として実装した。

本実装では、実デブリカタログの代わりに再現可能な乱数種 `seed=42` を用いた**50物体の模擬カタログ**を生成し、解析結果を `data/`、`results/`、`figures/`、`logs/` に保存した。

## 2. 使用した手法・アルゴリズムの概要

### 2.1 デブリカタログ生成

各物体について、軌道6要素（半長径、離心率、傾斜角、RAAN、AoP、真近点離角）と物理量（質量、代表寸法、断面積、RCS）を生成した。半長径は 6571–7371 km に制約し、LEO に相当する高度分布を与えた。軌道減衰は弾道係数

$$
B = \frac{m}{C_D A}
$$

と指数大気密度近似から評価し、日あたり減衰率と寿命を算出した。

### 2.2 ターゲット選定

衝突確率指標は、高度帯に対するガウス型の軌道密度関数、断面積、相対速度指標から

$$
P_{\mathrm{coll}} \propto \rho_{\mathrm{orb}}(h)\, A_{\mathrm{RCS}}\, v_{\mathrm{rel}}
$$

で定義した。除去効果指標は

$$
R_{\mathrm{eff}} = m A \left(\frac{1}{\tau_{\mathrm{decay}}}\right)
$$

とし、両者を 0–1 正規化した後、

$$
S = w_1 \hat{P}_{\mathrm{coll}} + w_2 \hat{R}_{\mathrm{eff}},\qquad (w_1=0.55,\; w_2=0.45)
$$

で統合した。

### 2.3 低推力軌道遷移

ターゲット間遷移には Edelbaum 近似を採用し、円軌道速度 $V_1,V_2$ と軌道傾斜角差 $\Delta i$ から

$$
\Delta V_{\mathrm{Edelbaum}} = \sqrt{V_1^2 + V_2^2 - 2V_1V_2\cos\left(\frac{\pi}{2}\Delta i\right)}
$$

で評価した。連続低推力加速度は $a=10^{-4}\,\mathrm{m/s^2}$ とし、遷移時間は

$$
T \approx \frac{\Delta V}{a}
$$

で見積もった。

### 2.4 ランデブー

近傍運動は Clohessy–Wiltshire (Hill) 方程式

$$
\ddot{x} - 2n\dot{y} - 3n^2x = 0,\quad
\ddot{y} + 2n\dot{x} = 0,\quad
\ddot{z} + n^2 z = 0
$$

の解析解を用いた。初期相対距離約 5 km の V-bar 接近を仮定し、状態遷移行列から 2 インパルス操舵量を求めた。

### 2.5 デブリ回転推定

円筒形ロケット胴体を仮定し、主慣性モーメントを

$$
I_{\mathrm{axial}} = \frac{1}{2}mr^2,\qquad
I_{\mathrm{trans}} = \frac{1}{12}m(3r^2 + L^2)
$$

で定義した。無外力オイラー方程式を `scipy.integrate.solve_ivp` で積分し、生成した擬似ライトカーブに対して FFT により回転周期を推定した。

### 2.6 捕獲機構

- **ロボットアーム**: 2関節平面アームの順運動学と接触力近似
- **ネット**: 展開半径の時間発展と包絡成功率近似
- **ハープーン**: 回転に伴う角度誤差を含む貫入成功率近似

を実装し、回転角速度に対する成功確率曲線を比較した。

### 2.7 ミッション最適化

上位 10 ターゲットを候補集合とし、燃料制約付き TSP 変種として扱った。初期解は nearest-neighbor、局所改善に 2-opt、全体探索に遺伝的アルゴリズムを用いた。評価関数は

$$
J = \Delta V + \lambda T - \alpha \sum S_k - \beta N
$$

とし、$\Delta V < 2000\,\mathrm{m/s}$ を満たす訪問系列のみを採用した。ここで $N$ は燃料制約内で実際に訪問できたターゲット数である。

## 3. 主要な結果と数値

### 3.1 カタログ統計

- 生成デブリ数: **50**
- 平均高度: **628.186 km**
- 高度標準偏差: **224.527 km**
- 質量中央値: **154.378 kg**
- 寸法中央値: **0.671 m**
- 平均軌道寿命: **9182.11 日**

### 3.2 ターゲット選定結果

上位候補は以下であった。

1. **DEBRIS-049**: 総合スコア **0.5669**
2. **DEBRIS-012**: 総合スコア **0.5197**
3. **DEBRIS-037**: 総合スコア **0.4502**

DEBRIS-049 は衝突確率スコアが **1.0000** と最大であり、DEBRIS-037 は除去効果スコアが **1.0000** と最大であった。すなわち、本スコアリングは「衝突危険度」と「除去便益」を異なる候補に対して切り分けて評価できている。

### 3.3 軌道遷移結果

- 上位 10 目標間の平均 pairwise ΔV: **3865.37 m/s**
- 最大 pairwise ΔV: **8854.43 m/s**
- 初期機から最も近い上位候補への遷移 ΔV 最小値: **129.32 m/s**

この結果は、上位スコア物体が必ずしも軌道力学的に近接していないことを示す。したがって、単純な「スコア順巡回」は ADR 初期設計では不適切であり、燃料制約付き系列最適化が必須である。

### 3.4 ランデブー結果

3シナリオのうち、最小 ΔV は **Scenario-A: 0.5577 m/s** であったが、残留距離は **568.94 m** であり厳密ドッキングには未到達であった。**Scenario-B** は総 ΔV **18.0989 m/s** を要した一方、最終距離は数値的に **0 m** であり、完全収束を達成した。**Scenario-C** は **0.5916 m/s**、残留距離 **1395.36 m** であった。

以上より、ドリフト条件下では「最小 ΔV」と「最終位置精度」は一致せず、実運用では終端拘束を明示した誘導則が必要である。

### 3.5 回転推定結果

- 真の回転周期: **30.0000 s**
- FFT 推定周期: **30.0125 s**
- 推定誤差: **0.0125 s**（約 **0.042%**）
- 平均スピンレート: **12.833 deg/s**
- `5 deg/s` 未満の自然捕獲ウィンドウ: **0 回**

したがって、本模擬物体は受動的に静穏化する時間帯を持たず、捕獲前に減速・拘束戦略が必要である。

### 3.6 捕獲機構比較

成功確率は回転角速度の増加とともに低下したが、低〜中程度スピン領域ではロボットアームが最も頑健であった。

- 0 deg/s: Arm **0.920**, Net **0.858**, Harpoon **0.583**
- 5 deg/s: Arm **0.917**, Net **0.656**, Harpoon **0.533**
- 10 deg/s: Arm **0.909**, Net **0.507**, Harpoon **0.486**

50% 成功確率に相当する回転閾値は、Net **10.25 deg/s**、Harpoon **8.5 deg/s** であった。ロボットアームは評価範囲 30 deg/s でも 50% を下回らず、最も安定である。

### 3.7 ミッション最適化結果

燃料制約 **2000 m/s** の下で、最適化器は上位 10 候補から **3 目標**を訪問する実行可能系列を選択した。

- 選択系列: **DEBRIS-001 → DEBRIS-049 → DEBRIS-012**
- 総 ΔV: **1533.78 m/s**
- 総遷移時間: **177.52 日**
- 累積スコア: **1.3196**
- 実行可能性: **fuel budget 内で feasible**

初期の全 10 物体巡回案は 9952.64 m/s と非実行であったが、制約付き最適化へ改良した結果、燃料内で意味のある 3 目標系列が抽出できた。

## 4. 考察と今後の展望

本システムは、ADR の概念設計に必要な主要構成要素を一貫して実装し、**ターゲット価値**と**軌道力学コスト**の緊張関係を定量化できた点に意義がある。特に、スコア上位目標をそのまま全件巡回する案が燃料制約を大きく逸脱すること、ならびに回転安定化なしには高精度捕獲が難しいことが明確になった。

一方で、以下の限界がある。

1. 軌道遷移は Edelbaum 近似であり、J2 摂動・昇交点整列・推力方向制約を含まない。
2. ランデブーは線形 CW モデルであり、非線形相対運動やセンサ雑音を考慮していない。
3. 捕獲成功率は工学的近似であり、高忠実度接触力学や柔軟体効果は未実装である。
4. 模擬カタログは再現性重視の合成データであり、実運用には実測 TLE/OD データとの同化が必要である。

今後は、(i) J2 を含む長期位相最適化、(ii) MPC/凸最適化による終端拘束付きランデブー、(iii) 姿勢制御と捕獲機構の協調設計、(iv) 実カタログ接続による運用レベル評価へ拡張するのが望ましい。

## 5. 生成したファイル一覧

### 主要コード
- `adr_mission/__init__.py`
- `adr_mission/debris_catalog.py`
- `adr_mission/target_selection.py`
- `adr_mission/orbit_transition.py`
- `adr_mission/rendezvous.py`
- `adr_mission/debris_rotation.py`
- `adr_mission/capture_mechanism.py`
- `adr_mission/mission_optimizer.py`
- `adr_mission/main.py`

### データ
- `data/debris_catalog.csv`
- `data/preprocessing-log.md`

### 数値結果
- `results/target_scores.csv`
- `results/delta_v_matrix.csv`
- `results/transfer_time_matrix_days.csv`
- `results/rendezvous_trajectories.csv`
- `results/rotation_analysis.csv`
- `results/capture_analysis.csv`
- `results/optimal_mission_sequence.json`
- `results/mission_summary.json`
- `results/statistical-summary.md`
- `results/target_selection_summary.json`
- `results/orbit_transition_summary.json`
- `results/rendezvous_summary.json`
- `results/rotation_summary.json`
- `results/capture_summary.json`

### 図
- `figures/target_selection.png`
- `figures/delta_v_heatmap.png`
- `figures/rendezvous_trajectory.png`
- `figures/debris_rotation.png`
- `figures/capture_mechanisms.png`
- `figures/mission_optimization.png`

### ログ
- `logs/process-log.jsonl`
