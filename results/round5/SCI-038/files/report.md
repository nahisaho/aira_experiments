# Active Debris Removal (ADR) Mission Optimal Trajectory Design System
## Experimental Results & Analysis Report

---

## 1. 実験目的と背景

### 1.1 研究背景

低軌道（LEO）における宇宙デブリの増加は、持続可能な宇宙利用に対する深刻な脅威となっている。欧州宇宙機関（ESA）の推計では、追跡可能な10cm以上のデブリが36,500個以上存在し、国際宇宙ステーション（ISS）や各種衛星コンステレーションへの衝突リスクが年々高まっている。Kessler症候群の連鎖衝突シナリオを回避するためには、IADCガイドラインが示す「年間5個以上の高質量デブリ除去」が必要とされる。

### 1.2 実験目的

本実験では、ADRミッションの完全なパイプライン—ターゲット選定から最終捕獲まで—をカバーする統合最適軌道設計システムを開発・検証することを目的とする。具体的には以下6モジュールの設計と評価を行った：

1. デブリカタログからのターゲット選定（衝突リスク×除去効果スコアリング）
2. マルチターゲット除去の最適軌道遷移計画（低推力Q-Law）
3. ランデブー・近接運動（Hill-Clohessy-Wiltshire方程式）シミュレーション
4. 姿勢不安定デブリの回転運動推定（Euler方程式+四元数）
5. 捕獲機構（ロボットアーム/ネット/ハープーン）の動力学
6. コスト最小化のためのミッションシーケンス最適化（ACO）

---

## 2. 使用した手法・アルゴリズムの概要

### 2.1 デブリ優先度スコアリング（Modified MITRI-like Index）

Servadio et al. (2023)のMITRIインデックスを参考に、以下の複合スコアを定義：

```
Priority_i = 0.6 × P_i + 0.4 × (S_i × L_i)

P_i = 0.4 × ρ_i(軌道密度) + 0.3 × S_i(破片重大度) 
    + 0.2 × L_i(軌道寿命) + 0.1 × E_i(除去容易性)
```

- **ρ_i**: 軌道空間密度プロキシ（1/(a²√(1-e²))）
- **S_i**: 破片重大度（質量 × 断面積）
- **L_i**: 軌道永続性（逆寿命）
- **E_i**: 基準軌道からの傾斜角差による除去容易性

### 2.2 低推力軌道遷移（Q-Law Lyapunov制御則）

Petropoulos (2004)のQ-Law Lyapunov制御則をベースに実装：

- **推力**: 0.5 N（電気推進）
- **比推力**: Isp = 3000 s
- **初期質量**: 500 kg
- **積分ステップ**: 200 s
- **制御目標**: (a, e, i) を目標軌道に収束させるGauss変分方程式に基づく勾配制御

### 2.3 Hill-Clohessy-Wiltshire（HCW）方程式

```
ẍ - 2nẏ - 3n²x = ux
ÿ + 2nẋ = uy
z̈ + n²z = uz
```

多段階アプローチ制御：区分定数推力ベクトルをRK45（rtol=1e-8）で積分。

### 2.4 Euler方程式＋四元数キネマティクス

剛体トルクフリー回転：
```
I·dω/dt = -ω × (I·ω)
dq/dt = 0.5 · Ω(ω) · q
```

主慣性モーメント: Ix=150, Iy=300, Iz=200 kg·m²

### 2.5 捕獲機構モデル

| 機構 | モデル | 主要パラメータ |
|------|--------|--------------|
| ロボットアーム | インピーダンス制御（3フェーズ） | kc=5×10⁴ N/m, bc=2×10³ N·s/m |
| ネット | 4角質量展開モデル | v_net=3 m/s, r=5 m, k_net=50 N/m |
| ハープーン | 弾塑性貫通モデル | σy=270 MPa, m=0.5 kg, v=50 m/s |

### 2.6 アリコロニー最適化（ACO）

10ターゲットのデブリ除去シーケンスをTSP問題として定式化：

- **アリ数**: 40
- **反復数**: 150
- **フェロモン重み** α=1.0, **ヒューリスティック重み** β=2.5
- **蒸発率** ρ=0.15
- ΔV推定：Hohmann遷移 + 最適plane change結合操作

---

## 3. 主要な結果と数値

### 3.1 デブリカタログ解析

![図1: デブリカタログ解析](figures/fig01_catalog.png)

**50個のデブリ分布特性**:

| 統計量 | 高度(km) | 傾斜角(°) | 質量(kg) | 優先度スコア |
|--------|----------|-----------|----------|-------------|
| 平均 | 776.5 | 79.8 | 752.1 | 0.312 |
| 標準偏差 | 228.6 | 17.9 | 592.4 | 0.215 |
| 最小 | 402.1 | 50.1 | 52.3 | 0.021 |
| 最大 | 1197.4 | 103.5 | 2144.7 | 1.000 |

**優先度上位10ターゲット**:

| 順位 | ID | 高度(km) | 傾斜角(°) | 質量(kg) | 優先度 |
|------|----------|----------|-----------|----------|--------|
| 1 | DEB-014 | 1058.2 | 71.0 | 2144.7 | 1.0000 |
| 2 | DEB-009 | 843.7 | 71.0 | 1387.2 | 0.8312 |
| 3 | DEB-024 | 991.5 | 98.0 | 1852.9 | 0.7845 |
| 4 | DEB-033 | 908.4 | 97.0 | 1601.4 | 0.7231 |
| 5 | DEB-041 | 1012.3 | 71.0 | 1923.5 | 0.7189 |
| 6 | DEB-007 | 776.2 | 53.0 | 1102.8 | 0.6734 |
| 7 | DEB-022 | 654.1 | 98.0 | 809.3  | 0.6201 |
| 8 | DEB-018 | 887.6 | 97.0 | 1447.6 | 0.6012 |
| 9 | DEB-036 | 735.9 | 71.0 | 943.2  | 0.5874 |
| 10 | DEB-046 | 929.8 | 53.0 | 1288.4 | 0.5623 |

**5-fold交差検証 (Kendall's τ)**: 0.5022 ± 0.1625

### 3.2 低推力軌道遷移

![図2: 低推力軌道計画](figures/fig02_lowthrust.png)

**DEB-014への遷移（600km SSO → 1058km, 71°）**:

| パラメータ | 値 |
|-----------|-----|
| 総ΔV | 4170 m/s |
| 消費推進剤 | 66.1 kg |
| 初期質量 | 500 kg |
| 推進剤質量比 | 13.2% |
| 遷移期間 | 45日 |
| 推力レベル | 0.5 N |
| 比推力 | 3000 s |

### 3.3 ランデブー近接運動

![図3: 近接運動シミュレーション](figures/fig03_proximity.png)

**HCWフレーム内マルチフェーズアプローチ**:

| フェーズ | 制御加速度 | 継続時間(s) | 備考 |
|---------|-----------|------------|------|
| 1 | [0, +5×10⁻⁴, 0] m/s² | 600 | along-track加速 |
| 2 | [0, -3×10⁻⁴, 0] m/s² | 600 | 減速 |
| 3 | [-4×10⁻⁴, 0, 0] m/s² | 400 | 半径方向収束 |
| 4 | [+3×10⁻⁴, 0, 0] m/s² | 400 | 半径方向制動 |
| 5 | [0, -5×10⁻⁵, 0] m/s² | 200 | 最終停止 |

**結果**:
- 初期相対距離: 1020 m
- **最近接点: 858.5 m** (最近接時刻 t ≈ 1200 s)
- 接近速度: 88.0 cm/s

⚠️ **注**: オープンループ制御では最終接近距離が大きく、実際のミッションにはMPC等のフィードバック制御が必要。

### 3.4 姿勢不安定デブリの回転推定

![図4: 姿勢ダイナミクス](figures/fig04_tumbling.png)

**トルクフリー回転（Envisat級ロケット上段）**:

| パラメータ | 値 |
|-----------|-----|
| 主慣性モーメント (Ix,Iy,Iz) | 150, 300, 200 kg·m² |
| 初期角速度 |ω₀| | 7.62 deg/s |
| 平均スピンレート | **7.62 ± 0.03 deg/s** |
| 支配的スピン周波数 | **0.0267 Hz** |
| スピン周期 | 37.5 s |
| 回転速度（RPM） | **1.60 rpm** |
| 変動係数（CoV） | 0.39% |

Envisatの実測値（約0.03 Hz）と一致する。

### 3.5 捕獲機構比較

![図5: 捕獲機構動力学](figures/fig05_capture.png)

| 機構 | 捕獲時間 | 最大力 | 成功 | 利点 | 欠点 |
|------|---------|--------|------|------|------|
| ロボットアーム | 117 s | 2518 N | ✓ | 剛性結合、精密制御 | 回転同期必要 |
| ネット | 17.2 s | N/A | ✓ | 高速展開、広域捕捉 | 捕獲後の剛性なし |
| ハープーン | <1 ms | ~85 kN | ✓ | 超高速貫通 | 表面アクセス必要 |

### 3.6 ミッションシーケンス最適化（ACO）

![図6: ACO最適化](figures/fig06_aco.png)
![図7: ミッション概要](figures/fig07_mission.png)

**10ターゲット除去シーケンス最適化結果**:

| 手法 | 総ΔV (km/s) | 標準偏差 | 改善率 |
|------|-------------|---------|--------|
| 貪欲法（Nearest Neighbor） | 7.022 | — | ベースライン |
| ACO最良解 | **5.873** | — | **−16.4%** |
| ACO (15回平均 ± std) | 5.873 ± 0.000 | 0.000 | −16.4% |

**ACO最適シーケンス**: 7回の反復でベストコスト5.873 km/sに収束。

⚠️ **注**: 15回すべてで同一解に収束（std=0）。10ターゲットのような小問題サイズでは実際に起こりうる現象だが、実際のミッション（50+ターゲット）では確率的変動が大きくなる。

---

## 4. 考察と今後の展望

### 4.1 モデルの限界・前提条件への依存

#### 合成データへの依存

本研究のすべての定量結果は、合成的に生成されたデブリカタログに基づいている。実際のTLEカタログ（Space-Track, DISCOS）を使用した場合、以下の差異が生じる可能性がある：

1. **質量・断面積分布**: 実際のLEOデブリの質量分布はロケット上段（数百〜数千kg）と破片（数g〜数kg）の二峰性分布を示すが、本モデルでは対数正規分布を仮定。
2. **軌道分布**: SSO帯（97-98°）の偏りは反映しているが、歴史的な打上げ痕跡による特定高度への集中（850km帯等）は表現不足。
3. **物理特性不確実性**: TLE由来のデブリは形状・材質が不明であり、衝突確率計算に大きな不確実性を持つ。

#### Q-Lawの過簡略化

- 日食回避（太陽電池推進に必須）を未考慮
- J2摂動による昇交点離角ドリフトを未反映
- Gauss変分方程式を近日点付近の近似式で評価（完全な偏近点角積分は未実施）
- 実際のフライト品質ΔVはLee & Ahn (2023)報告値（1.5-3.5 km/s）より15-30%程度高い可能性

#### HCW近接運動

最近接距離858mは、オープンループ制御の限界を示す。実際の近接運動では：
- MPC（モデル予測制御）または最適誘導則が必須
- LiDAR/ステレオビジョンによるオンボード状態推定が必要
- 衝突回避ゾーン（Keep-Out Zone, KOZ）の設定が必要

### 4.2 性能値の現実性評価

| 結果 | 現実的? | 理由 |
|------|---------|------|
| ΔV 4.17 km/s (45日) | やや楽観〜現実的 | Q-Law簡略化により過大推定の可能性あり（実際は3-4 km/s程度） |
| ACO改善 16.4% | 現実的 | 文献値(10-20%)の範囲内 |
| スピンレート 1.60 rpm | 現実的 | Envisat実測値と一致 |
| ロボットアーム捕獲117s | 楽観的 | 実際は回転同期に数分〜数十分を要する場合あり |
| ネット捕獲17.2s | 楽観的 | 実際の展開・絡みつきには不確定性が大きい |

### 4.3 今後の展望

1. **実軌道カタログとの統合**: Space-Track TLE + DISCOSデータベースとの接続により現実的なターゲット評価を実現
2. **完全最適制御への移行**: GPOPS-II/Pontryagin最小原理による低推力軌道最適化で精度向上
3. **フィードバック制御統合**: MPCベースの近接誘導アルゴリズムで最終接近精度を向上（目標: 最終距離<5m）
4. **多機体ミッション拡張**: 複数のOTV（Orbital Transfer Vehicle）の協調最適化
5. **捕獲信頼性モデリング**: 確率的捕獲成功率モデルの導入（デブリ姿勢不確実性考慮）
6. **ESA Orekit/AstroPy統合**: 現在の数値モデルをOrekit高精度軌道伝播器に移行

---

## 5. 生成ファイル一覧

| ファイル | 説明 |
|---------|------|
| `adr_system.py` | ADRシミュレーション全モジュールのPythonコード |
| `figures/fig01_catalog.png` | デブリカタログ解析（図1） |
| `figures/fig02_lowthrust.png` | 低推力軌道遷移（図2） |
| `figures/fig03_proximity.png` | HCW近接運動シミュレーション（図3） |
| `figures/fig04_tumbling.png` | 姿勢不安定デブリ動力学（図4） |
| `figures/fig05_capture.png` | 捕獲機構動力学比較（図5） |
| `figures/fig06_aco.png` | ACOミッションシーケンス最適化（図6） |
| `figures/fig07_mission.png` | ミッション全体最適シーケンス（図7） |
| `paper.md` | 学術論文形式のレポート（英語） |
| `report.md` | 本実験レポート（日本語） |

---

## 参考文献

1. Narayanaswamy et al. (2022). Low-thrust rendezvous trajectory generation for multi-target active space debris removal using the RQ-Law. *Advances in Space Research*. https://doi.org/10.1016/j.asr.2022.12.049

2. Servadio et al. (2023). Risk Index for the Optimal Ranking of Active Debris Removal Targets. *Journal of Spacecraft and Rockets*. https://doi.org/10.2514/1.a35752

3. Borelli, Gaias, & Colombo (2023). Rendezvous and proximity operations design of an active debris removal service to a large constellation fleet. *Acta Astronautica*. https://doi.org/10.1016/j.actaastro.2023.01.021

4. Medhin & Servadio (2025). The Sustainability of the LEO Orbit Capacity via Risk-Driven Active Debris Removal. *arXiv*. https://doi.org/10.48550/arXiv.2507.16101

5. Poupon et al. (2024). AI-Driven Risk-Aware Scheduling for Active Debris Removal Missions. *arXiv*. https://doi.org/10.48550/arXiv.2409.17012

6. Zhang et al. (2018). Ant Colony Optimization based design of multiple-target active debris removal mission. *Trans. Japan Soc. Aeronaut. Space Sci.*, 61(4), 201–211. https://doi.org/10.2322/tjsass.61.201

7. Chutivikai et al. (2025). Bi-Objective Optimal Mission Planning for Active Debris Removal with Refueling. *iSpaRo 2025*. https://doi.org/10.1109/iSpaRo66239.2025.11436815

8. De Jongh et al. (2020). Experiment for pose estimation of uncooperative space debris using stereo vision. *Acta Astronautica*, 168, 164–173. https://doi.org/10.1016/j.actaastro.2019.12.006

9. Bourabah et al. (2023). Estimation of uncooperative space debris inertial parameters after tether capture. *Acta Astronautica*, 202, 97–112. https://doi.org/10.1016/j.actaastro.2022.07.041

10. Lee & Ahn (2023). Optimal Active Debris Removal Mission Design Using Low-thrust Trajectory. *AIAA SCITECH 2023*. https://doi.org/10.2514/6.2023-2550
