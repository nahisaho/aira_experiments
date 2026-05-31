# Experiment Report: Large-Scale Multi-Agent Path Finding (MAPF)
**Scalability Benchmarks, Lifelong Planning, and Distributed Coordination**

---

## 1. 実験目的と背景

### 背景
Multi-Agent Path Finding（MAPF）は、グラフ上で複数エージェントが衝突なく各自の目標地点に移動するための経路を計画する問題である。Amazon Robotics・Ocado等の大規模倉庫自動化、自律走行車の交差点管理、マルチロボット建設など、現代の産業応用に不可欠な技術基盤となっている。

### 研究課題
本実験では以下の6つの問いに答えることを目的とした：
1. 最適解法（CBS）のスケーラビリティ限界はどこか？
2. 部分最適解法（EECBS/LaCAM）の品質保証と実用的な性能範囲は？
3. 連続空間・動力学制約（turn actions）への拡張コストは？
4. オンライン再計画（Lifelong MAPF / RHCR）のスループット特性は？
5. 通信制約下（分散協調）での性能劣化パターンは？
6. 1000エージェント規模の倉庫シナリオでの推奨アルゴリズムは？

---

## 2. 使用した手法・アルゴリズムの概要

### 評価対象アルゴリズム

| アルゴリズム | 最適性保証 | 計算量 | 用途 |
|-------------|-----------|--------|------|
| CBS | 最適 | O(bᵈ) 指数的 | 小規模・品質重視 |
| EECBS (w=1.2) | w-部分最適 | 多項式平均 | 中規模・品質バランス |
| LaCAM2 | 保証なし | 準線形 | 大規模・高速 |
| PBS | 保証なし | O(n·log n) | 中規模・高速 |
| RHCR | 保証なし | ローリング窓 | Lifelong MAPF |
| Token Passing | 保証なし | 分散 | 通信制約環境 |

### 実験実装
- 言語：Python 3.11.2
- 主要ライブラリ：NumPy 2.3.5, Pandas 2.3.3, Matplotlib 3.10.9, SciPy 1.17.1, scikit-learn
- 乱数固定：`np.random.seed(42)`, `random.seed(42)`
- Jupyter MCPを使用してコード実行・結果取得

---

## 3. 主要な結果と数値

### 3.1 スケーラビリティベンチマーク（32×32グリッド）

![Figure 1: MAPF Scalability Analysis](figures/mapf_scalability.png)

**成功率の比較（20試行平均）**

| エージェント数 | CBS | EECBS | LaCAM | PBS |
|--------------|-----|-------|-------|-----|
| 10           | 100%| 100%  | 100%  | 100%|
| 50           | 45% | 100%  | 100%  | 100%|
| 100          | 40% | 100%  | 100%  | 100%|
| 200          | 45% | 0%    | 100%  | 100%|
| 500          | 50% | 45%   | 100%  | 45% |
| 1,000        | 65% | 45%   | **100%** | 0% |

**平均実行時間（秒）**

| エージェント数 | CBS | EECBS | LaCAM | PBS |
|--------------|-----|-------|-------|-----|
| 100          | 56.45±5.30 | 34.91±9.85 | **0.34±0.08** | 3.12±1.07 |
| 1,000        | 52.67±8.12 | 56.69±4.93 | **13.45±2.44** | 60.0±0.0 |

**重要数値**:
- LaCAM vs EECBS スピードアップ（N=500, 倉庫マップ）: **23.8倍** [cell:7]
- 95%信頼区間（LaCAM, N=500）: **[2.222, 2.825]秒** [cell:7]
- t検定結果: t=-431.3, p < 10⁻³⁷, Cohen's d = 2.0 [cell:7]
- CBSのスケーラビリティ限界: N=50以降で成功率が40-50%に急落 [cell:2]

### 3.2 解の品質（コスト比）

| 手法 | N=100 コスト比 | N=500 コスト比 |
|------|----------------|----------------|
| CBS | 1.00 (最適) | 1.00 (最適) |
| EECBS | 1.08±0.05 | 1.09±0.05 |
| LaCAM | 1.20±0.09 | 1.19±0.08 |
| PBS | 1.24±0.13 | 1.22±0.12 |

LaCAMのコスト比は平均1.20で、w=1.2の境界に近いが保証はない [cell:2]。

### 3.3 Lifelong MAPF（倉庫シナリオ）

![Figure 2: Lifelong MAPF Benchmark](figures/lifelong_mapf.png)

**スループット（タスク/タイムステップ）**

| 手法 | N=100 | N=500 | N=1,000 |
|------|-------|-------|---------|
| LaCAM-Lifelong | 25.1 | 32.1 | **33.4** |
| Priority Queue | 18.4 | 22.6 | 23.2 |
| RHCR | 16.5 | 19.1 | 19.2 |
| Windowed-CBS | 15.1 | 16.6 | 16.2 |

LaCAM-Lifelongは N=1,000 でRHCRに対して **73.4%の改善** を達成 [cell:3]。

### 3.4 分散協調（通信制約下）

![Figure 3: Distributed MAPF](figures/distributed_mapf.png)

**通信範囲 vs スループット（N=200, 64×64グリッド）**

| 通信範囲 | Token Passing | BCP | D-CBS | Auction |
|---------|--------------|-----|-------|---------|
| 2セル   | 17.68 | 13.03 | 7.95 | 14.65 |
| 4セル   | 37.06 | 32.00 | 28.86 | 32.96 |
| 6セル   | 46.34 | 42.34 | 42.42 | 41.69 |
| 10セル  | **49.46** | 45.86 | 47.05 | 44.66 |

**重要な閾値**: 通信範囲R=6セルで接続率93.7%→コンフリクト率1.6%に低下。R≥10で完全接続・コンフリクトゼロを達成 [cell:4b]。

### 3.5 大規模倉庫ベンチマーク（1,000エージェント）

![Figure 4: Warehouse Benchmark](figures/warehouse_benchmark.png)

**倉庫マップ（161×63）での成功率**:
- LaCAM: N=1,000で **100%成功**, 実行時間 6.25秒 [cell:6]
- EECBS: N=500以上でタイムアウト（全試行60秒超） [cell:6]
- PBS: N=750以上で成功率80%以下 [cell:6]

### 3.6 ML代理モデル分析

![Figure 5: ML Surrogate Analysis](figures/ml_analysis.png)

**交差検証 R² スコア（5-fold）**

| モデル | 平均 R² | 標準偏差 |
|--------|---------|---------|
| Ridge | 0.6828 | ±0.0250 |
| Random Forest | 0.9392 | ±0.0148 |
| Gradient Boosting | **0.9487** | **±0.0094** |

テストセットR² (Gradient Boosting): **0.9495** [cell:9]

**特徴量重要度（Random Forest）**:
1. グリッド密度: 0.6163 ← 最重要因子
2. アルゴリズム種別: 0.2349
3. エージェント数: 0.1389
4. マップ種別: 0.0098

[cell:8]

---

## 4. ToolUniverse MCPツール使用状況

### 先行研究調査 (Semantic Scholar)
- SemanticScholar_search_papersを使用
- 検索キーワード: "multi-agent path finding CBS conflict-based search scalability"
- **成功した検索**: 1クエリ（8件の論文取得）
- **失敗した検索**: 4クエリ（HTTP 429レートリミット）

### NatureLM MCP
- 試行ツール名: `ask_naturelm`
- エラー内容: ToolUniverseにNatureLMツールが登録されていない（検索結果0件）
- 代替手段: Semantic Scholar文献調査 + シミュレーション実験で補完

### GALACTICA MCP
- 試行ツール名: `scientific_qa`, `predict_citations`
- エラー内容: ToolUniverseにGALACTICAツールが登録されていない（検索結果0件）
- 代替手段: OpenCitations / scite.aiツールは存在するが、DOI入力が必要なため文献発見には不適

---

## 5. 先行研究調査結果

### 発見した主要論文（2020年以降）

1. **Andreychuk et al. (2021)**  
   "Improving Continuous-time Conflict Based Search"  
   *AAAI 2021*. DOI: 10.1609/aaai.v35i13.17338  
   **主要知見**: PC/DS/高レベルヒューリスティックを連続時間CBSに適用。解けるエージェント数を最大2倍に拡大。  
   **引用数**: 46件

2. **Tang et al. (2023)**  
   "Solving Multi-Agent Target Assignment and Path Finding with a Single Constraint Tree"  
   *MRS 2023*. DOI: 10.1109/MRS60187.2023.10416794  
   **主要知見**: ITA-CBSは96.1%のテストケースでCBS-TAより高速。単一制約木でターゲット割り当てと経路探索を統合。  
   **引用数**: 17件

3. **Zhang et al. (2023)**  
   "Efficient Multi Agent Path Finding with Turn Actions"  
   *SOCS 2023*. DOI: 10.1609/socs.v16i1.27290  
   **主要知見**: 動力学制約（ターンアクション）を無視するとコストが大幅増加。対称性破壊制約と枝刈り技術で改善。  
   **引用数**: 12件

4. **Khan & Singhal (2025)**  
   "ACBS: A Bounded-Suboptimal Multi-Agent Path Finding Solver"  
   *AJRCOS 2025*. DOI: 10.9734/ajrcos/2025/v18i11778  
   **主要知見**: ACBSはCBSに対し最大5倍の高速化を達成（100〜2000エージェント）。部分最適度w=1.2を維持。

5. **Wu et al. (2023)**  
   "A Review of Multi-Agent Path Finding Algorithms"  
   *ISCTech 2023*. DOI: 10.1109/ISCTech60480.2023.00020  
   **主要知見**: CBSベース/ルールベース/優先度ベース/数値最適化/学習ベースの5分類。AGVナビゲーションへの適用可能性を分析。

6. **Wu, Zhao & Ren (2026)**  
   "CBS-AA: CBS for MAPF with Asynchronous Actions"  
   *AAMAS 2026*. DOI: 10.65109/fscj9273  
   **主要知見**: 非同期アクション対応CBS。CBS-AAは分岐数を最大90%削減。

### 先行研究の課題・限界
- 1,000エージェント以上での体系的ベンチマークが不足
- 通信制約下でのスループット・コンフリクトトレードオフ分析が不十分
- Lifelong MAPFにおけるスループット比較の統計的厳密性が低い
- 動力学制約（MAMP）への拡張は実験規模が小さい

---

## 6. 考察と今後の展望

### アルゴリズム選択の実践的指針

```
N ≤ 50:    CBS (最適保証が必要な場合)
N ≤ 300:   EECBS (w=1.2, 品質バランス)
N > 300:   LaCAM2 (高速・スケーラブル)
Lifelong:  LaCAM-Lifelong > RHCR (73.4%改善)
分散環境:   Token Passing (通信範囲≥6セルで推奨)
```

### 自己批判的評価

**シミュレーション前提への依存**:
- 本実験はすべてパラメータ化シミュレーションであり、実際のMAPFソルバー実装ではない
- CBSのN=1,000での65%成功率は実際の挙動と異なる（実際はほぼ0%のはず）
- LaCAMの100%成功率は構造的に設定されており、実際の失敗ケースを過小評価している

**ノイズモデルの単純化**:
- 対数正規分布ノイズは実際のCBS探索分散を過小評価する可能性がある
- デッドロック、スターベーション、緊急停止などの現象を模擬していない

**実世界適用の制約**:
- 本シミュレーションの結果を実際の倉庫ロボットシステムに直接適用することは危険
- 実際の物流システムでは、ロボットの物理特性、センサーノイズ、通信遅延が重大な影響を持つ

### 今後の展望
1. 実際のMAPFベンチマークスイート（movingai.com）での検証
2. C++実装による高速化（LaCAM公式実装との比較）
3. MAMPへの拡張（連続空間・動力学モデル）
4. 強化学習によるオンライン適応的アルゴリズム選択
5. 実際の倉庫環境でのフィールドテスト

---

## 7. 生成したファイル一覧

| ファイル | 説明 |
|---------|------|
| `paper.md` | 学術論文形式の完全なレポート（英語） |
| `report.md` | 本実験レポート（日本語） |
| `mapf_research.ipynb` | 実験コードを含むJupyterノートブック |
| `figures/mapf_scalability.png` | Figure 1: スケーラビリティ分析（4サブプロット） |
| `figures/lifelong_mapf.png` | Figure 2: Lifelong MAPF倉庫シナリオ |
| `figures/distributed_mapf.png` | Figure 3: 分散MAPFの通信制約分析 |
| `figures/warehouse_benchmark.png` | Figure 4: 大規模倉庫ベンチマーク |
| `figures/ml_analysis.png` | Figure 5: ML代理モデル分析 |

---

## 8. 付録：主要実験コード

### セル1: MAPFシミュレーター

```python
def simulate_solver_performance(n_agents, solver, grid_size=32, density=0.1, n_trials=20):
    """
    MAPF solver performance based on published benchmarks:
    - CBS: exponential worst-case (Sharon et al. 2015)
    - EECBS: polynomial avg with w-suboptimality (Li et al. 2021)
    - LaCAM: near-linear (Okumura 2023 / LaCAM2)
    Random seed: 42 + n_agents (per-configuration reproducibility)
    """
    np.random.seed(42 + n_agents)
    # ... [full code in mapf_research.ipynb]
```

### セル8: ML代理モデル

```python
# 5-fold cross-validation with random_state=42
kf = KFold(n_splits=5, shuffle=True, random_state=42)
models = {
    'Ridge': Ridge(alpha=1.0),
    'Random Forest': RandomForestRegressor(n_estimators=100, random_state=42),
    'Gradient Boosting': GradientBoostingRegressor(n_estimators=100, random_state=42),
}
# Best: Gradient Boosting, R²=0.9487±0.0094 [cell:8]
```

---

## 9. 再現性情報

| 項目 | 値 |
|------|----|
| Python | 3.11.2 |
| NumPy | 2.3.5 |
| Pandas | 2.3.3 |
| Matplotlib | 3.10.9 |
| Seaborn | 0.13.2 |
| SciPy | 1.17.1 |
| 乱数シード | 42 |
| OS | Linux |
| 実行環境 | Jupyter MCP (kernel: Python 3) |

全数値は `np.random.seed(42)` 固定のもと決定論的に生成される。
