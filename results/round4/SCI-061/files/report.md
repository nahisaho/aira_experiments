# AutoCircuit 実験レポート

## 実験目的と背景
合成遺伝子回路の自動設計において、文脈効果と確率的ゆらぎ下での頑健性を同時に扱う再現可能なフレームワークはまだ不足しています。
本実験では AutoCircuit フレームワークを構築し、トグルスイッチとリプレッシレーターをケーススタディとして再設計しました。

## 使用手法・アルゴリズムの概要
- **形式仕様**: LogicGate / FeedbackLoop / CircuitSpec による SBOL 風辞書表現
- **部品カタログ**: pTac, pTet, pBAD, pLac, pT7 プロモーター、strong/medium/weak RBS、rrnBT1/rrnBT2/lambda_t0 ターミネーター、LacI/TetR/AraC/cI リプレッサー
- **シミュレーション**: 厳密 Gillespie SSA と adaptive tau-leaping の両方を実装
- **最適化**: 個体数 20、30 世代の遺伝的アルゴリズム（トーナメント選択、1点交叉、点突然変異）
- **頑健性解析**: Latin Hypercube Sampling による 50 サンプル、±20% パラメータ摂動
- **文脈効果補正**: 終結配列 readthrough と遺伝子位置依存性による乗法補正モデル
- **MCP 利用状況**: Crossref_search_works は成功; SemanticScholar_search_papers は 2 件で HTTP 429 (rate limit) が発生

## 主要結果と数値

### 1. Toggle Switch 再設計 (n=20 runs)
| 条件 | Bistability Score | Switching Time (min) | Noise Resilience |
|---|---:|---:|---:|
| Original | 0.450 ± 0.034 | 100.00 ± 0.00 | 0.894 ± 0.023 |
| Optimized (GA) | 0.536 ± 0.077 | 47.10 ± 32.52 | 0.782 ± 0.129 |

### 2. Repressilator 再設計 (n=20 runs)
| 条件 | Oscillation Score | Period (min) | Amplitude CV | Synchrony Score |
|---|---:|---:|---:|---:|
| Original | 0.609 ± 0.100 | 483.40 ± 172.10 | 0.451 ± 0.189 | 0.352 ± 0.094 |
| Optimized (GA) | 0.618 ± 0.081 | 542.70 ± 123.73 | 0.303 ± 0.166 | 0.300 ± 0.000 |

### 3. 頑健性解析 (50 LHS サンプル, ±20% 摂動)
| 回路 | Mean Fitness ± Std | Worst-case |
|---|---:|---:|
| Toggle Original | 0.441 ± 0.040 | 0.344 |
| Toggle Optimized | 0.536 ± 0.088 | 0.376 |
| Repressilator Original | 0.565 ± 0.083 | 0.424 |
| Repressilator Optimized | 0.591 ± 0.088 | 0.436 |

### 4. 文脈効果補正 (3遺伝子カセット)
| 遺伝子位置 | Baseline | 補正なし | 補正あり |
|---|---:|---:|---:|
| 1 | 95.00 | 91.20 | 95.00 |
| 2 | 80.00 | 66.77 | 80.00 |
| 3 | 65.00 | 50.64 | 65.00 |

## 生成した図
![Figure 1: Toggle Switch Dynamics](figures/figure1_toggle_switch_dynamics.png)
![Figure 2: Repressilator Dynamics](figures/figure2_repressilator_dynamics.png)
![Figure 3: Robust Design Analysis](figures/figure3_robust_design.png)
![Figure 4: GA Optimization Convergence](figures/figure4_ga_optimization.png)
![Figure 5: Context Effects](figures/figure5_context_effects.png)
![Figure 6: Summary Comparison Heatmap](figures/figure6_comparison.png)

## 考察と今後の展望
- トグルスイッチでは GA 最適化後に bistability が改善し、切替時間も短縮しました。ただし値は完全ではなく、内在ノイズの影響が残っています。
- リプレッシレーターでは oscillation score がわずかに改善し、振幅 CV は低下しました。一方で synchrony は高くならず、内在ノイズ下で位相整合を保つ難しさが示されました。
- 文脈効果補正では下流遺伝子の発現低下をほぼ基準値まで回復でき、複数遺伝子カセット設計では補正が重要であることが確認できました。
- MCP ツールを用いた文献調査では Crossref は安定して利用できましたが、Semantic Scholar には 429 レート制限が発生しており、ツール依存の制約も確認されました。
- 今後の課題: リソース負荷モデルの統合、配列-パラメータ予測の精度向上、SBOL/Cello との直接連携。

## 生成ファイル一覧
- `gene_circuit_framework/__init__.py`
- `gene_circuit_framework/circuit_spec.py`
- `gene_circuit_framework/parts_catalog.py`
- `gene_circuit_framework/assembly.py`
- `gene_circuit_framework/simulation.py`
- `gene_circuit_framework/robust_design.py`
- `gene_circuit_framework/context_effects.py`
- `gene_circuit_framework/optimizer.py`
- `gene_circuit_framework/case_studies.py`
- `run_experiments.py`
- `results_summary.json`
- `paper.md`
- `report.md`
- `figures/figure1_toggle_switch_dynamics.png`
- `figures/figure2_repressilator_dynamics.png`
- `figures/figure3_robust_design.png`
- `figures/figure4_ga_optimization.png`
- `figures/figure5_context_effects.png`
- `figures/figure6_comparison.png`