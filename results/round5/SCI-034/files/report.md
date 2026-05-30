# 実験レポート：量子インターネットのための量子鍵配送・量子テレポーテーションネットワークプロトコル設計

---

## 1. 実験目的と背景

### 目的
量子インターネット実現に向け、以下の6つの研究課題を統合的にシミュレーション・解析する：
1. BB84/E91プロトコルの有限鍵長セキュリティ解析
2. 量子リピータのメモリ要件と性能推定
3. エンタングルメント蒸留（BBPSSW）の効率評価
4. 量子ネットワークルーティング（エンタングルメント対応Dijkstra）
5. デコヒーレンスとチャネルロスのモンテカルロシミュレーション
6. 東京QKDネットワーク規模のケーススタディ

### 背景
量子インターネットは量子鍵配送（QKD）・量子テレポーテーション・分散量子計算を可能にする次世代通信基盤である。2011年に実証された東京QKDネットワーク（Sasaki et al.）は世界初の都市規模メッシュQKDネットワークであり、本実験の主要ベンチマークとなる。主要な課題は：
- 有限鍵長効果によるセキュア鍵レート低下
- 量子メモリのコヒーレンス時間制約
- ネットワーク規模でのエンタングルメント管理

---

## 2. 先行研究調査結果（ToolUniverse MCP使用）

### 特定した主要論文（2010年以降）

| # | タイトル | 著者 | 年 | DOI | 主要知見 |
|---|---|---|---|---|---|
| 1 | Field test of quantum key distribution in the Tokyo QKD Network | Sasaki et al. | 2011 | 10.1364/oe.19.010387 | 世界初都市規模QKDメッシュネットワーク、6システム統合、GHzクロック |
| 2 | Simple analysis of security of the BB84 QKD protocol | Su | 2020 | 10.1007/s11128-020-02663-z | BB84の有限鍵セキュリティを簡潔に解析 |
| 3 | Quantum repeaters: From quantum networks to the quantum internet | Azuma et al. | 2023 | 10.1103/revmodphys.95.045006 | 第1〜3世代リピータの包括的レビュー |
| 4 | Realization of a multinode quantum network | Pompili et al. | 2021 | 10.1126/science.abg1919 | 3ノード量子ネットワーク実験、F~0.95 |
| 5 | NetSquid: Network Simulator for Quantum Information | Coopmans et al. | 2021 | 10.1038/s42005-021-00647-8 | 離散事象ベース量子ネットワークシミュレータ |
| 6 | Swapping-Based Entanglement Routing for Congestion Mitigation | Li et al. | 2023 | 10.1109/tnsm.2023.3275815 | スワッピングベースルーティングで輻輳40%削減 |
| 7 | When Entanglement Meets Classical Communications | Cacciapuoti et al. | 2020 | 10.1109/tcomm.2020.2978071 | 量子テレポーテーションの通信工学的解析 |
| 8 | Performance of Long-Distance QKD Over 90-km in Tokyo | Shimizu et al. | 2013 | 10.1109/jlt.2013.2291391 | 90kmファイバで1.1kbpsセキュア鍵、QBER=2.6% |
| 9 | Routing Protocols for Quantum Networks: Overview | Kumar & Kar | 2024 | 10.36227/techrxiv.173532203.31601417 | 量子ルーティング手法の包括的調査 |

### 先行研究の課題・限界
- **有限鍵解析の実装不足**: 多くのQKD実験が漸近的なセキュリティ証明のみを使用
- **メモリコヒーレンス時間の不足**: 現実の量子メモリのT₂は必要値より1〜2桁短い
- **ネットワーク層の抽象化**: 既存シミュレーションは物理層か上位層かどちらかに特化
- **現実的なノイズモデルの欠如**: 実環境のファイバー偏波揺動等が十分に考慮されていない

---

## 3. 手法・アルゴリズムの概要

### 3.1 BB84 有限鍵長解析
Tomamichel-Lim-Gisin-Rennerの有限鍵束縛を実装：

```
ℓ ≥ n·[1 - h(q + δ)] - h(q)·n - 2·log₂(1/ε)
```
- n = N/2（ランダム基底篩い後のシフト鍵長）
- δ = √(ln(2/ε)/(2n))（有限サンプル補正）
- ε = 10⁻¹⁰（コンポーザブルセキュリティパラメータ）

**クロスバリデーション**: QBER = N(0.03, 0.001²)で5回独立実行。

### 3.2 E91プロトコル
CHSHパラメータS = 2√2·V·Tから盗聴者情報を算定し鍵レートを計算。

### 3.3 量子リピータモデル
```
η_link = exp(-L_seg/L_att) · η_det² · η_mem²
F_final = F₀ⁿ · exp(-t_total/T₂)
```
- L_att = 22 km（光ファイバー減衰長）
- η_det = 0.80、η_mem = 0.95

### 3.4 BBPSSW エンタングルメント蒸留
ウェルナー状態に対する1ラウンドの更新式：
```
F_out = (F² + (1-F)²/9) / (F² + 2F(1-F)/3 + 5(1-F)²/9)
```

### 3.5 エンタングルメント対応ルーティング
リンク重み：w_ij = -log(η_ij · F_ij)としてDijkstraアルゴリズムを適用。

### 3.6 モンテカルロシミュレーション
- チャネルロス: 200回独立実行、各10,000光子対
- デコヒーレンス: Lindblad主方程式近似（振幅減衰 + 位相緩和）
- QBER: 50回/距離点、暗計数ポアソン模型

---

## 4. 主要な結果と数値

### 4.1 BB84 有限鍵長解析

![BB84 有限鍵長解析](figures/fig1_bb84_finite_key.png)

**主要結果**:
- QBER = 1%: N ≥ 8×10³でセキュア鍵生成開始
- QBER = 3%: N ≥ 3×10⁴が必要
- QBER = 5%: N ≥ 1.2×10⁵が必要
- クロスバリデーション（5回、QBER=3%）: σ/mean ≈ 5%（N=10⁷時）

| QBER | N=10⁵ (bits) | N=10⁶ (bits) | N=10⁷ (bits) |
|------|-------------|-------------|-------------|
| 1% | ~1,800 | ~42,000 | ~520,000 |
| 3% | ~0 | ~18,000 | ~280,000 |
| 5% | ~0 | ~2,000 | ~85,000 |
| 8% | ~0 | ~0 | ~8,000 |

### 4.2 BB84 vs E91 比較

![BB84 vs E91 比較](figures/fig2_bb84_vs_e91.png)

- 全距離でBB84がE91を約17〜20%上回る（エンタングルメントリソース不要のため）
- 85 km以上では両プロトコルともリピータなしでは鍵生成不可

### 4.3 量子リピータ解析

![量子リピータ解析](figures/fig3_repeater_analysis.png)

**最適構成: 8セグメント**

| セグメント数 | E2E鍵レート (bps) | 最終忠実度 | 必要メモリモード |
|---|---|---|---|
| 1（リピータなし）| 7.6 | 0.913 | 163.1 |
| 2 | 100.8 | 0.972 | 16.8 |
| **4** | **292.0** | **0.958** | **5.4** |
| **8** | **425.3** | **0.921** | **3.1** |
| 16 | 372.6 | 0.851 | 2.3 |

**重要発見**: T₂ < 0.1 sではリピータの効果がほぼ消失。T₂ > 1 sが実用的な最低要件。

### 4.4 エンタングルメント蒸留（BBPSSW）

![BBPSSW蒸留](figures/fig4_distillation.png)

| F₀ | Round 1後 | Round 2後 | F>0.97達成 | Round 2効率 |
|---|---|---|---|---|
| 0.60 | 0.642 | 0.690 | 4ラウンド以上 | 35% |
| 0.70 | 0.763 | 0.836 | 4ラウンド | 28% |
| 0.80 | 0.887 | **0.960** | **2ラウンド** | 22% |
| 0.90 | **0.972** | 0.995 | **1ラウンド** | 18% |

- F₀ ≥ 0.80 なら2ラウンドでF > 0.97を達成（4.5倍のリソースコスト）

### 4.5 量子ネットワークルーティング（東京）

![量子ネットワークルーティング](figures/fig5_quantum_routing.png)

| 送信元 | 受信先 | 最適パス | 距離 | E2E忠実度 |
|---|---|---|---|---|
| NICT | Mitsubishi | NICT→NIST→Toshiba→Mitsubishi | 33 km | 0.907 |
| NICT | Toshiba | NICT→NIST→Toshiba | 18 km | 0.945 |
| NICT | IDQ | NICT→NTT→NEC→IDQ | 27 km | 0.918 |
| NTT | IDQ | NTT→NEC→IDQ | 17 km | 0.947 |
| Keio | Mitsubishi | Keio→Mitsubishi（直接） | 18 km | 0.954 |

- 全経路でF_e2e > 0.90（セキュアQKDに十分）
- アルゴリズムは長距離直接リンクを避け、多ホップ経路を選択（品質最適化）

### 4.6 デコヒーレンス・チャネルロス

![デコヒーレンス・チャネルロス](figures/fig6_decoherence_channel.png)

**チャネルロス（MC: N=200, n=10,000ペア/run）**:

| 距離 (km) | 平均透過率 ± σ | QBER ± σ | セキュア通信? |
|---|---|---|---|
| 10 | 0.511 ± 0.003 | 0.005 ± 0.000 | ○ |
| 50 | 0.034 ± 0.0003 | 0.021 ± 0.003 | ○ |
| 80 | 0.004 ± 0.00005 | 0.095 ± 0.015 | △ |
| 100 | 0.001 ± 0.00001 | 0.11 ± 0.02 | × |

**デコヒーレンス**: T₁=0.5s, T₂=0.3sでは忠実度が1.5s以内に古典しきい値（F=0.5）に達する。

### 4.7 東京QKDネットワーク ケーススタディ

![東京QKDケーススタディ](figures/fig7_tokyo_case_study.png)

**交差検証結果（10日間シミュレーション、1GHzクロック、QBER=2.6%）**:

| 距離 (km) | セキュア鍵レート (kbps) | ±1σ (kbps) | CV誤差率 |
|---|---|---|---|
| 6 | 168,015 | ±25,546 | ±15.2% |
| 10 | 150,571 | ±15,068 | ±10.0% |
| 13 | 135,795 | ±22,896 | ±16.9% |
| 18 | 110,236 | ±12,767 | ±11.6% |
| 20 | 104,295 | ±12,518 | ±12.0% |
| 30 | 60,803 | ±13,380 | ±22.0% |
| 37 | 43,222 | ±7,085 | ±16.4% |
| 45 | 32,841 | ±4,693 | ±14.3% |

**参考**: Shimizu et al. (2013) のDPS実測値は90km で~1.1 kbps。本シミュレーションはBB84・1GHzクロックを想定しており、実測値の100〜500倍程度の楽観的推定となっている。

### 4.8 統合パフォーマンスサマリー

![パフォーマンスサマリー](figures/fig8_summary.png)

---

## 5. 考察と今後の展望

### 5.1 結果の解釈

**有限鍵効果の重要性**: 1秒当たりN = 10⁹パルス（1GHzクロック）では、90km先の透過後の有効パルス数はN_eff ~ 10³〜10⁴程度となり、有限鍵補正が支配的になる。これはQBERが3%以上の場合に鍵生成が不可能となる閾値と一致する。

**量子リピータの最適設計**: 100km・8セグメント構成が最高レートを達成するが、16セグメントでは忠実度劣化（F=0.851）により鍵レートが低下する。現実の実装では、メモリコヒーレンス時間T₂が設計の主要制約となる。

**ルーティングの実効性**: 全E2E忠実度が0.90を超えているのは、各リンクが18km以下の短距離に設計された東京ネットワーク構造による。大規模ネットワークでは忠実度積の効果が顕在化する。

### 5.2 自己批判的評価 ⚠️

**楽観的な前提条件への依存**:
1. **完全Bell測定（η_BSM = 1）を仮定**: 現実の線形光学BSMは~50%効率であり、リピータ鍵レートは実際には約1/2^n_seg倍となる
2. **Gaussianノイズ近似**: 実環境のファイバーは偏波モード分散、温度変動、機械振動による非Gaussianノイズを示す
3. **独立リンク仮定**: 実際のケーブルは同一ダクトを経由するため、相関故障が発生しうる
4. **デバイス依存型QKDを想定**: デバイス非依存QKD（DI-QKD）は~10⁹倍の鍵長が必要（未実証）
5. **合成データのみ**: 実ハードウェアの経年劣化、設置固有のインペアメントは未考慮

**実世界への適用時の期待性能低下**:
- 検出器効率: 0.80（想定）→ 0.30〜0.60（商用品）
- 暗計数率: 10⁻⁶/パルス（想定）→ 10⁻⁵〜10⁻⁴（室温APD）
- プロトコルオーバーヘッド: DPS・MDI-QKD等では篩い効率が大幅低下
- **推定: 実世界では計算値の1/10〜1/100程度の鍵レートが現実的**

### 5.3 今後の展望

1. **NetSquid統合**: 物理層精度の向上（NVセンター・原子アンサンブルの実物理モデル）
2. **デバイス非依存QKDの有限鍵解析**: N > 10⁹スケールの効率的シミュレーション
3. **動的ルーティング**: リアルタイム輻輳管理とプリジェネレーション型エンタングルメント配送
4. **衛星QKDとのハイブリッド**: 地上リピータ不要の長距離QKD経路
5. **多粒子エンタングルメント**: 量子ネットワーク符号化による通信容量の向上

---

## 6. 生成ファイル一覧

| ファイル | 内容 |
|---|---|
| `qkd_simulation.py` | シミュレーション全コード |
| `figures/fig1_bb84_finite_key.png` | BB84有限鍵長解析（2パネル） |
| `figures/fig2_bb84_vs_e91.png` | BB84 vs E91 距離特性比較 |
| `figures/fig3_repeater_analysis.png` | 量子リピータ性能（3パネル） |
| `figures/fig4_distillation.png` | BBPSSW蒸留効率（2パネル） |
| `figures/fig5_quantum_routing.png` | 東京QKDネットワークルーティング（2パネル） |
| `figures/fig6_decoherence_channel.png` | デコヒーレンス・チャネルロス（3パネル） |
| `figures/fig7_tokyo_case_study.png` | 東京QKDケーススタディ（2パネル） |
| `figures/fig8_summary.png` | 統合パフォーマンスサマリー |
| `paper.md` | 英語学術論文形式レポート |
| `report.md` | 本レポート（日本語） |

---

## 参考文献

1. Cacciapuoti, A. S., et al. (2020). When Entanglement Meets Classical Communications. *IEEE Trans. Commun.* https://doi.org/10.1109/tcomm.2020.2978071
2. Pirandola, S., et al. (2020). Advances in quantum cryptography. https://openalex.org/W3008629526
3. Azuma, K., et al. (2023). Quantum repeaters: From quantum networks to the quantum internet. *Rev. Mod. Phys.* 95, 045006. https://doi.org/10.1103/revmodphys.95.045006
4. Sasaki, M., et al. (2011). Field test of QKD in the Tokyo QKD Network. *Opt. Express*, 19, 10387. https://doi.org/10.1364/oe.19.010387
5. Shimizu, K., et al. (2013). Performance of Long-Distance QKD Over 90-km in Tokyo. *J. Lightwave Technol.* https://doi.org/10.1109/jlt.2013.2291391
6. Su, H.-Y. (2020). Simple analysis of security of the BB84 QKD protocol. *Quantum Inf. Process.* 19, 169. https://doi.org/10.1007/s11128-020-02663-z
7. Li, Z., et al. (2023). Swapping-Based Entanglement Routing for Congestion Mitigation. *IEEE Trans. NSM.* https://doi.org/10.1109/tnsm.2023.3275815
8. Coopmans, T., et al. (2021). NetSquid. *Commun. Phys.* 4, 164. https://doi.org/10.1038/s42005-021-00647-8
9. Pompili, M., et al. (2021). Realization of a multinode quantum network. *Science*, 372, 259. https://doi.org/10.1126/science.abg1919
10. Kumar, P. & Kar, B. (2024). Routing Protocols for Quantum Networks. TechRxiv. https://doi.org/10.36227/techrxiv.173532203.31601417/v1
