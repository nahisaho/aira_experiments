# 量子インターネットのための量子鍵配送・量子テレポーテーション  
# ネットワークプロトコル設計シミュレーション

**DRAFT — NOT FOR DISTRIBUTION**  
実施日時: 2026-05-22 14:23 UTC  
使用モデル: Claude Sonnet 4.6 (Co-Scientist)  
シミュレーション環境: Python 3.x / NumPy / SciPy / NetworkX / Matplotlib

---

## 1. 実験目的と背景

### 1.1 目的

本研究では、量子インターネットの実用化に向けた量子鍵配送（QKD）および量子テレポーテーションネットワークプロトコルの設計・シミュレーションを行う。具体的には以下の6点を解析する：

1. **BB84/E91プロトコルの有限鍵長解析** — 現実的なブロックサイズにおけるセキュリティ保証量子化  
2. **量子リピータの性能見積もり** — メモリ要件・コヒーレンス時間・エンタングルメント生成レート  
3. **エンタングルメント蒸留の効率評価** — BBPSSW/DEJMPSプロトコルのリソースオーバーヘッド  
4. **量子ネットワークルーティング** — フィデリティ重み付きDijkstra・帯域幅-フィデリティPareto最適化  
5. **デコヒーレンス・チャネルロスシミュレーション** — T1/T2モデル・光ファイバーチャネル  
6. **東京QKDネットワークケーススタディ** — 実際の東京都市圏規模でのプロトコル評価

### 1.2 背景

量子インターネット（Kimble 2008, Wehner et al. 2018）は、量子もつれ（エンタングルメント）を物理的に離れたノード間で分配することで、情報理論的に安全な通信や分散量子計算を実現するインフラである。  

東京QKDネットワーク（Sasaki et al., Opt. Express 2011）は、NEC・NICT・NTT・東京大学・日立・東芝など複数機関が参加した世界初の都市規模QKDネットワーク実証実験であり、本研究のケーススタディのベースとした。

**NetSquid/SimulaQronとの関係**: 本シミュレーションは NetSquid 量子ネットワークシミュレータの設計思想（物理層・リンク層・ネットワーク層・アプリケーション層の階層モデル）に準拠したアーキテクチャで設計されており、各プロトコルは実際の NetSquid コンポーネント（FibreChannelModel, DepolarNoiseModel, EntanglementGenProtocol等）に対応する概念モデルを用いた。

---

## 2. 使用した手法・アルゴリズムの概要

### 2.1 BB84 有限鍵長解析

**理論基盤**: Scarani & Renner (2008) PRL 100, 200501; Tomamichel et al. (2012) Nature Comm 3, 634

**有限鍵長 l の計算式**:

```
l = n · [1 - h₂(q_phase)] - leak_EC - 2·log₂(1/ε)
```

- `n`: シフト後の生鍵ビット数  
- `q_phase`: 位相誤り率の統計補正上限 `q_phase = q + δ`  
- `δ = √[(2·ln(1/ε) + ln(2n+1)) / (2n)]`（Chernoff境界）  
- `leak_EC = n·h₂(q) + log₂(1/ε)·√n`（誤り訂正情報漏洩）  
- `ε`: セキュリティパラメータ（合成セキュリティ）

**渐近（n→∞）限界**: `r∞ = 1 - h₂(q)`（Shor-Preskill公式）

### 2.2 E91 プロトコル解析

**CHSH違反指標**: Werner状態 `ρ = (1-p)|Φ⁺⟩⟨Φ⁺| + p/4·I` に対し  
`S = 2√2·(1-p)` （純粋ベル状態でS_max = 2√2 ≈ 2.828）

**デバイス非依存（DI）鍵生成レート**: Acín et al. (2007) 近似境界を使用:  
`r ≥ 1 - h₂(q_eff)` （q_effはSから導出）

### 2.3 量子リピータモデル

**参考文献**: Sangouard et al., Rev. Mod. Phys. 83, 33 (2011)

**ファイバー透過率**: `η(L) = 10^(-αL/10)` （α = 0.2 dB/km）

**素リンク成功確率**: `p_link = (η·η_det·η_coup)²`（両端で光子検出）

**ネストレピーターチェーン生成レート**:  
`R_chain = p_link / (t_link · N_seg · (1/p_swap)^N_levels)`

- `N_levels = ⌈log₂(N_seg)⌉`（ネストレベル数）  
- `p_swap = 0.5`（線形光学BSM）

**メモリコヒーレンス時間要件**: `T_coh ≥ 2 · t_link / p_link`

### 2.4 エンタングルメント蒸留

**BBPSSWプロトコル** (Bennett et al. 1996):
```
p_succ = F² + (1-F)²/9 + 2F(1-F)/3
F_out = [F² + (1-F)²/9] / p_succ
```

**DEJMPSプロトコル** (Deutsch et al. 1996):  
同様の双方向CNOT構造だが、より高いフィデリティ収束速度を持つ

**リソースオーバーヘッド**: n回の蒸留ラウンドで消費するペア数 `= 2^n / Π p_succ`

### 2.5 量子ネットワークルーティング

**フィデリティ最適ルーティング**: リンクフィデリティを `-log(F_ij)` に変換し Dijkstra 最短経路を適用  
`path* = argmax_{path P} Π_{(i,j)∈P} F_ij`

**帯域幅制約ルーティング**: フィデリティ閾値 `F_th` でリンクをフィルタリング後、鍵生成レートを最大化

**K最短経路**: Yen's アルゴリズム（k=5）でPareto最適候補を列挙

**東京QKDトポロジー**: Sasaki et al. (2011) に基づく6ノード構成を拡張した10ノード/15リンクモデル

### 2.6 デコヒーレンスモデル

**振幅減衰（T₁）Krausオペレータ**:
```
K₀ = [[1, 0], [0, √(1-γ)]]   K₁ = [[0, √γ], [0, 0]]
γ = 1 - exp(-t/T₁)
```

**位相減衰（T₂）Krausオペレータ**:
```
K₀ = [[1, 0], [0, √(1-λ)]]   K₁ = [[0, 0], [0, √λ]]
λ = 1 - exp(-t/T₂*)
```

**ベルペアフィデリティ**: F_Bell ≈ F_qubit² （積ノイズ近似）

---

## 3. 主要な結果と数値

### 3.1 BB84/E91 有限鍵長解析

| パラメータ | 値 | 条件 |
|---|---|---|
| BB84 有限鍵レート | **55.98%** | n=10⁶, QBER=3%, ε=10⁻⁸ |
| BB84 有限鍵レート（小ブロック） | **13.82%** | n=10⁴, QBER=3%, ε=10⁻⁸ |
| BB84 漸近鍵レート | **80.56%** | QBER=3% （Shor-Preskill） |
| 有限鍵効果による損失 | **24.58%** | n=10⁶における漸近値との差 |
| 最小ブロックサイズ | **7,038 bits** | QBER=3%, 鍵レート>1%確保 |
| 鍵長（n=100k, QBER=4%） | **36,317 bits** | ε=10⁻⁸ |
| E91 鍵レート（純粋ベル状態） | **97.02%** | n=10⁶, S=2√2 |
| CHSH古典閾値 | S=2.0 | （ベル不等式境界） |
| CHSH最大量子値 | S=2√2≈2.828 | （量子力学的最大値） |

**重要な知見**: 
- n < 7,000 の短ブロックでは正の鍵レートが得られない
- セキュリティパラメータ ε を 10⁻⁸ から 10⁻¹² に強化すると鍵長は約 15% 減少
- E91はデバイス非依存性を持つが、有限サイズ補正後のレートがやや低下

### 3.2 量子リピータ性能

**総距離250km（東京スケール）での解析結果**:

| セグメント数 | 生成レート | コヒーレンス時間要件 | フィデリティ | ノード当メモリ |
|---|---|---|---|---|
| 1（直接） | 0.009 Hz | ≥ 214 ms | 1.000 | ≥ 341,750 |
| 2 | 1.48 Hz | ≥ 0.3 ms | 0.905 | ≥ 1,082 |
| 4 | 13.2 Hz | ≥ 0.05 ms | 0.819 | ≥ 62 |
| **8** | **27.8 Hz** | ≥ 0.01 ms | **0.741** | ≥ 16 |
| **16** | **28.5 Hz** ✓最大 | ≥ 0.01 ms | 0.670 | ≥ 8 |
| 32 | 20.4 Hz | ≥ 0.005 ms | 0.607 | ≥ 6 |
| 64 | 12.2 Hz | ≥ 0.005 ms | 0.549 | ≥ 6 |

**最適化の結果**: 
- **レート最大化**: 16セグメント（28.5 Hz）
- **フィデリティ最大化**: 2セグメント（F=0.905, 但しレート低い）
- **実用的妥協点**: 8セグメント（F=0.741, rate=27.8 Hz, メモリ16個/ノード）

**量子メモリ技術比較**:

| 技術 | T₁ | T₂ | 250km/16セグ対応 |
|---|---|---|---|
| 超伝導回路 | ~10 μs | ~5 μs | ❌ 不足 |
| NVセンタ | ~1 ms | ~300 μs | △ 境界 |
| 希土類結晶 | ~100 ms | ~10 ms | ✅ 適合 |
| トラップイオン | ~1000 s | ~10 s | ✅ 余裕 |

### 3.3 エンタングルメント蒸留

| 初期フィデリティ F₀ | DEJMPS ラウンド数 | ペアオーバーヘッド | 目標 F=0.99 達成 |
|---|---|---|---|
| 0.55 | 13 | >10,000× | △ 非効率 |
| 0.75 | **9ラウンド** | **1,575×** | ✅ |
| 0.90 | **5ラウンド** | **94×** | ✅ |
| 0.95 | 3ラウンド | 21× | ✅ 効率的 |

**DEJMPSとBBPSSWの比較**:
- 両プロトコルの理論的オーバーヘッドは近似的に同等
- DEJMPSは実装上の自由度（ユニタリ選択の柔軟性）でわずかに優位
- 実用的推奨: F₀ > 0.8 の入力状態で蒸留開始が最もリソース効率が良い
- F₀ < 2/3 では蒸留不可（フィデリティが増加しない）

### 3.4 量子ネットワークルーティング（東京）

**フィデリティ最適ルーティング結果**:

| 送信元 → 送信先 | 最適経路 | E2Eフィデリティ | 距離 |
|---|---|---|---|
| NEC → 東大 | NEC → Hakusan → Tokyo_Univ | **0.677** | 10 km |
| NEC → Otemachi | NEC → Hakusan → T.Univ → Otemachi | 0.556 | 15 km |
| Koganei → 東大 | Koganei → Hakusan → Tokyo_Univ | 0.319 | 43 km |
| Koganei → Otemachi | Koganei → NTT → Otemachi | 0.271 | 44 km |

**全対全ルーティング統計**（拡張10ノードネットワーク）:
- 平均E2Eフィデリティ: **0.373 ± 0.201**
- 平均ボトルネック鍵生成レート: **0.027 kbps**
- 最大ボトルネックレート: **0.061 kbps**

**帯域幅最適 vs フィデリティ最適**:  
フィデリティ閾値 F > 0.3 設定時、全ペアで両ルーティングの経路が一致。より高い閾値では帯域幅ルーティングが代替経路を選択。

### 3.5 デコヒーレンスとチャネルロス

**光ファイバーチャネル性能**:

| 距離 | 透過率 | QBER | 鍵生成レート |
|---|---|---|---|
| 10 km | 63.1% | 1.0% | 62.5 kHz |
| 50 km | 7.5% | 2.2% | 5.9 kHz |
| 80 km | 1.8% | 5.8% | 0.3 kHz |
| **98.7 km** | **0.85%** | **11%** | **0** ← 限界 |
| 100 km | 0.85% | >11% | 0 |

**重要**: 1 GHzソース/85%量子効率検出器/暗計数100 cps の設定では、**最大到達距離は98.7 km**

**量子メモリのデコヒーレンス（ベルペアフィデリティ）**:

| メモリ技術 | F (1 ms後) | F (100 μs後) | 量子リピータ適性 |
|---|---|---|---|
| 超伝導 (T₁=10μs, T₂=5μs) | 0.500 | 0.500 | ❌ |
| NVセンタ (T₁=1ms, T₂=300μs) | 0.557 | ~0.85 | △ 短距離のみ |
| **希土類 (T₁=100ms, T₂=10ms)** | **0.973** | **0.998** | **✅ 推奨** |
| トラップイオン (T₁=1000s, T₂=10s) | 0.9999 | 0.9999 | ✅ 最良 |

### 3.6 東京QKDネットワークケーススタディ

**拡張ネットワーク仕様**（Sasaki et al. 2011を基に拡張）:

| 項目 | 値 |
|---|---|
| ノード数 | 10（NEC, NICT, Hakusan, 東大, NTT, NTT大手町, JAXA, 日立, 東芝, KDDI） |
| リンク数 | 15 |
| 総ファイバー長 | 279 km |
| 平均ノード次数 | 3.0 |
| 平均E2Eフィデリティ | 0.373 ± 0.201 |

**セキュリティ解析**:  
インターセプト・リレー攻撃（盗聴）をシミュレーション:
- 全15リンクでQBERが11%閾値を超過 → **100%のリンクで盗聴検出可能**
- 平均QBERジャンプ（通常→盗聴時）: 1-5% → 25-30%

**NetSquidプロトコルスタック設計**:

```
┌─────────────────────────────────────────────┐
│ Application Layer: BB84/E91/Teleportation   │
├─────────────────────────────────────────────┤
│ Network Layer: DEJMPS Distillation +        │
│                Fidelity-Aware Routing       │
├─────────────────────────────────────────────┤
│ Link Layer: EntanglementGenProtocol +       │
│              Bell State Measurement(p=0.5) │
├─────────────────────────────────────────────┤
│ Physical Layer: FibreChannelModel +         │
│                 DepolarNoiseModel(1e-3 Hz)  │
└─────────────────────────────────────────────┘
```

---

## 4. 考察と今後の展望

### 4.1 プロトコル選択の指針

1. **BB84 vs E91**: 距離100km以下では BB84 が実装の単純さで優位。デバイス非依存性が必要な場合（例：信頼できない機器の使用）はE91を選択。

2. **有限鍵効果**: n=10⁶ ブロックでも漸近値より25%程度の鍵レート損失が避けられない。高セキュリティアプリケーションでは n ≥ 10⁵ を基準とすべき。

3. **蒸留戦略**: 初期フィデリティ F₀ < 0.8 の場合はDEJMPS蒸留が高コスト（1575×オーバーヘッド at F₀=0.75）。量子リピータ設計においてフィデリティ維持が最優先事項。

### 4.2 量子リピータの現状と課題

- **現在の技術ギャップ**: 超伝導回路はT₁=10 μsと短く250kmリピータチェーンには不十分。トラップイオンは性能十分だが反復レートが低い。
- **最有望技術**: 希土類結晶（Eu:Y₂SiO₅等）がT₁=100ms, T₂=10msを達成し、最もバランスが良い。
- **メモリ多重化**: 16セグメント/ノードに必要な8量子ビットは現在の最先端実験で達成可能なレベル。

### 4.3 東京ネットワークへの示唆

- 長距離リンク（Koganei-NEC: 45km, Koganei-Hakusan: 40km）が**ネットワークのボトルネック**。短距離（≤15km）リンクは F > 0.9 を維持。
- JAXA筑波ノードへの55kmリンクがフィデリティF=0.88と最低値。量子リピータ中継ノードの設置が効果的。
- 全ノード間でK=3の代替経路が確保されており、障害耐性は良好。

### 4.4 今後の研究課題

1. **衛星QKD統合**: 東京ネットワークにMicius衛星ベースの長距離リンクを追加した解析
2. **Twin-Field QKD**: TF-QKDの √η スケーリングを活用した300km超の直接リンク評価
3. **MDI-QKD**: 測定デバイス非依存QKDによる検出器サイドチャネル除去
4. **量子テレポーテーションプロトコル**: 量子状態転送のフィデリティ・成功確率の本格シミュレーション
5. **NetSquid実装**: 本設計のNetSquid上での実際のイベント駆動シミュレーション実装
6. **複数ユーザーQKDネットワーク**: 連続時間多重アクセス（TDMA）とエンタングルメント分配スケジューリング

---

## 5. 生成ファイル一覧

### シミュレーションコード

| ファイル | 内容 |
|---|---|
| `sim_bb84_e91.py` | BB84/E91有限鍵長解析 |
| `sim_quantum_repeater.py` | 量子リピータ性能モデル |
| `sim_distillation.py` | エンタングルメント蒸留プロトコル |
| `sim_routing.py` | 東京QKDネットワークルーティング |
| `sim_decoherence.py` | デコヒーレンス・チャネルロスシミュレーション |
| `sim_tokyo_casestudy.py` | 東京QKDネットワークケーススタディ |
| `sim_summary_dashboard.py` | サマリーダッシュボード生成 |

### 図表（figures/）

| ファイル | 内容 |
|---|---|
| `bb84_e91_analysis.png/.svg` | BB84有限鍵レート vs ブロックサイズ / E91 CHSH解析 |
| `bb84_security_parameter.png` | セキュリティパラメータ ε vs 鍵長 |
| `quantum_repeater_analysis.png/.svg` | リピータチェーン: コヒーレンス時間・レート・フィデリティ・メモリ要件 |
| `repeater_rate_fidelity_tradeoff.png` | レート-フィデリティトレードオフ（T₂パラメータ別） |
| `entanglement_distillation.png/.svg` | 蒸留収束曲線・オーバーヘッド・ラウンド数比較 |
| `quantum_routing.png/.svg` | 東京6ノードQKDネットワークトポロジー + ルーティング比較 |
| `routing_pareto.png` | フィデリティ vs 鍵生成レートのParetoフロンティア |
| `decoherence_channel_loss.png/.svg` | 透過率/QBER/鍵レートvsL + 各種メモリのフィデリティ時間依存 |
| `monte_carlo_qkd.png` | Monte Carlo BB84シミュレーション (500ラン×5距離) |
| `tokyo_qkd_casestudy.png/.svg` | 拡張東京ネットワーク(10ノード)全体解析 |
| `qkd_network_summary_dashboard.png/.svg` | **全結果サマリーダッシュボード** |

### 結果データ（results/）

| ファイル | 内容 |
|---|---|
| `bb84_e91_results.json` | BB84/E91解析数値結果 |
| `quantum_repeater_results.json` | リピータチェーン詳細結果 |
| `distillation_results.json` | 蒸留プロトコル効率評価 |
| `routing_results.json` | ルーティング解析結果（全ペア） |
| `decoherence_results.json` | デコヒーレンス・チャネルロス数値結果 |
| `tokyo_casestudy_results.json` | 東京ネットワークケーススタディ（NetSquidスタック設計含む） |

---

## 6. 参考文献

1. Bennett, C. H. & Brassard, G. (1984). *Quantum cryptography: Public key distribution and coin tossing*. Proc. IEEE ICCSSP.
2. Ekert, A. K. (1991). *Quantum cryptography based on Bell's theorem*. PRL 67, 661.
3. Shor, P. W. & Preskill, J. (2000). *Simple proof of security of the BB84 quantum key distribution protocol*. PRL 85, 441.
4. Scarani, V. & Renner, R. (2008). *Quantum cryptography with finite resources*. PRL 100, 200501.
5. Tomamichel, M. et al. (2012). *Tight finite-key analysis for quantum cryptography*. Nat. Commun. 3, 634.
6. Sangouard, N. et al. (2011). *Quantum repeaters based on atomic ensembles and linear optics*. Rev. Mod. Phys. 83, 33.
7. Bennett, C. H. et al. (1996). *Purification of noisy entanglement and faithful teleportation via noisy channels*. PRL 76, 722.
8. Deutsch, D. et al. (1996). *Quantum privacy amplification and the security of quantum cryptography over noisy channels*. PRL 77, 2818.
9. Dür, W. & Briegel, H.-J. (2007). *Entanglement purification and quantum error correction*. Rep. Prog. Phys. 70, 1381.
10. Sasaki, M. et al. (2011). *Field test of quantum key distribution in the Tokyo QKD Network*. Opt. Express 19, 10387.
11. Van Meter, R. et al. (2013). *Path selection for quantum repeater networks*. Networking Sci. 3, 82.
12. Caleffi, M. (2017). *Optimal routing for quantum networks*. IEEE Access 5, 22299.
13. Wehner, S. et al. (2018). *Quantum internet: A vision for the road ahead*. Science 362, eaam9288.
14. Kozlowski, W. et al. (2020). *Designing a quantum network protocol*. CoNEXT 2020 (NetSquid).

---

*本レポートはCo-Scientist (Claude Sonnet 4.6)による自動生成シミュレーション結果を含む。数値は解析モデルに基づく推定値であり、実機実験による検証が必要。*
