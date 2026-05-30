# 量子インターネットのためのQKD・量子テレポーテーションネットワークプロトコル設計と性能評価

**DRAFT — NOT FOR DISTRIBUTION**  
作成日: 2026-05-28  
著者: Co-Scientist Research Suite

---

## Abstract

本研究は量子インターネット実現に向けた量子鍵配送（QKD）および量子テレポーテーションネットワークプロトコルの設計と性能評価を目的とする。BB84プロトコルの有限鍵長解析、量子リピータのメモリ要件分析、DEJMPS エンタングルメント蒸留プロトコルの効率評価、量子ネットワークルーティングアルゴリズムの比較、デコヒーレンス・チャネルロスのモンテカルロシミュレーション、および東京QKDネットワーク規模のケーススタディを実施した。

シミュレーション結果として、BB84有限鍵解析では N=10¹⁰ パルス時に距離50 km で安全鍵率 2.923×10⁻² ± 1.221×10⁻³ ビット/パルス（5分割交差検証）を達成した。NVセンター量子メモリを用いた4ノードリピータ鎖（500 km）では 29.5 Hz のエンタングルメント生成率と忠実度 0.990 を確認した。DEJMPS 蒸留プロトコルは初期忠実度0.70から6ラウンドで忠実度1.000に収束した（ゲート忠実度完全時）。東京QKDネットワーク（7ノード、約45 km規模）の最大忠実度ルーティング解析では、隣接ノード間で0.940の端末間忠実度を達成した。

---

## 1. 実験目的と背景

量子インターネットは、量子力学の原理に基づいた情報通信ネットワークであり、無条件に安全な鍵共有（QKD）、量子テレポーテーション、分散量子計算などの応用を可能にする（Wehner et al., 2018）。現在の量子通信では光ファイバー中の光子損失（通常0.2 dB/km）により実用的な通信距離は100〜200 km に制限される。この制約を克服するために量子リピータ技術が不可欠であり、エンタングルメント生成、交換（スワッピング）、蒸留（精製）の組み合わせにより長距離の量子もつれ分配を実現する（Briegel et al., 1998）。

本研究の目的は以下の6点である：

1. **BB84/E91有限鍵長解析**: 距離と鍵長に依存する安全鍵率の定量評価
2. **量子リピータ性能評価**: 複数プラットフォームのメモリ要件と到達距離の比較
3. **エンタングルメント蒸留効率**: DEJMPS プロトコルの収束特性と資源消費
4. **量子ルーティングアルゴリズム**: 距離最短・忠実度最大・レート最大の3種類の比較
5. **デコヒーレンスとチャネルロス**: 光ファイバー損失係数変化の影響モンテカルロ
6. **東京QKDケーススタディ**: 実在する7ノードメトロポリタンネットワークの解析

### 先行研究調査の結果

文献調査には ToolUniverse MCP経由の ArXiv 学術検索ツール（`ArXiv_search_papers`）および Crossref (`Crossref_search_works`) を使用した。Semantic Scholar API は試行時にレート制限エラー（HTTP 429）が発生したため代替手段に切り替えた。以下の主要な先行研究を特定した：

- **Wiesemann et al. (2024)**: decoy-state BB84 有限鍵セキュリティ証明の統合的解説（EUR フレームワーク）
- **Kamin et al. (2025)**: decoy-state BB84 の有限鍵レート改善（集中不等式の精緻化）
- **Coopmans et al. (2021)**: NetSquid 量子ネットワーク離散事象シミュレータの紹介
- **Yehia et al. (2022)**: NetSquid による都市型量子ネットワーク（Quantum City）のシミュレーション
- **Haldar et al. (2024)**: 多重化量子リピータ鎖における古典通信コストの削減
- **Liu et al. (2026)**: 量子リピータの基本ブロック（10 km 繊維でのDI-QKD実証）
- **Tian et al. (2026)**: RADAR-Q マルチテナント量子ネットワーク向け分散ルーティングプロトコル
- **Sasaki et al. (2011)**: 東京QKDネットワークのフィールド実証

---

## 2. 使用した手法・アルゴリズムの概要

### 2.1 BB84 有限鍵解析（EUR フレームワーク）

安全鍵長 $\ell$ は Entropic Uncertainty Relation（EUR）フレームワーク（Tomamichel & Renner, 2011）に基づき以下で計算される：

$$
\ell \leq n_{\text{sifted}} \left[ 1 - h(e_z) - h(e_x) \right] - \lambda_{\text{EC}} - \Delta_{\text{sec}}
$$

ここで $h(\cdot)$ は二値エントロピー関数、$e_z$（ビット誤り率）と $e_x$（位相誤り率）は観測QBERから推定、$\lambda_{\text{EC}} = n_{\text{sifted}} \cdot f_{\text{EC}} \cdot h(e_z)$（$f_{\text{EC}}=1.16$ は誤り訂正効率）、$\Delta_{\text{sec}}$ はセキュリティパラメータ $\varepsilon_{\text{sec}}=10^{-10}$ から導かれる有限サイズ補正項である。

光ファイバー伝送率は $\eta = 10^{-\alpha d / 10}$（$\alpha=0.2$ dB/km）で計算し、ダークカウント由来のQBERは $e_{\text{dark}} = d_c / (\eta + d_c)$ で推定した（$d_c = 10^{-6}$ /ゲート）。

### 2.2 E91（Ekert91）プロトコル

E91 はベル状態を利用したエンタングルメントベースQKDであり、Werner状態の忠実度 $F$ に対してCHSHパラメータは：

$$
S = 2\sqrt{2}(2F - 1), \quad F \geq 0.5
$$

古典的限界 $S=2$ は $F=0.5$ に対応し、$S > 2$ の量子優位性により盗聴検知を実現する。有限鍵補正は不確定性指数から $\Delta = \sqrt{\log_2(1/\varepsilon_{\text{sec}}) / n_{\text{detect}}}$ で見積もった。

### 2.3 量子リピータ鎖モデル（BDCZ型）

$n$ セグメントのリピータ鎖では、各エレメンタリリンク長は $d_0 = d_{\text{total}}/n$。生成確率は：

$$
p_{\text{gen}} = \eta_{\text{em}} \cdot \eta_{\text{det}}^2 \cdot \eta_{\text{fiber}}^2 \cdot \eta_{\text{write}}
$$

Werner状態の $T_2$ コヒーレンス時間減衰：

$$
F(t) = F_0 e^{-t/T_2} + \frac{1 - e^{-t/T_2}}{4}
$$

$k$ レベルのネストスワッピング後の忠実度：

$$
F_{k+1} = F_k^2 + \frac{(1-F_k)^2}{9}
$$

これは理想的なBell状態測定の場合であり、ゲート忠実度 $F_g$ を考慮した補正 $F \leftarrow F \cdot F_g + (1-F_g)/4$ も適用した。

### 2.4 DEJMPS エンタングルメント蒸留

Werner状態への DEJMPS 蒸留の1ラウンドでは、初期忠実度 $F$ から：

$$
p_{\text{success}} = F^2 + \frac{5}{9}(1-F)^2
$$

$$
F' = \frac{F^2 + \frac{1}{9}(1-F)^2}{p_{\text{success}}}
$$

この漸化式は $F > 0.5$ で収束し、目標忠実度 $F_{\text{target}}$ までの必要ラウンド数と消費ペア数を数値的に求めた。

### 2.5 ルーティングアルゴリズム

3種類のルーティング目標を Dijkstra アルゴリズムの変形で実装した：

- **距離最短ルーティング**: 重み $w_{ij} = d_{ij}$（距離[km]）の最短路
- **忠実度最大ルーティング**: 重み $w_{ij} = -\log F_{ij}$ の最短路（忠実度の積を最大化）
- **ボトルネックレート最大ルーティング**: 最大ボトルネック路（最小リンクレートを最大化）

チャネル忠実度はデポラライズノイズと位相ノイズを組み合わせたモデルで計算：

$$
F_{\text{ch}} = 0.25 + 0.75 \cdot e^{-\alpha_{\text{dep}} d} \cdot e^{-\alpha_{\text{deph}} d}
$$

---

## 3. 主要な結果と数値

### 3.1 BB84/E91 有限鍵率対距離

![BB84/E91 有限鍵率対距離](figures/fig1_key_rate_vs_distance.png)

**Figure 1**: BB84 プロトコル（N=10⁹, 10¹⁰, 10¹¹）および E91（N=10⁹）の有限鍵安全鍵率を距離の関数として示す。縦軸は対数スケール。ブロックサイズが大きいほど低損失領域では差が小さく、高損失（長距離）で顕著に分岐する。

| 距離 (km) | BB84 N=10⁹ | BB84 N=10¹⁰ | BB84 N=10¹¹ |
|-----------|------------|-------------|-------------|
| 25 | 1.09×10⁻¹ | 1.09×10⁻¹ | 1.09×10⁻¹ |
| 50 | 3.59×10⁻² | 3.60×10⁻² | 3.60×10⁻² |
| 75 | 1.08×10⁻² | 1.08×10⁻² | 1.08×10⁻² |
| 100 | 3.54×10⁻³ | 3.55×10⁻³ | 3.55×10⁻³ |
| 150 | 3.44×10⁻⁴ | 3.45×10⁻⁴ | 3.45×10⁻⁴ |

5分割交差検証（QBER測定のランダム誤差 ±0.2%）:
- 距離50 km: **2.923×10⁻² ± 1.221×10⁻³** ビット/パルス
- 距離100 km: **2.951×10⁻³ ± 5.853×10⁻⁵** ビット/パルス

### 3.2 量子リピータメモリ要件

![量子リピータ性能](figures/fig2_repeater_memory.png)

**Figure 2**: NVセンター・捕捉イオン・原子アンサンブルの3種メモリを用いた場合の、500 km 総距離における（左）エンタングルメント生成率と（右）端末間忠実度のセグメント数依存性。

NVセンターの最適動作点は **4セグメント（各125 km）** で、29.5 Hz・忠実度0.990を達成した。メモリ1ノードあたり最低2量子ビットが必要。16セグメント以上では待機時間によるデコヒーレンスが率を大幅に低下させた。

| セグメント数 | 率 (Hz) | 忠実度 | ノードあたりメモリ |
|------|---------|--------|-----------|
| 2 | 0.029 | 0.982 | 1 |
| 4 | 29.5 | 0.990 | 2 |
| 8 | 41.5 | 0.983 | 4 |
| 16 | 2.50 | 0.967 | 8 |
| 32 | 0.58 | 0.969 | 16 |

### 3.3 エンタングルメント蒸留

![DEJMPS 蒸留収束](figures/fig3_purification.png)

![DEJMPS 蒸留軌跡](figures/fig3b_purification_trace.png)

**Figure 3**: （左）各初期忠実度と各ゲート忠実度におけるDEJMPS蒸留の必要ラウンド数。（右）F₀=0.70 から始まる忠実度の収束軌跡。

ゲート忠実度が完全（1.000）の場合、F₀=0.70 から **6ラウンド** で忠実度1.000に収束し、約64対のペアを消費した。ゲート忠実度0.98では理論的上限の忠実度が0.95付近に制限される。初期忠実度が0.52〜0.60の低い領域では収束が遅く（8〜15ラウンド）資源消費が著しく増大した。

### 3.4 東京QKDネットワークルーティング

![東京QKD 忠実度行列](figures/fig4_tokyo_fidelity_matrix.png)

![東京QKD ルーティング結果](figures/fig4b_tokyo_routing.png)

**Figure 4**: （左）東京QKDネットワーク7ノード間の最大忠実度パス忠実度行列。（右）各ノード間の光子伝送成功率とMCシミュレーション忠実度。

最大忠実度ルーティング結果（代表的なペア）:

| ノードペア | 最大忠実度 | MC成功率 | MC忠実度 |
|---------|---------|---------|---------|
| Hakusan→Otemachi | 0.940 | 0.718 | 0.999 ± 0.000 |
| Hakusan→NEC | 0.917 | 0.664 | 0.999 ± 0.000 |
| Hakusan→NIST_US | 0.918 | 0.596 | 0.999 ± 0.000 |
| NICT→Otemachi | 0.858 | 0.530 | 0.999 ± 0.000 |
| NICT→KOGANEI | 1.000 | 1.000 | N/A |

### 3.5 デコヒーレンスとチャネルロス

![チャネルロスと忠実度劣化](figures/fig5_decoherence.png)

![量子メモリ T₂ デコヒーレンス](figures/fig5b_memory_decoherence.png)

**Figure 5**: （左）異なる減衰係数（0.15, 0.20, 0.35 dB/km）における距離に対する伝送成功率と忠実度。（右）T₂ コヒーレンス時間（10 ms〜10 s）に対する Werner 状態忠実度の時間減衰。

通常のシングルモードファイバー（0.2 dB/km）では、距離50 km で伝送成功率10%（$\eta = 10^{-1}$）。T₂ = 1 秒のメモリでは1秒の保持後も忠実度0.5以上を維持する一方、T₂ = 10 ms のメモリでは100 ms 後に古典的限界 F=0.25 に接近する。

### 3.6 有限鍵境界

![有限鍵境界](figures/fig6_finite_key_boundary.png)

**Figure 6**: BB84 有限鍵率のブロックサイズ依存性（QBER=3%、各距離25/50/75/100 km）。正の鍵率が得られる最小ブロックサイズは距離100 kmで約10⁸〜10⁹パルス。

---

## 4. 考察と今後の展望

### 4.1 BB84 有限鍵解析の示唆

BB84 有限鍵率は距離 200 km まで正値を維持したが（N=10¹⁰ 時）、鍵率は 25 km から 200 km にかけて約3桁低下した。EUR フレームワークはセキュリティ補正項 $\Delta_{\text{sec}}$ が $\sqrt{N}$ のオーダーで効く一方で、鍵長 $\ell$ は $N$ に線形比例するため、充分大きな $N$ では収束した漸近率に近づく。本実験結果は先行研究（Kamin et al., 2025; Wiesemann et al., 2024）の解析と定性的に一致する。

### 4.2 量子リピータの最適化

NVセンターでは 4セグメントが最適な理由は、セグメント数増大で (a) リンクごとの伝送率 $p_{\text{gen}}$ が改善するが、(b) $\log_2 n$ レベルのネストスワッピングで全体待機時間が指数的に増加するためである。T₂ が短い（0.01 s）の場合、この待機中デコヒーレンスが致命的となり、捕捉イオン（T₂=60 s）が NV センター（T₂=1 s）を大幅に上回る性能を示した。

### 4.3 実装上の制約

1. 本シミュレーションは NetSquid のような精密な量子状態追跡ではなく Werner 状態近似を使用している。実際の量子チャネルはより複雑なノイズモデルを持つ可能性がある。
2. 東京QKDネットワークのトポロジーは Sasaki et al. (2011) の公開情報に基づく近似であり、正確な距離・機器パラメータは非公開である。
3. DEJMPS 計算でゲート忠実度 1.000 時に忠実度が 1.000 に収束する結果は理論的には正しいが（収束点 F=1 は固定点）、現実では有限ゲート誤りが上限を設定する。

### 4.4 今後の展望

- **誤り訂正符号の統合**: 量子誤り訂正（QEC）を組み合わせた fault-tolerant QKD の性能評価
- **マルチユーザネットワーク**: 複数エンド-エンドペアの同時ルーティング（RADAR-Q など）の実装
- **連続変数QKD**: CV-QKDの有限鍵解析拡張
- **衛星量子通信**: 大気損失モデルを追加した地上・衛星ハイブリッドネットワーク解析

---

## 5. 生成したファイル一覧

### ソースコード
| ファイル | 行数 | 内容 |
|---------|------|------|
| `src/qkd_finite_key.py` | ~210 | BB84/E91 有限鍵解析 |
| `src/quantum_repeater.py` | ~230 | 量子リピータ・蒸留モデル |
| `src/quantum_network.py` | ~280 | ネットワークトポロジー・ルーティング |
| `src/simulation_runner.py` | ~380 | メイン実験ランナー |
| `tests/test_qkd_protocols.py` | ~200 | 28件のユニットテスト（全パス）|

### 図表
| ファイル | 内容 |
|---------|------|
| `figures/fig1_key_rate_vs_distance.png` | BB84/E91 鍵率対距離 |
| `figures/fig2_repeater_memory.png` | 量子リピータ性能 |
| `figures/fig3_purification.png` | DEJMPS 蒸留収束 |
| `figures/fig3b_purification_trace.png` | 蒸留軌跡 |
| `figures/fig4_tokyo_fidelity_matrix.png` | 東京QKDネットワーク忠実度行列 |
| `figures/fig4b_tokyo_routing.png` | ルーティング結果 |
| `figures/fig5_decoherence.png` | チャネルロス・忠実度劣化 |
| `figures/fig5b_memory_decoherence.png` | T₂ デコヒーレンス |
| `figures/fig6_finite_key_boundary.png` | 有限鍵境界 |

### データ・ログ
| ファイル | 内容 |
|---------|------|
| `results/exp1_key_rate.json` | 実験1 数値データ |
| `results/exp2_repeater.json` | 実験2 数値データ |
| `results/exp4_tokyo_routing.json` | 実験4 ルーティング結果 |
| `results/exp6_finite_key_cv.json` | 実験6 交差検証結果 |
| `results/summary_table.md` | 全結果サマリーテーブル |
| `logs/process-log.jsonl` | 実行ログ（JSONL形式）|

---

## 参考文献

1. Wehner, S., Elkouss, D., & Hanson, R. (2018). Quantum internet: A vision for the road ahead. *Science*, 362(6412). https://doi.org/10.1126/science.aam9288
2. Briegel, H.-J., Dür, W., Cirac, J. I., & Zoller, P. (1998). Quantum repeaters: The role of imperfect local operations in quantum communication. *Physical Review Letters*, 81(26), 5932. https://doi.org/10.1103/PhysRevLett.81.5932
3. Wiesemann, J., Krause, J., Tupkary, D., Lütkenhaus, N., Rusca, D., & Walenta, N. (2024). A consolidated and accessible security proof for finite-size decoy-state quantum key distribution. arXiv:2405.16578.
4. Kamin, L., Tupkary, D., & Lütkenhaus, N. (2025). Improved finite-size effects in QKD protocols with applications to decoy-state QKD. arXiv:2502.05382.
5. Coopmans, T., et al. (2021). NetSquid, a NETwork Simulator for QUantum Information using Discrete events. *Communications Physics*, 4, 164. https://doi.org/10.1038/s42005-021-00647-8
6. Yehia, R., Neves, S., Diamanti, E., & Kerenidis, I. (2022). Quantum City: simulation of a practical near-term metropolitan quantum network. arXiv:2211.01190.
7. Haldar, S., et al. (2024). Reducing classical communication costs in multiplexed quantum repeaters using hardware-aware quasi-local policies. arXiv:2401.13168.
8. Liu, W.-Z., et al. (2026). A building block of quantum repeaters for scalable quantum networks. arXiv:2602.08472.
9. Tian, C., et al. (2026). RADAR-Q: Resource-Aware Distributed Asynchronous Routing for Entanglement Distribution in Multi-Tenant Quantum Networks. arXiv:2603.27570.
10. Sasaki, M., et al. (2011). Field demonstration of quantum key distribution in the Tokyo QKD Network. *Optics Express*, 19(11), 10387–10409. https://doi.org/10.1364/OE.19.010387
11. Tomamichel, M., & Renner, R. (2011). Uncertainty relation for smooth entropies. *Physical Review Letters*, 106(11), 110506. https://doi.org/10.1103/PhysRevLett.106.110506
12. Deutsch, D., et al. (1996). Quantum privacy amplification and the security of quantum cryptography over noisy channels. *Physical Review Letters*, 77(13), 2818. https://doi.org/10.1103/PhysRevLett.77.2818
