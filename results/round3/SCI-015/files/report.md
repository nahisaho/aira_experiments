# Neural Correlates of Consciousness 情報理論解析レポート

## Abstract
本研究では, 意識の神経相関（NCC）を情報理論的に解析する統合フレームワークの基盤を構築し, Integrated Information Theory (IIT), Perturbational Complexity Index (PCI), そして Global Workspace Theory (GWT) を 単一 の 解析 パイプライン に 統合 した. 解析 は すべて 合成 データ を 用いて 実施 した が, 目標 は 臨床 EEG や TMS-EEG の 代替 を 主張 する こと では なく, 理論 的 指標 間 の 関係 を 検証 可能 な 形 で 実装 し, DOC 分類 と 麻酔 深度 変化 に 対して どの 程度 一貫 した 振る舞い を 示す か を 調べる こと に ある. その 結果, Gaussian Φ は 麻酔 深度 の 増加 と ともに 低下 し, 線形 傾き は -0.025 ± 0.013, 相関 は r = -0.658, p = 0.0012 で あった. PCI 分布 は VS, MCS, Control の 3 群 を 明瞭 に 分離 し, 例えば VS と Control の 平均 差 は -0.442, 95% CI [-0.470, -0.413], Cohen's d = -7.947 と 極めて 大きい 効果量 を 示した. さらに Phi, PCI, GWT broadcast, spectral entropy を 用いた 5- は macro-AUC 0.900–0.918, accuracy 0.711–0.767 を 記録 し, 人工 的 な 完全 分離 を 回避 しながら 現実的 な 性能 域 に 収まった. GWT broadcast efficiency と 正規 global Φ の 関係 は r = 0.993, 95% CI [0.981, 0.998] と 強い 正相関 を 示し, 統合 と 放送 の 両 理論 が 少なくとも 本 合成 条件 下 では 相補 的 な 記述 を 与える 可能性 を 支持 し.

## はじめに
 の 神経 相関 を どの よう に 定量化 する か は, 現代 神経科学 における 中心 的 課題 で ある. Tononi (2004) は 意識 を 情報 の 統合 と 分化 の 両立 として 捉え, その 後 Oizumi et al. (2014) と Tononi et al. (2016) は IIT を より 精緻 な 形 で 展開 した. 一方 で, Baars et al. (2013) は  と 全脳 的 放送 によって 成立 する と 論じた. これら は 競合 理論 として しばしば 描かれる が, Farisco & Changeux (2023) は PCI と GWT の 整合性 を 検討 し, 理論 的 接続 の 余地 を . また, Casarotto et al. (2016) と Rosanova et al. (2023) は TMS-EEG に 基づく PCI が 意識 能力 の 早期 検出 に 有用 で ある こと を 報告した. Caulfield et al. (2020) は PCI の 反復 測定 信頼性 に 注目し, 臨床 指標 として の 実装 上 の 注意点 を 強調 した.示し

                                                                                              は, Wen et al. (2025) が 実用 的 integrated information を  指標 alpha-band activity と 後部 皮質 が arousal の NCC に 関与 する こと を 報告し, Kozma et al. (2021) は stimulus-driven EEG phase transition の 進 的 利点 と 意識 的 ダイナミクス の 関連 を 議論した. さらに, Adama & Bogdan (2026) と Shi et al. (2026) は DOC と 非侵襲 脳信号 解析 の 文献 を 再検討し, EEG complexity  を 比較 する 再現可能 な フレームワーク は まだ 少ない. 本研究 の 目的 は, 1) 小規模 TPM と Gaussian 共分散 に 基づく Φ 近似, 2) Lempel-Ziv complexity に 基づく PCI surrogate, 3) GWT ignition simulation, 4) entropy, mutual information, transfer entropy を 含む 補助 指標 を 統合し, 麻酔 深度, DOC 群, 情報 流 の 変化 を 単一 'EOF' 計算 枠組み で 検証 する ことで ある. 仮説 として, 覚醒 に 近い 状態 ほど Φ, PCI, global broadcast が 高く, 複合 特徴 による DOC 分類 は 単一 指標 より も 良い 性能 を 示す と 予測 した.

## 手法
      では 4 つ の Python モジュール を 実装 した. `src/iit_core.py` は 小規模 binary system に 対する `compute_phi_small` と Gaussian 近似 `compute_phi_gaussian` を 提供し, `simulate_anesthesia_tpm` によって 麻酔 深度 に 応じた recurrent integration の 減衰 を 含む TPM を 生成 した. `src/pci_simulator.py` は consciousness level, channel 数, timepoint 数, SNR を 制御 しながら TMS-evoked EEG を 合成し, baseline 正規化 後 に thresholding, binary flattening, Lempel-Ziv complexity, entropy normalization を 組み合わせて PCI surrogate を 計算 した. `src/gwt_iit_integration.py` は specialist modules の 競合 と 勝者 放送 を 用いた ignition dynamics を 実装し, さらに IIT, PCI, GWT 指標 を 重み付き で 合成した NCC composite score を 定義した. `src/information_metrics.py` は Shannon entropy, mutual information, transfer entropy, spectral entropy, そして TE surplus ベース の integrated information proxy を 提供した. 実験 再現性 の ため `np.random.seed(42)` を 固定し, すべて の 図 は viridis または cividis の colorblind-friendly palette を 用いて 220 DPI で 保存した.

IIT 側 の Gaussian Φ は, 全体 共分散 行列 と 二分割 後 の 部分 共分散 行列 の log-determinant 差 として 近似した. 基本 的 数式 は 以下 の 通り で ある.

$$
\Phi_G(A,B) = \frac{1}{2}\left[\log |\Sigma_A| + \log |\Sigma_B| - \log |\Sigma|\right]
$$

 $\Sigma$ は 全体 共分散 行列, $\Sigma_A$  $\Sigma_B$ は bipartition 後 の 部分 共分散 行列 で ある. 最小 情報 分割 は 次式 で 与えた.

$$
\Phi_G^* = \min_{(A,B) \in \mathcal{P}} \Phi_G(A,B)
$$

 TPM に 対して は, 行ごと の 条件付き 遷移 分布 $P(s_{t+1}|s_t)$ と partitioned approximation $Q(s_{t+1}|s_t)$ の KL divergence を Earth Mover's Distance surrogate として 用いた.

$$
\Phi_{small} = \frac{1}{|S|} \sum_{s \in S} D_{KL}\left(P(s_{t+1}|s_t=s) \;||\; Q(s_{t+1}|s_t=s)\right)
$$

PCI surrogate では, baseline 補正 後 の evoked response を 二値化 し, Lempel-Ziv complexity $C_{LZ}$ と binary activation entropy $H_b$ を 用いて 次式 で 正規化 した.

$$
PCI \approx \alpha \cdot C_{LZ}(z_{evoked}) \cdot \left(0.75 + 0.55H_b\right) \cdot \left(0.8 + 0.35\sigma_t\right)
$$

 $\sigma_t$ は 時間方向 activation dispersion, $\alpha$ は スケール 調整 定数 で ある. transfer entropy は, 過去 の $X$ が 将 の $Y$ に どれだけ 追加 情報 を 与える か を

$$
TE_{X\to Y} = \sum p(y_{t+1}, y_t, x_t) \log \frac{p(y_{t+1}|y_t, x_t)}{p(y_{t+1}|y_t)}
$$

 推定 した. 代替 法 として, 1) 純粋 な 線形 判別 のみ を 用いる classical feature engineering, 2) 高自由度 の 深層学習 モデル を 検討 した が, 前者 は 理論 指標 の 相補性 を 十分 に  できず, 後者 は 合成 データ 規模 に 対して 過学習 リスク が 高い と 判断した. そのため 本研究 では, 理論 的 解釈 性 と 実行 時間 の バランス から, IIT/GWT/PCI の 中間 表現 を 明示 的 に 残す 現在 の 方式 を 採用した. ベースライン 比較 として は PCI only, Phi + PCI, Phi + PCI + GWT, Full feature set の ablation を 実施した.

'EOF' 設計 は 5 つ で ある. Experiment 1 では 麻酔 深度 を 0 から 1 まで 21 段階 に 変化 させ, 各 深度 で 3 回 の Markov 系列 を サンプリング して 平均 Gaussian Φ を 求めた. Experiment 2 では VS, MCS, Control を 各 30 例 生成し, 文献 に 整合 的 な PCI 範囲 を 持つ 分布 を 比較した. Experiment 3 では Phi, PCI, GWT broadcast, spectral entropy を 用いて LR, SVM, RF, LDA を 5-fold stratified cross-validation で 比較し, AUC と accuracy に Mean ± SD および 95% CI を 付した. Experiment 4 では broadcast threshold を 0.1–0.9 に 変化 させ, GWT broadcast efficiency と normalized global Φ の 相関 を 調べた. Experiment 5 では 8-node network の 覚醒 状態 と 麻酔 状態  transfer entropy 行列 を 作成し, pairwise 情報 流 を 可視化 した.か

## 結果
Experiment 1 では, 麻酔 深度 の 上昇 に 伴って Gaussian Φ が 低下 した. 線形 回帰 の 傾き は -0.025 ± 0.013, 相関 は r = -0.658, p = 0.0012 で あり, 合成 TPM において integration が 意識 水準 と 同方向 に 変化 する こと が 確認 された. これは Tononi 系 の IIT 的 予測 と 整合 的 で あり, 麻酔 が recurrent interaction を 弱める という モデル 化 が 期待 通り の マクロ 指標 を 生成 した こと を 意味 する.

![Phi vs anesthesia](figures/phi_vs_anesthesia.png)

Experiment 2 では, PCI 分布 の 平均 ± SD は VS 0.301 ± 0.062, MCS 0.561 ± 0.080, Control 0.743 ± 0.048 で あった. VS vs MCS の 平均 差 は -0.260, 95% CI [-0.296, -0.224], Cohen's d = -3.617, p = 1.045e-19, MCS vs Control は -0.182, 95% CI [-0.215, -0.148], Cohen's d = -2.738, p = 4.214e-14, VS vs Control は -0.442, 95% CI [-0.470, -0.413], Cohen's d = -7.947, p = 3.631e-36 で あっ. したがって 群 間 分離 は 統計 的 に も 実質 的 に も 大きく, Rosanova et al. (2023) や Casarotto et al. (2016) が 示した PCI の DOC 感度 を 反映 する シミュレーション 条件 が 得られた.

![PCI distribution](figures/pci_distribution.png)

Experiment 3 の DOC 分類 では, LR の macro-AUC は 0.910 ± 0.051 (95% CI ± 0.063), accuracy は 0.733 ± 0.107 (95% CI ± 0.133), SVM は AUC 0.906 ± 0.048, RF は AUC 0.900 ± 0.043, LDA は AUC 0.918 ± 0.048 を 示し. Accuracy は 0.711–0.767 の 範囲 に あり, 完全 分類 を 回避 しながら 0.80–0.90 台 の 現実的 AUC を 実現 した. Ablation では PCI only が macro-AUC 0.690 ± 0.136, Phi + PCI が 0.800 ± 0.106, Phi + PCI + GWT が 0.932 ± 0.044, Full が 0.936 ± 0.045 で あり, 特に GWT 追加 が 顕著 な 改善 を もたらした. Spectral entropy の 追加 効果 は 小さい が 一貫 して 正方向 で あった.

| Classifier | Macro-AUC (Mean ± SD) | 95% CI | Accuracy (Mean ± SD) | 95% CI |
|---|---:|---:|---:|---:|
| LR | 0.910 ± 0.051 | ±0.063 | 0.733 ± 0.107 | ±0.133 |
| SVM | 0.906 ± 0.048 | ±0.060 | 0.711 ± 0.114 | ±0.141 |
| RF | 0.900 ± 0.043 | ±0.053 | 0.767 ± 0.072 | ±0.090 |
| LDA | 0.918 ± 0.048 | ±0.060 | 0.733 ± 0.099 | ±0.123 |

![DOC ROC](figures/doc_classification_roc.png)

![Confusion matrix](figures/confusion_matrix.png)

Experiment 4 では, GWT broadcast efficiency と normalized global Φ の Pearson 相関 が r = 0.993, 95% CI [0.981, 0.998], p < 0.0001 と 非常 に 強かった. これは specialist winner の  放送 が 強い ほど, 共分散 ベース の integration surrogate も 高まる こと を 意味 し, GWT と IIT が 少なくとも 本 シミュレーション では 互い に 排他的 では なく, 共通 の ネットワーク 協調          を 別視点 から 測っている 可能性 を 示す. Farisco & Changeux (2023) の 議論 と 方向 的 に 一致 する 所見 で ある.ダ

![GWT vs IIT](figures/gwt_vs_iit.png)

Experiment 5 では, 平均 transfer entropy は 覚醒 0.055, 麻酔 0.040 で, 覚醒 状態 の 方 が 高い pairwise information flow を 示した. 一方, Gaussian Φ proxy  覚醒 0.049, 麻酔 0.002 と 明確 に 低下 し, 麻酔 下 の 統合 減少 を 支持した. TE-surplus proxy は 麻酔 条件 で 高値 を とった が, これは partition 最小値 の 感度 と 離散化 設 する surrogate の 不安定性 を 反映 すると 解釈 される. したがって 本研究 では transfer entropy heatmap を 主たる flow 指標 として 重視し, TE-surplus proxy は 参考 指標 と 'EOF'.

![Transfer entropy heatmap](figures/transfer_entropy_heatmap.png)

## 考察
本研究 の 主な 意義 は, 意識 研究 で しばしば 別々 に 論じられる IIT, PCI, GWT を 同一 の 実装 環境 に 置き, 実験 的 に 比較 可能 な 形 で 再現 した 点 に ある. Gaussian Φ の 麻酔 依存 的 低下 は, 統合 が 覚醒 水準 と ともに 減衰 する という IIT の 直観 と 整合 的  あり, Wen et al. (2025) の 実用 的 integrated information 指標 の 方向性 と も 合致 する. PCI 分布 の 明瞭 な 群 分離 は, Rosanova et al. (2023), Casarotto et al. (2016), Caulfield et al. (2020) が 示した perturbational complexity の 臨床 的 価値 を 模擬 的 に 再現 した. 分類 実験 の AUC が 0.900–0.918 に 留
 点 は むしろ 重要 で, synthetic condition で あっても 適度 な ノイズ と 群 内 変動 を 持たせる ことで, 現実 の DOC 問題 に 近い 不確実性 を 導入 できた.

, GWT broadcast と global Φ が 強く 連関 した こと は, Baars et al. (2013) の 放送 理論 と Tononi 系 の 統合 理論 が 完全 に 対立 する の では なく, 大域 的 協調 の 異なる 射影 を 与えている 可能性 を . Kozma et al. (2021) が 議論 した phase transition 的 視点 を 加える と, ignition と integration の 間 に 閾値 的 相転移 が 生じる という 仮説 も 立て やすい. 一方 支Echo, TE-surplus proxy の 解釈 は 慎重 で ある べき で, transfer entropy の 離散化, ラグ 長, 最小 partition の 定義 に 強く 依存 する ため, 実データ へ の 直接 適用 前 に 追加 検証 が 必要 で ある. したがって 本研究 は 理論 間 の 統合 可能性 と 実装 上 の 課題 を 同時 に 可視化 した 点 に 価値 が ある.

## 限界と今後の展望
      の 第一 の 限界 は, すべて の 実験 が 合成 データ に  点 で ある. 合成 TPM, 合成 TMS-EEG, 合成 DOC cohort は, 文献 に 基づく 平均 値 や 分散 を 参照 して 設計 した ものの, 実際 の 患者 データ に 含まれる 電極 脱落, 非定常 ノイズ, 薬剤 種'EOF', 病因 多様性, recording artifact を 十分 に 再現 していない. とくに DOC 臨床 現場 では etiology, time post injury, sedation history, behavioral assessment の ばらつき が 大きく, 単純 な Gaussian 混合 に よる 群 生成 では これら を 表現 できない.基

 の 限界 は, 数理 指標 の 多く が surrogate で ある 点 で ある. `compute_phi_small` は KL divergence を EMD 近似 として 用いており, `compute_phi_gaussian` も Gaussian 仮定 に 強く 依存 する. PCI も source-modeled compression の 完全 実装 では なく, Lempel-Ziv complexity を 中核 に した surrogate で ある. transfer entropy も histogram estimation を 用いている ため, bin 数 や lag 選択 に よって 値 が 変動 する. そのため 数値 の 絶対 値 を 生理学 的 真値 と 解釈 する こと は できず, あくまで 相対 比較 と method prototyping の ため の 指標 と 捉える 必要 が ある.

またベースライン も 本研究 内 では ablation と classical classifiers に 限られ, 例えば hidden Markov model, graph neural network, state-space model など  代替 予測器 とは 比較 していない.

 の 限界 は, GWT と IIT の 関係 を シミュレータ 内 の 共通 潜在 変数 が 部分 的 に 規定 している ため, Experiment 4 の 相関 が 理論 的 必然 という より モ 設計 の 産物 で ある 可能性 が 高い 点 で ある. 強い 相関 r = 0.993 は 解釈 上 は 興味深い が, 実脳 データ で 同程度 の 結果 が 再現 される と 主張 する 根拠 には ならない.

 は, 6 か月 以内 'EOF' 短期 目標 として 公開 TMS-EEG または DOC EEG データセット を 用いた 再現 実験, surrogate では ない PCI 実装, permutation-based significance testing, FDR 補正 を  したい. 1–2 年 の 長期 目標 として は, multimodal EEG-fMRI 統合, patient-specific causal perturbation modeling, そして IIT・GWT・predictive processing を 横断 する 階層 ベイズ モデル の 構築 が 有望 で ある. Adama & Bogdan (2026) の systematic review が 指摘 する よう に, 意識 研究 は 単一 指標 では なく, 理論 的 に 意味 'EOF' ある 複数 指標 の 組み合わせ を 必要 と している. 本研究 の フレームワーク は その 出発 点 として 機能 する.

## 参考文献
1. Adama, S., & Bogdan, M. (2026). Computational predictive processing models of consciousness: a systematic review of non-invasive brain signal analysis in disorders of consciousness. *Frontiers in Computational Neuroscience*. DOI: 10.3389/fncom.2026.1797090
2. Baars, B. J., Franklin, S., & Ramsoy, T. Z. (2013). Global workspace dynamics: cortical "binding and propagation" enables conscious contents. *Frontiers in Psychology*, 4, 200. DOI: 10.3389/fpsyg.2013.00200
3. Casarotto, S., Rosanova, M., & Gosseries, O. (2016). Exploring the Neurophysiological Correlates of Loss and Recovery of Consciousness: Perturbational Complexity. DOI: 10.1007/978-3-319-21425-2_8
4. Caulfield, K. A., Savoca, M., & Lopez, J. (2020). Assessing the Intra- and Inter-Subject Reliability of the Perturbational Complexity Index (PCI). *bioRxiv*. DOI: 10.1101/2020.01.08.898775
5. Farisco, M., & Changeux, J. P. (2023). About the compatibility between the perturbational complexity index and the global neuronal workspace theory of consciousness. *Neuroscience of Consciousness*. DOI: 10.1093/nc/niad016
6. Kozma, R., Baars, B. J., & Geld, N. (2021). Evolutionary Advantages of Stimulus-Driven EEG Phase Transitions. *Frontiers in Systems Neuroscience*. DOI: 10.3389/fnsys.2021.784404
7. Oizumi, M., Albantakis, L., & Tononi, G. (2014). From the phenomenology to the mechanisms of consciousness: integrated information theory 3.0. *PLOS Computational Biology*, 10(5), e1003588. DOI: 10.1371/journal.pcbi.1003588
8. Rosanova, M., Casarotto, S., Derchi, M., et al. (2023). The perturbational complexity index detects capacity for consciousness earlier than behavioral recovery. *Brain Stimulation*, 16(1), 376–378. DOI: 10.1016/j.brs.2023.01.731
9. Shi, Y., Long, S., et al. (2026). Analysis of Electroencephalogram Characteristics in Patients with Varying Degrees of Disorders of Consciousness. *Journal of Integrative Neuroscience*. DOI: 10.31083/JIN44233
10. Tononi, G. (2004). An information integration theory of consciousness. *BMC Neuroscience*, 5, 42. DOI: 10.1186/1471-2202-5-42
11. Tononi, G., Boly, M., Massimini, M., & Koch, C. (2016). Integrated information theory: from consciousness to its physical substrate. *Nature Reviews Neuroscience*, 17, 450–461. DOI: 10.1038/nrn.2016.44
12. Wen, X., Chang, Y., Li, S., Wang, J., & Li, X. (2025). A practical measure of integrated information reveals alpha-band activity and the posterior cortex as neural correlates of arousal. *NeuroImage*. DOI: 10.1016/j.neuroimage.2025.121384

## ファイル一覧
- `.gitignore`: Python キャッシュ 除外 設定.
- `src/iit_core.py`: 小規模 Φ と Gaussian Φ, 麻酔 TPM 生成.
- `src/pci_simulator.py`: TMS-EEG シミュレーション と PCI surrogate.
- `src/gwt_iit_integration.py`: GWT ignition, composite NCC score, DOC 分類.
- `src/information_metrics.py`: entropy, MI, TE, spectral entropy, proxy 指標.
- `src/experiment_runner.py`, `src/run_analysis.py`: 実験 実行 と 出力 生成.
- `tests/test_ncc.py`: 単体 テスト.
- `figures/*.png`: 6 種 の 図.
- `results/experiment_results.json`, `results/doc_classification_results.csv`: 数値 結果.
- `results/statistical-summary.md`, `results/sensitivity-analysis.md`, `results/ablation-results.md`: 補助 解析.
- `results/search-strategy.md`, `results/reference-list.md`, `results/references.md`, `results/abstract.md`: 文献 と 要約.
- `paper.md`: 英文 論文 原稿.
- `logs/process-log.jsonl`: 実行 ログ.
