# 実験レポート: 意識の神経相関の情報理論的解析フレームワーク

## 1. 実験目的と背景
本実験の目的は、意識の神経相関（NCC）を情報理論的に解析する統合フレームワークを構築し、IITの統合情報量 \(\Phi\)、PCI（Perturbational Complexity Index）、および Global Workspace Theory (GWT) の ignition 指標を同一パイプラインで扱えるようにすることである。対象は覚醒、NREM睡眠、プロポフォール麻酔、ケタミン、植物状態（VS）、最小意識状態（MCS）の6状態であり、EEGシミュレーション、特徴抽出、分類、統計解析、図生成までを一括実装した。

## 2. 先行研究調査結果
- **Sinitsyn et al. (2020)**: PCI は MCS と UWS/VS の識別に高い感度を示した（DOI: 10.3390/brainsci10120917）。
- **Wang et al. (2023)**: PCIst は CRS-R と相関し、自然周波数変化も DOC 評価に有効だった（DOI: 10.1016/j.compbiomed.2023.107547）。
- **Farnes et al. (2020)**: ケタミンでは PCI が大きく変化しない一方、自発脳活動の Lempel-Ziv complexity は増加した（DOI: 10.1371/journal.pone.0242056）。
- **Edlow et al. (2020)**: DOC では covert consciousness / CMD の見逃しが臨床的課題である（DOI: 10.1038/s41582-020-00428-x）。
- **Comanducci et al. (2020)**: EEG と高度神経生理指標の併用が重要（DOI: 10.1016/j.clinph.2020.07.015）。
- **Frohlich et al. (2021)**: 低周波優位でも複雑性の観点から意識を再評価すべき場合がある（DOI: 10.1093/brain/awab095）。
- **Safron (2020)**: IIT と GWT を統合する IWMT を提案（DOI: 10.3389/frai.2020.00030）。
- **Storm et al. (2024)**: 意識理論の統合的・マルチスケール理解の必要性を示した（DOI: 10.1016/j.neuron.2024.02.004）。
- **Butlin et al. (2023)**: AI における意識評価でも情報統合や全体放送の観点が重要（DOI: 10.48550/arxiv.2308.08708）。
- **Caulfield et al. (2020)**: PCI の部位間再現性が検討された（DOI: 10.1101/2020.01.08.898775）。

## 3. NatureLM MCPによる科学的知見
NatureLM により以下の知見を整理した。
- IIT \(\Phi\): 意識ありで概ね **0.2–0.6**、意識低下でより低値。
- PCI: 意識ありで高く、無意識状態では低下。
- EEG複雑性（LZc）: 麻酔下で **0.2–0.4**、覚醒で **0.8 前後**。
- Spectral entropy: 無意識で低く、意識ありで高い。
- 覚醒では alpha/beta 優位、VS では delta 優位。

本フレームワークでは、これらの知見を完全な厳密再現ではなく、**シミュレーションの事前分布・キャリブレーション目標**として利用した。

## 4. 使用した手法・アルゴリズムの概要
### EEG生成
6チャネル EEG を状態別パラメータで生成した。覚醒では alpha/beta を強くし、NREM と VS では delta を強め、プロポフォールでは burst-suppression、ケタミンでは高周波ゆらぎ、MCS では断続的 alpha を導入した。

### IIT \(\Phi^*\) 近似
6ノードの二値化状態列を作り、
\[
EI = I(X_t; X_{t+1})
\]
を全体の有効情報量とした。さらに二分割 \(A|B\) ごとに
\[
EI_{A|B}=I(A_t;A_{t+1})+I(B_t;B_{t+1})
\]
を計算し、MIP を探索して \(\Phi^*\) を近似した。

### PCIシミュレーション
TMS様パルスを小規模再帰ネットワークに入力し、時空間応答を二値化して Lempel-Ziv complexity を計算し、正規化した PCI-like 指標を得た。

### GWT ignition
刺激後に多ノードへ分散しつつ一定時間持続する高活動状態を ignition と定義し、試行間の発生確率を GWT ignition probability とした。

### 分類器
5クラス（Awake, NREM, Propofol, VS, MCS）に対し、Random Forest + SVM のアンサンブルを 5-fold stratified cross-validation で評価した。

## 5. 主要な結果と数値
### 状態別 \(\Phi\) と PCI
- Awake: \(\Phi = 0.314 \pm 0.044\), PCI = \(0.915 \pm 0.051\)
- NREM: \(\Phi = 0.164 \pm 0.030\), PCI = \(0.605 \pm 0.045\)
- Propofol: \(\Phi = 0.106 \pm 0.033\), PCI = \(0.525 \pm 0.056\)
- Ketamine: \(\Phi = 0.305 \pm 0.054\), PCI = \(0.896 \pm 0.057\)
- VS: \(\Phi = 0.112 \pm 0.029\), PCI = \(0.484 \pm 0.045\)
- MCS: \(\Phi = 0.186 \pm 0.037\), PCI = \(0.650 \pm 0.035\)

![Figure 1](figures/fig1_phi_comparison.png)

![Figure 2](figures/fig2_pci_simulation.png)

### EEGシミュレーション例
![Figure 3](figures/fig3_eeg_signals.png)

### 特徴空間とクラス分離
\(\Phi\)–PCI 平面では Awake / Ketamine が高複雑性側、VS / Propofol が低複雑性側、MCS が中間に位置した。

![Figure 4](figures/fig4_feature_space.png)

### 分類性能
- Accuracy: **0.786 ± 0.039**
- Weighted F1: **0.780 ± 0.040**
- Macro AUC: **0.954**
- クラス別 AUC: Awake 1.000, NREM 0.973, Propofol 0.899, VS 0.918, MCS 0.978

主な誤分類は、**Propofol ↔ VS** の相互混同、MCS の一部が VS / NREM に寄るパターンであり、完全分離ではなく臨床的に近い状態の重なりを再現できた。

![Figure 5](figures/fig5_classifier_performance.png)

### GWT ignition
Awake で最も強い ignition、MCS で中等度、VS で弱い spread が観察された。

![Figure 6](figures/fig6_gwt_ignition.png)

### VS と MCS の統計比較
有意差が大きかった特徴は以下の通り。
- \(\Phi\): \(p = 5.55 \times 10^{-10}\), Cohen's d = 2.83
- PCI: \(p = 2.19 \times 10^{-8}\), d = 2.16
- Spectral entropy: \(p = 1.40 \times 10^{-10}\), d = 6.53
- Sample entropy: \(p = 1.40 \times 10^{-10}\), d = 6.07
- LZ complexity: \(p = 1.40 \times 10^{-10}\), d = 4.75
- Mean coherence: \(p = 1.56 \times 10^{-9}\), d = 2.40

Random Forest の feature importance 上位5項目は、**spectral entropy, sample entropy, LZ complexity, GWT ignition probability, mean coherence** だった。

## 6. 考察と今後の展望
本実装により、IIT・PCI・GWT を別々の理論ではなく、同一データ上の相補的な特徴群として扱えることが示された。特に、VS と MCS を単一指標ではなく多特徴量で分ける設計は、DOC 臨床に近い発想である。また、ケタミンを Awake に近い高複雑性側へ配置した点は、PCI と自発複雑性の乖離を議論する上でも有用である。

今後は以下を進める余地がある。
1. 実EEG/TMS-EEG データへの適用
2. PyPhi 互換のより厳密な \(\Phi\) 計算
3. グラフ理論指標や source-space connectivity の追加
4. 患者縦断データへの回復予測モデル拡張

## 生成したファイル一覧
- `ncc_framework.py` — NCC 解析フレームワーク本体
- `ncc_results.json` — 数値結果の保存
- `figures/fig1_phi_comparison.png`
- `figures/fig2_pci_simulation.png`
- `figures/fig3_eeg_signals.png`
- `figures/fig4_feature_space.png`
- `figures/fig5_classifier_performance.png`
- `figures/fig6_gwt_ignition.png`
- `paper.md`
- `report.md`
