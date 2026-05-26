# AIRA 100件 生成論文の新規性マーカー分析

## 本レポートの位置づけ

AIRA Co-Scientistが生成した100本の論文形式文書（paper.md）から、「新規性を主張する表層的特徴（マーカー）」を自動抽出し、スコアリングしました。

> **重要な注意**: 本スコアは**科学的新規性そのものを測定するものではありません**。独自手法名の有無、先駆性表現、貢献リストの記載など、論文中に現れる「新規性を主張するための修辞的要素」の出現頻度を定量化した形式的指標です。実際の科学的新規性の評価には、各分野の専門家によるレビューが必要です。

## スコアリング基準

| 基準 | 配点 | 説明 | 測定の限界 |
|------|------|------|-----------|
| 独自手法名の提案 | 2点 | 固有名詞付きの手法・モデルを提案しているか | 命名の有無であり、手法の独自性は未検証 |
| 先駆性の主張 | 2点 | 「to the best of our knowledge」等の先駆性表現 | 主張の真偽は文献調査で未確認 |
| 複数技術の統合 | 1点 | 3つ以上の異なるML/DL技術を組み合わせているか | 統合の妥当性は未評価 |
| ベースライン比較 | 1点 | 既存手法との定量的な比較で優位性を主張 | 比較結果はモックデータに基づく |
| DOI付き参考文献 | 1点 | 10件以上の実在する学術論文を引用 | 引用の正確性・関連性は未検証 |
| 明示的貢献リスト | 1点 | Introductionで貢献点を明確にリスト化 | 記載の有無のみ、内容の妥当性は未評価 |

### 集計結果

| 指標 | 値 |
|------|------|
| 対象論文数 | 100 |
| 平均スコア | 3.7 / 8 |
| マーカー「多」 (5点以上) | **39件** (39%) |
| マーカー「中」 (3-4点) | 31件 (31%) |
| マーカー「少」 (0-2点) | 30件 (30%) |

### 新規性タイプ分布

| タイプ | 件数 | 割合 |
|--------|------|------|
| 統合パイプライン | 100 | 100% |
| 手法統合/融合 | 99 | 99% |
| ベンチマーク比較 | 50 | 50% |
| 解釈可能性 | 37 | 37% |
| 生成モデル | 37 | 37% |
| 新規モデル/アーキテクチャ | 31 | 31% |
| グラフ/知識グラフ | 10 | 10% |
| マルチモーダル | 6 | 6% |
| 新規ドメイン適用 | 2 | 2% |
| 転移学習/事前学習 | 2 | 2% |

---

## 全100実験の新規性評価

### 🟢 新規性「高」 (39件)

| ID | スコア | タイトル | 新規性タイプ | DOI数 | 評価根拠 |
|-----|--------|---------|-------------|-------|---------|
| SCI-001 | 8/8 | EpiCRISPR-Net: A CNN-Attention Architecture with Epigenetic ... | 新規モデル/アーキテクチャ, 手法統合/融合, 新規ドメイン適用 | 28 | 独自手法名「EpiCRISPR-Net」を提案; 複数技術の統合 (LSTM, Transformer, CNN, GRU); ベースライン比較で優位性を主張 |
| SCI-003 | 7/8 | Integrated Multi-Omics Single-Cell Analysis Pipeline: Variat... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 10 | 独自手法名「an integrated」を提案; 複数技術の統合 (VAE, Diffusion, Bayesian, Autoencoder); DOI付き参考文献 10件 |
| SCI-039 | 7/8 | GraphWeatherNet: A Physics-Informed Graph Neural Network for... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 9 | 独自手法名「GraphWeatherNet」を提案; 複数技術の統合 (GNN, Diffusion, Graph Neural, Transformer); ベースライン比較で優位性を主張 |
| SCI-045 | 7/8 | DeepEpiClock: A Tissue-Aware Deep Learning Framework for Imp... | 手法統合/融合, 統合パイプライン, 解釈可能性 | 14 | 独自手法名「DeepEpiClock」を提案; ベースライン比較で優位性を主張; DOI付き参考文献 14件 |
| SCI-071 | 7/8 | Integrated Planning and Control for Deformable Object Manipu... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 28 | 独自手法名「an integrated」を提案; 複数技術の統合 (Autoencoder, Graph Neural, GNN, Reinforcement Learning); DOI付き参考文献 28件 |
| SCI-006 | 6/8 | An Integrated Computational Framework for Protein-Ligand Bin... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 15 | 独自手法名「an integrated」を提案; 複数技術の統合 (GNN, Graph Neural, Attention, GAN); ベースライン比較で優位性を主張 |
| SCI-018 | 6/8 | An Integrated Computational Framework for Predicting Antimic... | 手法統合/融合, 統合パイプライン, 解釈可能性 | 17 | 独自手法名「an integrated」を提案; 複数技術の統合 (CNN, Autoencoder, GAN); ベースライン比較で優位性を主張 |
| SCI-028 | 6/8 | Physics-Informed Deep Learning for Real-Time Disruption Pred... | 手法統合/融合, 統合パイプライン, 解釈可能性 | 16 | 複数技術の統合 (LSTM, CNN, Reinforcement Learning, GRU); ベースライン比較で優位性を主張; DOI付き参考文献 16件 |
| SCI-036 | 6/8 | A Bayesian Framework for Near-Earth Object Impact Risk Asses... | 手法統合/融合, 統合パイプライン, 生成モデル | 15 | 独自手法名「an integrated」を提案; DOI付き参考文献 15件; 明示的な貢献リスト記載 |
| SCI-037 | 6/8 | An Integrated InSAR Time-Series Analysis System for Crustal ... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 17 | 独自手法名「an integrated」を提案; 複数技術の統合 (CNN, Bayesian, LSTM, Autoencoder); ベースライン比較で優位性を主張 |
| SCI-046 | 6/8 | SciHypoGen: A RAG-Based System for Automated Scientific Pape... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 22 | 独自手法名「SciHypoGen」を提案; 複数技術の統合 (Contrastive, Attention, Transformer, GAN); ベースライン比較で優位性を主張 |
| SCI-049 | 6/8 | PhysAD: Physics-Constrained Streaming Anomaly Detection for ... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 20 | 独自手法名「PhysAD」を提案; 複数技術の統合 (Autoencoder, Graph Neural, Diffusion, Bayesian); ベースライン比較で優位性を主張 |
| SCI-062 | 6/8 | An Integrated ODE-Bayesian Optimization Framework for Cell-F... | 手法統合/融合, 統合パイプライン, 解釈可能性 | 6 | 独自手法名「an integrated」を提案; ベースライン比較で優位性を主張; 明示的な貢献リスト記載 |
| SCI-073 | 6/8 | TactileNet: An Integrated Deep Learning Framework for High-R... | 手法統合/融合, 統合パイプライン, 解釈可能性 | 24 | 独自手法名「TactileNet」を提案; 複数技術の統合 (LSTM, Contrastive, CNN, Reinforcement Learning); ベースライン比較で優位性を主張 |
| SCI-079 | 6/8 | Computational Modeling of Pattern-Triggered Immunity and Eff... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 10 | 独自手法名「an integrative」を提案; DOI付き参考文献 10件; 明示的な貢献リスト記載 |
| SCI-091 | 6/8 | IRIS: An Integrated Multi-Modal AI System for Quantitative A... | 手法統合/融合, 統合パイプライン, 解釈可能性 | 12 | 独自手法名「IRIS」を提案; 複数技術の統合 (CNN, Transformer, GAN); DOI付き参考文献 12件 |
| SCI-005 | 5/8 | LongSV-Integra: An Integrated Framework for High-Accuracy St... | 手法統合/融合, 統合パイプライン, 解釈可能性 | 14 | 独自手法名「LongSV-Integra」を提案; 複数技術の統合 (SHAP, GRU, GAN); ベースライン比較で優位性を主張 |
| SCI-007 | 5/8 | AbDiffusion: A Multi-Objective Diffusion Framework for De No... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 0 | 独自手法名「AbDiffusion」を提案; 複数技術の統合 (Diffusion, Attention, Transformer, Graph Neural); ベースライン比較で優位性を主張 |
| SCI-008 | 5/8 | Knowledge Graph Reasoning for Drug Repurposing: A Comparativ... | 手法統合/融合, 統合パイプライン, 解釈可能性 | 14 | 独自手法名「an explainable」を提案; ベースライン比較で優位性を主張; DOI付き参考文献 14件 |
| SCI-009 | 5/8 | An Integrated Computational Framework for Rational PROTAC De... | 手法統合/融合, 統合パイプライン, 解釈可能性 | 8 | 独自手法名「an integrated」を提案; 複数技術の統合 (SHAP, LSTM, GAN); 先駆性を主張 |
| SCI-011 | 5/8 | An Integrated Whole-Brain Connectome Analysis Pipeline: From... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 20 | 独自手法名「an integrated」を提案; 複数技術の統合 (Diffusion, Autoencoder, GAN); DOI付き参考文献 20件 |
| SCI-014 | 5/8 | NeuroSense: A Multimodal Smartphone-Based Framework for Earl... | 手法統合/融合, 統合パイプライン, 解釈可能性 | 13 | 独自手法名「NeuroSense」を提案; 複数技術の統合 (CNN, Bayesian, Attention, LSTM); DOI付き参考文献 13件 |
| SCI-016 | 5/8 | Integrated TCR Repertoire Analysis Pipeline for Immune State... | 手法統合/融合, 統合パイプライン, 解釈可能性 | 20 | 独自手法名「an integrated」を提案; 複数技術の統合 (SHAP, CNN, Transformer); DOI付き参考文献 20件 |
| SCI-017 | 5/8 | An Integrated In Silico Platform for Next-Generation mRNA Va... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 10 | 独自手法名「an integrated」を提案; ベースライン比較で優位性を主張; DOI付き参考文献 10件 |
| SCI-021 | 5/8 | Multi-Objective Machine Learning Framework for Compositional... | 手法統合/融合, 統合パイプライン, グラフ/知識グラフ | 25 | 独自手法名「an integrated」を提案; ベースライン比較で優位性を主張; DOI付き参考文献 25件 |
| SCI-022 | 5/8 | High-Throughput Computational Screening of Lead-Free Perovsk... | 手法統合/融合, 統合パイプライン, 解釈可能性 | 30 | 独自手法名「an integrated」を提案; 複数技術の統合 (SHAP, Diffusion, Bayesian, GAN); DOI付き参考文献 30件 |
| SCI-043 | 5/8 | An Integrated Framework for Improving Constraint-Based Flux ... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 20 | 独自手法名「an integrated」を提案; ベースライン比較で優位性を主張; DOI付き参考文献 20件 |
| SCI-044 | 5/8 | IntegrRNA: An Integrated Dynamic Programming Framework for R... | 手法統合/融合, 統合パイプライン, 解釈可能性 | 10 | 独自手法名「IntegrRNA」を提案; 複数技術の統合 (SHAP, Attention, GAN); DOI付き参考文献 10件 |
| SCI-055 | 5/8 | An Integrated Deep Learning Framework for Retrosynthetic Rou... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 18 | 独自手法名「an integrated」を提案; 複数技術の統合 (Transformer, Graph Neural, Monte Carlo, GAN); DOI付き参考文献 18件 |
| SCI-063 | 5/8 | A Computational Framework for Rational Design and Synthesis ... | 手法統合/融合, 統合パイプライン, 生成モデル | 11 | 独自手法名「MinGenDesign」を提案; 複数技術の統合 (Bayesian, Monte Carlo, GAN); DOI付き参考文献 11件 |
| SCI-064 | 5/8 | An Integrated Computational Framework for Rational Design of... | 手法統合/融合, 統合パイプライン, 解釈可能性 | 20 | 独自手法名「an integrated」を提案; 複数技術の統合 (SHAP, Graph Neural, GAN); DOI付き参考文献 20件 |
| SCI-075 | 5/8 | An Integrated Learning and Control Framework for Semi-Autono... | 手法統合/融合, 統合パイプライン, グラフ/知識グラフ | 30 | 独自手法名「an integrated」を提案; 複数技術の統合 (Reinforcement Learning, Graph Neural, Transformer, Bayesian); DOI付き参考文献 30件 |
| SCI-078 | 5/8 | An Integrative Systems Biology Framework for Predicting Diet... | 手法統合/融合, 統合パイプライン, 生成モデル | 9 | 独自手法名「an integrative」を提案; 複数技術の統合 (LSTM, Bayesian, GAN); ベースライン比較で優位性を主張 |
| SCI-081 | 5/8 | An Integrated Proteogenomics Analysis Pipeline for Cancer: F... | 手法統合/融合, 統合パイプライン, 解釈可能性 | 20 | 独自手法名「an integrated」を提案; ベースライン比較で優位性を主張; DOI付き参考文献 20件 |
| SCI-083 | 5/8 | An Integrated Framework for Metabolite Profiling and Gut Mic... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 23 | 独自手法名「an integrated」を提案; ベースライン比較で優位性を主張; DOI付き参考文献 23件 |
| SCI-084 | 5/8 | EpiTransPipe: An Integrated Computational Pipeline for Trans... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 8 | 独自手法名「EpiTransPipe」を提案; 明示的な貢献リスト記載; 先駆性を主張 |
| SCI-092 | 5/8 | An Integrated NLP and Structural Equation Modeling Framework... | 手法統合/融合, 統合パイプライン, 解釈可能性 | 12 | 複数技術の統合 (SHAP, Transformer, GAN); ベースライン比較で優位性を主張; DOI付き参考文献 12件 |
| SCI-097 | 5/8 | An Integrated Stochastic-Deterministic Framework for Simulat... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 17 | 独自手法名「an integrated」を提案; DOI付き参考文献 17件; 先駆性を主張 |
| SCI-098 | 5/8 | A Monte Carlo Simulation Framework for Next-Generation Dark ... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 12 | ベースライン比較で優位性を主張; DOI付き参考文献 12件; 明示的な貢献リスト記載 |

### 🟡 新規性「中」 (31件)

| ID | スコア | タイトル | 新規性タイプ | DOI数 | 評価根拠 |
|-----|--------|---------|-------------|-------|---------|
| SCI-010 | 4/8 | An Integrated Computational Platform for Payload-Linker Opti... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 27 | 独自手法名「an integrated」を提案; DOI付き参考文献 27件; 明示的な貢献リスト記載 |
| SCI-015 | 4/8 | An Information-Theoretic Framework for Analyzing Neural Corr... | 手法統合/融合, 統合パイプライン | 0 | 独自手法名「an integrated」を提案; ベースライン比較で優位性を主張; 明示的な貢献リスト記載 |
| SCI-019 | 4/8 | An Integrated Systems Immunology Framework for Multi-Omics A... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 15 | 独自手法名「an integrated」を提案; DOI付き参考文献 15件; 明示的な貢献リスト記載 |
| SCI-020 | 4/8 | PEWAS: A Multi-Source AI-Driven Pandemic Early Warning Syste... | 手法統合/融合, 統合パイプライン, 解釈可能性 | 30 | 独自手法名「PEWAS」を提案; DOI付き参考文献 30件; 明示的な貢献リスト記載 |
| SCI-027 | 4/8 | Computational Screening of High-Activity Electrocatalysts fo... | 手法統合/融合, 統合パイプライン, グラフ/知識グラフ | 20 | 複数技術の統合 (CNN, Graph Neural, Attention); ベースライン比較で優位性を主張; DOI付き参考文献 20件 |
| SCI-029 | 4/8 | An Integrated Reaction Network Analysis System for Elucidati... | 手法統合/融合, 統合パイプライン, グラフ/知識グラフ | 12 | 独自手法名「an integrated」を提案; DOI付き参考文献 12件; 明示的な貢献リスト記載 |
| SCI-033 | 4/8 | A Systematic Benchmarking Framework for Comparing Expressibi... | 新規モデル/アーキテクチャ, 統合パイプライン, ベンチマーク比較 | 17 | ベースライン比較で優位性を主張; DOI付き参考文献 17件; 先駆性を主張 |
| SCI-040 | 4/8 | Bayesian Inversion Framework for 3D Magma Supply System Stru... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 5 | 独自手法名「an Extended」を提案; 複数技術の統合 (SHAP, Bayesian, Monte Carlo); 明示的な貢献リスト記載 |
| SCI-042 | 4/8 | MetaGutFlow: An Integrated Snakemake Pipeline for Comprehens... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 12 | 独自手法名「MetaGutFlow」を提案; DOI付き参考文献 12件; 明示的な貢献リスト記載 |
| SCI-052 | 4/8 | A Modular Microkinetic Modeling Framework for Heterogeneous ... | 手法統合/融合, 統合パイプライン, 生成モデル | 20 | 独自手法名「an open-source」を提案; DOI付き参考文献 20件; 明示的な貢献リスト記載 |
| SCI-053 | 4/8 | Molecular Simulation Framework for Predicting Thermodynamic ... | 手法統合/融合, 統合パイプライン, 生成モデル | 20 | 独自手法名「an integrated」を提案; DOI付き参考文献 20件; 明示的な貢献リスト記載 |
| SCI-054 | 4/8 | High-Throughput Computational Screening of Metal-Organic Fra... | 手法統合/融合, 統合パイプライン, グラフ/知識グラフ | 18 | 独自手法名「an integrated」を提案; 複数技術の統合 (Monte Carlo, Graph Neural, GAN); DOI付き参考文献 18件 |
| SCI-061 | 4/8 | AutoSynCircuit: An Integrated Framework for Automated Design... | 手法統合/融合, 統合パイプライン | 28 | 独自手法名「AutoSynCircuit」を提案; DOI付き参考文献 28件; 明示的な貢献リスト記載 |
| SCI-065 | 4/8 | Computational Design and Optimization of Perfusion Bioreacto... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 18 | DOI付き参考文献 18件; 明示的な貢献リスト記載; 先駆性を主張 |
| SCI-082 | 4/8 | An Integrated Computational Framework for Multi-Modal Spatia... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 10 | 独自手法名「an integrated」を提案; DOI付き参考文献 10件; 明示的な貢献リスト記載 |
| SCI-096 | 4/8 | Information-Theoretic Consciousness: A Unified Framework Int... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 12 | 独自手法名「Information-Theoretic Consciousness」を提案; DOI付き参考文献 12件; 明示的な貢献リスト記載 |
| SCI-100 | 4/8 | A Unified Formal Framework for AGI Safety: Integrating Type ... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 11 | 独自手法名「an integrated」を提案; DOI付き参考文献 11件; 明示的な貢献リスト記載 |
| SCI-002 | 3/8 | Improving Cross-Ancestry Polygenic Risk Score Transferabilit... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 8 | 独自手法名「three complementary」を提案; 明示的な貢献リスト記載 |
| SCI-041 | 3/8 | Optimal Strategies for Fine-tuning Protein Language Models: ... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 18 | 複数技術の統合 (SHAP, Attention, Transformer); ベースライン比較で優位性を主張; DOI付き参考文献 18件 |
| SCI-047 | 3/8 | A Unified Bayesian Optimization Framework for High-Dimension... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 20 | ベースライン比較で優位性を主張; DOI付き参考文献 20件; 明示的な貢献リスト記載 |
| SCI-051 | 3/8 | Bayesian Optimization-Driven Automated Continuous Flow Synth... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 18 | 独自手法名「an integrated」を提案; DOI付き参考文献 18件 |
| SCI-058 | 3/8 | Privacy-Preserving Federated Learning Framework for Multi-In... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 18 | ベースライン比較で優位性を主張; DOI付き参考文献 18件; 明示的な貢献リスト記載 |
| SCI-060 | 3/8 | A Comprehensive Methodological Framework for Estimating Vacc... | 手法統合/融合, 統合パイプライン | 24 | 独自手法名「and evaluate」を提案; DOI付き参考文献 24件 |
| SCI-066 | 3/8 | Physics-Constrained Deep Learning Emulators for Earth System... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 10 | 複数技術の統合 (LSTM, Attention, Transformer); DOI付き参考文献 10件; 明示的な貢献リスト記載 |
| SCI-067 | 3/8 | AutoLCA: An AI-Driven Pipeline for Automated Life Cycle Asse... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 11 | 複数技術の統合 (Transformer, Graph Neural, Monte Carlo, GAN); ベースライン比較で優位性を主張; DOI付き参考文献 11件 |
| SCI-068 | 3/8 | An Integrated Modeling Framework for Predicting Ocean Acidif... | 手法統合/融合, 統合パイプライン, 生成モデル | 10 | 独自手法名「an integrated」を提案; DOI付き参考文献 10件 |
| SCI-076 | 3/8 | Multimodal Deep Learning for Paddy Rice Yield Prediction: In... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 11 | 複数技術の統合 (CNN, Attention, LSTM, GAN); ベースライン比較で優位性を主張; DOI付き参考文献 11件 |
| SCI-080 | 3/8 | An Integrated AI System for Food Supply Chain Safety Risk Pr... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 17 | 複数技術の統合 (SHAP, LSTM, Transformer, GAN); ベースライン比較で優位性を主張; DOI付き参考文献 17件 |
| SCI-086 | 3/8 | A Patient-Specific Cardiac Digital Twin Framework Integratin... | 手法統合/融合, 統合パイプライン, マルチモーダル | 12 | 独自手法名「an integrated」を提案; DOI付き参考文献 12件 |
| SCI-094 | 3/8 | EthicAI-Bench: A Unified Quantitative Framework for Multi-Di... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 10 | 独自手法名「EthicAI-Bench」を提案; DOI付き参考文献 10件 |
| SCI-099 | 3/8 | An Integrated ODE Framework for Multi-Hallmark Aging Dynamic... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 12 | 独自手法名「an integrated」を提案; DOI付き参考文献 12件 |

### 🔴 新規性「低」 (30件)

| ID | スコア | タイトル | 新規性タイプ | DOI数 | 評価根拠 |
|-----|--------|---------|-------------|-------|---------|
| SCI-004 | 2/8 | An Integrated Pharmacogenomics Framework for Drug Response P... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 24 | 複数技術の統合 (Graph Neural, Transformer, GAN); DOI付き参考文献 24件 |
| SCI-012 | 2/8 | An Efficient Simulation Framework for Large-Scale Spiking Ne... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 12 | DOI付き参考文献 12件; 明示的な貢献リスト記載 |
| SCI-023 | 2/8 | Multiscale Coarse-Grained Molecular Dynamics Prediction of B... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 19 | DOI付き参考文献 19件; 明示的な貢献リスト記載 |
| SCI-024 | 2/8 | An Integrated Computational Framework for Theoretical Design... | 手法統合/融合, 統合パイプライン | 24 | DOI付き参考文献 24件; 明示的な貢献リスト記載 |
| SCI-032 | 2/8 | Efficient Simulation Framework for Logical Error Rate Estima... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 28 | DOI付き参考文献 28件; 明示的な貢献リスト記載 |
| SCI-034 | 2/8 | Design and Performance Evaluation of Quantum Key Distributio... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 13 | ベースライン比較で優位性を主張; DOI付き参考文献 13件 |
| SCI-038 | 2/8 | Integrated Optimal Trajectory Design System for Multi-Target... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 22 | 複数技術の統合 (Reinforcement Learning, Attention, GAN); DOI付き参考文献 22件 |
| SCI-048 | 2/8 | Extending the Applicability of Physics-Informed Neural Netwo... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 16 | 複数技術の統合 (Diffusion, Attention, GAN); DOI付き参考文献 16件 |
| SCI-050 | 2/8 | A Systematic Comparison Framework for Causal Effect Estimati... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 8 | 先駆性を主張 |
| SCI-056 | 2/8 | A Structural Selection Framework for Infectious Disease Math... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 10 | ベースライン比較で優位性を主張; DOI付き参考文献 10件 |
| SCI-069 | 2/8 | An Integrated WRF-UCM Framework for Quantitative Prediction ... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 26 | ベースライン比較で優位性を主張; DOI付き参考文献 26件 |
| SCI-072 | 2/8 | Scalable Multi-Agent Path Finding: A Comparative Study of Op... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 18 | DOI付き参考文献 18件; 明示的な貢献リスト記載 |
| SCI-074 | 2/8 | Integrated VSLAM and Obstacle Avoidance System for Autonomou... | 手法統合/融合, 統合パイプライン | 38 | ベースライン比較で優位性を主張; DOI付き参考文献 38件 |
| SCI-085 | 2/8 | A Modular Computational Framework for Perturb-seq Data Analy... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 20 | DOI付き参考文献 20件; 明示的な貢献リスト記載 |
| SCI-087 | 2/8 | Digital Twin Framework for Injection Molding Quality Predict... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 7 | 複数技術の統合 (LSTM, Attention, GRU); ベースライン比較で優位性を主張 |
| SCI-088 | 2/8 | Integrated Urban Traffic Microsimulation and Real-Time Adapt... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 26 | 複数技術の統合 (SHAP, Reinforcement Learning, Graph Neural); DOI付き参考文献 26件 |
| SCI-089 | 2/8 | An Integrated Real-Time Simulation Framework for Power Grids... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 26 | ベースライン比較で優位性を主張; DOI付き参考文献 26件 |
| SCI-013 | 1/8 | N/A | 手法統合/融合, 統合パイプライン, 解釈可能性 | 8 | 複数技術の統合 (LSTM, Attention, Transformer) |
| SCI-025 | 1/8 | A Computational Framework for Molecular Design of Biodegrada... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 16 | DOI付き参考文献 16件 |
| SCI-026 | 1/8 | First-Principles Computational Framework for Elucidating Int... | 手法統合/融合, 新規ドメイン適用, 統合パイプライン | 10 | DOI付き参考文献 10件 |
| SCI-030 | 1/8 | An Integrated Simulation Framework for Supercritical Enhance... | 手法統合/融合, 統合パイプライン, 生成モデル | 24 | DOI付き参考文献 24件 |
| SCI-031 | 1/8 | Enhancing Noise Resilience of the Variational Quantum Eigens... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 24 | DOI付き参考文献 24件 |
| SCI-035 | 1/8 | A Systematic Performance Evaluation Framework for Quantum An... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 0 | ベースライン比較で優位性を主張 |
| SCI-057 | 1/8 | A Comprehensive Analytical Framework for Estimating Causal R... | 手法統合/融合, 統合パイプライン, 解釈可能性 | 12 | DOI付き参考文献 12件 |
| SCI-059 | 1/8 | A Geostatistical Framework for Disease Risk Spatial Pattern ... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 32 | DOI付き参考文献 32件 |
| SCI-070 | 1/8 | An Integrated Framework for Economic Valuation of Ecosystem ... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 8 | 明示的な貢献リスト記載 |
| SCI-077 | 1/8 | An Integrated Multiscale Modeling Framework for Predicting F... | 手法統合/融合, 統合パイプライン, ベンチマーク比較 | 30 | DOI付き参考文献 30件 |
| SCI-095 | 1/8 | A Comprehensive Quantitative Framework for Assessing the Imp... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 10 | DOI付き参考文献 10件 |
| SCI-090 | 0/8 | An Integrated BIM-Based Environmental Performance Simulation... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 7 | - |
| SCI-093 | 0/8 | Optimizing Efficiency and Equity in Research Funding Allocat... | 新規モデル/アーキテクチャ, 手法統合/融合, 統合パイプライン | 8 | - |

---

## 分析と考察

### 新規性の傾向

1. **統合パイプラインが主流** — 100件中100件が複数の既存技術を統合した「組み合わせ的新規性」を持つ。AIは個別技術を新たに発明するのではなく、既存技術を新しいパターンで組み合わせることで研究アプローチを構築する傾向が見られた。

2. **手法融合が高頻度** — 99%の論文形式文書が複数の手法を融合するアプローチを採用。CNN+Attention、GNN+Transformer等の組み合わせが典型的。

3. **独自命名は限定的** — 固有の手法名（EpiCRISPR-Net, DeepEpiClock等）を提案したのは約31%。多くは「An Integrated Framework for...」のような一般的な命名。なお、自動抽出では「an integrated」自体が手法名として誤検出されるケースがあり（SCI-003, SCI-006等）、スコアリングにはこの系統的バイアスが含まれる。

4. **DOI形式の参考文献の引用** — 97%の論文形式文書がDOI形式の参考文献を含む。平均16.3件のDOI形式引用（DOIの実在性および引用内容の正確性は未検証）。

5. **ベンチマーク比較は半数** — 50%の論文形式文書が既存手法との比較で優位性を主張。残り半数は比較なしまたは定性的な議論に留まる。

### 新規性の限界

1. **実験未実施** — 全論文形式文書のコード・結果はシミュレーションまたはモックデータに基づく。実データでの検証は行われていない。

2. **組み合わせ的新規性の位置づけ** — 100件すべてが「既存手法A + 既存手法B + ドメインC」という組み合わせ的アプローチに基づく。これは科学研究における正当な新規性の形態（多くの査読付き論文も同様の構造を持つ）であるが、既存のパラダイムを覆すような根本的に新しいアルゴリズムの提案は本分析の範囲では確認されていない。

3. **先駆性主張の信頼性** — 「to the best of our knowledge, this is the first...」等の表現が使われるが、候補文献の収集に基づく部分的な調査であり、網羅的な文献調査による検証ではない。

4. **再現性の課題** — 生成されたコードは構造的には妥当だが、実行して同等の結果が得られるかの検証は未実施。

---

## 結論: AI for Scienceにおける観察と示唆

### AIは「組み合わせ的新規性の生成」に長けるが、「問い」は設定していない

本ベンチマークの結果は、現時点のAIの能力と限界について以下の観察を示しています。

100件すべてにおいて、AIRAは与えられた研究テーマに対して、先行研究調査・コード実装・論文形式の文書作成を一貫して完了しました。ワークフロー完了率100%、平均12.5分という速度は、研究成果物の生成におけるAIの実行力を示しています。

さらに注目すべきは、**100件中100件が複数の既存技術を組み合わせた新しい研究アプローチを構築している**点です。例えば、CNNとAttention機構とエピジェネティクスデータを統合した「EpiCRISPR-Net」のように、個々の技術は既存でも、その組み合わせとドメイン適用は独自のものです。このような「組み合わせ的新規性」は、多くの査読付き科学論文が持つ新規性の形態と同種であり、AIがこれを高速かつ大量に生成できることは、AI for Scienceにおける重要な能力です。

一方で、新規性マーカー分析の結果は、AIの出力パターンの現時点での境界も示唆しています。

自動分類の範囲では、100件すべてが「組み合わせ」型のアプローチであり、既存のパラダイムを覆すような根本的に新しいアルゴリズムの提案と判断できる事例は確認されていません（ただし、専門家による個別評価は未実施です）。

### 観察から得られた示唆

> **注**: 以下は100件の観察から得られた示唆であり、実験的に証明された結論ではありません。AIに研究テーマ自体を設定させる対照実験を行っていないため、「AIが問いを立てられない」ことを本実験から断定することはできません。

本ベンチマークにおいて、以下の活動はすべて人間によって行われ、AIは関与していませんでした。

1. **問いの設定（What to solve）** — 100件の研究テーマはすべて人間が設計したプロンプトに由来する。
2. **研究意義の判断（Why it matters）** — テーマの選定理由や社会的重要性の判断は人間が行った。
3. **実験設計全体の方向づけ（How to approach）** — プロンプトの3段階構成やイテレーション改善の方針は人間が決定した。

この観察は、**AI for Scienceにおいて、AIは組み合わせ的新規性の生成に長けるが、「何を問うか」を決める研究者の役割は依然として重要**という仮説と整合的です。ただし、これを実験的に検証するには、AIに自律的に研究テーマを設定させた場合の出力との比較が必要です。

### AI for Scienceにおける現時点の分業モデル（仮説）

| 役割 | 担い手 | 本ベンチマークでの観察 |
|------|--------|---------------------|
| 課題設定・仮説構築 | 🧑‍🔬 研究者 | 100件すべて人間がテーマ設計。AIによるテーマ設定は未検証 |
| 文献調査・知識統合 | 🤖 AI + 🧑‍🔬 研究者 | AIは平均16.3件のDOI付き文献を収集。引用の正確性は未検証 |
| 実験設計・実装 | 🤖 AI（研究者監督下） | パイプライン構築を自動化。コードの実行検証は未実施 |
| 結果の解釈・考察 | 🧑‍🔬 研究者 | AI生成の考察は表層的パターンの記述にとどまる傾向 |
| 論文執筆 | 🤖 AI + 🧑‍🔬 研究者 | 論文形式の文書を自動生成。科学的推敲は研究者が必要 |

### 今後の検証課題

本ベンチマークの結論をより堅固にするためには、以下の追加実験・検証が必要です。

1. **AIによる問い設定実験** — AIに研究テーマ自体を設定させ、人間設定テーマとの比較を行う
2. **専門家によるブラインド評価** — 生成論文のうち層化抽出した20件について、各分野の研究者が科学的妥当性・新規性・実行可能性を評価する
3. **引用正確性の監査** — DOIの実在確認、引用内容と本文記述の整合性検証
4. **生成コードの実行検証** — 生成されたPythonコードが実際に動作し、論文中の結果と一致するかの確認
5. **再現性分析** — 同一プロンプトの複数回実行によるLLM出力の安定性評価
6. **コントロール群の設定** — 人間の研究者が同じテーマで書いた論文との比較

### まとめ

AI for Scienceの進展において、本ベンチマークが示唆するのは、**AIの能力を正確に理解し、適切な期待値を持つこと**の重要性です。AIは課題が明確に定義されれば、研究成果物の生成において高い実行力を発揮します。しかし、生成物の科学的価値——仮説の妥当性、結果の信頼性、知見の新規性——は、研究者による批判的評価なしには判断できません。

本ベンチマークが、AI支援研究の可能性と限界を議論するための一つの材料となれば幸いです。
