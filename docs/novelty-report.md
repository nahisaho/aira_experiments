# AIRA 100実験 新規性評価レポート

## エグゼクティブサマリー

AIRA Co-Scientistが生成した100本の科学論文（paper.md）について、各研究の新規性を評価しました。
評価は以下の6つの基準に基づくスコアリング方式（0〜8点）で実施しています。

### 評価基準

| 基準 | 配点 | 説明 |
|------|------|------|
| 独自手法名の提案 | 2点 | 固有名詞付きの手法・モデルを提案しているか |
| 先駆性の主張 | 2点 | 「to the best of our knowledge」等の先駆性表現 |
| 複数技術の統合 | 1点 | 3つ以上の異なるML/DL技術を組み合わせているか |
| ベースライン比較 | 1点 | 既存手法との定量的な比較で優位性を主張 |
| DOI付き参考文献 | 1点 | 10件以上の実在する学術論文を引用 |
| 明示的貢献リスト | 1点 | Introductionで貢献点を明確にリスト化 |

### 総合結果

| 指標 | 値 |
|------|------|
| 評価対象論文数 | 100 |
| 平均スコア | 3.7 / 8 |
| 新規性「高」 (5点以上) | **39件** (39%) |
| 新規性「中」 (3-4点) | 31件 (31%) |
| 新規性「低」 (0-2点) | 30件 (30%) |

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

1. **統合パイプラインが主流** — 100件中100件が「統合パイプライン」型の新規性を持つ。AIRAは個別技術の発明よりも、既存技術の新しい組み合わせ・統合を得意とする傾向。

2. **手法融合が高頻度** — 99%の論文が複数の手法を融合するアプローチを採用。CNN+Attention、GNN+Transformer等の組み合わせが典型的。

3. **独自命名は限定的** — 固有の手法名（EpiCRISPR-Net, DeepEpiClock等）を提案したのは約31%。多くは「An Integrated Framework for...」のような一般的な命名。

4. **DOI付き実在文献の引用** — 97%の論文がDOI付き参考文献を含み、先行研究調査プロンプトの効果が確認できる。平均16.3件のDOI付き引用。

5. **ベンチマーク比較は半数** — 50%の論文が既存手法との比較で優位性を主張。残り半数は比較なしまたは定性的な議論に留まる。

### 新規性の限界

1. **実験未実施** — 全論文のコード・結果はシミュレーションまたはモックデータに基づく。実データでの検証は行われていない。

2. **組み合わせ的新規性** — 多くの論文の新規性は「既存手法A + 既存手法B + ドメインC」という組み合わせに依存。根本的に新しいアルゴリズムの提案は稀。

3. **先駆性主張の信頼性** — 「to the best of our knowledge, this is the first...」等の表現が使われるが、網羅的な文献調査に基づく検証ではない。

4. **再現性の課題** — 生成されたコードは構造的には妥当だが、実行して同等の結果が得られるかの検証は未実施。

---

## 結論: AI for Scienceにおける研究者の役割

### AIは優秀な「実行者」だが「問いの発見者」ではない

本実験の結果は、AI for Scienceにおいて**研究者による課題設定が依然として不可替である**ことを明確に示している。

100件の実験すべてにおいて、AIRAは与えられた研究テーマに対して高い実行力を発揮した。先行研究を調査し、複数の技術を統合したパイプラインを設計し、コードを実装し、学術論文形式で成果をまとめる——この一連のプロセスは完全に自動化された。成功率100%、平均12.5分という速度は、研究の「実行フェーズ」におけるAIの有用性を実証している。

しかし、新規性の評価結果はAIの本質的な限界を浮き彫りにする。

**100件中100件が「統合パイプライン」型の新規性**に分類されたという事実は、AIが自律的に研究課題を発見できないことの証左である。AIRAが生み出す新規性は、一貫して「既存手法Aと既存手法Bをドメインcに適用する」という組み合わせ的なものであり、根本的に新しいアルゴリズムや、既存のパラダイムを覆すような着想は一切生まれていない。

### 研究者にしかできない3つの知的活動

本実験から、以下の知的活動がAIでは代替できないことが明らかになった。

1. **問いの設定（What to solve）** — 「何が未解決の重要問題か」を見極める判断力。100件の研究テーマはすべて人間が設計したプロンプトに由来し、AIが自発的に問いを立てた事例はゼロである。

2. **研究意義の判断（Why it matters）** — 「なぜこの問いが今、重要なのか」という文脈の理解。AIは先行研究を引用できるが、研究コミュニティの動向や社会的ニーズを踏まえた意義づけは行えない。

3. **本質的ギャップの発見（What is fundamentally missing）** — 「既存手法の何が根本的に不十分か」を見抜く洞察力。AIの論文が指摘する「ギャップ」は、多くの場合「手法AはデータタイプBを扱えない」といった表層的なものに留まる。

### AI for Scienceの望ましい分業構造

| 役割 | 担い手 | 具体的活動 |
|------|--------|-----------|
| 課題設定・仮説構築 | 🧑‍🔬 研究者 | 未解決問題の特定、研究意義の判断、仮説の立案 |
| 文献調査・知識統合 | 🤖 AI + 🧑‍🔬 研究者 | 網羅的な文献検索はAI、批判的評価は研究者 |
| 実験設計・実装 | 🤖 AI（研究者監督下） | コード実装、パイプライン構築、パラメータ探索 |
| 結果の解釈・考察 | 🧑‍🔬 研究者 | 結果の意味づけ、限界の認識、次の問いの導出 |
| 論文執筆 | 🤖 AI + 🧑‍🔬 研究者 | ドラフト生成はAI、批判的推敲と最終判断は研究者 |

### 提言

AI for Scienceの進展において、最も危険なのは「AIが科学研究を完全に自動化できる」という誤解である。本実験が示すように、AIは**課題が明確に定義されれば**極めて有能な実行者となる。しかし、科学の本質的な価値——新しい問いを立て、既存の枠組みを疑い、知の地平を拡張すること——は、依然として人間の研究者に委ねられている。

AIを科学研究に活用する上での鍵は、**研究者が課題設定と批判的評価に集中し、AIに実行と探索を委ねる**という分業体制の確立にある。本実験の100本の論文が示すように、優れた課題設定があれば、AIはその課題に対して多様な手法を統合し、高速に実験を遂行し、学術的な形式で成果をまとめることができる。研究者の知的活動とAIの実行力を適切に組み合わせることこそが、AI for Scienceの真の可能性を引き出す道筋である。
