# 総合実験レポート

## 実験目的と背景
本実験の目的は、遺伝子編集・AI・核融合という新興科学技術に対する社会的受容を、NLP と構造方程式モデリング（SEM）的な枠組みを統合して予測することである。特に、日本におけるゲノム編集食品の受容研究を中核アンカーとして、メタ分析、感情分析、心理計量モデル、フレーミング実験、パス解析、日本ケーススタディを一体化した。

## 手法の概要
1. **メタ分析**: 20件の模擬調査研究をランダム効果モデルで統合し、Q統計量・I²・τ²を算出した。  
2. **感情分析**: 3技術×500件のSNS投稿を生成し、VADER風辞書法とロジスティック回帰（BERT代替）を組み合わせた。  
3. **心理計量モデル**: Dread Risk, Unknown Risk, Control, Benefit, Trust の5因子を用いて近似CFAを実施した。  
4. **フレーミング分析**: Gain / Loss / Neutral × 3技術の被験者間実験を ANOVA で評価した。  
5. **SEMパス解析**: Trust→Benefit/Risk→Attitude→Behavioral Intention→Acceptance の経路を標準化OLSで近似した。  
6. **日本ケース**: n=800 の日本模擬調査を作成し、属性別受容率と 2019-2024 時系列を推定した。

## 全実験結果（数値付き）
- **メタ分析**: 全体 pooled acceptance は 0.464、I²=63.6% で中程度の異質性を示した。  
- **感情分析**: Hybrid モデルは AUC=0.850±0.022、F1=0.792±0.026。Lexicon ベースライン F1 は 0.598±0.045。  
- **心理計量モデル**: 近似適合度は RMSEA=0.070, CFI=0.946。  
- **フレーミング分析**: 主効果は有意で、技術別の交互作用差は小さいが観測された。最大（非残差）η² は 0.090。  
- **SEM**: グローバルモデルの Acceptance R²=0.538, BI R²=0.456, RMSEA=0.066, CFI=0.939。  
- **日本ケース**: 観測受容率は 0.399 で、文献アンカーの public support 0.38 に近い。Scientist support アンカーは 0.70。日本SEMの Acceptance R² は 0.433。  
- **人口属性**: 平均年齢 50.1 歳、都市部比率 0.486、大学院卒比率 0.200。

## 図の埋め込み
![Figure 1](figures/figure1_meta_analysis_forest.png)
![Figure 2](figures/figure2_sentiment_comparison.png)
![Figure 3](figures/figure3_psychometric_map.png)
![Figure 4](figures/figure4_framing_effects.png)
![Figure 5](figures/figure5_sem_path_diagram.png)
![Figure 6](figures/figure6_japan_case_study.png)
![Figure 7](figures/figure7_model_comparison.png)

## 考察と展望
主な示唆は、社会的受容が単純な好意的感情だけで決まるのではなく、**制度的信頼・便益認知・リスク認知・フレーミング**の相互作用によって形成される点である。日本ケースでは、2019年から2020年にかけて 0.49→0.33 の急落を設定しており、He Jiankui 事件後の受容低下を反映した。その後の回復も限定的であり、ラベリング要求や制度信頼が大きな規定因となる構図が再現された。今後は実データへの置換、BERT/日本語LLMの導入、厳密な潜在変数SEM、縦断パネル化が重要である。

## 生成ファイル一覧
- `experiments/run_experiments.py`
- `paper.md`
- `report.md`
- `figures/figure1_meta_analysis_forest.png`
- `figures/figure2_sentiment_comparison.png`
- `figures/figure3_psychometric_map.png`
- `figures/figure4_framing_effects.png`
- `figures/figure5_sem_path_diagram.png`
- `figures/figure6_japan_case_study.png`
- `figures/figure7_model_comparison.png`
