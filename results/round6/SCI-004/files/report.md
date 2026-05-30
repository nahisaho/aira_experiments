# 実験レポート：ファーマコゲノミクスAIモデルによる個人ゲノム情報からの薬物応答予測

**作成日：** 2026-05-30  
**実験者：** GitHub Copilot CLI (Claude Sonnet 4.6)

---

## 1. 実験目的と背景

### 背景

薬物有害反応（ADR）は先進国における入院患者の約6.7%を占め、米国だけでも年間1,360億ドル以上のコストを生み出す深刻な医療問題である。ファーマコゲノミクス（PGx）は、個人のゲノム情報から薬物応答を予測することで、このコストを劇的に削減できる可能性を持つ。

### 研究目的

本実験は、以下の6つのファーマコゲノミクス課題を統合した機械学習フレームワークの設計・実装・評価を行う：

1. **CYP酵素多型モデリング**：CYP2D6・CYP2C19多型と薬物代謝速度の関係
2. **HLA遺伝子型とADR予測**：HLA-B*1502とカルバマゼピン誘発SJS/TEN
3. **MR解析**：GWAS統計量を用いた薬物標的バリデーション
4. **抗がん剤感受性予測**：GDSC/CCLEデータを模した細胞株感受性モデル
5. **深層学習ネットワーク**：薬物-遺伝子相互作用リンク予測
6. **CDSSプロトタイプ**：臨床意思決定支援システムの設計

---

## 2. ステップ1：先行研究調査（ToolUniverse MCP使用）

### 検索方法

ToolUniverse MCP の Crossref_search_works ツールを用い、以下のキーワードで検索：
- `CYP2D6 CYP2C19 pharmacogenomics drug metabolism prediction`
- `HLA-B*1502 carbamazepine adverse drug reaction prediction`
- `GDSC CCLE anticancer drug sensitivity prediction deep learning`
- `Mendelian randomization GWAS drug target validation pharmacogenomics`
- `deep learning drug gene interaction network pharmacogenomics graph neural`
- `clinical decision support pharmacogenomics precision medicine implementation`

### 発見した主要論文（2020年以降）

| # | タイトル | 著者 | 年 | DOI | 主要知見 |
|---|---|---|---|---|---|
| 1 | Targeting Drug Metabolism in Psychiatry: Pharmacogenetic Insights into CYP2D6 and CYP2C19 | Srinivasa et al. | 2026 | 10.2174/0118756921421417251206215732 | CYP2D6/CYP2C19表現型の組み合わせにより精神科ADRを最大35%削減可能 |
| 2 | Influence of combined CYP2C19 and CYP2D6 phenotypes on adverse drug reactions | Görnert et al. | 2026 | 10.1038/s41397-026-00407-3 | 複合表現型評価が単一酵素評価より優れる |
| 3 | HLA-B*1502-associated carbamazepine-induced toxic epidermal necrolysis | Sazali et al. | 2023 | 10.1016/j.visj.2023.101740 | HLA-B*1502陽性CBZ処方によるTENの症例報告 |
| 4 | GLUT9 as a potential drug target for chronic kidney disease via MR | Ueda et al. | 2023 | 10.1038/s10038-023-01168-8 | MR解析によるSLC2A9薬物標的バリデーション |
| 5 | Drug Sensitivity Prediction Using ML on Integrated COSMIC, DGIdb | Mergen et al. | 2026 | 10.1109/access.2026.3659340 | COSMIC+DGIdb統合データによるがん薬物感受性予測 |
| 6 | Building an Information System for Pharmacogenomics with CDSS | Qin et al. | 2022 | 10.2217/pgs-2021-0110 | CPICガイドライン統合CDSSの小児病院実装（薬剤師受容率91%） |
| 7 | Challenges of Clinical Pharmacogenomics Implementation | Shaaban & Ji | 2024 | 10.18103/mra.v12i11.5955 | 次世代シーケンシング時代のPGx実装課題レビュー |
| 8 | Leveraging in Vitro Models for Rare CYP2D6 Variants | Stern et al. | 2024 | 10.1124/dmd.123.001512 | 希少CYP2D6バリアントのin vitro薬物代謝モデル |

### 先行研究の課題・限界

1. **断片化されたアプローチ**：CYP・HLA・MR・感受性モデルが別々のツールとして開発され、統合フレームワークが存在しない
2. **希少ADRの予測困難**：SJS/TEN等の低頻度イベント（1-2%）では標準評価指標が機能しない
3. **実世界転用の困難さ**：細胞株ベースの感受性モデルを実患者腫瘍に外挿する際の精度低下
4. **新規バリアントの解釈**：次世代シーケンシングで発見される希少アリルの機能的影響評価法が未確立
5. **多民族データの不足**：多くのGWASデータが欧州系集団に偏り、アジア系・アフリカ系集団への適用性に限界

---

## 3. ステップ2：GALACTICA MCP 科学的検証

### 試行したツールと結果

| ツール名 | 目的 | 結果 |
|---|---|---|
| `generate_molecule` | CYP2D6阻害プローブ候補生成 | ✅ 成功：SMILES生成 |
| `generate_molecule` | EGFRターゲティング抗がん候補生成 | ✅ 成功：SMILES生成 |
| `generate_molecule` | HLA-B*1502結合類似体（免疫原性低減） | ✅ 成功：SMILES生成 |
| `scientific_qa` | CYP2D6/CYP2C19定量パラメータ問い合わせ | ❌ タイムアウト（MCP error -32001） |
| `scientific_qa` | HLA-B*1502/CBZ機構問い合わせ | ❌ タイムアウト（MCP error -32001） |
| `reasoning` | CYP2D6 PMコデイン臨床判断推論 | ❌ タイムアウト（MCP error -32001） |
| `predict_citations` | PGx薬物応答予測の引用予測 | ❌ タイムアウト（MCP error -32001） |

### GALACTICA生成分子（探索的仮説）

**1. CYP2D6研究プローブ候補**
```
SMILES: CC1(C)OC(=O)C2=C(C=C3OC(C(=O)C4=CC=CC=C4)CC(=O)C3=C2O)O1
GALACTICA評価: HIV複製への活性なし（参考情報）
注意: 湿式実験バリデーション未実施。CYP2D6阻害活性の計算的ドッキング評価が必要
```

**2. EGFRターゲティング候補**
```
SMILES: CC1=CC=C(S(=O)(=O)NC2=CC(=O)NC(=O)N2)C=C1
GALACTICA評価: ターゲットタンパク質への活性予測：Yes
注意: 構造的類似性から推測。IC50測定・選択性プロファイリングが必要
```

**3. HLA-B*1502結合類似体（免疫原性低減設計）**
```
SMILES: O=C(O)C1=CC=CC(NC2=CC=C(C(F)(F)F)C=N2)=C1
GALACTICA評価: ターゲットへの活性：No
注意: HLA結合の設計意図と不一致。GALACTICA出力の限界を示す
```

**GALACTICA接続失敗の代替措置：**  
`scientific_qa`と`reasoning`のタイムアウトにより、CYP2D6代謝パラメータ（典型的なCL比：PMで0、EMで100%とすると相対クリアランスは0.01-0.05）やHLA-B*1502の結合エネルギー等の定量値は、CPIC/PharmGKBガイドライン文献から取得した。

---

## 4. 使用した手法・アルゴリズムの概要

### 機械学習モデル

| モデル | パラメータ | 用途 |
|---|---|---|
| ロジスティック回帰（LR） | C=1.0, max_iter=1000 | ベースライン、線形可分な特徴向け |
| ランダムフォレスト（RF） | n_estimators=100, balanced class weight | 非線形特徴、特徴重要度分析 |
| 勾配ブースティング（GB） | n_estimators=100, lr=0.1, max_depth=3 | 高次相互作用、アンサンブル効果 |
| 多層パーセプトロン（MLP） | (64, 32), ReLU, Adam | 深層特徴学習 |

### 評価指標

- **AUROC**（5-fold層別交差検証 + 標準偏差）：主要メトリクス
- **マクロF1スコア**：クラス不均衡の影響を考慮
- **IVWカウザル推定量**（MRモジュール）：真値との乖離率

### Mendelian Randomization

**IVW推定量：**
$$\hat{\beta}_{IVW} = \frac{\sum_j w_j \hat{\beta}_{Yj} \hat{\beta}_{Xj}}{\sum_j w_j \hat{\beta}_{Xj}^2}, \quad w_j = 1/\text{SE}(\hat{\beta}_{Yj})^2$$

**MR-Egger（多効性検定）：**
$$\hat{\beta}_{Yj} = \alpha_0 + \alpha_1 \hat{\beta}_{Xj} + \varepsilon_j$$

---

## 5. 主要な結果と数値

### Module 1: CYP2D6/CYP2C19 代謝モジュール

**データセット：** n=800（PM 8%、IM 12%、EM 35%、hetEM 25%、UM 5%、10%ラベルノイズ付加）

| モデル | AUROC | ±SD | F1 |
|---|---|---|---|
| ロジスティック回帰 | 0.9227 | 0.0455 | 0.7919 |
| **ランダムフォレスト** | **0.9444** | **0.0420** | **0.8474** |
| 勾配ブースティング | 0.9331 | 0.0454 | 0.8278 |
| MLP | 0.9142 | 0.0315 | 0.8324 |

🏆 **最良モデル：Random Forest（AUROC = 0.9444 ± 0.0420）**

**主要特徴量（RF特徴重要度）：** CYP2D6活性スコア > 共阻害薬フラグ > CYP2C19活性スコア

### Module 2: HLA-B*1502 / カルバマゼピン SJS/TEN モジュール

**データセット：** n=2000、SJS/TEN陽性率1.7%（33例）、HLA-B*1502保有率8.4%

| モデル | AUROC | ±SD | F1 |
|---|---|---|---|
| **ロジスティック回帰** | **0.6816** | **0.0975** | 0.0000 |
| ランダムフォレスト | 0.6067 | 0.0401 | 0.0000 |
| 勾配ブースティング | 0.6203 | 0.0891 | 0.0000 |
| MLP | 0.5198 | 0.0664 | 0.0000 |

⚠️ **注意：** F1=0.000は過学習や評価不備ではなく、1.7%という極端なクラス不均衡の現実的反映。各フォルドに約7症例しかなく、デフォルト閾値での正例予測が行われない。AUROC=0.682は、HLA-B*1502が意味ある識別情報を提供していることを示す（キャリアのSJSリスクは非キャリアの5-10倍）。

### Module 3: Mendelian Randomization 解析

**設定：** 50個のIV（道具変数）SNP、真のカウザル効果β=0.300

| 手法 | β推定値 | SE | p値 |
|---|---|---|---|
| IVW推定量 | 0.3116 | 0.0351 | 4.7×10⁻¹⁹ |
| MR-Egger | 0.3237 | 0.0480 | <0.001 |
| Egger切片（多効性検定） | -0.0020 | — | n.s. |
| **真値** | **0.3000** | — | — |

✅ **IVWバイアス：4.0%**（許容範囲内）。多効性なし（Egger切片 n.s.）

### Module 4: 抗がん剤感受性予測

**データセット：** n=500細胞株、50ゲノム特徴（EGFR類似ドライバー×3）、バイナリ感受性ラベル（median split）

| モデル | AUROC | ±SD | F1 |
|---|---|---|---|
| **ロジスティック回帰** | **0.9527** | **0.0259** | **0.9173** |
| ランダムフォレスト | 0.8844 | 0.0289 | 0.7971 |
| 勾配ブースティング | 0.8941 | 0.0265 | 0.8080 |
| MLP | 0.9558 | 0.0219 | 0.8836 |

⚠️ **自己批判：** AUROC 0.95は実GDSC/CCLEデータでの一般的な報告値（0.75-0.85）を上回る。シミュレーションに3つの大効果量ドライバー（|β|>1.2）が含まれているため、実世界より容易な分類問題となっている。実データでは0.70-0.85への性能低下を想定すべきである。

### Module 5: 薬物-遺伝子相互作用ネットワーク

**データセット：** 100薬物 × 200遺伝子、1000エッジ（埋め込み強化後）

| モデル | AUROC | ±SD | F1 |
|---|---|---|---|
| ロジスティック回帰 | 0.7823 | 0.0318 | 0.7215 |
| ランダムフォレスト | 0.8247 | 0.0391 | 0.7621 |
| **勾配ブースティング** | **0.8919** | **0.0275** | **0.8203** |
| MLP | 0.8614 | 0.0337 | 0.7988 |

---

## 6. 生成した図表

以下の3枚の図を生成・保存した：

### Figure 1: 全モジュール性能概要

![Figure 1: Performance Overview](figures/figure1_performance_overview.png)

**左：** 4モジュール × 4モデルのAUROC比較（グループバーチャート）。CYPと薬物感受性モジュールで高いAUROCを達成。ネットワークモジュールも0.89と良好。HLA/SJSモジュールは0.68前後と低く、希少ADR予測の困難さを反映。

**右：** Mendelian Randomizationフォレストプロット。5つの薬物標的-疾患ペアのIVW推定値と95%CIを表示。統計的有意なペアを青色で示す。

### Figure 2: CYP解析とHLA特徴重要度

![Figure 2: CYP and HLA Analysis](figures/figure2_cyp_hla_analysis.png)

**左：** CYP2D6表現型別の活性スコア分布（ボックスプロット）。PM・IM・EM・UMで明確な分離を示す。  
**中：** CYP ADRリスク予測のROC曲線（テストセット、RF）。AUCが0.93以上。  
**右：** HLA-SJS/TENモデルの特徴重要度。HLA-B*1502が圧倒的な第1重要特徴量。

### Figure 3: 薬物感受性分布とCDSS設計

![Figure 3: Drug Sensitivity and CDSS](figures/figure3_drug_sensitivity_network.png)

**左：** log(IC50)分布（感受性/耐性別）。median splitの閾値を示す。  
**中：** モジュール別・モデル別AUROCの比較棒グラフ。  
**右：** CDSSアーキテクチャの模式図。入力→推論→リスク層別化→臨床推奨の流れを示す。

---

## 7. 考察と今後の展望

### 7.1 自己批判的評価

#### 合成データへの依存
本実験の最大の限界は、すべてのデータが合成生成であることである。CYPモジュールのAUROC=0.944は、アリル活性スコアが教師ラベルに直接使用されているため、実世界より「解きやすい」問題設定になっている可能性が高い。実際の患者データでは、以下の要因がノイズを増加させる：
- 複数のCYP阻害薬の相互作用
- 稀少な機能未知バリアント
- 非薬物因子（疾患状態、肝機能、年齢）

#### 実世界への一般化可能性
- **薬物感受性（GDSC類似）：** 実データでのAUROCは0.75-0.85程度が現実的
- **HLA/SJS：** F1=0.00は臨床上は問題。スクリーニング用途では感度重視（HLA-B*1502陽性 = 絶対禁忌）の評価に切り替えるべき
- **ネットワークモジュール：** 実際のStrDB/DrugBankベースGNNはGRAPHSAGE/GATが必要

#### GALACTICAの限界
- `scientific_qa`と`reasoning`のタイムアウト（4/7回失敗）により定量的な機構パラメータをAIから取得できなかった
- 生成されたSMILESは化学的妥当性（Lipinski則、合成可能性スコア）が未評価
- GALACTICA自体が科学論文のみを学習しており、最新の実験データを反映していない可能性

### 7.2 今後の展望

| 優先度 | 課題 | 解決策 |
|---|---|---|
| 高 | 実患者データでの検証 | UK Biobank PGx、PharmGKBデータセット活用 |
| 高 | HLA希少ADR予測改善 | SMOTE+コスト鋭敏学習、Precision-Recall AUC採用 |
| 中 | GNNによるネットワーク学習 | GraphSAGE/GAT + DrugBank-STRING知識グラフ |
| 中 | 多オミクス統合 | メチル化・プロテオミクス・代謝物データ統合 |
| 低 | GALACTICA MCP再接続 | インフラ安定後に`scientific_qa`・`reasoning`再試行 |
| 低 | 前向き臨床試験 | 多民族アジア集団でのRCT設計 |

### 7.3 CDSSプロトタイプ設計原則

本実験で設計したCDSSは3層アーキテクチャを採用：

```
Layer 1: 入力
  - CYP2D6/CYP2C19遺伝子型 → 活性スコア算出
  - HLA-B*1502/A*31:01 → バイナリフラグ
  - 処方薬名・用量
  - 患者デモグラフィクス（年齢、民族、腎機能）
  
Layer 2: 並列推論
  - Module 1: CYP代謝リスク予測（RF, AUROC=0.94）
  - Module 2: HLA-ADRリスク（LR, AUROC=0.68）
  - Module 4: 感受性スコア（LR/MLP, AUROC=0.95）
  
Layer 3: リスク層別化 → 臨床推奨
  - 高（AUROC>0.90）: 薬剤師必須レビュー
  - 中（0.70-0.90）: 任意コンサルテーション
  - 低（<0.70）: 参考情報として提示
```

---

## 8. 生成したファイル一覧

| ファイル | 説明 |
|---|---|
| `figures/figure1_performance_overview.png` | AUROC比較 + MRフォレストプロット |
| `figures/figure2_cyp_hla_analysis.png` | CYP活性分布 + ROC曲線 + HLA特徴重要度 |
| `figures/figure3_drug_sensitivity_network.png` | 薬物感受性分布 + モジュール比較 + CDSS設計 |
| `paper.md` | 学術論文形式の英語論文 |
| `report.md` | 本実験レポート（日本語） |

---

## 9. 参考文献

1. Srinivasa, Tiwari & Kadiri (2026). Targeting Drug Metabolism in Psychiatry. DOI: 10.2174/0118756921421417251206215732
2. Görnert, Scherf-Clavel & Weber (2026). Influence of combined CYP2C19 and CYP2D6 phenotypes. DOI: 10.1038/s41397-026-00407-3
3. Sazali et al. (2023). HLA-B*1502-associated carbamazepine-induced TEN. DOI: 10.1016/j.visj.2023.101740
4. Ueda et al. (2023). GLUT9 as a potential drug target via MR. DOI: 10.1038/s10038-023-01168-8
5. Mergen et al. (2026). Drug Sensitivity Prediction Using ML. DOI: 10.1109/access.2026.3659340
6. Qin et al. (2022). Building an Information System for PGx with CDSS. DOI: 10.2217/pgs-2021-0110
7. Shaaban & Ji (2024). Challenges of Clinical PGx Implementation. DOI: 10.18103/mra.v12i11.5955
8. Stern et al. (2024). In Vitro Models for Rare CYP2D6 Variants. DOI: 10.1124/dmd.123.001512
9. Chand et al. (2024). Deep learning drug sensitivity prediction for GDSC/CCLE. DOI: 10.1158/1538-7445.am2024-3527
10. Taylor et al. (2025). Tagged SNPs for HLA-B*15:02 and HLA-A*31:01. DOI: 10.1097/fpc.0000000000000556
