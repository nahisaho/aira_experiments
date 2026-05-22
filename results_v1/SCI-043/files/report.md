# ゲノムスケール代謝モデル（GEM）制約条件ベースフラックス解析フレームワーク

**DRAFT — NOT FOR DISTRIBUTION**

**実行日時**: 2026-05-23  
**モデル**: E. coli core model (e_coli_core) — COBRApy textbook model  
**ツール**: COBRApy 0.31.1, Python 3.12

---

## 1. 実験目的と背景

ゲノムスケール代謝モデル（GEM）の制約条件ベースフラックス解析（Constraint-Based Flux Analysis）を改善する統合フレームワークを設計・実装した。大腸菌コアモデルを用いて、以下の6つのモジュールを体系的に検証した：

1. **FBA制約条件設定最適化** — 標準FBA、pFBA、FVAの比較と感度解析
2. **13C-MFA統合** — 13C代謝フラックス解析データによるFBA制約の強化
3. **動的FBA（dFBA）** — バッチ培養の時間変化動態シミュレーション
4. **酵素容量制約（GECKO/sMOMENT）** — プロテオーム制約がフラックス分布に与える影響
5. **条件特異的モデル構築** — RNA-seqデータ統合によるGIMMEアルゴリズム
6. **リシン生産最適化** — 代謝工学のケーススタディ（生産エンベロープ、OptKnock）

---

## 2. 使用した手法・アルゴリズムの概要

### 2.1 FBA制約条件最適化（Module 1）

| 手法 | 概要 |
|------|------|
| **Standard FBA** | 線形計画法によるバイオマス最大化 |
| **pFBA** (Parsimonious FBA) | 最適成長を維持しつつ総フラックスを最小化 |
| **FVA** (Flux Variability Analysis) | 最適解の90%以上で各反応のフラックス範囲を探索 |
| **感度解析** | グルコース・酸素取り込み速度の系統的変化 |

### 2.2 13C-MFA統合（Module 2）

- **EMU法**（簡略版）による中心炭素代謝12反応の13C標識パターンシミュレーション
- **合成測定データ**生成（5%ノイズ）と2σ制約によるモデル拘束
- **χ²適合度検定**によるモデル—実験データ整合性評価
- **フラックス分配比**解析（G6Pノード、PEP/PYRノード）

### 2.3 動的FBA（Module 3）

- **SOA**（Static Optimization Approach）：各タイムステップでFBAを解くオイラー積分法
- **Michaelis-Menten動力学**：グルコース（Km=0.5 mM）、酸素（Km=0.003 mM）
- **ジオーキシック成長**検出：グルコース枯渇後のアセテート再利用
- 初期グルコース濃度（5, 10, 20, 40 mM）の比較

### 2.4 酵素容量制約（Module 4）

- **sMOMENT風**アプローチ：BRENDAデータベースのkcat値（20酵素）と分子量を使用
- 酵素量 E = v × MW / kcat（二分探索 + pFBA）
- 酵素予算を0.2x〜1.5xに変化させた成長率への影響解析

### 2.5 条件特異的モデル（Module 5）

- **GIMMEアルゴリズム**：低発現遺伝子の反応を制約し、代謝活性反応を特定
- 3条件（好気、嫌気、ストレス）の合成RNA-seqデータ生成
- **Jaccard類似度**による条件間の活性反応パターン比較

### 2.6 リシン生産最適化（Module 6）

- **DAP経路**（lumped）をコアモデルに追加: OAA + PYR + 2NADPH + Glu + SuccCoA → Lys
- **生産エンベロープ**：成長率vs最大リシン生産のトレードオフ解析
- **OptKnock**探索：単一・二重遺伝子ノックアウトの網羅的スクリーニング
- **過剰発現ターゲット**解析：主要反応の上限2倍化による効果

---

## 3. 主要な結果と数値

### 3.1 FBA制約条件最適化

| 指標 | 値 |
|------|-----|
| 標準FBA成長率 | **0.8739 h⁻¹** |
| pFBA総フラックス | 518.42 mmol/gDW/h |
| FVA平均フラックス範囲 | 31.84 mmol/gDW/h |
| 最大変動反応 | FRD7, SUCDi, FORt2 |
| 好気条件成長率 | 0.8739 h⁻¹ |
| 嫌気条件成長率 | 0.2117 h⁻¹ |

- pFBAと標準FBAのフラックス分布は本コアモデルでは同一であった（最適解が一意に近い）
- FVAにより、FRD7/SUCDi間の可逆的フラックスが最大変動（range=1000）を示した

### 3.2 13C-MFA統合

| 指標 | 値 |
|------|-----|
| 13C制約後の成長率 | **0.8739 h⁻¹**（変化なし） |
| 制約反応数 | 9反応 |
| χ²値 | 6.42 (dof=9) |
| Reduced χ² | **0.713**（良好な適合） |
| G6P→解糖系分配比 | 0.601 |
| G6P→PPP分配比 | 0.399 |
| PEP→PYR分配比 | 0.413 |
| PEP→OAA分配比 | 0.588 |

- Reduced χ² < 1 は、合成データのノイズレベルに対してモデルが良好に適合していることを示す
- フラックス分配比は13C制約前後で安定（コアモデルの解空間が狭い）

### 3.3 動的FBA

| 指標 | 値 |
|------|-----|
| 最大バイオマス | **1.078 gDW/L** |
| グルコース枯渇時間 | **6.2 h** |
| 最大アセテート蓄積 | **9.28 mM** |
| 初期グルコース40mM時の最大バイオマス | 2.093 gDW/L |

- バッチ培養シミュレーションでは、6.2時間でグルコース枯渇 → ジオーキシック移行を再現
- 初期グルコース濃度に対してバイオマスは概ね線形に増加

### 3.4 酵素容量制約

| 酵素予算 | 成長率 (h⁻¹) | 成長率低下 |
|----------|-------------|-----------|
| 0.2x | 0.2945 | **66.3%** |
| 0.4x | 0.6629 | **24.1%** |
| 0.6x | 0.7598 | **13.1%** |
| 0.8x | 0.8178 | **6.4%** |
| 1.0x | 0.8739 | 0.0% |
| 1.5x | 0.8739 | 0.0% |

**酵素コストTop 3（1.0x予算時）**:
1. FBA (F-bisP aldolase): 19.06 mg/gDW — kcat=17 s⁻¹（最も遅い酵素）
2. ACONTa (Aconitase): 17.52 mg/gDW — kcat=18 s⁻¹
3. AKGDH (α-KG DH): 11.72 mg/gDW — kcat=48 s⁻¹, MW=400 kDa

- 低kcatかつ高分子量の酵素がボトルネック（FBA, ACONTa/b, AKGDH）
- 酵素予算20%では成長率が66%低下し、プロテオーム制約の重要性を示した

### 3.5 条件特異的モデル

| 条件 | 成長率 (h⁻¹) | 活性反応数 |
|------|-------------|-----------|
| 好気 | **0.8739** | 48/95 |
| 嫌気 | **0.2117** | 47/95 |
| ストレス | **0.8739** | 48/95 |

- **嫌気特異的反応（12反応）**: ACKr, PTAr, PFL, ACALD, ALCD2x, ETOHt2r, ACt2r, FORti, THD2 等
  → 混合酸発酵経路（酢酸、エタノール、ギ酸の生成経路）が活性化
- 好気—嫌気間のJaccard類似度: **0.583**（活性反応パターンに明確な差異）
- 好気—ストレス間のJaccard類似度: **1.000**（コアモデルではストレス応答の差異が限定的）

### 3.6 リシン生産最適化

| 指標 | 値 |
|------|-----|
| 野生型成長率 | 0.8739 h⁻¹ |
| 最大理論リシン生産量 | **7.994 mmol/gDW/h** |
| 成長共役リシン生産（≥10%成長時） | **7.209 mmol/gDW/h** |
| 理論炭素収率 | **0.799 mol_C_lys/mol_C_glc** |
| 50%成長時の最大リシン | 3.87 mmol/gDW/h |

- **生産エンベロープ**: 成長率とリシン生産は強い負の相関。成長率0でリシンは最大7.99
- **OptKnock**: コアモデルの小規模（137遺伝子）のため、成長必須条件と両立する有効なKO戦略は検出されず
  → ゲノムスケールモデル（iML1515等）では有効戦略が期待される
- **過剰発現ターゲット**: コアモデルでは上限拡大の効果がなく（既にフラックス制約がボトルネックでない）、全酵素で生産量変化0.0

---

## 4. 考察と今後の展望

### 4.1 フレームワークの有用性

本フレームワークは、FBAの6つの拡張手法を統合パイプラインとして実装し、COBRApy/Cameoベースの再現可能な解析環境を提供する。各モジュールは独立に実行可能であり、研究目的に応じた柔軟な組み合わせが可能である。

### 4.2 コアモデルの制限と全ゲノムモデルへの展開

- **コアモデル**（95反応, 72代謝物, 137遺伝子）は教育・デモ用途には適しているが、代謝工学の実用的予測にはゲノムスケールモデル（例：iML1515, 2719反応）が必要
- OptKnockやGIMMEの有効性は、モデルの規模と解空間の広さに依存する
- 酵素容量制約はGECKO（ecYeast8等）やsMOMENTの本格的実装で、より生物学的に妥当な予測が期待される

### 4.3 今後の展望

1. **ゲノムスケールモデル（iML1515）への適用**: 全6モジュールをフルスケールで実行
2. **13C-MFA**: 実測データ（GC-MS / LC-MS MDV）の直接入力インターフェース構築
3. **dFBA高度化**: RK4積分法、代謝遺伝子制御（動的酵素発現）の導入
4. **GECKO完全実装**: ECモデル自動構築ツール（geckopy）との連携
5. **マルチオミクス統合**: RNA-seq + プロテオミクス + メタボロミクスの同時制約
6. **メタ学習**: 複数条件のモデル解を統合したアンサンブルフラックス推定
7. **産業応用**: リシンに加え、コハク酸、1,4-BDO、イソプレンなど高付加価値化学品への拡張

### 4.4 制限事項

- 合成データ（13C-MFA, RNA-seq）を使用しており、実験データでの検証が必要
- コアモデルの小規模さにより、OptKnock等の組合せ最適化の効果が限定的
- dFBAのオイラー積分は数値安定性に課題あり（特にグルコース枯渇付近）
- 酵素kcat値はBRENDA由来の代表値であり、条件依存的な変動は未考慮

---

## 5. 生成ファイル一覧

### スクリプト

| ファイル | 内容 |
|---------|------|
| `scripts/s01_fba_optimization.py` | FBA/pFBA/FVA/感度解析 |
| `scripts/s02_13c_mfa_integration.py` | 13C-MFA統合フレームワーク |
| `scripts/s03_dynamic_fba.py` | 動的FBA (SOA) |
| `scripts/s04_enzyme_constraints.py` | 酵素容量制約 (sMOMENT) |
| `scripts/s05_condition_specific.py` | 条件特異的モデル (GIMME) |
| `scripts/s06_lysine_optimization.py` | リシン生産最適化 |
| `scripts/run_all.py` | 全モジュール一括実行 |

### 図表 (`figures/`)

| ファイル | 内容 |
|---------|------|
| `01_fba_optimization.png/svg` | グルコース/O₂感度、FVA、pFBA比較 |
| `02_13c_mfa_integration.png/svg` | FBA vs 13C-MFA、分配比、残差解析 |
| `03_dynamic_fba.png/svg` | バイオマス/グルコース動態、初期条件比較 |
| `04_enzyme_constraints.png/svg` | 酵素予算vs成長率、酵素配分、kcat分布 |
| `05_condition_specific.png/svg` | 条件別成長率、活性反応、Jaccard類似度 |
| `06_lysine_optimization.png/svg` | 生産エンベロープ、KO戦略、フラックスマップ |

### 結果データ (`results/`)

| ファイル | 内容 |
|---------|------|
| `01_fba_optimization.json` | FBA/pFBA/FVA数値結果 |
| `01_fva_results.csv` | 全反応FVA結果テーブル |
| `02_13c_measured_data.csv` | 13C合成測定データ |
| `02_13c_mfa_results.json` | 13C統合解析結果 |
| `03_dfba_trajectory.csv` | dFBA時系列データ |
| `03_dfba_results.json` | dFBA要約統計量 |
| `04_enzyme_constraints.json` | 酵素制約解析結果 |
| `05_condition_specific.json` | 条件特異的モデル結果 |
| `06_lysine_optimization.json` | リシン最適化結果 |

### データ (`data/`)

| ファイル | 内容 |
|---------|------|
| `05_rnaseq_aerobic.csv` | 好気条件合成RNA-seqデータ |
| `05_rnaseq_anaerobic.csv` | 嫌気条件合成RNA-seqデータ |
| `05_rnaseq_stress.csv` | ストレス条件合成RNA-seqデータ |

### ログ (`logs/`)

| ファイル | 内容 |
|---------|------|
| `process-log.jsonl` | 実行トレース（タイムスタンプ、フェーズ、生成ファイル） |

---

## 参考文献

1. Orth JD et al. (2010) "What is flux balance analysis?" *Nat Biotechnol* 28:245-248
2. Lewis NE et al. (2012) "Constraining the metabolic genotype-phenotype relationship using a phylogeny of in silico methods" *Nat Rev Microbiol* 10:291-305
3. Sánchez BJ et al. (2017) "Improving the phenotype predictions of a yeast genome-scale metabolic model by incorporating enzymatic constraints" *Mol Syst Biol* 13:935 (GECKO)
4. Becker SA & Palsson BØ (2008) "Context-specific metabolic networks are consistent with experiments" *PLoS Comput Biol* 4:e1000082 (GIMME)
5. Mahadevan R & Schilling CH (2003) "The effects of alternate optimal solutions in constraint-based genome-scale metabolic models" *Metab Eng* 5:264-276 (dFBA)
6. Burgard AP et al. (2003) "OptKnock: a bilevel programming framework for identifying gene knockout strategies" *Biotechnol Bioeng* 84:647-657
7. Wiechert W (2001) "13C metabolic flux analysis" *Metab Eng* 3:195-206
