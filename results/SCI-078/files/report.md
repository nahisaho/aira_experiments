# 実験レポート: 食事成分と腸内細菌叢の相互作用を予測するシステムバイオロジーフレームワーク

## 1. 実験目的と背景

本研究では、食事成分と腸内細菌叢（gut microbiota）の相互作用を包括的に予測するシステムバイオロジーフレームワークを設計・実装した。腸内細菌叢は宿主の健康状態に深く関与しており、食事パターンがその組成と代謝機能に与える影響の理解は、個別化栄養学やプロバイオティクス治療戦略の開発において極めて重要である。

本フレームワークは以下の6つのモジュールから構成される：

1. SHIME模擬に基づく食品成分の消化・吸収動態モデル
2. 一般化Lotka-Volterra（gLV）方程式による腸内細菌群集の資源競争モデル
3. 短鎖脂肪酸（SCFA）生成のフラックス予測
4. 食事パターンと菌叢組成の長期動態シミュレーション
5. プロバイオティクス/プレバイオティクスの効果予測
6. 発酵食品摂取の菌叢多様性への影響ケーススタディ
7. MICOM/gapseqベースのコミュニティ代謝モデリング

## 2. 使用した手法・アルゴリズムの概要

### 2.1 SHIME模擬消化モデル
Simulator of the Human Intestinal Microbial Ecosystem（SHIME）を参考に、胃・小腸・上行結腸・横行結腸・下行結腸の5コンパートメントモデルを構築した。各栄養素（デンプン、食物繊維、タンパク質、脂肪、単純糖質）の通過・消化・発酵過程を常微分方程式で記述した。

### 2.2 一般化Lotka-Volterra（gLV）モデル
8種の主要腸内細菌（*Bacteroides*, *Faecalibacterium*, *Bifidobacterium*, *Roseburia*, *Lactobacillus*, *Prevotella*, *Clostridium*, *Akkermansia*）の群集動態を、gLV方程式で記述した：

$$\frac{dx_i}{dt} = x_i \left( \mu_i + \sum_{j=1}^{N} A_{ij} x_j \right)$$

ここで、$x_i$は種$i$の存在量、$\mu_i$は固有増殖速度、$A_{ij}$は種間相互作用係数である。

### 2.3 SCFA フラックス予測
各細菌種のSCFA（酢酸・プロピオン酸・酪酸）生産収率係数行列を定義し、菌叢組成からSCFA生産量を予測した。

### 2.4 コミュニティ代謝モデリング（MICOM/gapseq準拠）
ゲノムスケール代謝モデルに基づく簡易コミュニティ代謝モデルを構築し、協調的トレードオフ解析を実施した。

## 3. 主要な結果と数値

### 3.1 SHIME消化動態

5つの栄養素の消化管通過・消化動態をシミュレーションした。食物繊維は胃・小腸でほぼ消化されず、結腸に到達して微生物発酵を受けることが確認された。

![SHIME消化動態](figures/shime_digestion.png)

### 3.2 gLV群集動態

8種の腸内細菌の動態シミュレーションにおいて、*Bacteroides*が最も高い存在量を示し、酪酸産生菌（*Faecalibacterium*, *Roseburia*）が交差栄養（cross-feeding）を通じて安定的に共存することが示された。

![gLV群集動態](figures/glv_dynamics.png)

種間相互作用マトリクスの可視化：

![種間相互作用マトリクス](figures/interaction_matrix.png)

### 3.3 SCFA生産動態

ベースライン条件下でのSCFA生産は定常状態に達し、酢酸が最大割合（約50%）を占めた。

![SCFA生産動態](figures/scfa_dynamics.png)

### 3.4 食事パターン別長期動態（90日間）

4つの食事パターンの比較結果：

| 食事パターン | 優勢種 | Shannon多様性 | 総SCFA (mmol/L) | 酪酸 (mmol/L) |
|---|---|---|---|---|
| 西洋食 | Bifidobacterium (22.1%) | 1.831 | 444.0 | 149.1 |
| 高食物繊維食 | Bifidobacterium (19.7%) | 1.903 | 547.3 | 179.3 |
| 地中海食 | Bifidobacterium (19.0%) | 1.926 | 531.0 | 169.0 |
| 低FODMAP食 | Bacteroides (20.9%) | 1.858 | 398.3 | 129.5 |

高食物繊維食は最大のSCFA生産量（+23.3% vs 西洋食）を示し、地中海食は最高のShannon多様性（1.926）を達成した。

![食事パターン別菌叢組成](figures/dietary_patterns.png)

![食事パターン別SCFA比較](figures/diet_scfa_comparison.png)

### 3.5 プロバイオティクス/プレバイオティクス効果

プレバイオティクス（イヌリン）介入は酪酸生産を24.6%増加させた（142.5 → 177.6 mmol/L）。シンバイオティクス（プロバイオティクス+プレバイオティクス）は同等の効果を示した。

![プロバイオティクス/プレバイオティクス効果](figures/probiotic_prebiotic.png)

![多様性の変化](figures/diversity_interventions.png)

### 3.6 発酵食品ケーススタディ

発酵食品（ヨーグルト/キムチ/ケフィア）摂取の3相試験（ベースライン30日→介入30日→追跡30日）では、介入期間中に*Lactobacillus*と*Bifidobacterium*の増加、SCFA生産の一時的上昇が観察された。

![発酵食品ケーススタディ](figures/fermented_food.png)

### 3.7 コミュニティ代謝モデリング

MICOM/gapseqベースの代謝フラックス解析により、高食物繊維食条件下で酪酸産生菌の交差栄養が最も活発であることが示された。

![コミュニティ代謝モデリング](figures/community_metabolic.png)

## 4. 考察と今後の展望

### 4.1 主要な知見
- 高食物繊維食および地中海食が腸内環境に最も好ましい影響を与える
- 食物繊維は結腸微生物発酵を通じてSCFA生産を促進し、特に酪酸生産が健康に重要
- プレバイオティクス介入はプロバイオティクス単独より効果的な場合がある
- 発酵食品は一時的な菌叢変化を誘導するが、長期的な定着には継続的摂取が必要

### 4.2 限界
- 本モデルは8種のみを対象としており、実際の腸内細菌叢（>1000種）の複雑さを完全には再現できない
- 宿主免疫応答、腸管上皮バリア機能、胆汁酸代謝などの要因は含まれていない
- gLVパラメータはin vitroデータに基づく推定値であり、in vivoとの差異がありうる

### 4.3 今後の展望
- ゲノムスケール代謝モデル（GEM）の完全統合（AGORA2データベース）
- 個人差を反映したパーソナライズドモデルの開発
- 機械学習（RNN/Transformer）による動態予測精度の向上
- 臨床試験データとの検証

## 5. 生成ファイル一覧

| ファイル | 説明 |
|---|---|
| `simulation.py` | メインシミュレーションスクリプト |
| `report.md` | 本レポート |
| `paper.md` | 学術論文形式文書 |
| `figures/shime_digestion.png` | SHIME消化動態図 |
| `figures/glv_dynamics.png` | gLV群集動態図 |
| `figures/interaction_matrix.png` | 種間相互作用マトリクス |
| `figures/scfa_dynamics.png` | SCFA生産動態図 |
| `figures/dietary_patterns.png` | 食事パターン別菌叢組成図 |
| `figures/diet_scfa_comparison.png` | 食事パターン別SCFA比較図 |
| `figures/probiotic_prebiotic.png` | プロバイオティクス/プレバイオティクス効果図 |
| `figures/diversity_interventions.png` | 多様性介入効果図 |
| `figures/fermented_food.png` | 発酵食品ケーススタディ図 |
| `figures/community_metabolic.png` | コミュニティ代謝モデリング図 |
