---
name: scaffold-config
description: >-
  scRNA テンプレの層3（config.yaml の事実 + AGENTS.md「このプロジェクト固有」節）を、着手時に
  /grill-me で対話的に埋める薄い scaffold。loader 種別・mito pattern・batch key・解像度・marker
  hint 等の scRNA 固有論点を順に確定する。共有データ由来なら README.yaml から下書きする。
  Use when: 新規 scRNA プロジェクトを開始し config.yaml をまだ埋めていないとき、
  says「解析を始める」「config を埋める」「セットアップ」「scaffold」。
---

# 層3 scaffold — scRNA プロジェクト固有設定を grill-me で埋める

このテンプレの層3（プロジェクト固有）は 2 つに分かれる（① 指示 3 層モデル）:
- **`config.yaml`** … 機械可読な「事実」（loader / QC 閾値 / batch_key / 次元 / marker_hint 等）
- **`AGENTS.md`「このプロジェクト固有」節** … 非自明な narrative 規約（この検体特有の注意等）

研究員が手書きせずに済むよう、着手時にこの scaffold が `/grill-me` を **scRNA 固有の論点**で走らせ、
対話で両者を確定する。

## 手順

1. **共有データ由来なら下書きを作る**: 入力が `/data/shared/datasets/<dataset>/` 由来なら、その
   `README.yaml`（provenance マニフェスト）を読み、organism（→ `mito_pattern`）・サンプル数・
   modality から `config.yaml` の下書きを作る（手入力を減らし来歴を継承する）。

2. **`/grill-me` を下記の論点で起動**し、1 つずつ確定する（各問いは平易に説明してから問う）:
   - **loader**: データ形式は 10x Cell Ranger 出力（`tenx`）か、per-sample の密行列 gz（`dense`・GEO
     processed 等）か。`dense` なら `dense_pattern` / `dense_sample_regex` をファイル名に合わせる。
   - **mito_pattern**: 生物種。ヒト=`^MT-` / マウス=`^mt-`。**種を間違えると percent.mt が全部 0** になる。
   - **batch_key**: 統合（Harmony）でバッチ効果を除く軸。通常 `sample`。**癌種などの生物学的変数を
     batch key にしない**（本物の差まで消える）。
   - **qc**: `nFeature_min/max`・`percent_mt_max`。まず violin plot で分布を見てから決めるのが定石。
   - **resolution / resolutions / n_pcs / dims / n_variable**: primary 解像度と走査解像度群・次元。既定は汎用値。
     標準は単一固定でなく複数解像度スキャンを求める（`resolutions`）。
   - **marker_hint**: アノテーションの当たり付けに使う「遺伝子 → 細胞種」。組織に合わせて増減する。
   - **sample_class**（アノテ標準の前段軸）: `mixed`（腫瘍+正常/免疫/間質の混在組織・v0 の確定操作 scope）か
     `pure_cell_line`（純細胞株・CNV ゲート省略/composition 無効/反復単位=culture replicate）か。方法論の適用が分岐する。
   - **doublet**: 除去の有無・手法を**記録**する（`method: none`＋`rationale`、または `scDblFinder`）。除去可否・手法・
     閾値の科学判断はここで grill する（形の記録は必須）。
   - **cell_cycle**: `score`（S/G2M を直交軸で保持・既定 true）/ `regress`（無条件除去はしない・既定 false）。
   - **cnv_gate**（悪性判定）: 対象が aneuploidy を持ちうる固形上皮 carcinoma なら `enabled: true` + `tool` + 非悪性
     `normal_reference` を検討（near-diploid 腫瘍種・純細胞株は非適用）。ゲートは非対称（CNV 陰性→保留）。
   - **reference**（two-track）: 適合参照アトラスが在れば `enabled: true` + `method`（SingleR/Azimuth 等）で必須
     complementary クロスチェックにする（無条件 primary にはしない）。

3. **確定値を書き込む**: `config.yaml` を更新し、非自明な固有規約があれば `AGENTS.md` の
   「このプロジェクト固有」節に短く追記する（共通・自明は書かず上位層＝AGENTS-base / 本 AGENTS.md
   に委ねる＝重複させない）。

4. **確認**: `Rscript -e 'source("analysis/00_config.R"); str(CONFIG)'` で config.yaml が読めることを確かめる。
