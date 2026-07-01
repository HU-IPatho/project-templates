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
   - **resolution / n_pcs / dims / n_variable**: クラスタ解像度と次元。既定は汎用値。
   - **marker_hint**: アノテーションの当たり付けに使う「遺伝子 → 細胞種」。組織に合わせて増減する。

3. **確定値を書き込む**: `config.yaml` を更新し、非自明な固有規約があれば `AGENTS.md` の
   「このプロジェクト固有」節に短く追記する（共通・自明は書かず上位層＝AGENTS-base / 本 AGENTS.md
   に委ねる＝重複させない）。

4. **確認**: `Rscript -e 'source("analysis/00_config.R"); str(CONFIG)'` で config.yaml が読めることを確かめる。
