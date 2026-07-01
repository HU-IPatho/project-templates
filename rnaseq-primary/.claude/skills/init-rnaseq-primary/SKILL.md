---
name: init-rnaseq-primary
description: >-
  rnaseq-primary（bulk RNA-seq 1 次解析: FASTQ→Salmon→tximeta→gene-level SE）テンプレートの
  層 3（config.yaml の事実 と AGENTS.md「このプロジェクト固有」節）を、着手時に /grill-me で
  対話的に埋める scaffold。研究員が config を手書きせず、対話で参照・サンプル・設計を確定する。
  Use when: プロジェクト着手直後 / 「このプロジェクトを初期化」「config を埋めて」「セットアップ」。
---

# rnaseq-primary 層 3 初期化（grill-me 駆動）

このテンプレートから始めた新規プロジェクトの層 3（プロジェクト固有の設定と規約）を確定する。
層 1（workspace 共通 `AGENTS-base.md`）と層 2（テンプレ汎用規約 = `AGENTS.md` 上部）は既に与えられている。
ここで埋めるのは **`config.yaml`（事実）** と **`AGENTS.md`「このプロジェクト固有」節（narrative）** だけ。

## 手順

1. **共有データ由来なら README.yaml から下書きする。**
   入力が `/data/shared/datasets/<dataset>/` 由来なら、その `README.yaml`（provenance マニフェスト）を読み、
   organism / 参照 release / サンプル一覧を `config.yaml` の下書きに反映する（手入力を減らし来歴を継承）。

2. **`/grill-me` を以下の rnaseq-primary 固有論点で起動し、`config.yaml` を確定する。**
   一問ずつ、平易な説明を添えて詰める:
   - **organism / paired_end**: ヒト/マウス、ペアエンドかシングルエンドか（`project.*`）。
   - **参照（reference.*）**: 管理者が `/work/shared` に構築済みの **decoy-aware salmon index** のパス、
     `tx2gene`、`gtf`。tximeta 来歴用 `txome`（source/organism/release/genome/fasta/gtf）は
     **index を作った参照と一致**させる（不一致だと来歴・集約が壊れる）。
   - **サンプル表（samples）**: 各 `id` と `data/raw` 配下の `fq1`/`fq2`。
     公開データを取るなら `fetch.accessions`（SRR/ENA）。
   - **定量オプション（salmon.*）**: libtype（既定 A=自動）、threads、`extra_flags`
     （**`--gcBias` は契約上必須**）。fastp の threads / extra_flags。

3. **`AGENTS.md`「このプロジェクト固有」節を確定する。**
   非自明な固有規約のみ（例: 非標準の参照を使う理由、特殊な libtype、実験デザイン上の注意）。
   workspace 共通・テンプレ汎用・自明な内容は書かず、上位層へのポインタに留める（drift 防止）。

4. 確定後の動作確認は README のクイックスタート（install_deps.R → 01→02→03→04）に従う。
   長時間の salmon 定量は必ず `job-run` で回す（AGENTS.md「長時間ジョブ」参照）。
