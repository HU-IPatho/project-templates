---
name: scaffold-rnaseq-dtu
description: >-
  transcript-level DTE/DTU 解析テンプレ（rnaseq-dtu）の着手時 scaffold。/grill-me を
  テンプレ固有の論点で走らせ、config.yaml の事実（reference/design/contrast/dte/dtu/samples）と
  AGENTS.md「このプロジェクト固有」節（非自明な narrative）を対話で埋める。
  入力が共有領域（/data/shared）由来なら README.yaml から config を下書きする。
  Use when: このテンプレで新規プロジェクトを始める / config を埋める / DTU 解析設計を確定する。
---

# scaffold: rnaseq-dtu の層 3 を埋める

このテンプレ（analysis-template-standard 準拠）の**層 3＝プロジェクト固有**を確定する薄い
scaffold。層 1（`AGENTS-base.md`）・層 2（本テンプレ `AGENTS.md` 上部）は固定物。ここで埋めるのは
**事実（`config.yaml`）**と**非自明な固有規約（`AGENTS.md`「このプロジェクト固有」節）**のみ。

## 手順

1. **共有由来データの下書き（あれば）**: 入力が `/data/shared/datasets/<dataset>/` 由来なら、
   その `README.yaml`（provenance マニフェスト・`specs/shared-data-curation`）を読み、
   `project.organism` / `reference`（release）/ サンプル一覧（`metadata/samples.tsv` の group 割当）を
   `config.yaml` の下書きに反映する（手入力を減らし来歴を継承する）。

2. **`/grill-me` を下記の DTU 固有論点で起動**し、ユーザーと対話して各分岐を 1 つずつ解消する。
   各質問の前に論点を平易に説明してから問う（grill-me の作法）。

   - **問いの種別（最重要）**: 見たいのは **DTE（transcript の発現量差）** か **DTU（isoform 使用比の
     変化＝ switching）** か、両方か。DTU が主眼なら DRIMSeq+stageR（既定）で回す。DTE の色分け FDR は `dte.fdr`。
   - **複製（必須ゲート）**: **各群のサンプル数**。DTU/DTE は各群 n>=2 が必須。**各群 n=1 なら本テンプレは
     使えない**（01/05 が fail-fast）→ 量差だけなら rnaseq-secondary（② edgeR n=1 HK）へ誘導する。
   - **実験デザイン**: 群（`group_col` の水準）・対照（`reference_level`）・共変量（`covariates`・バッチ等・
     無ければ `[]`）・対比（`contrasts` を `名前: [num, den]` で列挙）。
   - **参照（reference.*）**: decoy-aware salmon index / `tx2gene`（**DRIMSeq の gene グルーピングに必須**）/
     `gtf`。index と tx2gene の由来を一致させる（不一致だと gene 割当が壊れる）。
   - **DTU 手法と閾値（dtu.*）**: `methods`（既定 `[drimseq]`。上級で `swish`/`dexseq` を追加）・`fdr`（OFDR）・
     `n_top`（主役図の遺伝子数）・dmFilter 閾値（`min_gene_expr`/`min_feature_expr`/`min_feature_prop`）。
   - **swish を使うか**: 使うなら **`salmon.num_bootstraps` を >0** に設定する必要（inferential replicates が必須。
     未取得だと `run_swish` が fail-fast）。使わないなら 0 のまま（軽い salmon）。
   - **非自明な固有規約**: このデータ固有の落とし穴（特殊なアノテーション整合・除外サンプル・
     参照の非標準運用など）があれば `AGENTS.md`「このプロジェクト固有」節に narrative で残す。

3. **書き込み**: 合意した事実を `config.yaml` に、非自明な規約のみを `AGENTS.md`
   「このプロジェクト固有」節に反映する。共通・自明な内容は書かず上位層を指すに留める
   （重複＝drift 防止）。事実を `AGENTS.md` に、規約を `config.yaml` に混ぜないこと。

4. **確認**: `python3 analysis/01_fetch.py --dry-run`（複製検査が通る）と、salmon quant 後に
   `Rscript analysis/04_import_tx.R` が `se_tx.rds` を作れる最小状態（reference と samples の整合）まで導く。
   長時間の salmon 定量は必ず `job-run` で回す（AGENTS.md / 03 スクリプト参照）。
