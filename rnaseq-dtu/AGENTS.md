# transcript-level DTE/DTU 解析プロジェクト規約（Salmon → tximport → DRIMSeq+stageR）

このテンプレートは、既知 isoform の **DTE（発現量差）/ DTU（使用比変化）** を FASTQ から
**再現性高く** 回すための雛形。② rnaseq-secondary（gene-level DE）の姉妹＝**transcript-level 2 次解析**。
**まず本ファイルを最後まで読むこと。** 環境・共通作法は上位層（`AGENTS-base.md`）を参照し、
ここには本テンプレ固有の非自明な規約のみを置く（重複記載しない＝drift 防止）。

## この解析の骨格（analysis-template-standard 準拠）

- **正準 object = transcript-level SummarizedExperiment**（`data/processed/se_tx.rds`）。
  `analysis/04_import_tx.R` が保存する。② の gene-level `se.rds` とは別 object（名前で区別）。
  検定結果は `data/processed/dtu.rds`（DTE + DTU をまとめた list）。
- **パイプラインは番号順**: `01_fetch` → `02_qc_fastp` → `03_salmon_quant` → `04_import_tx` →
  `05_dtu` → `06_figures`。各スクリプトはルートから実行し、パスは `here::here()` でルート相対解決。
- **図表は必ず出力ハーネス経由**（`R/helpers.R`）: `save_fig`（PNG+PDF・固定 DPI・統一テーマ）/
  `save_table`（`outputs/tables` へ TSV）/ `outputs/captions.tsv`（図表→説明→由来 script）。
- **レポートは再計算しない**: `reports/report.qmd` は `se_tx.rds`/`dtu.rds` と `outputs/` を読むだけ。

## 編集するのは `config.yaml` だけ（事実層）

スクリプト本体（`analysis/`・`R/`・`python/`）は触らない。着手時に **`/scaffold-rnaseq-dtu`** を実行すると
`/grill-me` が固有論点を対話で引き出し、`config.yaml` と本ファイル「このプロジェクト固有」節を埋める。主なキー:

- **`reference.*`**: 管理者が `/work/shared` に構築済みの decoy-aware salmon index / `tx2gene` / `gtf`
  （③ rnaseq-primary と同形）。**`tx2gene` は DRIMSeq の gene グルーピングに必須**（index と整合させる）。
- **`design` / `contrasts`**: 主比較変数（`group_col`）・対照（`reference_level`）・多対比（`名前: [num, den]`）。
- **`dte`**: DTE エンジン（`deseq2` / `edger` / `both`）と FDR。
- **`dtu`**: DTU 手法（`methods`）・OFDR（`fdr`）・主役図の遺伝子数（`n_top`）・DRIMSeq dmFilter 閾値。
- **`salmon.num_bootstraps`**: `swish` を使うときだけ >0（後述）。既定 0＝軽い salmon。

## ★DTE と DTU は別の問い（取り違え厳禁）

- **DTE（differential transcript expression）** = 各 transcript の発現「量」が群間で違うか。
  ② の DESeq2/edgeR エンジンを **transcript-level counts** に適用する（`R/dte.R`）。出力=volcano・表。
- **DTU（differential transcript usage）** = 遺伝子内の isoform「使用比」が群間で変わるか
  （＝ isoform switching）。総発現量は同じでも比率が変われば DTU。**主役**。出力=isoform 使用比プロット・表。
- 例: 遺伝子 G の総発現は不変だが isoform A→B に切替 → **DTE では検出されにくいが DTU で検出**。

## ★DTU は rnaseqDTU（Love/Soneson/Patro F1000）に準拠

- **import**: `04_import_tx.R` が `tximport(txOut=TRUE)` を 2 通りの `countsFromAbundance` で呼び、
  se_tx.rds に 2 つの count assay を持たせる（gene 集約はしない・DTE/DTU とも transcript 単位で検定）。
  - assay `counts` = **`dtuScaledTPM`**（DTU 用）: TPM を「遺伝子内 isoform 間で比較可能な尺度」に
    整えた擬似 count で、DRIMSeq/DEXSeq が期待する isoform 使用比向け入力。
  - assay `counts_dte` = **`lengthScaledTPM`**（DTE 用）: 各転写産物を自身の長さで scale した count。
    DTE（発現「量」差）の妥当な count 表現。**dtuScaledTPM を DTE に流用すると p 値/FDR が
    miscalibrated になる**ため、DTE は必ず `counts_dte` を使う（`R/dte.R`）。
- **DRIMSeq**: `dmFilter → dmPrecision → dmFit → dmTest`。dmFilter の発現閾値は `config.dtu`
  （`min_gene_expr`/`min_feature_expr`/`min_feature_prop`）。標本数依存の `min_samps_*` は
  実行時に n（総標本）/ n.small（小さい群）から自動決定する。
- **stageR 二段（OFDR 制御）**: screening=遺伝子（DRIMSeq gene-level p）→ confirmation=transcript、
  `method="dtu"`, `alpha=config.dtu.fdr`。screening を通った遺伝子内でのみ transcript を確認する。
  結果表の `gene_padj`=screening OFDR、`tx_padj`=confirmation OFDR。

## ★複製は必須（各群 n>=2）

DTU/DTE は各群に複製が要る。`01_fetch.py`（Python）と `05_dtu.R` 冒頭（`check_replicates`）で
群別複製数を検査し、**n=1 の群があれば fail-fast で停止**する。量の差だけを見たい n=1 実験は
**rnaseq-secondary（② edgeR n=1 HK 分散）** を使う（本テンプレは誘導メッセージで案内する）。

## 任意手法（上級・既定では走らない）

`config.dtu.methods` に足すと有効化（既定は `[drimseq]`）:

- **`swish`（fishpond）**: inferential replicates を使う不確実性込みの **DTU** 検定。fishpond の DTU
  レシピに忠実（`scaleInfReps → labelKeep → keep フィルタ → isoformProportions → swish(x=condition)`）で、
  `isoformProportions` により「isoform 使用比」を検定対象にする（発現量そのものではない）。
  **salmon bootstrap が必須**: `config.salmon.num_bootstraps` を >0 にして `03_salmon_quant.sh` を
  再実行してから使う（未取得のまま `swish` を指定すると `run_swish` が fail-fast する）。
- **`dexseq`（DEXSeq）**: DRIMSeq と同じ dmFilter を共有し、`perGeneQValue`+stageR 二段で OFDR 制御。
- **IsoformSwitchAnalyzeR（ISA）**: switch 図・機能帰結（配列解析）まで見る上級ツール。**依存が重いため
  `install_deps.R` には同梱しない**。使う場合のみ各自で `renv::install("IsoformSwitchAnalyzeR")` し、
  `renv::snapshot()` で lock を更新する。

## 再現性（版ピン 3 系統）

- **R + Bioconductor** = `renv.lock`（**Bioconductor 3.23** を明示・DRIMSeq/DEXSeq/stageR/fishpond/
  tximport/DESeq2/edgeR ほか）。`install_deps.R` が `renv::restore()`、不完全なら fresh install + snapshot。
- **Python** = `pyproject.toml`（`uv lock && uv sync` で `uv.lock` を生成・追跡）。取得系の最小構成。
- **システムツール**（salmon/fastp/fastqc/multiqc/quarto）= omics-dev イメージのタグで固定（config でなくイメージが担体）。
- `renv::update()` を安易に実行しない（版が動くと結果が変わる）。更新したら `renv::snapshot()`。
  ※ 同梱の `renv.lock`/`renv/activate.R` は管理者の dogfood e2e（実 index+workspace）で確定する seed。

## データ操作（ワークスペース共通ヘルパー・詳細は AGENTS-base.md）

```bash
work-fetch /data/shared/datasets/<dataset>/   # HDD → data/raw（不変入力層）
work-save  data/processed outputs reports      # SSD → HDD（永続・必須）
work-share data/processed/se_tx.rds -m "説明"  # 正準 object を共有 _inbox へ昇格
```

- `data/` は git 非追跡（`raw`=不変入力 / `interim`=再生成可能 / `processed`=正準 object）。
- `outputs/`・`reports/` は git 追跡（PR/diff でレビュー）。生データ（`/data/shared`）は read-only。

---

## このプロジェクト固有

<!-- 層 3（narrative）: 着手時に /scaffold-rnaseq-dtu（/grill-me 駆動）が
     非自明な固有規約のみを追記する。事実（design/contrast/samples 等）は config.yaml へ。
     共通・自明な内容はここに書かず上位層（AGENTS-base.md・本ファイル上部）を指す。 -->

_（未設定。`/scaffold-rnaseq-dtu` を実行して埋める。）_
