# bulk RNA-seq 1 次解析プロジェクト規約（Salmon → tximeta → gene-level SE）

このテンプレートは bulk RNA-seq を FASTQ から **再現性高く** 定量し、来歴付きの gene-level
`SummarizedExperiment`（`data/processed/se.rds`）を作る雛形。**着手時にまず本ファイルを最後まで読むこと。**

> **指示 3 層モデル**（この規約の位置づけ）
> - **層 1 = workspace 共通**: `AGENTS-base.md`（全プロジェクト共通の SSoT）。ここには重複記載せず層 1 を参照する。
> - **層 2 = テンプレ汎用規約**: 本ファイルの以下の節（rnaseq-primary 共通の解析規約）。
> - **層 3 = プロジェクト固有**: 事実は `config.yaml`、narrative は本ファイル末尾の「このプロジェクト固有」節。
>   着手時は `/init-rnaseq-primary`（`.claude/skills/init-rnaseq-primary`）が `/grill-me` で層 3 を対話的に埋める。

## クイックスタート（この順で実行する。すべてプロジェクトルートから叩く）

```bash
# 1) 解析環境を復元（初回だけ・renv.lock で R+Bioconductor の版を固定）
Rscript install_deps.R                        # → INSTALL_DONE / tximeta / Bioconductor 3.23
uv sync                                        # → Python 取得系（pysradb/ffq）を uv.lock で固定

# 2) 層 3 を埋める（config.yaml と本 AGENTS.md 固有節）— 対話で確定
#    Claude: /init-rnaseq-primary   （手書きより grill-me 駆動を推奨）

# 3) データ取得（公開データなら 01・共有領域由来なら work-fetch）
python3 analysis/01_fetch.py                   # config.fetch.accessions を data/raw へ

# 4) QC トリム → FastQC
bash analysis/02_qc_fastp.sh                    # → QC_FASTP_DONE

# 5) salmon 定量（★長時間 → job-run。03 は内部で job-run/job-wait を使う）→ MultiQC 集約
bash analysis/03_salmon_quant.sh               # → SALMON_QUANT_DONE / reports/multiqc/

# 6) tximeta で来歴付 gene-level SE を生成（= 下流 rnaseq-secondary の入力）
Rscript analysis/04_tximeta_se.R               # → data/processed/se.rds / TXIMETA_SE_DONE

# 7) QC レポート（保存済み object を読むだけ）
quarto render reports/qc_report.qmd

# 8) 成果を HDD へ退避（/work は揮発的・セッション終了前に必須）
work-save data/processed/ outputs/ reports/
```

## 自分のデータに合わせる — 編集するのは `config.yaml` だけ

`analysis/`・`R/`・`python/` のスクリプト本体は触らない。`config.yaml` の値を変えるだけで動く。

- **`reference.*`**: 管理者が `/work/shared` に構築済みの **decoy-aware salmon index** / `tx2gene` / `gtf` のパス。
  `txome`（tximeta 来歴用）は **index を作った参照 FASTA/GTF と一致**させる（不一致だと来歴・gene 集約が壊れる）。
- **`samples`**: 各 `id` と `data/raw` 配下の `fq1`/`fq2`。`project.paired_end` でペア/シングルを切替。
- **`salmon.extra_flags`**: selective alignment の補正。**`--gcBias` は必須**（GC バイアス補正）。

## 長時間の salmon 定量は job-run で回す（必須）

salmon 定量はサンプルあたり長時間になる。前景直実行は (a) ツールがタイムアウトし (b) 通信断でジョブごと死ぬ。
`analysis/03_salmon_quant.sh` は各サンプルを `job-run --exclusive salmon-quant` で detach 実行し `job-wait` で待つ
（同一リソースを直列化）。`job-wait` の終了コード: `0`=成功 / `124`=未完（もう一度待つ）/ その他=失敗。

## パイプラインの中身と順序（ベストプラクティス）

1. **01_fetch.py**: 公開 FASTQ を `data/raw`（不変入力層）へ。共有由来は `work-fetch` で `data/raw` へ。
2. **02_qc_fastp.sh**: fastp トリム（`data/raw`→`data/interim/trimmed`）→ FastQC。fastp/FastQC ログは `data/interim`。
3. **03_salmon_quant.sh**: **salmon selective alignment（decoy-aware index・`-l A`・`--gcBias`）** で定量
   → 末尾で **MultiQC が fastp/FastQC/salmon を集約し `reports/multiqc/` に HTML**。
4. **04_tximeta_se.R**: **tximeta（linkedTxome で来歴担保）→ `summarizeToGene` で gene-level SE** →
   `data/processed/se.rds`。あわせて QC（library size/mapping rate/PCA）を計算し出力ハーネスで図表保存。

順序の勘所（崩すと壊れる）:
- **decoy-aware index が前提**。decoy 無し index では selective alignment のマッピング特異性が落ちる。
- **`txome` は index と一致**。tximeta は来歴ハッシュで参照を照合するため、config の `txome` が実 index とずれると
  gene 集約（`summarizeToGene`）に必要な tx2gene が引けず失敗する（その場合は tx2gene での tximport にフォールバック）。
- **重い計算は `analysis/` が済ませ `data/processed` に保存**。Quarto レポートは読むだけ（再計算しない）。

## 出力の作法（出力ハーネス経由）

図表は生の `ggsave`/`write.csv` を直呼びせず **`R/helpers.R` の `save_fig` / `save_table`** を使う。
`save_fig` は 1 呼出しで **PNG+PDF を固定 DPI・統一テーマ**で `outputs/figures/figNN_<desc>` に出し、
`outputs/captions.tsv`（図表番号→説明→由来スクリプト）に索引する。表は `outputs/tables/tabNN_<desc>.tsv`。

## ディレクトリ規約

| パス | 用途 | git |
|---|---|---|
| `data/raw/` | 生 FASTQ（不変入力層。`01_fetch.py`/`work-fetch` で配置） | 非追跡 |
| `data/interim/` | トリム済み FASTQ・FastQC・salmon quant（再生成可能） | 非追跡 |
| `data/processed/` | 正準 object（`se.rds`・`qc.rds`）＝共有昇格の単位 | 非追跡 |
| `analysis/` | 番号付き実行スクリプト（01-04・R/Python/shell 混在） | 追跡 |
| `R/` `python/` | 再利用ロジック（`helpers.R`/`qc_se.R`・`config.py`/`fetch.py`） | 追跡 |
| `outputs/{tables,figures}` | 表・図（図は PNG+PDF）・`captions.tsv` | 追跡 |
| `reports/` | Quarto レポート・MultiQC 集約 HTML | 追跡 |

## 共有領域との連携（昇格と取り込み・詳細は層 1）

- **取り込み**: `work-fetch /data/shared/datasets/<dataset>/` → `data/raw`（不変入力層）。
- **昇格**: `data/processed` の正準 object を `work-share` で `_inbox` へ提出 → 管理者キュレーション → 共有 `processed/`。
- 共有 `curated` データは read-only（直接書き換えない）。詳細な helper 仕様は層 1（`AGENTS-base.md`）を参照。

## 再現性（版ピン 3 系統）

- **R + Bioconductor** = `renv.lock`（**Bioconductor 3.23** を明示）。`install_deps.R` が `renv::restore()` で復元。
- **Python** = `uv.lock`（`pyproject.toml` の依存を `uv lock`/`uv sync` で固定）。
- **システムツール**（salmon/fastp/fastqc/multiqc/quarto）= omics-dev イメージのタグで固定（config でなくイメージが担体）。
- `renv::update()` を安易に走らせない（版が動くと結果が変わる）。更新時は `renv::snapshot()` して commit。
  ※ 同梱の `renv.lock` は Bioc 3.23 の seed。管理者の dogfood e2e（実 index+workspace）で完全な版に確定する。
  ※ 同梱の `uv.lock` も手書き seed（host はネットワーク不可で `uv lock` を回せないため最小 seed を固定）。
    管理者の dogfood e2e（omics-dev イメージ内で `uv lock`）で hash 付き完全 lock に確定する。

---

## このプロジェクト固有（層 3・narrative）

<!--
ここには **非自明なプロジェクト固有規約のみ** を書く（workspace 共通・テンプレ汎用・自明は書かない）。
着手時に /init-rnaseq-primary（/grill-me）が対話で埋める。例:
- 非標準の参照を使う理由（例: 独自アセンブリ・特定 GENCODE release 固定）
- 実験デザイン上の注意（バッチ・条件・technical replicate の扱い）
- libtype を A 以外に固定する根拠
-->

（未確定。`/init-rnaseq-primary` を実行して埋める）
