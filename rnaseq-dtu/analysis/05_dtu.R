# 05: transcript-level DTE + DTU 検定。
#   DTE = 転写産物の発現「量」差（② の DESeq2/edgeR エンジンを tx-level に適用・R/dte.R）。
#   DTU = 遺伝子内 isoform「使用比」の変化（DRIMSeq: dmFilter→dmPrecision→dmFit→dmTest →
#         stageR 二段確認で OFDR 制御・R/dtu.R）。config.dtu.methods で分岐。
#   任意 = swish（fishpond・bootstrap 必須）/ dexseq。
# 結果を data/processed/dtu.rds に保存し、表を出力ハーネスで outputs/tables に出す。
# 実行（ルートから）:  Rscript analysis/05_dtu.R
suppressPackageStartupMessages({
  library(here); library(yaml); library(SummarizedExperiment)
})
source(here::here("R", "helpers.R"))   # save_table / %||%
source(here::here("R", "design.R"))    # check_replicates / deseq2_formula / edger_design ...
source(here::here("R", "dte.R"))       # run_dte / tidy_dte
source(here::here("R", "dtu.R"))       # run_dtu(DRIMSeq+stageR) / run_swish / run_dexseq

cfg <- yaml::read_yaml(here::here("config.yaml"))
se  <- readRDS(here::here("data", "processed", "se_tx.rds"))
coldata <- as.data.frame(SummarizedExperiment::colData(se))

# --- 複製ガード（DTU/DTE は各群 n>=2 が必須。n=1 は fail-fast で ② へ誘導）---
check_replicates(coldata, cfg)

fdr     <- cfg$dtu$fdr %||% 0.05
methods <- cfg$dtu$methods %||% c("drimseq")

# --- DTE（transcript 発現差）---
dte      <- run_dte(se, cfg)
dte_tidy <- lapply(dte, tidy_dte)
# 選択した DTE エンジンごとに per-contrast 表を保存（dte.method: both なら DESeq2/edgeR 両方）。
# DTE 表は num=k 帯（下の DTU 表 num=10+k 帯と分離し、1 contrast での tabNN 衝突を避ける）。
for (eng in names(dte_tidy)) {
  tid <- dte_tidy[[eng]]
  for (k in seq_along(tid)) {
    save_table(tid[[k]], num = k,
               desc = sprintf("dte_%s_%s", tolower(eng), names(tid)[k]),
               script = "analysis/05_dtu.R")
  }
}

# --- DTU（isoform 使用比）: config.dtu.methods で分岐（複数エンジン並走可）---
dtu <- list()
if ("drimseq" %in% methods) dtu$drimseq <- run_dtu(se, cfg)
if ("swish"  %in% methods) dtu$swish  <- run_swish(cfg)      # bootstrap 必須（内部で fail-fast）
if ("dexseq" %in% methods) dtu$dexseq <- run_dexseq(se, cfg)

# 選択した各 DTU エンジンの per-contrast 表を engine 修飾 desc で保存（drimseq/dexseq は $table、
# swish は data.frame 直）。DTU 表は num=10+k 帯（DTE 表の num=k 帯と衝突しない）。
.dtu_table <- function(x) if (is.data.frame(x)) x else x$table
for (eng in names(dtu)) {
  per <- dtu[[eng]]
  for (k in seq_along(per)) {
    save_table(.dtu_table(per[[k]]), num = 10 + k,
               desc = sprintf("dtu_%s_%s", eng, names(per)[k]),
               script = "analysis/05_dtu.R")
  }
}

result <- list(
  dte = dte, dte_tidy = dte_tidy, dtu = dtu,
  meta = list(fdr = fdr, methods = methods,
              contrasts  = names(cfg$contrasts),
              dte_method = cfg$dte$method %||% "deseq2"))
saveRDS(result, here::here("data", "processed", "dtu.rds"))
cat("DTU_DONE: DTE(", paste(names(dte), collapse = ","), ") ",
    "DTU(", paste(names(dtu), collapse = ","), ")\n", sep = "")
