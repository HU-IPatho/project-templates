# 02: SE → DEG。config$method（deseq2 / edger / both）で選択。両走時は一致度も出力。
# 実行: プロジェクトルートから  Rscript analysis/02_de.R
suppressPackageStartupMessages({
  library(here); library(yaml); library(SummarizedExperiment)
})
source(here::here("R", "helpers.R"))       # save_table / %||%
source(here::here("R", "de_common.R"))     # tidy_de / compare_methods
source(here::here("R", "design.R"))        # edger_design / edger_contrasts / deseq2_formula
source(here::here("R", "housekeeping.R"))  # load_housekeeping
source(here::here("R", "de_deseq2.R"))     # run_deseq2
source(here::here("R", "de_edger.R"))      # run_edger

cfg <- yaml::read_yaml(here::here("config.yaml"))
se  <- readRDS(here::here("data", "processed", "se.rds"))
method <- cfg$method %||% "edger"
fdr    <- cfg$fdr %||% 0.05

de_list <- list()
if (method %in% c("deseq2", "both")) {
  de_list$DESeq2 <- run_deseq2(se, cfg)
}
if (method %in% c("edger", "both")) {
  hk <- load_housekeeping(cfg)
  de_list$edgeR <- run_edger(se, cfg, hk_genes = hk)
}
if (length(de_list) == 0) stop("02_de: 未知の method: ", method, "（deseq2 / edger / both）")

# DEG 表を統一列で保存（エンジン×対比ごと）
tidy_all <- list()
for (eng in names(de_list)) {
  tid <- tidy_de(de_list[[eng]])
  tidy_all[[eng]] <- tid
  for (nm in names(tid)) {
    save_table(tid[[nm]],
               tbl_id = sprintf("deg_%s_%s", tolower(eng), nm),
               desc   = sprintf("DEG %s %s", eng, nm),
               script = "analysis/02_de.R")
  }
}

# 両走時は DESeq2 vs edgeR 一致度
if (method == "both") {
  agree <- compare_methods(tidy_all$DESeq2, tidy_all$edgeR, fdr = fdr)
  save_table(agree, "method_agreement", "DESeq2 vs edgeR agreement", "analysis/02_de.R")
  cat("AGREEMENT:\n"); print(agree)
}

saveRDS(de_list, here::here("data", "processed", "de.rds"))
cat("DE_DONE:", paste(names(de_list), collapse = ","),
    "| method =", method, "\n")
