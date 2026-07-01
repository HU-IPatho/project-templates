# 01: counts 行列 + metadata(config) → SummarizedExperiment → data/processed/se.rds
# 実行: プロジェクトルートから  Rscript analysis/01_build_se.R
# パスは here::here() でルート相対解決（実行位置非依存・spec 準拠）。
suppressPackageStartupMessages({
  library(here); library(yaml); library(SummarizedExperiment)
})
source(here::here("R", "helpers.R"))     # %||%
source(here::here("R", "build_se.R"))

cfg <- yaml::read_yaml(here::here("config.yaml"))

# config$samples（list of {sample_id, group, batch, ...}）を data.frame 化
samples <- do.call(rbind, lapply(cfg$samples, function(s) {
  as.data.frame(s, stringsAsFactors = FALSE)
}))
rownames(samples) <- samples$sample_id

counts_file <- here::here(cfg$counts_file)
se <- build_se(counts_file, samples)

out <- here::here("data", "processed", "se.rds")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
saveRDS(se, out)
cat(sprintf("SE_DONE: %d genes x %d samples -> %s\n",
            nrow(se), ncol(se), out))
