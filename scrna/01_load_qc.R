# 01: データ読込 → merge → QC。設定=00_config.R、ローダ=R/load_data.R。
# 実行: Rscript 01_load_qc.R  （数分以内なら前景で可。大量サンプルなら job-run 推奨）
suppressPackageStartupMessages({ library(Seurat); library(dplyr) })
source("00_config.R")
source("R/load_data.R")

objs <- load_samples(CONFIG)
if (length(objs) == 1) {
  merged <- objs[[1]]
} else {
  merged <- merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs))
}
merged <- JoinLayers(merged)   # Seurat v5: sample ごとの layer を 1 つに結合
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = CONFIG$mito_pattern)
cat(sprintf("pre-QC: %d cells\n", ncol(merged)))

q <- CONFIG$qc
merged <- subset(merged, subset = nFeature_RNA > q$nFeature_min &
                                  nFeature_RNA < q$nFeature_max &
                                  percent.mt   < q$percent_mt_max)
cat(sprintf("post-QC: %d cells, %d %s\n",
            ncol(merged), length(unique(merged@meta.data[[CONFIG$batch_key]])), CONFIG$batch_key))

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
saveRDS(merged, "data/processed/seurat_qc.rds")
cat("QC_DONE\n")
