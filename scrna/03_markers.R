# 03: クラスタマーカー同定 + CONFIG のマーカーでアノテーション目安。
# 実行: Rscript 03_markers.R （クラスタ数が多いと時間がかかる場合は job-run 推奨）
suppressPackageStartupMessages({ library(Seurat); library(dplyr) })
source("00_config.R")

o <- readRDS("data/processed/seurat_integrated.rds")
o <- JoinLayers(o)          # FindAllMarkers の前に layer を結合（Seurat v5）
Idents(o) <- "seurat_clusters"

mk <- FindAllMarkers(o, only.pos = TRUE, min.pct = 0.25,
                     logfc.threshold = 0.25, verbose = FALSE)
dir.create("results", showWarnings = FALSE)
write.csv(mk, "results/markers_all.csv", row.names = FALSE)
top <- mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10) %>% ungroup()
write.csv(top, "results/markers_top10.csv", row.names = FALSE)

avail <- intersect(names(CONFIG$marker_hint), rownames(o))
cat("MARKERS_DONE clusters:", length(unique(mk$cluster)),
    "| hint markers found:", length(avail), "/", length(CONFIG$marker_hint), "\n")
