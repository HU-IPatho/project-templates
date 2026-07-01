# 03: クラスタマーカー同定 + CONFIG のマーカーでアノテーション目安。
# 実行: プロジェクトルートから  Rscript analysis/03_markers.R
# 入力: data/processed/seurat_integrated.rds（02 の出力・正準 object）
# 出力: outputs/tables/tbl01_markers_all.csv・tbl02_markers_top10.csv（ハーネス経由）
suppressPackageStartupMessages({ library(here); library(Seurat); library(dplyr) })
source(here::here("analysis", "00_config.R"))
source(here::here("R", "helpers.R"))

o <- readRDS(here::here("data", "processed", "seurat_integrated.rds"))
o <- JoinLayers(o)          # FindAllMarkers の前に layer を結合（Seurat v5）
Idents(o) <- "seurat_clusters"

mk <- FindAllMarkers(o, only.pos = TRUE, min.pct = 0.25,
                     logfc.threshold = 0.25, verbose = FALSE)
save_table(mk, "tbl01", "markers_all", "03_markers.R")
top <- mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10) %>% ungroup()
save_table(top, "tbl02", "markers_top10", "03_markers.R")

avail <- intersect(names(CONFIG$marker_hint), rownames(o))
cat("MARKERS_DONE clusters:", length(unique(mk$cluster)),
    "| hint markers found:", length(avail), "/", length(CONFIG$marker_hint), "\n")
