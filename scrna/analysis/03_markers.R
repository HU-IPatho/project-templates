# 03: クラスタマーカー同定（cluster-then-marker の marker 段）。アノテーション確定は 04_annotate.R。
# 実行: プロジェクトルートから  Rscript analysis/03_markers.R
# 入力: data/processed/seurat_integrated.rds（02 の出力・正準 object）
# 出力: outputs/tables/tbl01_markers_all.csv・tbl02_markers_top10.csv（ハーネス経由）
#
# ★double-dipping の注意（標準・必須注記）: 同一データで clustering と marker 検定を行うと p 値は過大化する
#   （二度漬け）。クラスタマーカーの p 値を絶対視せず、effect size（avg_log2FC）・発現特異性（pct.1/pct.2）・
#   文献照合・標本横断の再現性で総合判断する。p 値でクラスタの「有意性」を主張しない。
suppressPackageStartupMessages({ library(here); library(Seurat); library(dplyr) })
source(here::here("analysis", "00_config.R"))
source(here::here("R", "helpers.R"))

o <- readRDS(here::here("data", "processed", "seurat_integrated.rds"))
o <- JoinLayers(o)          # FindAllMarkers の前に layer を結合（Seurat v5）
Idents(o) <- "seurat_clusters"

# positive-only の cluster-marker 検定。p 値は double-dipping で過大化する点に留意（上記注記）。
mk <- FindAllMarkers(o, only.pos = TRUE, min.pct = 0.25,
                     logfc.threshold = 0.25, verbose = FALSE)
save_table(mk, "tbl01", "markers_all", "03_markers.R")
top <- mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10) %>% ungroup()
save_table(top, "tbl02", "markers_top10", "03_markers.R")

avail <- intersect(names(CONFIG$marker_hint), rownames(o))
cat("MARKERS_DONE clusters:", length(unique(mk$cluster)),
    "| hint markers found:", length(avail), "/", length(CONFIG$marker_hint), "\n")
cat("NOTE: marker p 値は double-dipping で過大化。次段 04_annotate.R で直交レイヤ注釈 + confidence を付す。\n")
