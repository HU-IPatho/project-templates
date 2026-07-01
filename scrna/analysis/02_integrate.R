# 02: 正規化 → PCA → Harmony 統合 → クラスタ → UMAP。
# ★ここが長時間ジョブ（Harmony）。前景で直接実行せず job-run で回すこと（AGENTS.md 参照）:
#     JID=$(job-run --label integrate -- Rscript analysis/02_integrate.R)
#     job-wait "$JID" --timeout 570   # 未完(124)なら間を空けず再度
# 入力: data/interim/seurat_qc.rds（01 の出力・中間物）
# 出力: data/processed/seurat_integrated.rds（＝正準 object・共有昇格の単位）+ UMAP 図
suppressPackageStartupMessages({ library(here); library(Seurat); library(ggplot2) })
source(here::here("analysis", "00_config.R"))
source(here::here("R", "helpers.R"))

o <- readRDS(here::here("data", "interim", "seurat_qc.rds"))
o <- NormalizeData(o, verbose = FALSE)
o <- FindVariableFeatures(o, nfeatures = CONFIG$n_variable, verbose = FALSE)
o <- ScaleData(o, verbose = FALSE)
o <- RunPCA(o, npcs = CONFIG$n_pcs, verbose = FALSE)

# バッチが複数あるときだけ Harmony 統合。単一バッチなら PCA をそのまま使う。
n_batch <- length(unique(o@meta.data[[CONFIG$batch_key]]))
if (n_batch > 1) {
  suppressPackageStartupMessages(library(harmony))
  cat(sprintf("Harmony 統合 (batch_key=%s, %d バッチ) ...\n", CONFIG$batch_key, n_batch))
  o <- RunHarmony(o, group.by.vars = CONFIG$batch_key, verbose = TRUE)
  reduction <- "harmony"
} else {
  cat("単一バッチ → Harmony スキップ（PCA を使用）\n")
  reduction <- "pca"
}

o <- FindNeighbors(o, reduction = reduction, dims = 1:CONFIG$dims, verbose = FALSE)
o <- FindClusters(o, resolution = CONFIG$resolution, verbose = FALSE)
o <- RunUMAP(o, reduction = reduction, dims = 1:CONFIG$dims, verbose = FALSE)

# 正準 object を data/processed へ保存（下流 03 とレポートの開始点・共有昇格候補）。
dir.create(here::here("data", "processed"), showWarnings = FALSE, recursive = TRUE)
saveRDS(o, here::here("data", "processed", "seurat_integrated.rds"))

# UMAP 図は出力ハーネス経由（PNG+PDF・統一テーマ・captions.tsv 記録）。
save_fig(DimPlot(o, group.by = "seurat_clusters", label = TRUE) + ggtitle("clusters"),
         "fig01", "umap_clusters", "02_integrate.R", width = 7, height = 6)
save_fig(DimPlot(o, group.by = CONFIG$batch_key) + ggtitle(paste("by", CONFIG$batch_key)),
         "fig02", "umap_batch", "02_integrate.R", width = 8, height = 6)
cat(sprintf("INTEGRATE_DONE: %d clusters, %d cells\n",
            length(levels(o$seurat_clusters)), ncol(o)))
