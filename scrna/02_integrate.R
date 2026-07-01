# 02: 正規化 → PCA → Harmony 統合 → クラスタ → UMAP
# ★ここが長時間ジョブ（Harmony）。前景で直接実行せず job-run で回すこと（AGENTS.md 参照）:
#     JID=$(job-run --label integrate -- Rscript 02_integrate.R)
#     job-wait "$JID" --timeout 570   # 未完(124)なら間を空けず再度
suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
source("00_config.R")

o <- readRDS("data/processed/seurat_qc.rds")
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

saveRDS(o, "data/processed/seurat_integrated.rds")
dir.create("reports", showWarnings = FALSE)
ggsave("reports/umap_clusters.png",
       DimPlot(o, group.by = "seurat_clusters", label = TRUE) + ggtitle("clusters"),
       width = 7, height = 6, dpi = 120)
ggsave("reports/umap_batch.png",
       DimPlot(o, group.by = CONFIG$batch_key) + ggtitle(paste("by", CONFIG$batch_key)),
       width = 8, height = 6, dpi = 120)
cat(sprintf("INTEGRATE_DONE: %d clusters, %d cells\n",
            length(levels(o$seurat_clusters)), ncol(o)))
