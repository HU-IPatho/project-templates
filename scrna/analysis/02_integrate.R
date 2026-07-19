# 02: 正規化 → (cell-cycle スコア) → PCA → Harmony 統合 → 複数解像度クラスタ → UMAP。
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

# --- cell-cycle スコア（identity と直交な別軸・標準）---
# S/G2M スコアと Phase を付す。無条件の regression 除去はしない（cell_cycle.regress=TRUE のときだけ回帰で除去）。
# cc.genes.updated.2019 はヒト遺伝子記号。マウス等では ortholog 供給が要る（該当遺伝子が無ければ scoring は skip）。
regress_vars <- NULL
if (isTRUE(CONFIG$cell_cycle$score)) {
  ok <- tryCatch({
    o <- CellCycleScoring(o,
                          s.features   = cc.genes.updated.2019$s.genes,
                          g2m.features = cc.genes.updated.2019$g2m.genes,
                          set.ident = FALSE)
    TRUE
  }, error = function(e) { warning("cell-cycle scoring を skip（該当遺伝子不足の可能性）: ", conditionMessage(e)); FALSE })
  if (ok && isTRUE(CONFIG$cell_cycle$regress)) {
    regress_vars <- c("S.Score", "G2M.Score")   # identity 歪みが実証された時のみ（既定 FALSE）
    cat("cell-cycle: S.Score/G2M.Score を ScaleData で回帰除去します（regress=TRUE）\n")
  }
}

o <- FindVariableFeatures(o, nfeatures = CONFIG$n_variable, verbose = FALSE)
o <- ScaleData(o, vars.to.regress = regress_vars, verbose = FALSE)
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

# --- 複数解像度スキャン（標準: 単一固定でなく安定性・over/under-clustering を評価）---
# CONFIG$resolutions 各々でクラスタし <assay>_snn_res.<r> 列に保持。primary（CONFIG$resolution）を
# seurat_clusters に採用して下流 03/04 が使う。cluster 数を解像度別に報告し over-clustering の目安にする。
o <- FindClusters(o, resolution = CONFIG$resolutions, verbose = FALSE)
res_cols <- paste0(DefaultAssay(o), "_snn_res.", CONFIG$resolutions)
for (i in seq_along(CONFIG$resolutions)) {
  rc <- res_cols[i]
  if (rc %in% colnames(o@meta.data))
    cat(sprintf("  resolution %g → %d clusters\n", CONFIG$resolutions[i], length(unique(o@meta.data[[rc]]))))
}
primary_col <- paste0(DefaultAssay(o), "_snn_res.", CONFIG$resolution)
if (!primary_col %in% colnames(o@meta.data)) {
  o <- FindClusters(o, resolution = CONFIG$resolution, verbose = FALSE)  # フォールバック
  primary_col <- "seurat_clusters"
}
o$seurat_clusters <- factor(o@meta.data[[primary_col]])
Idents(o) <- "seurat_clusters"
cat(sprintf("primary resolution=%g → %d clusters（seurat_clusters）\n",
            CONFIG$resolution, length(levels(o$seurat_clusters))))

o <- RunUMAP(o, reduction = reduction, dims = 1:CONFIG$dims, verbose = FALSE)

# 正準 object を data/processed へ保存（下流 03/04 とレポートの開始点・共有昇格候補）。
dir.create(here::here("data", "processed"), showWarnings = FALSE, recursive = TRUE)
saveRDS(o, here::here("data", "processed", "seurat_integrated.rds"))

# UMAP 図は出力ハーネス経由（PNG+PDF・統一テーマ・captions.tsv 記録）。
save_fig(DimPlot(o, group.by = "seurat_clusters", label = TRUE) + ggtitle("clusters (primary)"),
         "fig01", "umap_clusters", "02_integrate.R", width = 7, height = 6)
save_fig(DimPlot(o, group.by = CONFIG$batch_key) + ggtitle(paste("by", CONFIG$batch_key)),
         "fig02", "umap_batch", "02_integrate.R", width = 8, height = 6)
# 解像度比較（over/under-clustering を目で確認する・標準の複数解像度評価）。
res_cols_present <- res_cols[res_cols %in% colnames(o@meta.data)]
if (length(res_cols_present) > 1) {
  save_fig(DimPlot(o, group.by = res_cols_present, label = TRUE, combine = TRUE),
           "fig03", "umap_resolution_scan", "02_integrate.R",
           width = 6 * length(res_cols_present), height = 6)
}
# cell-cycle（直交軸）の可視化。
if ("Phase" %in% colnames(o@meta.data)) {
  save_fig(DimPlot(o, group.by = "Phase") + ggtitle("cell-cycle phase"),
           "fig04", "umap_cellcycle", "02_integrate.R", width = 7, height = 6)
}
cat(sprintf("INTEGRATE_DONE: %d clusters, %d cells\n",
            length(levels(o$seurat_clusters)), ncol(o)))
