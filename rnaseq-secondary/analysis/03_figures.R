# 03: 図。PCA / MA / volcano / top-DEG heatmap を save_fig（PNG+PDF 両形式）で保存。
# 実行: プロジェクトルートから  Rscript analysis/03_figures.R
# 重い計算はしない（01/02 が生成した se.rds / de.rds を読むだけ）。
suppressPackageStartupMessages({
  library(here); library(yaml); library(ggplot2)
  library(SummarizedExperiment); library(edgeR)
})
source(here::here("R", "helpers.R"))     # save_fig / %||%
source(here::here("R", "de_common.R"))   # tidy_de

cfg <- yaml::read_yaml(here::here("config.yaml"))
se  <- readRDS(here::here("data", "processed", "se.rds"))
de  <- readRDS(here::here("data", "processed", "de.rds"))
fdr <- cfg$fdr %||% 0.05
gcol <- cfg$design$group_col
grp_vec <- as.character(as.data.frame(SummarizedExperiment::colData(se))[[gcol]])

logcpm <- edgeR::cpm(SummarizedExperiment::assay(se, "counts"),
                     log = TRUE, prior.count = 2)

# ---- fig01: PCA（上位変動遺伝子の log-CPM）------------------------------------
vars <- apply(logcpm, 1, stats::var)
topv <- utils::head(order(vars, decreasing = TRUE), min(500L, nrow(logcpm)))
pc   <- stats::prcomp(t(logcpm[topv, , drop = FALSE]), scale. = FALSE)
pv   <- (pc$sdev^2 / sum(pc$sdev^2)) * 100
pca_df <- data.frame(
  sample = colnames(logcpm),
  PC1 = pc$x[, 1], PC2 = pc$x[, 2],
  group = as.data.frame(SummarizedExperiment::colData(se))[[gcol]])
p_pca <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  geom_text(aes(label = sample), vjust = -0.8, size = 3, show.legend = FALSE) +
  labs(x = sprintf("PC1 (%.1f%%)", pv[1]),
       y = sprintf("PC2 (%.1f%%)", pv[2]),
       title = "PCA (top variable genes, log-CPM)")
save_fig(p_pca, "fig01", "pca", "analysis/03_figures.R", width = 6, height = 5)

# ---- fig02: MDS（記述解析の併走・n=1 screening standard で必須・edgeR plotMDS）----
# 標準: n=1 KD screening lane では記述解析（FC ランキング + MDS・有意性主張なし）を必ず併走出力する
# （edgeR §2.13 option 1）。MDS は TMM 正規化後の log-CPM 距離で試料配置を見る（複製ありでも診断に有用）。
dge <- edgeR::calcNormFactors(
  edgeR::DGEList(counts = SummarizedExperiment::assay(se, "counts"), group = factor(grp_vec)),
  method = "TMM")
save_fig(function() {
  edgeR::plotMDS(dge, labels = colnames(dge), col = as.integer(factor(grp_vec)),
                 main = "MDS (TMM log-CPM distances)")
}, "fig02", "mds", "analysis/03_figures.R", width = 6, height = 5)

# ---- MA / volcano（主エンジンの各対比）---------------------------------------
primary <- if (!is.null(de$edgeR)) "edgeR" else "DESeq2"
tid     <- tidy_de(de[[primary]])
for (k in seq_along(tid)) {
  nm <- names(tid)[k]
  d  <- tid[[nm]]
  d$sig <- !is.na(d$padj) & d$padj < fdr

  # fig id は対比数に依存しない別バンドに分離（MA=fig1NN / volcano=fig2NN）。
  # 同一 id だと captions.tsv の upsert で上書きが起きるため、対比 >10 でも衝突させない。
  p_ma <- ggplot(d, aes(expr, logFC, color = sig)) +
    geom_point(alpha = 0.5, size = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick")) +
    labs(x = "mean expression (log)", y = "log2 fold-change",
         title = sprintf("MA: %s (%s)", nm, primary))
  save_fig(p_ma, sprintf("fig1%02d", k), paste0("ma_", nm),
           "analysis/03_figures.R", width = 6, height = 5)

  p_vol <- ggplot(d, aes(logFC, -log10(pvalue), color = sig)) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_color_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick")) +
    labs(x = "log2 fold-change", y = "-log10 p-value",
         title = sprintf("Volcano: %s (%s)", nm, primary))
  save_fig(p_vol, sprintf("fig2%02d", k), paste0("volcano_", nm),
           "analysis/03_figures.R", width = 6, height = 5)
}

# ---- fig40: top-DEG heatmap（第 1 対比の上位遺伝子・z-score 化 log-CPM）--------
d1    <- tid[[1]]
d1    <- d1[!is.na(d1$padj), ]
d1    <- d1[order(d1$padj), ]
n_top <- min(cfg$heatmap$n_top %||% 30L, nrow(d1))
topg  <- utils::head(d1$gene, n_top)
mat   <- logcpm[topg, , drop = FALSE]
matz  <- t(scale(t(mat)))                 # 遺伝子ごとに z-score
ann   <- as.data.frame(SummarizedExperiment::colData(se))[, gcol, drop = FALSE]
save_fig(function() {
  pheatmap::pheatmap(matz, annotation_col = ann, silent = FALSE,
                     main = sprintf("Top %d DE genes (%s, z log-CPM)", n_top, names(tid)[1]))
}, "fig40", "top_deg_heatmap", "analysis/03_figures.R", width = 7, height = 8)

cat("FIGURES_DONE: PCA + MA/volcano x", length(tid), "+ heatmap\n")
