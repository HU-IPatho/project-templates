# 03: 図。PCA / MA / volcano / top-DEG heatmap を save_fig（PNG+PDF 両形式）で保存。
# 実行: プロジェクトルートから  Rscript analysis/03_figures.R
# 重い計算はしない（01/02 が生成した se.rds / de.rds を読むだけ）。
suppressPackageStartupMessages({
  library(here); library(yaml); library(ggplot2)
  library(SummarizedExperiment); library(edgeR); library(limma)
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
  # plotMDS は limma の generic（edgeR は DGEList メソッドを提供）。edgeR:: では export されない。
  limma::plotMDS(dge, labels = colnames(dge), col = as.integer(factor(grp_vec)),
                 main = "MDS (TMM log-CPM distances)")
}, "fig02", "mds", "analysis/03_figures.R", width = 6, height = 5)

# ---- MA / volcano（主エンジンの各対比・lane 別）------------------------------
# 標準: n=1 screening lane の出力は screening-grade。較正済み FDR/adjusted-p を確定的有意として
# 図に提示してはならない（MUST NOT）。ヒットは FC/収縮 LFC ランキングに限る（SHALL）。
# → screening 対比は padj 有意色分け・volcano(p 値提示)を出さず、|logFC| 上位候補で着色し表題に明示する。
# 複製あり lane は較正済み FDR ゆえ従来どおり padj<fdr 有意色分け・volcano を出す。
primary <- if (!is.null(de$edgeR)) "edgeR" else "DESeq2"
tid     <- tidy_de(de[[primary]])
lane_of <- de[[primary]]$lane                       # 名前付き contrast→lane（DESeq2 primary では NULL）
get_lane <- function(nm) if (!is.null(lane_of)) (lane_of[[nm]] %||% "replicate") else "replicate"
cand_top <- cfg$screening$candidate_rank_top %||% 100L

for (k in seq_along(tid)) {
  nm   <- names(tid)[k]
  d    <- tid[[nm]]
  lane <- get_lane(nm)

  if (lane == "screening") {
    # screening-grade: |logFC| 上位を「候補」として着色（padj 非依存）。volcano は出さない。
    ord  <- order(-abs(d$logFC))
    ncand <- min(cand_top, nrow(d))
    d$flag <- FALSE; d$flag[ord[seq_len(ncand)]] <- TRUE
    ttl  <- sprintf("MA: %s [screening-grade: |logFC| 上位%d 候補・確定FDRでない]", nm, ncand)
    lgd  <- "FC候補"
  } else {
    d$flag <- !is.na(d$padj) & d$padj < fdr          # 複製あり: 較正済み FDR 有意
    ttl  <- sprintf("MA: %s [replicate, FDR<%.2g] (%s)", nm, fdr, primary)
    lgd  <- sprintf("padj<%.2g", fdr)
  }

  # fig id は対比数に依存しない別バンドに分離（MA=fig1NN / volcano=fig2NN）。
  p_ma <- ggplot(d, aes(expr, logFC, color = flag)) +
    geom_point(alpha = 0.5, size = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick"),
                       labels = c(`FALSE` = "other", `TRUE` = lgd), name = NULL) +
    labs(x = "mean expression (log)", y = "log2 fold-change", title = ttl)
  save_fig(p_ma, sprintf("fig1%02d", k), paste0("ma_", nm),
           "analysis/03_figures.R", width = 6, height = 5)

  if (lane != "screening") {
    # volcano は p 値提示ゆえ複製あり lane のみ（screening は確定 FDR を提示しない）
    p_vol <- ggplot(d, aes(logFC, -log10(pvalue), color = flag)) +
      geom_point(alpha = 0.5, size = 0.8) +
      scale_color_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick"),
                         labels = c(`FALSE` = "other", `TRUE` = lgd), name = NULL) +
      labs(x = "log2 fold-change", y = "-log10 p-value",
           title = sprintf("Volcano: %s [replicate] (%s)", nm, primary))
    save_fig(p_vol, sprintf("fig2%02d", k), paste0("volcano_", nm),
             "analysis/03_figures.R", width = 6, height = 5)
  }
}

# ---- fig40: top-gene heatmap（第 1 対比の上位遺伝子・z-score 化 log-CPM）--------
# lane 別のランキング基準で上位を選ぶ: 複製あり=padj 昇順（較正済み有意）/ screening=|logFC| 降順
# （screening-grade ゆえ FC 候補・padj 非依存）。全 NA/空/1 行なら heatmap は skip（pheatmap は
# クラスタリングに >=2 行要・G1 記述降格で padj 全 NA でも落ちない）。
nm1   <- names(tid)[1]
lane1 <- get_lane(nm1)
d1    <- tid[[1]]
if (lane1 == "screening") {
  d1  <- d1[is.finite(d1$logFC), ]
  d1  <- d1[order(-abs(d1$logFC)), ]
  basis <- sprintf("|logFC| top (screening-grade: %s)", nm1)
} else {
  d1  <- d1[is.finite(d1$padj), ]
  d1  <- d1[order(d1$padj), ]
  basis <- sprintf("FDR top (%s)", nm1)
}
n_top <- min(cfg$heatmap$n_top %||% 30L, nrow(d1))
if (n_top >= 2) {
  topg  <- utils::head(d1$gene, n_top)
  mat   <- logcpm[topg, , drop = FALSE]
  matz  <- t(scale(t(mat)))                 # 遺伝子ごとに z-score
  ann   <- as.data.frame(SummarizedExperiment::colData(se))[, gcol, drop = FALSE]
  save_fig(function() {
    pheatmap::pheatmap(matz, annotation_col = ann, silent = FALSE,
                       main = sprintf("Top %d genes [%s, z log-CPM]", n_top, basis))
  }, "fig40", "top_gene_heatmap", "analysis/03_figures.R", width = 7, height = 8)
} else {
  message("fig40 heatmap: 上位遺伝子が 2 未満（記述降格 等）→ heatmap を skip。")
}

cat("FIGURES_DONE: PCA + MDS + MA/volcano(lane別) x", length(tid),
    "+ heatmap(", if (n_top >= 2) "描画" else "skip", ")\n")
