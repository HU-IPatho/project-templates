# ============================================================================
# DTE/DTU 可視化 — 図は ggplot オブジェクトを返すだけ（保存は R/helpers.R の save_fig）。
# ----------------------------------------------------------------------------
# 主役 = isoform 使用比プロット（群 × transcript proportion）。DTU 有意遺伝子について、
#   各 transcript の「遺伝子内 count 比率」を群ごとに示す。比率の群間シフト = isoform switching。
#   比率は se_tx（dtuScaledTPM）の観測 counts から算出する（05 の検定とは独立の可視化用集計）。
# 併走 = DTE volcano（transcript 発現差）。
# ============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(SummarizedExperiment)
})

# se_tx の counts から、指定 gene 群の transcript 比率（sample→群平均）を長形式で返す。
#   proportion(tx, group) = mean over group of ( count_tx / 遺伝子内 count 合計 )。
#   単一 isoform の遺伝子は DTU 概念が無いためスキップ（proportion は常に 1）。
isoform_proportion_df <- function(se, genes, gcol) {
  rd  <- as.data.frame(SummarizedExperiment::rowData(se))
  grp <- as.data.frame(SummarizedExperiment::colData(se))[[gcol]]
  cts <- as.matrix(SummarizedExperiment::assay(se, "counts"))
  rows <- lapply(genes, function(g) {
    idx <- which(rd$gene_id == g)
    if (length(idx) < 2) return(NULL)
    sub  <- cts[idx, , drop = FALSE]
    prop <- sweep(sub, 2, pmax(colSums(sub), 1), "/")   # tx x sample の遺伝子内比率
    do.call(rbind, lapply(unique(grp), function(gg) {
      cols <- which(grp == gg)
      data.frame(gene_id    = g,
                 feature_id = rd$tx_id[idx],
                 group      = gg,
                 proportion = rowMeans(prop[, cols, drop = FALSE]),
                 stringsAsFactors = FALSE)
    }))
  })
  do.call(rbind, rows)
}

# isoform 使用比プロット（遺伝子ごとに facet・群で色分けした dodge 棒）。
plot_isoform_proportions <- function(se, genes, gcol,
                                     title = "Isoform usage (top DTU genes)") {
  df <- isoform_proportion_df(se, genes, gcol)
  if (is.null(df) || nrow(df) == 0) {
    return(ggplot2::ggplot() +
             ggplot2::labs(title = paste0(title, " — 対象 DTU 遺伝子なし")))
  }
  ggplot2::ggplot(df, ggplot2::aes(x = feature_id, y = proportion, fill = group)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8)) +
    ggplot2::facet_wrap(~ gene_id, scales = "free_x") +
    ggplot2::labs(title = title, x = "transcript", y = "isoform proportion (群平均)") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6))
}

# DTE volcano（transcript 発現差）。tidy_dte() の統一列（feature/logFC/pvalue/padj）を受ける。
plot_dte_volcano <- function(d, fdr = 0.05, title = "DTE volcano") {
  d$sig <- !is.na(d$padj) & d$padj < fdr
  ggplot2::ggplot(d, ggplot2::aes(x = logFC, y = -log10(pvalue), color = sig)) +
    ggplot2::geom_point(alpha = 0.5, size = 0.8) +
    ggplot2::scale_color_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick")) +
    ggplot2::labs(title = title, x = "log2 fold-change (transcript)",
                  y = "-log10 p-value", color = sprintf("padj<%.2g", fdr))
}
