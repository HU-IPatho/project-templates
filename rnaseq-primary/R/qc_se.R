# ============================================================================
# gene-level SE の QC ロジック（library size / mapping rate / PCA）。
# 重い計算はここ（analysis から source）で行い data/processed に object 保存する。
# Quarto レポートは保存済み object を読むだけ（spec「Quarto 計算分離」）。
# ============================================================================
suppressPackageStartupMessages({
  library(SummarizedExperiment)
})

# NULL 合体（jsonlite が欠損キーで NULL を返す場合の保護）
`%||%` <- function(a, b) if (is.null(a)) b else a

#' salmon の各サンプル出力から mapping 率を読む（aux_info/meta_info.json）。
#' @param salmon_dirs 各サンプルの salmon 出力ディレクトリ（quant.sf の親）
#' @return data.frame(sample, num_processed, num_mapped, percent_mapped)
read_salmon_mapping <- function(salmon_dirs) {
  rows <- lapply(names(salmon_dirs), function(s) {
    meta_path <- file.path(salmon_dirs[[s]], "aux_info", "meta_info.json")
    if (!file.exists(meta_path)) {
      return(data.frame(sample = s, num_processed = NA_real_,
                        num_mapped = NA_real_, percent_mapped = NA_real_))
    }
    m <- jsonlite::fromJSON(meta_path)
    data.frame(sample = s,
               num_processed  = as.numeric(m$num_processed %||% NA),
               num_mapped     = as.numeric(m$num_mapped %||% NA),
               percent_mapped = as.numeric(m$percent_mapped %||% NA))
  })
  do.call(rbind, rows)
}

#' SE から QC 指標を計算する（library size / PCA）。mapping はオプションで結合。
#' @param se          gene-level SummarizedExperiment（assay "counts" を持つ）
#' @param mapping_df  read_salmon_mapping() の返り値（任意）
#' @param ntop        PCA に使う高分散遺伝子数
#' @return list(lib_size, mapping, pca, summary)
compute_se_qc <- function(se, mapping_df = NULL, ntop = 500) {
  counts <- SummarizedExperiment::assay(se, "counts")
  lib_size <- data.frame(sample = colnames(se),
                         library_size = colSums(counts))

  # log-CPM 変換 → 高分散遺伝子で PCA（DESeq2 非依存の軽量版）
  cpm  <- t(t(counts) / pmax(colSums(counts), 1)) * 1e6
  logc <- log2(cpm + 1)
  vars <- matrixStats::rowVars(as.matrix(logc))
  sel  <- order(vars, decreasing = TRUE)[seq_len(min(ntop, nrow(logc)))]
  pc   <- stats::prcomp(t(logc[sel, , drop = FALSE]), scale. = FALSE)
  pct  <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  pca  <- data.frame(sample = colnames(se),
                     PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                     pc1_var = pct[1], pc2_var = pct[2])

  summary <- merge(lib_size,
                   if (!is.null(mapping_df)) mapping_df else
                     data.frame(sample = colnames(se)),
                   by = "sample", all.x = TRUE)

  list(lib_size = lib_size, mapping = mapping_df, pca = pca, summary = summary)
}

# ---- 描画（ggplot オブジェクトを返すだけ。保存は save_fig が担う）----
plot_lib_size <- function(qc) {
  ggplot2::ggplot(qc$lib_size,
                  ggplot2::aes(x = stats::reorder(sample, -library_size),
                               y = library_size)) +
    ggplot2::geom_col() +
    ggplot2::labs(title = "Library size", x = "sample", y = "assigned reads")
}

plot_mapping_rate <- function(qc) {
  stopifnot(!is.null(qc$mapping))
  ggplot2::ggplot(qc$mapping,
                  ggplot2::aes(x = stats::reorder(sample, -percent_mapped),
                               y = percent_mapped)) +
    ggplot2::geom_col() +
    ggplot2::labs(title = "Salmon mapping rate", x = "sample", y = "% mapped")
}

plot_pca <- function(qc) {
  ggplot2::ggplot(qc$pca, ggplot2::aes(x = PC1, y = PC2, label = sample)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(vjust = -0.8, size = 3) +
    ggplot2::labs(title = "PCA (log-CPM, top-variance genes)",
                  x = sprintf("PC1 (%.1f%%)", qc$pca$pc1_var[1]),
                  y = sprintf("PC2 (%.1f%%)", qc$pca$pc2_var[1]))
}
