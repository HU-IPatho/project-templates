# ============================================================================
# DE 結果の共通整形と 2 エンジンの一致度
# ----------------------------------------------------------------------------
# DESeq2 / edgeR は列名・統計量が異なるため、下流（図・一致度）は統一列で扱う。
#   統一列: gene, logFC, expr(平均発現の log 尺度), pvalue, padj, contrast
# ============================================================================

# --- run_deseq2 / run_edger の結果を統一列の per-contrast data.frame へ ----------
tidy_de <- function(de) {
  lapply(de$results, function(df) {
    if (de$method == "DESeq2") {
      data.frame(
        gene     = df$gene,
        logFC    = df$log2FoldChange,
        expr     = log2(df$baseMean + 1),       # MA 図の A 軸（平均発現）
        pvalue   = df$pvalue,
        padj     = df$padj,
        contrast = df$contrast,
        stringsAsFactors = FALSE)
    } else {
      data.frame(
        gene     = df$gene,
        logFC    = df$logFC,
        expr     = df$logCPM,                    # edgeR は logCPM が平均発現尺度
        pvalue   = df$PValue,
        padj     = df$FDR,
        contrast = df$contrast,
        stringsAsFactors = FALSE)
    }
  })
}

# --- DESeq2 vs edgeR の一致度（両走時のみ）--------------------------------------
# 対比ごとに: 有意遺伝子集合の Jaccard・logFC の相関・logFC 符号一致率を返す。
compare_methods <- function(tidy_a, tidy_b, fdr = 0.05) {
  contrasts <- intersect(names(tidy_a), names(tidy_b))
  do.call(rbind, lapply(contrasts, function(nm) {
    a <- tidy_a[[nm]]
    b <- tidy_b[[nm]]
    m <- merge(a, b, by = "gene", suffixes = c("_a", "_b"))
    sig_a <- m$gene[!is.na(m$padj_a) & m$padj_a < fdr]
    sig_b <- m$gene[!is.na(m$padj_b) & m$padj_b < fdr]
    inter <- length(intersect(sig_a, sig_b))
    uni   <- length(union(sig_a, sig_b))
    data.frame(
      contrast         = nm,
      n_sig_a          = length(sig_a),
      n_sig_b          = length(sig_b),
      jaccard          = if (uni > 0) inter / uni else NA_real_,
      logFC_cor        = stats::cor(m$logFC_a, m$logFC_b, use = "complete.obs"),
      sign_concordance = mean(sign(m$logFC_a) == sign(m$logFC_b), na.rm = TRUE),
      stringsAsFactors = FALSE)
  }))
}
