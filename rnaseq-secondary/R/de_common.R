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

# --- cross-engine concordance（DESeq2 vs edgeR・補助 sanity check）----------------
# 標準（specs/bulk-secondary-deg-standard/spec.md）:
#   - cross-engine concordance は「頑健性の保証ではなく補助的な sanity check」（MUST）。
#     「一致=真/高信頼」と読める含意を持たせない（MUST NOT。両エンジンは NB-GLM 族・「大多数は非DE」
#     前提・類似 shrinkage を共有し failure mode が相関する＝共有系統誤差には concordant に誤りうる）。
#   - 主指標は Pearson でなく Spearman 順位相関（MUST）。Jaccard・符号一致率・非対称 overlap を
#     effect-size bin と方向で層別報告する（SHALL）。一致水準は n・effect size・方向に依存し小 n で
#     低く非対称化する。
#   - n=1 KD 経路には適用しない（MUST NOT・DESeq2 は複製なしで null）。呼出側は複製あり lane で両走した
#     対比のみを渡す（run_deseq2 は replicate lane のみ返すため intersect で自然に除外される）。
# 返り値: long data.frame（contrast × stratum・stratum ∈ {all, large, medium, small}＝|logFC_a| 三分位）。
compare_methods <- function(tidy_a, tidy_b, fdr = 0.05) {
  contrasts <- intersect(names(tidy_a), names(tidy_b))
  spearman <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3) return(NA_real_)
    suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
  }
  one_stratum <- function(m, nm, stratum) {
    sig_a <- m$gene[!is.na(m$padj_a) & m$padj_a < fdr]
    sig_b <- m$gene[!is.na(m$padj_b) & m$padj_b < fdr]
    inter <- length(intersect(sig_a, sig_b)); uni <- length(union(sig_a, sig_b))
    data.frame(
      contrast          = nm,
      stratum           = stratum,
      n_genes           = nrow(m),
      n_sig_a           = length(sig_a),
      n_sig_b           = length(sig_b),
      jaccard           = if (uni > 0) inter / uni else NA_real_,
      overlap_a_in_b    = if (length(sig_a) > 0) inter / length(sig_a) else NA_real_,  # 非対称 overlap
      overlap_b_in_a    = if (length(sig_b) > 0) inter / length(sig_b) else NA_real_,
      spearman_logFC    = spearman(m$logFC_a, m$logFC_b),
      sign_concordance  = mean(sign(m$logFC_a) == sign(m$logFC_b), na.rm = TRUE),
      stringsAsFactors  = FALSE)
  }
  do.call(rbind, lapply(contrasts, function(nm) {
    m <- merge(tidy_a[[nm]], tidy_b[[nm]], by = "gene", suffixes = c("_a", "_b"))
    rows <- list(one_stratum(m, nm, "all"))
    # effect-size 層別（|logFC_a| 三分位）: 小 n で低く非対称化する依存性を surface する
    absfc <- abs(m$logFC_a)
    if (sum(is.finite(absfc)) >= 6) {
      qs <- stats::quantile(absfc, probs = c(1/3, 2/3), na.rm = TRUE)
      bin <- ifelse(absfc <= qs[1], "small", ifelse(absfc <= qs[2], "medium", "large"))
      for (s in c("large", "medium", "small")) {
        ms <- m[!is.na(bin) & bin == s, , drop = FALSE]
        if (nrow(ms) > 0) rows[[length(rows) + 1]] <- one_stratum(ms, nm, s)
      }
    }
    do.call(rbind, rows)
  }))
}
