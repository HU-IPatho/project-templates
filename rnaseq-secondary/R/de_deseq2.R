# ============================================================================
# DESeq2 パス — DESeqDataSet で DEG（median-of-ratios 正規化・複製ありの標準経路）
# ----------------------------------------------------------------------------
# design/contrast は config.yaml（design.R が formula を組む）。バッチは共変量で補正。
# 前処理: 低発現 prefilter（DESeq2 の正規化・分散推定は内部で median-of-ratios）。
# 注: DESeq2 は各群に複製が要る。n=1（各群 1 サンプル）では残差自由度が無く走らないため、
#     n=1 は edgeR の HK 分散パス（de_edger.R）を使う（config の method で選択）。
# ============================================================================
suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(S4Vectors)
})

run_deseq2 <- function(se, cfg) {
  coldata <- as.data.frame(SummarizedExperiment::colData(se))
  gcol <- cfg$design$group_col
  covs <- (cfg$design$covariates %||% character(0))
  covs <- covs[nzchar(covs)]
  for (v in c(gcol, covs)) coldata[[v]] <- factor(coldata[[v]])
  if (!is.null(cfg$design$reference_level)) {
    coldata[[gcol]] <- stats::relevel(coldata[[gcol]], ref = cfg$design$reference_level)
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = SummarizedExperiment::assay(se, "counts"),
    colData   = coldata,
    design    = deseq2_formula(coldata, cfg))

  # 軽い prefilter（全サンプル合計 count が閾値未満の遺伝子を落とす）
  min_sum <- cfg$prefilter$min_count_sum %||% 10
  dds <- dds[rowSums(DESeq2::counts(dds)) >= min_sum, ]

  dds <- DESeq2::DESeq(dds, quiet = TRUE)

  shrink <- isTRUE(cfg$deseq2$shrink_lfc)
  res <- lapply(names(cfg$contrasts), function(nm) {
    spec <- cfg$contrasts[[nm]]
    ct   <- c(gcol, spec[[1]], spec[[2]])
    r <- DESeq2::results(dds, contrast = ct)
    if (shrink) {
      # apeglm は coef 指定が要るため、shrink 時は coef ベースへ切替（reference に対する対比）
      coef <- paste0(gcol, "_", spec[[1]], "_vs_", spec[[2]])
      r <- tryCatch(
        DESeq2::lfcShrink(dds, coef = coef, type = "apeglm", res = r, quiet = TRUE),
        error = function(e) { warning("lfcShrink 失敗（生 LFC を使用）: ",
                                      conditionMessage(e)); r })
    }
    df <- as.data.frame(r)
    df$gene     <- rownames(df)
    df$contrast <- nm
    df
  })
  names(res) <- names(cfg$contrasts)
  list(results = res, method = "DESeq2", dds = dds)
}
