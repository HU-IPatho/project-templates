# ============================================================================
# edgeR パス — DGEList で DEG。★n=1 DEG（各群 1 サンプル）が本テンプレの主眼★
# ----------------------------------------------------------------------------
# 標準（複製あり）: filterByExpr → TMM → estimateDisp(design) → glmFit/glmLRT。
#
# n=1 DEG（NC/sh1/sh2 各 n=1 等・最重要）:
#   複製が無いと分散を各遺伝子から推定できない。そこで housekeeping 遺伝子（群間で
#   発現不変が前提）の群間ばらつきを分散の代理として使い common.dispersion を推定する。
#   手順（契約 SSoT）:
#     1) group を全 1 に潰した DGEList の HK 部分集合に対し
#        estimateDisp(HK, trend.method='none', tagwise=FALSE) → common.dispersion
#     2) その common.dispersion を全遺伝子の DGEList へ移植
#     3) glmFit/glmLRT で多対比（NC vs sh1, NC vs sh2 …）
#   HK が不足/推定失敗なら BCV 固定へフォールバック（human 0.4 / 細胞株 0.1 / technical 0.01）。
#   探索的であり FDR の解釈は限定的である旨を warning で明示する。
# ============================================================================
suppressPackageStartupMessages({
  library(edgeR)
  library(SummarizedExperiment)
})

run_edger <- function(se, cfg, hk_genes = NULL) {
  counts  <- SummarizedExperiment::assay(se, "counts")
  coldata <- as.data.frame(SummarizedExperiment::colData(se))
  gcol    <- cfg$design$group_col
  group   <- factor(coldata[[gcol]])

  y      <- edgeR::DGEList(counts = counts, group = group)
  design <- edger_design(coldata, cfg)

  keep <- edgeR::filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y, method = "TMM")

  is_n1 <- all(table(group) == 1)                 # 各群 1 サンプルか
  if (is_n1) {
    est <- estimate_hk_dispersion(y, hk_genes, cfg)
    y$common.dispersion  <- est$common.dispersion # 全遺伝子へ移植（trended/tagwise は付けない）
    y$trended.dispersion <- NULL
    y$tagwise.dispersion <- NULL
    disp_source <- est$source
    warning("n=1 DEG: ", disp_source,
            " による探索的解析。FDR/有意性の解釈は限定的（複製が無いため）。")
  } else {
    y <- edgeR::estimateDisp(y, design)
    disp_source <- "estimateDisp（複製あり）"
  }

  fit <- edgeR::glmFit(y, design)
  contrasts <- edger_contrasts(design, cfg)
  res <- lapply(names(contrasts), function(nm) {
    lrt <- edgeR::glmLRT(fit, contrast = contrasts[[nm]])
    tt  <- edgeR::topTags(lrt, n = Inf, sort.by = "none")$table
    tt$gene     <- rownames(tt)
    tt$contrast <- nm
    tt
  })
  names(res) <- names(contrasts)

  list(results          = res,
       method           = "edgeR",
       n1               = is_n1,
       common.dispersion = y$common.dispersion,
       bcv              = sqrt(y$common.dispersion),
       disp_source      = disp_source)
}

# --- HK 遺伝子で common.dispersion を推定（group を全 1 に潰す）--------------------
# 返り値 list(common.dispersion, source)。失敗時は BCV 固定にフォールバック。
estimate_hk_dispersion <- function(y, hk_genes, cfg) {
  min_hk   <- cfg$edger$min_hk_genes %||% 10
  bcv_fb   <- cfg$edger$bcv_fallback %||% 0.4

  fallback <- function(reason) {
    list(common.dispersion = bcv_fb^2,
         source = sprintf("BCV 固定=%.2f（%s）", bcv_fb, reason))
  }
  if (is.null(hk_genes) || length(hk_genes) == 0) return(fallback("HK リスト未指定"))

  present <- intersect(hk_genes, rownames(y))
  if (length(present) < min_hk) {
    return(fallback(sprintf("HK 遺伝子 %d < 必要 %d", length(present), min_hk)))
  }

  y_hk <- y[present, , keep.lib.sizes = FALSE]
  y_hk$samples$group <- factor(rep(1L, ncol(y_hk)))     # 全サンプルを 1 群扱い
  est <- tryCatch(
    edgeR::estimateDisp(y_hk, design = NULL,
                        trend.method = "none", tagwise = FALSE),
    error = function(e) NULL)
  if (is.null(est) || !is.finite(est$common.dispersion) || est$common.dispersion <= 0) {
    return(fallback("estimateDisp（HK）失敗"))
  }
  list(common.dispersion = est$common.dispersion,
       source = sprintf("HK common.dispersion（%d 遺伝子・BCV=%.3f）",
                        length(present), sqrt(est$common.dispersion)))
}
