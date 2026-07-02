# ============================================================================
# DTE（differential transcript expression）— transcript「発現量」の群間差
# ----------------------------------------------------------------------------
# ② rnaseq-secondary の DESeq2/edgeR エンジンを **transcript-level** counts に適用する。
# 入力は se_tx の assay(counts_dte)＝04 が countsFromAbundance="lengthScaledTPM" で作った count。
#   lengthScaledTPM は各転写産物を自身の長さで scale した count で、DTE（発現「量」差）の妥当な
#   count 表現。dtuScaledTPM（assay counts）は遺伝子内 isoform 比率専用ゆえ DTE には使わない
#   （dtuScaledTPM を DTE に流用すると p 値/FDR が miscalibrated になる）。
# design/contrast は R/design.R（config.yaml から formula/contrast を組む）。
#
# DTE と DTU の違い（AGENTS.md 参照）:
#   DTE = 各 transcript の発現「量」が群間で違うか（本ファイル・DESeq2/edgeR）。
#   DTU = 遺伝子内の isoform「使用比」が群間で違うか（R/dtu.R・DRIMSeq+stageR）。
# 注: 複製は 05 冒頭の check_replicates() で保証済み（各群 n>=2）。
# ============================================================================
suppressPackageStartupMessages({
  library(SummarizedExperiment)
})

# --- ディスパッチ: config.dte.method（deseq2 / edger / both）--------------------
run_dte <- function(se, cfg) {
  method <- cfg$dte$method %||% "deseq2"
  res <- list()
  if (method %in% c("deseq2", "both")) res$DESeq2 <- .dte_deseq2(se, cfg)
  if (method %in% c("edger",  "both")) res$edgeR  <- .dte_edger(se, cfg)
  if (length(res) == 0) stop("run_dte: 未知の dte.method: ", method, "（deseq2 / edger / both）")
  res
}

# --- DESeq2 パス（tx-level・median-of-ratios）------------------------------------
# DTE は lengthScaledTPM の assay(counts_dte) を使う（dtuScaledTPM は DTU 専用のため DTE には
# 使わない）。lengthScaledTPM は非整数 → DESeq2 は整数 count を要するため round する。
.dte_deseq2 <- function(se, cfg) {
  suppressPackageStartupMessages({ library(DESeq2); library(S4Vectors) })
  coldata <- as.data.frame(SummarizedExperiment::colData(se))
  gcol <- cfg$design$group_col
  covs <- (cfg$design$covariates %||% character(0)); covs <- covs[nzchar(covs)]
  for (v in c(gcol, covs)) coldata[[v]] <- factor(coldata[[v]])
  if (!is.null(cfg$design$reference_level)) {
    coldata[[gcol]] <- stats::relevel(coldata[[gcol]], ref = cfg$design$reference_level)
  }
  cts <- round(as.matrix(SummarizedExperiment::assay(se, "counts_dte")))
  storage.mode(cts) <- "integer"

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = cts, colData = coldata, design = deseq2_formula(coldata, cfg))
  # 低発現 transcript の軽い prefilter（全サンプル合計 count が閾値未満を落とす）。
  # DTE 専用キー cfg$dte$min_count_sum（② rnaseq-secondary の prefilter$min_count_sum に相当）。
  min_sum <- cfg$dte$min_count_sum %||% 10
  dds <- dds[rowSums(DESeq2::counts(dds)) >= min_sum, ]
  dds <- DESeq2::DESeq(dds, quiet = TRUE)

  res <- lapply(names(cfg$contrasts), function(nm) {
    spec <- cfg$contrasts[[nm]]
    r <- DESeq2::results(dds, contrast = c(gcol, spec[[1]], spec[[2]]))
    df <- as.data.frame(r)
    df$feature  <- rownames(df)
    df$contrast <- nm
    df
  })
  names(res) <- names(cfg$contrasts)
  list(results = res, method = "DESeq2")
}

# --- edgeR パス（tx-level・複製あり標準経路: filterByExpr→TMM→estimateDisp→glmLRT）--
.dte_edger <- function(se, cfg) {
  suppressPackageStartupMessages({ library(edgeR) })
  coldata <- as.data.frame(SummarizedExperiment::colData(se))
  gcol    <- cfg$design$group_col
  group   <- factor(coldata[[gcol]])
  cts     <- round(as.matrix(SummarizedExperiment::assay(se, "counts_dte")))  # DTE=lengthScaledTPM・非整数の警告回避

  y      <- edgeR::DGEList(counts = cts, group = group)
  design <- edger_design(coldata, cfg)
  keep <- edgeR::filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y, method = "TMM")
  y <- edgeR::estimateDisp(y, design)     # 複製ありなので通常の分散推定（HK 分散は使わない）
  fit <- edgeR::glmFit(y, design)

  contrasts <- edger_contrasts(design, cfg)
  res <- lapply(names(contrasts), function(nm) {
    lrt <- edgeR::glmLRT(fit, contrast = contrasts[[nm]])
    tt  <- edgeR::topTags(lrt, n = Inf, sort.by = "none")$table
    tt$feature  <- rownames(tt)
    tt$contrast <- nm
    tt
  })
  names(res) <- names(contrasts)
  list(results = res, method = "edgeR")
}

# --- 2 エンジンの結果を統一列（feature/logFC/expr/pvalue/padj/contrast）へ整形 -----
# 下流（volcano・表）は統一列で扱う。gene 単位ではなく transcript(feature)単位である点が
# ② tidy_de との違い。
tidy_dte <- function(de) {
  lapply(de$results, function(df) {
    if (de$method == "DESeq2") {
      data.frame(
        feature  = df$feature,
        logFC    = df$log2FoldChange,
        expr     = log2(df$baseMean + 1),   # volcano/MA の平均発現尺度
        pvalue   = df$pvalue,
        padj     = df$padj,
        contrast = df$contrast,
        stringsAsFactors = FALSE)
    } else {
      data.frame(
        feature  = df$feature,
        logFC    = df$logFC,
        expr     = df$logCPM,               # edgeR は logCPM が平均発現尺度
        pvalue   = df$PValue,
        padj     = df$FDR,
        contrast = df$contrast,
        stringsAsFactors = FALSE)
    }
  })
}
