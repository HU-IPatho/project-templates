# ============================================================================
# DESeq2 パス — median-of-ratios 正規化・複製ありの較正済み FDR lane
# ----------------------------------------------------------------------------
# 標準（specs/bulk-secondary-deg-standard/spec.md）に conform:
#   - 対照（NC）を reference level に明示 relevel し関心対比の coef 存在を保証する。
#   - 低発現 prefilter は「必須でない」（速度/メモリ目的）。FDR 統制は results() の
#     independent filtering が担う（別概念・用語区別）。
#   - LFC 収縮: apeglm は単一係数（参照との対比）専用。非参照間対比（sh1_vs_sh2 等）は
#     type="ashr"（contrast 対応）へ分岐する。apeglm 適用不可で生 LFC に落ちるときは
#     「収縮なし」と明示ラベルし silent 化しない（MUST NOT）。
#   - DESeq2 は複製なしで走らない（n=1 対比は edgeR screening lane が担う）。ここは routing で
#     複製あり（min(実効群サイズ)>=2）に振られた対比のみ処理する。
# ============================================================================
suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(S4Vectors)
})

# --- LFC 収縮の実行（収縮種別を明示ラベルで返す・生 LFC fallback を silent 化しない）------
# 返り値: list(res, shrink_used)。shrink_used ∈ {"none","apeglm","ashr"}。
.deseq2_shrink <- function(dds, r, ct, num, den, ref, gcol, mode) {
  if (identical(mode, "none")) return(list(res = r, shrink_used = "none"))
  try_shrink <- function(type, ...) {
    tryCatch(list(res = DESeq2::lfcShrink(dds, ..., type = type, res = r, quiet = TRUE),
                  shrink_used = type),
             error = function(e) {
               warning(sprintf("lfcShrink(type=%s) 失敗 → 収縮なし（生 LFC・明示ラベル shrink=none）: %s",
                               type, conditionMessage(e)))
               list(res = r, shrink_used = "none")
             })
  }
  if (identical(mode, "apeglm")) {
    if (identical(den, ref)) {                       # 単一係数（参照との対比）でのみ apeglm 可
      coef <- paste0(gcol, "_", num, "_vs_", den)
      if (coef %in% DESeq2::resultsNames(dds)) return(try_shrink("apeglm", coef = coef))
      warning(sprintf("apeglm: coef %s が resultsNames に不在 → ashr（contrast 対応）へ切替。", coef))
      return(try_shrink("ashr", contrast = ct))
    }
    warning(sprintf("apeglm は単一係数専用。非参照間対比（den=%s≠ref=%s）は ashr へ分岐。", den, ref))
    return(try_shrink("ashr", contrast = ct))
  }
  if (identical(mode, "ashr")) return(try_shrink("ashr", contrast = ct))
  warning("deseq2.shrink は none/apeglm/ashr のいずれか（指定=", mode, "）→ 収縮なし。")
  list(res = r, shrink_used = "none")
}

# routing: route_contrasts() の返り値。複製あり lane に振られた対比のみ処理する。
run_deseq2 <- function(se, cfg, routing = NULL) {
  coldata <- as.data.frame(SummarizedExperiment::colData(se))
  gcol <- cfg$design$group_col
  covs <- (cfg$design$covariates %||% character(0)); covs <- covs[nzchar(covs)]
  for (v in c(gcol, covs)) coldata[[v]] <- factor(coldata[[v]])
  ref <- cfg$design$reference_level
  if (!is.null(ref)) coldata[[gcol]] <- stats::relevel(coldata[[gcol]], ref = ref)

  # 処理対象の対比: routing があれば replicate lane のみ、無ければ全対比（後方互換）
  if (!is.null(routing)) {
    rep_names <- routing$contrast[routing$lane == "replicate"]
  } else {
    rep_names <- names(cfg$contrasts)
  }
  if (length(rep_names) == 0) {
    return(list(results = list(), method = "DESeq2", dds = NULL, note = "複製あり対比なし（DESeq2 は非適用）"))
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = SummarizedExperiment::assay(se, "counts"),
    colData   = coldata,
    design    = deseq2_formula(coldata, cfg))

  # 低発現 prefilter（任意・速度/メモリ目的。FDR 統制は independent filtering が担う）
  min_sum <- cfg$prefilter$min_count_sum %||% 10
  dds <- dds[rowSums(DESeq2::counts(dds)) >= min_sum, ]

  dds <- DESeq2::DESeq(dds, quiet = TRUE)

  mode <- tolower(cfg$deseq2$shrink %||% "ashr")
  res <- lapply(rep_names, function(nm) {
    spec <- cfg$contrasts[[nm]]; num <- spec[[1]]; den <- spec[[2]]
    ct   <- c(gcol, num, den)
    r    <- DESeq2::results(dds, contrast = ct)                # FDR は independent filtering
    sh   <- .deseq2_shrink(dds, r, ct, num, den, ref, gcol, mode)
    df <- as.data.frame(sh$res)
    df$gene     <- rownames(df)
    df$contrast <- nm
    df$shrink   <- sh$shrink_used                              # 明示ラベル（none/apeglm/ashr）
    df
  })
  names(res) <- rep_names
  list(results = res, method = "DESeq2", dds = dds, contrasts_run = rep_names)
}
