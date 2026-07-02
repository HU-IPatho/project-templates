# ============================================================================
# design / contrast の解決 — config.yaml の宣言を各エンジンの入力に変換する
# ----------------------------------------------------------------------------
# config.yaml（事実層）:
#   design:
#     group_col: group            # 主たる比較変数（列名）
#     covariates: [batch]         # DTE モデル共変量（バッチ等・空なら []）
#     reference_level: ctrl       # group の基準レベル
#   contrasts:
#     trt_vs_ctrl: [trt, ctrl]    # [numerator, denominator]（group のレベル名）
# ここは DTE（② の DESeq2/edgeR エンジン）の design/contrast を組む層。DRIMSeq（DTU）は
# 各 contrast の 2 群を subset し ~condition の 2 群比較に落とす（R/dtu.R が担う）。
# バッチは DTE モデル共変量として補正する。
# ============================================================================

# --- 複製ガード（DTU/DTE は各群 n>=2 が必須。n=1 は fail-fast で ② へ誘導）---------
# 05_dtu.R の冒頭で呼ぶ。量の差だけを見たい n=1 は rnaseq-secondary（edgeR n=1 HK 分散）へ。
check_replicates <- function(coldata, cfg) {
  gcol <- cfg$design$group_col
  if (!gcol %in% colnames(coldata)) stop("check_replicates: colData に group 列が無い: ", gcol)
  tab <- table(coldata[[gcol]])
  singletons <- names(tab)[tab < 2]
  if (length(singletons) > 0) {
    stop("DTU/DTE は各群に複製が必要（各群 n>=2）。n=1 の群: ",
         paste(sprintf("%s(n=%d)", singletons, tab[singletons]), collapse = ", "),
         "。量の差だけを見たいなら rnaseq-secondary（② edgeR n=1 HK 分散）を使う。")
  }
  invisible(TRUE)
}

# 使用可能な共変量に絞る（単一レベルの因子は model.matrix でエラーになるため除く）。
resolve_covariates <- function(coldata, cfg) {
  covs <- cfg$design$covariates %||% character(0)
  covs <- covs[nzchar(covs)]
  keep <- vapply(covs, function(v) {
    if (!v %in% colnames(coldata)) stop("design: colData に共変量が無い: ", v)
    nlevels(factor(coldata[[v]])) >= 2
  }, logical(1))
  dropped <- covs[!keep]
  if (length(dropped) > 0) {
    warning("design: 単一レベルの共変量を除外（補正不能）: ",
            paste(dropped, collapse = ", "))
  }
  covs[keep]
}

# --- edgeR 用 design 行列: model.matrix(~ 0 + group [+ covariates]) --------------
# group を means-model（各群 1 係数）にすると、多対比 contrast が「群係数の差」で表せる。
edger_design <- function(coldata, cfg) {
  gcol <- cfg$design$group_col
  covs <- resolve_covariates(coldata, cfg)
  for (v in c(gcol, covs)) {
    if (!v %in% colnames(coldata)) stop("edger_design: colData に列が無い: ", v)
    coldata[[v]] <- factor(coldata[[v]])
  }
  terms <- c(paste0("0 + ", gcol), covs)
  form  <- stats::as.formula(paste("~", paste(terms, collapse = " + ")))
  mm <- stats::model.matrix(form, data = coldata)
  attr(mm, "group_col") <- gcol
  mm
}

# --- edgeR 用 contrast ベクトル: groupNUM - groupDEN（共変量列は 0）---------------
edger_contrasts <- function(design, cfg) {
  gcol <- cfg$design$group_col
  cols <- colnames(design)
  lapply(cfg$contrasts, function(spec) {
    num <- paste0(gcol, spec[[1]])
    den <- paste0(gcol, spec[[2]])
    if (!num %in% cols || !den %in% cols) {
      stop("edger_contrasts: design 列に ", num, " / ", den,
           " が無い（group レベル名と contrast 指定を確認）")
    }
    v <- stats::setNames(numeric(length(cols)), cols)
    v[num] <- 1
    v[den] <- -1
    v
  })
}

# --- DESeq2 用 design formula: ~ covariates + group（対象を最後に置く慣習）-------
deseq2_formula <- function(coldata, cfg) {
  gcol <- cfg$design$group_col
  covs <- resolve_covariates(coldata, cfg)
  stats::as.formula(paste("~", paste(c(covs, gcol), collapse = " + ")))
}
