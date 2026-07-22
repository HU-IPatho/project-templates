# ============================================================================
# edgeR パス — 対比を lane（複製あり / n=1 screening）別に routing して DEG
# ----------------------------------------------------------------------------
# 標準（specs/bulk-secondary-deg-standard/spec.md）に conform:
#   - 複製あり lane（min(実効群サイズ)>=2）: filterByExpr → TMM → estimateDisp(design) →
#     glmQLFit + glmQLFTest を既定（MUST）。glmLRT は意図的逸脱（理由必須）のみ。
#   - n=1 screening lane（min<2）: HK 検証（G3）→ HK common.dispersion（BCV バンドより優先）or
#     BCV バンド代表値 → glmFit/glmLRT で screening-grade。G1 大域シフトで破棄/記述降格。
#     記述解析（FC ランキング・dispersion-free）を必ず併走。BCV 感度スイープで順位安定性を出す。
#   - バッチ加法補正は複製あり lane 限定（screening lane は素の ~0+group）。
# 各 lane 固有の gate/診断ロジックは R/screening.R。design/contrast 解決は R/design.R。
# ============================================================================
suppressPackageStartupMessages({
  library(edgeR)
  library(SummarizedExperiment)
})

# routing: route_contrasts() の返り値（contrast × lane）。hk_genes: 供給 HK 候補（要データ内検証）。
# run_replicate=FALSE で複製あり lane を edgeR で走らせない（method=deseq2 で複製対比を DESeq2 が担い
#   edgeR は screening 対比のみ担うケース。screening lane は n=1 の唯一経路ゆえ常に走る）。
run_edger <- function(se, cfg, routing, hk_genes = NULL, run_replicate = TRUE) {
  counts  <- SummarizedExperiment::assay(se, "counts")
  coldata <- as.data.frame(SummarizedExperiment::colData(se))
  gcol    <- cfg$design$group_col
  group   <- factor(coldata[[gcol]])

  y <- edgeR::DGEList(counts = counts, group = group)
  design_rep <- edger_design(coldata, cfg)              # 複製あり lane（共変量=batch 込み）
  design_scr <- edger_design_group_only(coldata, cfg)   # screening lane（素の ~0+group）

  keep <- edgeR::filterByExpr(y, design = design_scr)    # 群のみの design で共通 filter（保守側）
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y, method = "TMM")

  scr_names <- routing$contrast[routing$lane == "screening"]
  rep_names <- if (isTRUE(run_replicate)) routing$contrast[routing$lane == "replicate"] else character(0)

  results <- list(); diagnostics <- list(); disp_source <- NULL; screening_grade <- NULL

  if (length(rep_names) > 0) {
    rr <- .edger_replicate_lane(y, design_rep, cfg, rep_names)
    results <- c(results, rr$results); diagnostics$replicate <- rr$diag
  }
  if (length(scr_names) > 0) {
    ss <- .edger_screening_lane(y, design_scr, cfg, scr_names, hk_genes)
    results <- c(results, ss$results); diagnostics$screening <- ss$diag
    disp_source <- ss$disp_source; screening_grade <- ss$screening_grade
  }

  results <- results[intersect(names(cfg$contrasts), names(results))]
  list(results         = results,
       method          = "edgeR",
       lane            = stats::setNames(routing$lane, routing$contrast),
       screening_contrasts = scr_names,
       replicate_contrasts = rep_names,
       screening_grade = screening_grade,
       disp_source     = disp_source,
       diagnostics     = diagnostics)
}

# --- 複製あり lane（較正済み FDR・glmQLFTest 既定）------------------------------------
.edger_replicate_lane <- function(y, design, cfg, contrasts_names) {
  y <- edgeR::estimateDisp(y, design)
  test <- toupper(cfg$edger$test %||% "QL")
  contrasts <- edger_contrasts(design, cfg)[contrasts_names]

  if (test == "QL") {
    fit <- edgeR::glmQLFit(y, design)
    res <- lapply(names(contrasts), function(nm) {
      qlf <- edgeR::glmQLFTest(fit, contrast = contrasts[[nm]])
      tt  <- edgeR::topTags(qlf, n = Inf, sort.by = "none")$table
      tt$gene <- rownames(tt); tt$contrast <- nm; tt
    })
    disp <- "estimateDisp + glmQLFit/glmQLFTest（複製あり既定・very reliable FDR control）"
  } else if (test == "LRT") {
    reason <- cfg$edger$test_deviation_reason %||% ""
    if (!nzchar(reason)) {
      stop("edger.test=LRT は複製ありでは非推奨の意図的逸脱です。edger.test_deviation_reason に ",
           "「なぜ QL でなく LRT か」の理由と適用条件を記載してください（標準 SHALL）。")
    }
    warning("複製あり lane で glmLRT を使用（traditional・非推奨の意図的逸脱）: ", reason)
    fit <- edgeR::glmFit(y, design)
    res <- lapply(names(contrasts), function(nm) {
      lrt <- edgeR::glmLRT(fit, contrast = contrasts[[nm]])
      tt  <- edgeR::topTags(lrt, n = Inf, sort.by = "none")$table
      tt$gene <- rownames(tt); tt$contrast <- nm; tt
    })
    disp <- sprintf("estimateDisp + glmFit/glmLRT（意図的逸脱: %s）", reason)
  } else {
    stop("edger.test は QL / LRT のいずれか（指定=", test, "）")
  }
  names(res) <- names(contrasts)
  list(results = res, diag = list(test = test, disp_source = disp,
                                  common.dispersion = y$common.dispersion))
}

# --- n=1 screening lane（screening-grade・ゲート G1–G5・BCV スイープ・記述併走）--------
.edger_screening_lane <- function(y, design_scr, cfg, contrasts_names, hk_genes) {
  fdr <- cfg$fdr %||% 0.05

  # G2: 既知 master regulator 短リスト（適用禁止 + リスト外 warning）
  mr <- master_regulator_status(cfg)
  if (length(mr$forbidden) > 0) stop(mr$message)      # HK-dispersion 経路 適用禁止
  if (!mr$list_present) warning(mr$message)

  # G4: seed off-target スクリーン（enabled + 未実装は fail-closed stop）
  message(seed_offtarget_gate(cfg))

  # G3: 供給 HK のデータ内検証 + CNV アーム除外（無検証使用禁止）
  hkv <- validate_hk(hk_genes, y, cfg)
  validated_hk <- hkv$validated
  if (length(hkv$flagged) > 0) {
    warning(sprintf("G3: 供給 HK %d 個が当該データで応答的（経験的 control 基準外）→ 除外。検証済 HK=%d。",
                    length(hkv$flagged), length(validated_hk)))
  }
  if (length(hkv$cnv_dropped) > 0) {
    warning(sprintf("G3: CNV アーム上の HK 候補 %d 個を除外。", length(hkv$cnv_dropped)))
  }

  # G1: 大域シフトゲート（hard）
  gs <- global_shift_gate(y, validated_hk, cfg)

  band <- resolve_bcv_band(cfg)
  disp <- estimate_screening_dispersion(y, validated_hk, cfg, band)
  contrasts <- edger_contrasts(design_scr, cfg)[contrasts_names]

  # 記述解析（dispersion-free FC ランキング）— 必ず併走
  desc <- descriptive_tables(y, cfg, contrasts_names)

  if (isTRUE(gs$flag)) {
    # G1 破綻 → screening DEG を破棄し記述解析へ降格（有意性主張なし）
    warning("G1 大域シフト検出 → screening DEG を破棄し記述解析（FC ランキング）へ降格。", gs$detail)
    res <- desc
    grade <- "descriptive_degraded"
    sweep <- NULL
  } else {
    y2 <- y
    y2$common.dispersion  <- disp$common.dispersion   # 全遺伝子へ移植（trended/tagwise は付けない）
    y2$trended.dispersion <- NULL
    y2$tagwise.dispersion <- NULL
    fit <- edgeR::glmFit(y2, design_scr)
    res <- lapply(names(contrasts), function(nm) {
      lrt <- edgeR::glmLRT(fit, contrast = contrasts[[nm]])
      tt  <- edgeR::topTags(lrt, n = Inf, sort.by = "none")$table
      tt$gene <- rownames(tt); tt$contrast <- nm; tt
    })
    names(res) <- names(contrasts)
    grade <- "screening_grade"
    sweep <- bcv_sensitivity_sweep(y, design_scr, contrasts, band, cfg, fdr = fdr)
  }

  warning(sprintf(paste0(
    "n=1 screening（%s）: 出力は screening-grade（候補絞り・仮説生成）。較正済み FDR/有意性を",
    " 確定的主張として扱わない。記述解析（FC ランキング + MDS）を併走出力。anti-conservative 方向性",
    " バイアス（HK の低分散を高分散遺伝子へ移植 → 偽陽性が興味深い変動遺伝子に集中）に注意。"),
    disp$source))

  anti <- anticonservative_diagnostic(y, res, validated_hk)

  list(results = res, disp_source = disp$source, screening_grade = grade,
       diag = list(global_shift = gs, master_regulator = mr,
                   hk_validation = list(n_supplied = hkv$n_supplied,
                                        n_validated = length(validated_hk),
                                        n_flagged = length(hkv$flagged),
                                        n_cnv_dropped = length(hkv$cnv_dropped),
                                        metrics = hkv$metrics),
                   bcv_band = band, bcv_sweep = sweep, descriptive = desc,
                   anticonservative = anti))
}
