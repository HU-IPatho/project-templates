# ============================================================================
# n=1 KD スクリーニングスタンダード — screening-grade 経路とその必須ゲート G1–G6
# ----------------------------------------------------------------------------
# 標準（specs/bulk-secondary-deg-standard/spec.md）に conform する scaffold。
# 位置づけ（標準）: 複製取得不能時の documented degraded path。出力は screening-grade（候補絞り・
#   仮説生成）であり較正済み FDR を確定的主張として出さない。ヒットは fold-change / 収縮 LFC の
#   ランキングに限る。記述解析（FC ランキング + MDS）を必ず併走出力する（edgeR §2.13 option 1）。
#
# ここは「方法論の構造」（routing 済 screening lane・ゲート G1–G6・BCV バンドスイープ）を用意し、
# domain/dataset 固有の numeric（HK membership・BCV バンド値・各種閾値・G2 リスト）は config + /grill-me
# に委ねる（数値をここに固定しない）。★n=1 HK/固定 BCV 経路の妥当性は biology-conditional であり
# （大域シフト無し・当該摂動下で HK 非DE の precondition 下でのみ成立）、無条件の組織/遺伝子非依存では
# ない——ゆえに下記ゲートで前提破綻を検出し破棄/降格する。
suppressPackageStartupMessages({
  library(edgeR)
})

# --- 小さな遺伝子リストローダ（1 行 1 シンボル・# と空行は無視・空パスは空ベクトル）--------
load_gene_list <- function(path) {
  path <- path %||% ""
  if (!nzchar(path)) return(character(0))
  full <- here::here(path)
  if (!file.exists(full)) stop("screening: 遺伝子リストが無い: ", full)
  g <- trimws(readLines(full, warn = FALSE))
  unique(g[nzchar(g) & !startsWith(g, "#")])
}

# --- 群平均 logCPM と群間 |logFC|・発現・CV（分散に依存しない記述量）--------------
# 経験的 control 導出（G3）と大域シフトゲート（G1）の共通材料。dispersion を要しない点が要点
# （logFC の点推定は正規化済み群平均の差＝dispersion 非依存）。
.gene_descriptives <- function(y) {
  logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 2)
  grp    <- as.character(y$samples$group)
  levs   <- unique(grp)
  gmean  <- sapply(levs, function(g) rowMeans(logcpm[, grp == g, drop = FALSE]))
  colnames(gmean) <- levs
  mean_expr <- rowMeans(logcpm)
  # CV は線形 CPM 尺度で（低分散 = 安定発現）。
  cpm_lin <- edgeR::cpm(y, log = FALSE)
  mu  <- rowMeans(cpm_lin)
  sdv <- apply(cpm_lin, 1, stats::sd)
  cv  <- ifelse(mu > 0, sdv / mu, NA_real_)
  list(logcpm = logcpm, group_mean = gmean, mean_expr = mean_expr, cv = cv)
}

# --- G3: データ内経験的 control 集合の導出（低 |群間 logFC|・高発現・低 CV）-----------
# 汎用 HK リストの無検証使用を禁じ（MUST NOT）、当該データから経験的 control 集合を導出する（MUST）。
# 閾値は分位（config 例示・TBD）。返り値: list(controls, metrics)。
derive_empirical_controls <- function(y, cfg) {
  # 下の %||% fallback は config.yaml の「(例)」既定を鏡写しにした非規範のセーフガード（規範値でない・
  # 確定は config 事実層 / grill ゲート）。config キーが在れば常にそちらが権威。
  ec <- cfg$screening$empirical_control %||% list()
  q_logfc <- ec$max_abs_logfc_quantile %||% 0.5   # 非規範 fallback（config 例示の鏡写し）
  q_expr  <- ec$min_expr_quantile      %||% 0.5
  q_cv    <- ec$max_cv_quantile        %||% 0.5

  d <- .gene_descriptives(y)
  # 全群ペアの最大 |logFC|（どの群間でも動かない遺伝子を control とする）
  gm <- d$group_mean
  max_abs_logfc <- if (ncol(gm) >= 2) {
    apply(gm, 1, function(r) max(abs(outer(r, r, "-")), na.rm = TRUE))
  } else rep(0, nrow(gm))

  thr_logfc <- stats::quantile(max_abs_logfc, probs = q_logfc, na.rm = TRUE)
  thr_expr  <- stats::quantile(d$mean_expr,   probs = q_expr,  na.rm = TRUE)
  thr_cv    <- stats::quantile(d$cv,          probs = q_cv,    na.rm = TRUE)

  is_ctrl <- !is.na(max_abs_logfc) & max_abs_logfc <= thr_logfc &
             !is.na(d$mean_expr)   & d$mean_expr   >= thr_expr &
             !is.na(d$cv)          & d$cv          <= thr_cv
  genes <- rownames(y)[is_ctrl]
  metrics <- data.frame(gene = rownames(y), max_abs_logfc = max_abs_logfc,
                        mean_expr = d$mean_expr, cv = d$cv,
                        empirical_control = is_ctrl, stringsAsFactors = FALSE)
  list(controls = genes, metrics = metrics)
}

# --- G3: 供給 HK リストのデータ内検証（無検証使用禁止）+ CNV アーム除外 --------------
# 供給 HK（同梱 or config 指定）のうち、当該データで経験的 control 基準を満たすものだけを「検証済 HK」
# として採用する。基準を満たさない HK（当該データで応答的）は flag して除外する（silent に使わない）。
# 既知 CNV アーム上の遺伝子（config$screening$cnv_excluded_arms_file）も HK 候補から除外する。
validate_hk <- function(hk_genes, y, cfg) {
  present <- intersect(hk_genes %||% character(0), rownames(y))
  cnv_excl <- load_gene_list(cfg$screening$cnv_excluded_arms_file)
  present_no_cnv <- setdiff(present, cnv_excl)
  cnv_dropped <- intersect(present, cnv_excl)

  emp <- derive_empirical_controls(y, cfg)
  validated <- intersect(present_no_cnv, emp$controls)   # データ内で非応答が確認された HK のみ採用
  flagged   <- setdiff(present_no_cnv, emp$controls)      # 供給 HK だが当該データで応答的 → 除外

  list(validated = validated, flagged = flagged, cnv_dropped = cnv_dropped,
       n_supplied = length(present), emp_controls = emp$controls, metrics = emp$metrics)
}

# --- G1: 大域シフトゲート（hard）------------------------------------------------
# HK 群の群間 logFC 分布を診断する。HK 中央値が 0 から有意に外れる／HK 分散が全遺伝子中央に匹敵する
# なら「HK 前提破綻・dispersion 信頼不能」を flag し、DEG 出力を破棄 or 記述解析へ降格する。
# 閾値は config 例示（TBD）。返り値: list(flag, hk_median_abs_logfc, hk_var_ratio, action, detail)。
global_shift_gate <- function(y, hk_present, cfg) {
  gs <- cfg$screening$global_shift %||% list()
  thr_med <- gs$hk_median_abs_logfc_max %||% 0.5
  thr_var <- gs$hk_var_ratio_max        %||% 1.0

  d  <- .gene_descriptives(y)
  gm <- d$group_mean
  hk_present <- intersect(hk_present %||% character(0), rownames(gm))
  if (length(hk_present) < 2 || ncol(gm) < 2) {
    return(list(flag = TRUE, hk_median_abs_logfc = NA_real_, hk_var_ratio = NA_real_,
                action = "degrade_to_descriptive",
                detail = "HK 群が不足（<2）で大域シフトを評価できない → 記述解析へ降格（安全側）。"))
  }
  # (a) 系統的シフト: HK 群の全群ペア logFC（符号付き）の中央値が 0 から外れるか
  pair_idx <- utils::combn(ncol(gm), 2)
  hk_logfc <- unlist(lapply(seq_len(ncol(pair_idx)), function(k) {
    gm[hk_present, pair_idx[1, k]] - gm[hk_present, pair_idx[2, k]]
  }))
  med_signed <- stats::median(hk_logfc, na.rm = TRUE)
  med_abs    <- abs(med_signed)                    # 系統的シフトの大きさ
  # (b) 分散比: per-gene 群間分散の HK 中央値 / 全遺伝子中央値（HK が全遺伝子中央に匹敵するか）
  per_gene_bg_var <- apply(gm, 1, stats::var)
  hk_var  <- stats::median(per_gene_bg_var[hk_present], na.rm = TRUE)
  all_var <- stats::median(per_gene_bg_var, na.rm = TRUE)
  var_ratio <- if (is.finite(all_var) && all_var > 0) hk_var / all_var else NA_real_

  breached <- (is.finite(med_abs) && med_abs > thr_med) ||
              (is.finite(var_ratio) && var_ratio > thr_var)
  list(flag = breached,
       hk_median_abs_logfc = med_abs, hk_var_ratio = var_ratio,
       action = if (breached) "discard_or_degrade" else "ok",
       detail = if (breached)
         "HK 前提破綻の疑い（大域シフト）。HK-dispersion と TMM/median-of-ratios 正規化は同一前提を共有し独立でないため、DEG 出力を破棄するか記述解析（FC+MDS）へ降格する。"
       else "HK 大域シフトの兆候なし（前提は当該データで棄却されず）。")
}

# --- G2: 既知グローバル制御因子（master regulator）の除外（適用禁止 + リスト外 warning）----
# KD 標的が既知 master regulator 短リストに載れば HK 経路 適用禁止（spike-in/直交検証要求）。
# リスト外は warning。リスト未整備（空）は「常に warning」（網羅不能ゆえ短リスト + warning の二段）。
# 標的宣言は hairpin_map のキー（target→hairpins）から取る。返り値: list(list_present, targets,
#   forbidden, message)。forbidden が非空なら HK 経路を適用してはならない（呼出側が停止）。
master_regulator_status <- function(cfg) {
  mr <- load_gene_list(cfg$screening$master_regulator_file)
  targets <- names(cfg$screening$hairpin_map %||% list())   # KD 標的は hairpin_map のキーから申告
  if (length(mr) == 0) {
    return(list(list_present = FALSE, undeclared = (length(targets) == 0), targets = targets,
                forbidden = character(0),
                message = "G2: master regulator 短リスト未整備（config$screening$master_regulator_file 空）。既知グローバル制御因子の KD では HK-dispersion 経路が破綻しうる（例 MYC）。リスト整備までは warning に留め機械強制しない。"))
  }
  if (length(targets) == 0) {
    # リスト在るが KD 標的が hairpin_map で未申告 → G2 を評価できない（silent 通過させない）
    return(list(list_present = TRUE, undeclared = TRUE, targets = targets, forbidden = character(0),
                message = "G2: master regulator 短リストは在りますが KD 標的が未申告（hairpin_map 空）で G2 を評価できません。各 KD 標的を hairpin_map（target→hairpins）に宣言してください。"))
  }
  forbidden <- intersect(targets, mr)
  msg <- if (length(forbidden) > 0)
    sprintf("G2: KD 標的 %s は master regulator 短リストに載る → HK-dispersion 経路 適用禁止。スパイクイン（ERCC 等）正規化 and/or 直交検証を要求する。", paste(forbidden, collapse = ", "))
  else
    sprintf("G2: 宣言された KD 標的（%s）は master regulator 短リストに無い（リスト外は warning・網羅は不可能）。", paste(targets, collapse = ", "))
  list(list_present = TRUE, undeclared = FALSE, targets = targets, forbidden = forbidden, message = msg)
}

# --- G4: shRNA seed off-target スクリーン（既定 ON 推奨・要ツール）------------------
# enabled + 未実装フックは fail-closed で停止（silent に off-target を見逃さない）。既定 disabled では note。
seed_offtarget_gate <- function(cfg) {
  if (!isTRUE(cfg$screening$seed_offtarget$enabled)) {
    return("G4: seed off-target スクリーン off（既定）。shRNA seed 介在 off-target は数百遺伝子規模の広域変動を起こしうる。ON 推奨（逸脱時は理由記録）。")
  }
  # tool-agnostic 未実装フック: SeedMatchR 等での seed off-target 分離は dataset/hairpin-specific ゆえ
  # 研究員が /grill-me で本フックを完成させる。scaffold は捏造せず明示停止する。
  stop("G4: seed_offtarget.enabled=true ですが未実装フックです。研究員が R/screening.R の seed_offtarget_gate() を完成させてください（SeedMatchR 等で on-target 下流とハーピン seed off-target を分離する診断・要 renv 追加）。手順は /grill-me（grill ゲート）で確定する。")
}

# --- G5: cross-hairpin concordance（第一信頼フィルタ）-------------------------------
# hairpin_map（target→hairpins）が在れば、target ごとに複数 hairpin の対比で「符号一致かつ両者候補閾値超」
# を第一信頼に使う。単一ハーピンのみのヒットは低信頼ダウンランク。tidy 済 per-contrast 表（logFC・padj を持つ）
# のリストを受ける。返り値: data.frame(gene, target, confidence, n_hairpins, concordant)。
# 候補閾値は screening-grade ゆえ較正済み FDR でなく |logFC| ランキング（上位 candidate_rank_top）を使う
# （標準 SHALL: n=1 のヒットは FC/収縮 LFC ランキングに限定）。fdr 引数は後方互換で受けるが選抜には使わない。
cross_hairpin_concordance <- function(tidy_list, cfg, fdr = 0.05) {
  hmap <- cfg$screening$hairpin_map %||% list()
  if (length(hmap) == 0) return(NULL)
  cand_top <- cfg$screening$candidate_rank_top %||% 100L
  ref <- cfg$design$reference_level

  # numerator(group) → その群を numerator とする対比名（複数可・後勝ち上書きを避けリスト集約）。
  num_to_contrast <- list()
  for (nm in names(cfg$contrasts)) {
    num <- cfg$contrasts[[nm]][[1]]
    num_to_contrast[[num]] <- c(num_to_contrast[[num]], nm)
  }
  # hairpin(group) の代表対比を選ぶ: KD-vs-対照（den==reference）を優先し、無ければ最初。
  pick_contrast <- function(h) {
    cs <- num_to_contrast[[h]]
    if (is.null(cs)) return(NULL)
    ref_c <- cs[vapply(cs, function(c) identical(cfg$contrasts[[c]][[2]], ref), logical(1))]
    if (length(ref_c) > 0) ref_c[1] else cs[1]
  }
  # 対比の候補集合（|logFC| 上位 cand_top）
  candidates_of <- function(cn) {
    d <- tidy_list[[cn]]; d <- d[is.finite(d$logFC), ]; d <- d[order(-abs(d$logFC)), ]
    utils::head(d$gene, min(cand_top, nrow(d)))
  }

  out <- list()
  for (tgt in names(hmap)) {
    cnames <- unique(unlist(lapply(hmap[[tgt]], pick_contrast)))
    cnames <- intersect(cnames, names(tidy_list))
    if (length(cnames) == 0) next
    if (length(cnames) < 2) {
      # 単一ハーピン → 低信頼ダウンランク（直交検証 G6 を要求）。候補は |logFC| 上位。
      cand <- candidates_of(cnames[1])
      if (length(cand) > 0)
        out[[tgt]] <- data.frame(gene = cand, target = tgt, confidence = "low_single_hairpin",
                                 n_hairpins = 1L, concordant = NA, stringsAsFactors = FALSE)
      next
    }
    # 複数 hairpin: 符号一致 かつ 全対比で |logFC| 候補入り を高信頼 concordant に。
    mats <- lapply(cnames, function(cn) {
      d <- tidy_list[[cn]]; cand <- candidates_of(cn)
      data.frame(gene = d$gene, logFC = d$logFC, is_cand = d$gene %in% cand,
                 stringsAsFactors = FALSE)
    })
    m <- Reduce(function(a, b) merge(a, b, by = "gene"), mats)
    fc_cols <- grep("^logFC",  names(m), value = TRUE)
    cd_cols <- grep("^is_cand", names(m), value = TRUE)
    sign_ok  <- apply(sign(m[, fc_cols, drop = FALSE]), 1,
                      function(s) all(!is.na(s)) && length(unique(s)) == 1)
    all_cand <- apply(m[, cd_cols, drop = FALSE], 1, all)   # 両者で候補閾値超（|logFC| 上位）
    any_cand <- apply(m[, cd_cols, drop = FALSE], 1, any)
    concordant <- sign_ok & all_cand
    conf <- ifelse(concordant, "high_cross_hairpin",
                   ifelse(any_cand, "low_discordant", "not_candidate"))
    keep <- conf != "not_candidate"
    if (any(keep))
      out[[tgt]] <- data.frame(gene = m$gene[keep], target = tgt, confidence = conf[keep],
                               n_hairpins = length(cnames), concordant = concordant[keep],
                               stringsAsFactors = FALSE)
  }
  if (length(out) == 0) return(NULL)
  do.call(rbind, out)
}

# --- BCV 感度スイープ（バンド全域でヒット順位安定性を出す）--------------------------
# 固定単一 BCV を確定運用しない（最感度パラメータ）。バンドの各 BCV で fit→test し、screening 対比ごとに
# 上位 K 集合の重なり（Jaccard）と有意集合サイズの変動を出す。返り値: data.frame（contrast × bcv の要約）。
bcv_sensitivity_sweep <- function(y, design, contrasts, band, cfg, fdr = 0.05, top_k = 100L) {
  top_k <- cfg$screening$candidate_rank_top %||% top_k
  rows <- list()
  # 各 BCV での上位 K 集合を貯める（対比別）
  topsets <- lapply(names(contrasts), function(x) list())
  names(topsets) <- names(contrasts)
  for (b in band) {
    y2 <- y
    y2$common.dispersion  <- b^2
    y2$trended.dispersion <- NULL
    y2$tagwise.dispersion <- NULL
    fit <- edgeR::glmFit(y2, design)
    for (nm in names(contrasts)) {
      lrt <- edgeR::glmLRT(fit, contrast = contrasts[[nm]])
      tt  <- edgeR::topTags(lrt, n = Inf, sort.by = "PValue")$table
      topg <- utils::head(rownames(tt), min(top_k, nrow(tt)))
      topsets[[nm]][[as.character(b)]] <- topg
      n_sig <- sum(tt$FDR < fdr, na.rm = TRUE)
      rows[[length(rows) + 1]] <- data.frame(contrast = nm, bcv = b,
                                             n_sig = n_sig, top_k = length(topg),
                                             stringsAsFactors = FALSE)
    }
  }
  summ <- do.call(rbind, rows)
  # 対比ごとにバンド全域の top-K 平均ペア Jaccard（順位安定性）
  stab <- do.call(rbind, lapply(names(topsets), function(nm) {
    sets <- topsets[[nm]]
    if (length(sets) < 2) return(data.frame(contrast = nm, mean_topk_jaccard = NA_real_,
                                            stringsAsFactors = FALSE))
    js <- c()
    for (i in seq_len(length(sets) - 1)) for (j in (i + 1):length(sets)) {
      a <- sets[[i]]; bb <- sets[[j]]
      u <- length(union(a, bb)); js <- c(js, if (u > 0) length(intersect(a, bb)) / u else NA_real_)
    }
    data.frame(contrast = nm, mean_topk_jaccard = mean(js, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  merge(summ, stab, by = "contrast", all.x = TRUE)
}

# --- screening 用 dispersion 推定（HK 優先・不足時は BCV バンド代表値）----------------
# HK 由来 common.dispersion 推定が成功したら固定 BCV バンドより優先する（標準・降格順）。
# 検証済 HK が min_hk_genes 未満なら BCV バンドの代表値（中央値）へフォールバック（感度は sweep が示す）。
estimate_screening_dispersion <- function(y, validated_hk, cfg, band) {
  min_hk  <- cfg$screening$min_hk_genes %||% 10
  present <- intersect(validated_hk, rownames(y))
  if (length(present) >= min_hk) {
    y_hk <- y[present, , keep.lib.sizes = FALSE]
    y_hk$samples$group <- factor(rep(1L, ncol(y_hk)))       # 全サンプルを 1 群扱い
    est <- tryCatch(
      edgeR::estimateDisp(y_hk, design = NULL, trend.method = "none", tagwise = FALSE),
      error = function(e) NULL)
    if (!is.null(est) && is.finite(est$common.dispersion) && est$common.dispersion > 0) {
      return(list(common.dispersion = est$common.dispersion,
                  source = sprintf("HK common.dispersion（検証済 HK %d・BCV=%.3f・固定 BCV より優先）",
                                   length(present), sqrt(est$common.dispersion))))
    }
  }
  bcv_rep <- stats::median(band)
  list(common.dispersion = bcv_rep^2,
       source = sprintf("BCV バンド代表値=%.3f（検証済 HK %d < 必要 %d ゆえフォールバック・感度は BCV sweep 参照）",
                        bcv_rep, length(present), min_hk))
}

# --- 記述解析テーブル（dispersion-free の FC ランキング・n=1 で必ず併走）---------------
# 群平均 logCPM 差＝dispersion 非依存の logFC。有意性主張なし（PValue/FDR=NA）。edgeR §2.13 option 1。
# 返り値: named list（contrast → edgeR 互換列を持つ df）。
descriptive_tables <- function(y, cfg, contrasts_names) {
  d  <- .gene_descriptives(y)
  gm <- d$group_mean
  res <- lapply(contrasts_names, function(nm) {
    spec <- cfg$contrasts[[nm]]; num <- spec[[1]]; den <- spec[[2]]
    lfc  <- gm[, num] - gm[, den]
    data.frame(logFC = lfc, logCPM = d$mean_expr, PValue = NA_real_, FDR = NA_real_,
               gene = rownames(y), contrast = nm, stringsAsFactors = FALSE)
  })
  names(res) <- contrasts_names
  res
}

# --- anti-conservative 方向性バイアス診断 -------------------------------------------
# HK common.dispersion を全遺伝子へ移植すると本来高分散の遺伝子の分散が過小評価され偽陽性が
# 生物学的に興味深い変動遺伝子に集中する。上位ヒットが低分散 HK 近傍に偏っていないか（分散-発現プロファイル）
# を surface する。返り値: data.frame(gene, contrast, mean_expr, cv, near_hk_low_var)。
anticonservative_diagnostic <- function(y, tidy_list, validated_hk, top_k = 50L) {
  d <- .gene_descriptives(y)
  hk_cv <- d$cv[intersect(validated_hk, rownames(y))]
  hk_cv_med <- if (length(hk_cv) > 0) stats::median(hk_cv, na.rm = TRUE) else NA_real_
  out <- lapply(names(tidy_list), function(nm) {
    tt <- tidy_list[[nm]]
    tt <- tt[order(-abs(tt$logFC)), ]
    topg <- utils::head(tt$gene, min(top_k, nrow(tt)))
    data.frame(gene = topg, contrast = nm,
               mean_expr = d$mean_expr[topg], cv = d$cv[topg],
               near_hk_low_var = !is.na(d$cv[topg]) & !is.na(hk_cv_med) & d$cv[topg] <= hk_cv_med,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}
