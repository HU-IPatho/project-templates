# 02: SE → DEG。bulk-secondary-deg-standard（v0・方法論）に conform。
#   (1) scope/provenance ゲート（fail-closed）→ (2) 対比を lane（複製あり / n=1 screening）へ routing →
#   (3) 各 lane で DEG（複製あり=glmQLFTest/DESeq2・screening=HK 分散/BCV バンド・screening-grade）→
#   (4) cross-engine concordance（補助 sanity check・複製あり両走のみ）+ cross-hairpin concordance（G5）。
# 実行: プロジェクトルートから  Rscript analysis/02_de.R
suppressPackageStartupMessages({
  library(here); library(yaml); library(SummarizedExperiment)
})
source(here::here("R", "helpers.R"))       # save_table / %||%
source(here::here("R", "scope.R"))         # check_scope / resolve_bcv_band
source(here::here("R", "routing.R"))       # route_contrasts / effective_group_sizes
source(here::here("R", "de_common.R"))     # tidy_de / compare_methods
source(here::here("R", "design.R"))        # edger_design / edger_design_group_only / contrasts / formula
source(here::here("R", "housekeeping.R"))  # load_housekeeping
source(here::here("R", "screening.R"))     # G1–G5 ゲート・BCV sweep・cross-hairpin・記述解析
source(here::here("R", "de_deseq2.R"))     # run_deseq2
source(here::here("R", "de_edger.R"))      # run_edger

cfg <- yaml::read_yaml(here::here("config.yaml"))
se  <- readRDS(here::here("data", "processed", "se.rds"))
method <- cfg$method %||% "edger"
fdr    <- cfg$fdr %||% 0.05

# --- (1) scope / provenance ゲート（範囲外は fail-closed で停止）---------------------
check_scope(cfg)

# --- (2) 複製構造 routing（対比単位 min(実効群サイズ)）------------------------------
group   <- factor(as.data.frame(SummarizedExperiment::colData(se))[[cfg$design$group_col]])
routing <- route_contrasts(group, cfg)
save_table(routing, "routing", "contrast lane routing (min group)", "analysis/02_de.R")
cat("ROUTING:\n"); print(routing)
has_screening <- any(routing$lane == "screening")
has_replicate <- any(routing$lane == "replicate")

# --- (3) lane 別 DEG --------------------------------------------------------------
# screening 対比は engine 選択に依らず必ず edgeR（HK 分散・n=1 の唯一経路）。
# 複製あり対比は method の engine（edger/deseq2/both）。
edger_replicate <- method %in% c("edger", "both")
de_list <- list()

if (method %in% c("deseq2", "both") && has_replicate) {
  de_list$DESeq2 <- run_deseq2(se, cfg, routing)
}
if (has_screening || edger_replicate) {
  hk <- if (has_screening) load_housekeeping(cfg) else NULL
  de_list$edgeR <- run_edger(se, cfg, routing, hk_genes = hk, run_replicate = edger_replicate)
}
if (length(de_list) == 0) stop("02_de: 実行対象の対比がありません（method/routing を確認）。")

# --- DEG 表を統一列で保存（エンジン×対比ごと）------------------------------------
tidy_all <- list()
for (eng in names(de_list)) {
  if (length(de_list[[eng]]$results) == 0) next
  tid <- tidy_de(de_list[[eng]])
  tidy_all[[eng]] <- tid
  for (nm in names(tid)) {
    lane <- de_list[[eng]]$lane[[nm]] %||% "replicate"
    save_table(tid[[nm]],
               tbl_id = sprintf("deg_%s_%s", tolower(eng), nm),
               desc   = sprintf("DEG %s %s [%s]", eng, nm, lane),
               script = "analysis/02_de.R")
  }
}

# --- screening lane 診断・provenance の保存（screening-grade の来歴と必須ゲート）------
if (!is.null(de_list$edgeR) && length(de_list$edgeR$screening_contrasts) > 0) {
  scr_diag <- de_list$edgeR$diagnostics$screening
  edg <- de_list$edgeR

  # provenance 記録（screening-grade・disp_source・G1/G2 status・複製取得不能理由）
  prov <- data.frame(
    key = c("screening_grade", "disp_source", "global_shift_flag", "global_shift_action",
            "master_regulator_list_present", "master_regulator_forbidden",
            "replication_unavailable_reason", "orthogonal_validation_required",
            "hk_supplied", "hk_validated", "hk_flagged", "hk_cnv_dropped"),
    value = c(edg$screening_grade %||% NA,
              edg$disp_source %||% NA,
              as.character(scr_diag$global_shift$flag),
              scr_diag$global_shift$action,
              as.character(scr_diag$master_regulator$list_present),
              paste(scr_diag$master_regulator$forbidden, collapse = ";"),
              cfg$provenance$replication_unavailable_reason %||% "",
              as.character(cfg$screening$orthogonal_validation$required %||% TRUE),
              scr_diag$hk_validation$n_supplied, scr_diag$hk_validation$n_validated,
              scr_diag$hk_validation$n_flagged, scr_diag$hk_validation$n_cnv_dropped),
    stringsAsFactors = FALSE)
  save_table(prov, "screening_provenance", "n=1 screening provenance & gates", "analysis/02_de.R")

  # BCV 感度スイープ（順位安定性）
  if (!is.null(scr_diag$bcv_sweep)) {
    save_table(scr_diag$bcv_sweep, "bcv_sensitivity", "BCV band sensitivity sweep (rank stability)", "analysis/02_de.R")
  }
  # anti-conservative 方向性バイアス診断（上位ヒットの分散-発現プロファイル）
  if (!is.null(scr_diag$anticonservative)) {
    save_table(scr_diag$anticonservative, "anticonservative_diag",
               "top-hit variance-expression profile vs HK", "analysis/02_de.R")
  }
  # G3: HK 検証メトリクス（経験的 control 集合）
  if (!is.null(scr_diag$hk_validation$metrics)) {
    save_table(scr_diag$hk_validation$metrics, "hk_validation",
               "empirical control / HK validation metrics (G3)", "analysis/02_de.R")
  }
  # G5: cross-hairpin concordance（第一信頼フィルタ・hairpin_map があるとき）
  chc <- cross_hairpin_concordance(tidy_all$edgeR, cfg, fdr = fdr)
  if (!is.null(chc)) {
    save_table(chc, "cross_hairpin", "cross-hairpin concordance confidence (G5)", "analysis/02_de.R")
    cat("CROSS_HAIRPIN (G5):\n"); print(utils::head(chc))
  } else {
    cat("CROSS_HAIRPIN (G5): hairpin_map 未宣言 → concordance フィルタ無効（単一ハーピンは低信頼扱い）。\n")
  }
}

# --- (4) cross-engine concordance（補助 sanity check・複製あり両走のみ・n=1 非適用）----
if (method == "both" && !is.null(tidy_all$DESeq2) && !is.null(tidy_all$edgeR)) {
  agree <- compare_methods(tidy_all$DESeq2, tidy_all$edgeR, fdr = fdr)
  if (!is.null(agree) && nrow(agree) > 0) {
    save_table(agree, "concordance", "cross-engine concordance (sanity check, Spearman/strata)", "analysis/02_de.R")
    cat("CONCORDANCE (sanity check・頑健性保証でない):\n"); print(agree)
  } else {
    cat("CONCORDANCE: 複製あり両走の共通対比なし（n=1 screening は DESeq2 非適用ゆえ除外）。\n")
  }
}

saveRDS(de_list, here::here("data", "processed", "de.rds"))
cat("DE_DONE:", paste(names(de_list), collapse = ","),
    "| method =", method,
    "| screening =", sum(routing$lane == "screening"),
    "| replicate =", sum(routing$lane == "replicate"), "\n")
