# ============================================================================
# scope / provenance ゲート — bulk-secondary-deg-standard v0 の適用範囲を fail-closed で守る
# ----------------------------------------------------------------------------
# 標準（specs/bulk-secondary-deg-standard/spec.md）:
#   - 適用範囲は human/mouse ＋ full-length polyA bulk-RNA-seq に明示限定する。範囲外を silent に
#     処理しない（MUST NOT）。HK バンド・技術分散が未較正のため範囲外は fail-closed で停止する。
#   - config の BCV バンドとデータ provenance（cell_line / patient / technical）の整合を機械チェックし、
#     不整合・欠落は fail-closed で停止する（MUST）。
#
# ここは「形（scope 宣言・provenance 整合）」を機械強制するゲート。範囲外プロトコルへ適用したいときの
# 再較正の科学判断そのものは規定しない（TBD・grill ゲート）。数値（バンドの中身）は config 例示。
# ============================================================================

# v0 の scope（規範）。organism は HK リストが在る種、library_prep は full-length polyA bulk に限る。
.SCOPE_ORGANISMS   <- c("human", "mouse")
.SCOPE_LIBRARY_PREP <- c("polya_fulllength")  # full-length polyA bulk（total/3'-tag/UMI/FFPE は範囲外）

# --- 適用範囲チェック（範囲外は fail-closed で停止）------------------------------
# organism / library_prep が v0 scope 内か検証する。範囲外は明示 warning ＋ provenance を促し stop。
check_scope <- function(cfg) {
  organism <- tolower(cfg$organism %||% "")
  libprep  <- tolower(cfg$library_prep %||% "")

  if (!nzchar(organism)) {
    stop("scope: config$organism 未指定。v0 scope は human/mouse（HK リストが在る種）。")
  }
  if (!organism %in% .SCOPE_ORGANISMS) {
    stop(sprintf(paste0(
      "scope: organism=%s は v0 の適用範囲外です（scope = %s）。HK リスト不在の種は BCV バンドが",
      " 未較正のため silent フォールバックせず fail-closed で停止します。適用するには当該種の HK パネルと",
      " BCV バンドを再較正し（grill ゲートで確定）scope を拡張してください。"),
      organism, paste(.SCOPE_ORGANISMS, collapse = "/")))
  }

  if (!nzchar(libprep)) {
    stop("scope: config$library_prep 未指定（必須メタ）。v0 scope は full-length polyA bulk = polya_fulllength。")
  }
  if (!libprep %in% .SCOPE_LIBRARY_PREP) {
    stop(sprintf(paste0(
      "scope: library_prep=%s は v0 の適用範囲外です（scope = %s = full-length polyA bulk）。",
      " 3'-tag / UMI / 低入力 / FFPE（劣化 RNA）は count 分布・技術分散が異なり BCV バンドと HK 特性が",
      " 未較正のため「要 recalibration」で fail-closed です。silent に適用しません。"),
      libprep, paste(.SCOPE_LIBRARY_PREP, collapse = "/")))
  }
  invisible(TRUE)
}

# --- BCV バンドと provenance の整合チェック（不整合・欠落は fail-closed）-----------
# config$bcv$band に config$provenance$source（cell_line/patient/technical）に対応するバンドが
# 在るかを検証する。無ければ停止（どの分散帯で感度スイープすべきか未確定なため）。
# 返り値: 該当 provenance の BCV バンド（数値ベクトル・昇順）。
resolve_bcv_band <- function(cfg) {
  source <- cfg$provenance$source %||% ""
  if (!nzchar(source)) {
    stop("bcv: config$provenance$source 未指定（cell_line/patient/technical）。BCV バンド整合に必須。")
  }
  bands <- cfg$bcv$band
  if (is.null(bands) || is.null(bands[[source]])) {
    stop(sprintf(paste0(
      "bcv: provenance.source=%s に対応する BCV バンド（config$bcv$band$%s）が未定義です。",
      " 単一固定 BCV の確定運用は禁止（最感度パラメータ）。当該 provenance の分散帯をバンドで宣言してください",
      "（数値は自分のデータ由来で確定・grill ゲート）。"),
      source, source))
  }
  band <- suppressWarnings(as.numeric(bands[[source]]))
  band <- band[is.finite(band) & band > 0]
  if (length(band) < 1) {
    stop(sprintf("bcv: config$bcv$band$%s に有効な BCV 値がありません（正の数値を列挙）。", source))
  }
  sort(unique(band))
}
