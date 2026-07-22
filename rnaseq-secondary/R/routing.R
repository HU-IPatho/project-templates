# ============================================================================
# 複製構造による routing — 対比単位の min(group) で screening / replicate lane を決める
# ----------------------------------------------------------------------------
# 標準（specs/bulk-secondary-deg-standard/spec.md）:
#   - 経路分岐は対比単位の最小群サイズ min(table(group)) を判定キーにする（MUST）。
#     二値の all(table(group)==1) を分岐キーにしてはならない（MUST NOT・部分複製デザインの silent 漏れ）。
#   - min(table(group)) < 2 の対比を n=1 スクリーニング経路へ降格する（MUST）。
#   - 複製の独立性（biological/technical/pseudo）を必須メタとし、technical のみの群は biological n=1 扱い
#     （technical 分散は生物変動を計上せず Type I error を膨張させる）。
# ============================================================================

# --- 独立性を反映した「実効」群サイズ ------------------------------------------
# group ごとに、その群が technical のみ（生物学的複製でない）なら実効サイズ 1、そうでなければ実本数。
# config$replicate_independence$default（全群既定）+ per_group（群ごと上書き）で申告する。
#   biological : 真の生物学的複製 → 実本数を計上
#   technical  : 同一検体の再測定のみ → biological n=1 扱い（実効 1）
#   pseudo     : 疑似複製（過度な近接）→ biological n=1 扱い（実効 1・保守側）
effective_group_sizes <- function(group, cfg) {
  group   <- factor(group)
  n_raw   <- table(group)
  levs    <- names(n_raw)
  ind_def <- cfg$replicate_independence$default %||% "biological"
  ind_pg  <- cfg$replicate_independence$per_group %||% list()

  eff <- vapply(levs, function(g) {
    ind <- ind_pg[[g]] %||% ind_def
    if (identical(ind, "biological")) as.integer(n_raw[[g]]) else 1L  # technical/pseudo は実効 1
  }, integer(1))
  names(eff) <- levs
  eff
}

# --- 対比ごとに lane を決める ---------------------------------------------------
# 返り値: data.frame(contrast, num, den, min_group, lane)。lane ∈ {"screening","replicate"}。
#   lane="screening": min(実効群サイズ) < 2（n=1 スクリーニングスタンダードへ降格）
#   lane="replicate": min(実効群サイズ) >= 2（較正済み FDR lane）
route_contrasts <- function(group, cfg) {
  eff <- effective_group_sizes(group, cfg)
  do.call(rbind, lapply(names(cfg$contrasts), function(nm) {
    spec <- cfg$contrasts[[nm]]
    num  <- spec[[1]]; den <- spec[[2]]
    for (g in c(num, den)) {
      if (!g %in% names(eff)) {
        stop(sprintf("routing: contrast %s の群 %s が samples に存在しません。", nm, g))
      }
    }
    mn <- min(eff[[num]], eff[[den]])
    data.frame(contrast = nm, num = num, den = den, min_group = as.integer(mn),
               lane = if (mn < 2) "screening" else "replicate",
               stringsAsFactors = FALSE)
  }))
}
