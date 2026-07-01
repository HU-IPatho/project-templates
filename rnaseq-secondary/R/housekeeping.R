# ============================================================================
# housekeeping 遺伝子リストの読込（n=1 DEG の分散推定に使う・config で差替可）
# ----------------------------------------------------------------------------
# 既定は resources/housekeeping/ の同梱リスト:
#   - hk_human.txt : ヒト（Eisenberg & Levanon 2013 の安定 HK パネル + 古典 HK の curated 部分集合）
#   - hk_mouse.txt : マウス（古典 HK の ortholog / 慣用シンボル）
# config.yaml の edger.hk_gene_file を指定すれば任意リスト（例 Eisenberg 全 ~3800 遺伝子）に差替可。
# ファイル形式: 1 行 1 遺伝子シンボル。'#' 始まりの行と空行は無視する。
# ============================================================================

load_housekeeping <- function(cfg) {
  path <- cfg$edger$hk_gene_file %||% ""
  if (!nzchar(path)) {
    path <- switch(
      cfg$organism %||% "human",
      human = here::here("resources", "housekeeping", "hk_human.txt"),
      mouse = here::here("resources", "housekeeping", "hk_mouse.txt"),
      stop("load_housekeeping: 未知の organism: ", cfg$organism,
           "（human / mouse、または edger.hk_gene_file を指定）"))
  } else {
    path <- here::here(path)
  }
  if (!file.exists(path)) stop("load_housekeeping: HK リストが無い: ", path)
  genes <- trimws(readLines(path, warn = FALSE))
  genes <- genes[nzchar(genes) & !startsWith(genes, "#")]
  unique(genes)
}
