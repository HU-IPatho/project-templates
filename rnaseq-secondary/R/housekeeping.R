# ============================================================================
# housekeeping 候補リストの読込（n=1 screening lane の分散推定「候補」・要データ内検証）
# ----------------------------------------------------------------------------
# ★標準（G3）: 汎用 HK リストの無検証使用は禁止。ここで読むのは「候補」であり、当該データで
#   非応答性を経験的に検証（低 |群間 logFC|・高発現・低 CV）してから採用する（検証は R/screening.R の
#   validate_hk）。樹立細胞株は異数体/CNV で増幅・欠失アーム上の HK が無効化しうる。
# 既定候補は resources/housekeeping/ の同梱リスト:
#   - hk_human.txt : ヒト（Eisenberg & Levanon 2013 の安定 HK パネル + 古典 HK の curated 部分集合）
#   - hk_mouse.txt : マウス（古典 HK の ortholog / 慣用シンボル）
# config.yaml の screening.hk_gene_file を指定すれば任意リストに差替可。
# ファイル形式: 1 行 1 遺伝子シンボル。'#' 始まりの行と空行は無視する。
# ============================================================================

load_housekeeping <- function(cfg) {
  path <- cfg$screening$hk_gene_file %||% ""
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
