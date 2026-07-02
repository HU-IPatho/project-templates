# ============================================================================
# 出力ハーネス — 図表を一元化して保存する（spec「出力ハーネス」要件）。
#   save_fig()  : 1 呼出しで PNG+PDF を固定 DPI・統一テーマで outputs/figures へ
#   save_table(): 表を outputs/tables へ規定形式（TSV）で
#   どちらも outputs/captions.tsv に「図表番号 → 説明 → 由来スクリプト」を索引する
# 研究員は生の ggsave/write.csv を直呼びせず本ハーネス経由で出力する（体裁・来歴の一貫性）。
# ============================================================================
suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
})

# NULL 合体（config の任意キーに既定値を与える）。DTE/DTU の各 R 関数が参照する。
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

# 全図で統一する体裁（テーマ）と DPI。ここを変えると全図に一貫して反映される。
.FIG_DPI   <- 300
.FIG_THEME <- ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    plot.title       = ggplot2::element_text(face = "bold")
  )

.captions_path <- function() here::here("outputs", "captions.tsv")

# captions.tsv に 1 行追記（無ければヘッダ付きで作る）。列: id / description / script
.append_caption <- function(id, description, script) {
  path <- .captions_path()
  if (!file.exists(path)) {
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    writeLines("id\tdescription\tscript", path)
  }
  cat(sprintf("%s\t%s\t%s\n", id, description, script), file = path, append = TRUE)
}

#' 図を PNG+PDF 両形式で統一体裁・固定 DPI で保存し captions に索引する。
#' @param plot   ggplot オブジェクト
#' @param num    図番号（figNN の NN）
#' @param desc   短い説明（ファイル名 figNN_<desc> と captions に使う）
#' @param script 由来スクリプト名（例 "analysis/06_figures.R"）
#' @param width,height インチ
save_fig <- function(plot, num, desc, script, width = 6, height = 4) {
  stopifnot(inherits(plot, "ggplot"))
  fig_id <- sprintf("fig%02d_%s", num, desc)
  dir    <- here::here("outputs", "figures")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  p <- plot + .FIG_THEME
  # PNG（閲覧用ラスタ）と PDF（出版用ベクタ）を 1 呼出しで両方出す
  ggplot2::ggsave(file.path(dir, paste0(fig_id, ".png")), p,
                  width = width, height = height, dpi = .FIG_DPI)
  ggplot2::ggsave(file.path(dir, paste0(fig_id, ".pdf")), p,
                  width = width, height = height)
  .append_caption(fig_id, desc, script)
  invisible(fig_id)
}

#' 表を outputs/tables に TSV で保存し captions に索引する。
#' @param df     data.frame
#' @param num    表番号（tabNN の NN）
#' @param desc   短い説明
#' @param script 由来スクリプト名
save_table <- function(df, num, desc, script) {
  stopifnot(is.data.frame(df))
  tab_id <- sprintf("tab%02d_%s", num, desc)
  dir    <- here::here("outputs", "tables")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  utils::write.table(df, file.path(dir, paste0(tab_id, ".tsv")),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  .append_caption(tab_id, desc, script)
  invisible(tab_id)
}
