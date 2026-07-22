# ============================================================================
# 出力ハーネス — save_fig / save_table / captions.tsv（analysis-template-standard 準拠）
# ----------------------------------------------------------------------------
# spec（specs/analysis-template-standard/spec.md「出力ハーネス」）:
#   - save_fig は 1 呼出しで PNG+PDF 両形式を、固定 DPI・統一テーマで出力する。
#   - save_table は表を outputs/tables へ規定形式（TSV）で保存する。
#   - captions.tsv は 図表番号 → 説明 → 由来スクリプト を索引し provenance の乖離を防ぐ。
# 研究員は生の描画 API を直呼びせず、必ずこのハーネスを介して図表を出力する。
# ============================================================================
suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
})

# 統一体裁（スクリプト間で図の見た目を揃える単一の定義点）
.OUT_DPI   <- 300                       # 固定 DPI（出版品質のラスタ）
.FIG_THEME <- ggplot2::theme_bw(base_size = 11)

# NULL 合体（config の任意キーに既定値を与える）
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

# desc をファイル名の一部に使う前にパス危険文字（/ \）を無害化する（キャプション本文には原文を使う）
.fs_safe <- function(s) gsub("[/\\\\]", "-", s)

# --- 図の保存: figNN_<desc>.png と .pdf を outputs/figures に両形式で出す ---------
# plot は ggplot オブジェクト、または「現在のデバイスへ描画する関数」（base/grid 図・
# pheatmap 等）のいずれか。関数を渡せば png/pdf デバイスを開いて captured 描画する。
save_fig <- function(plot, fig_id, desc, script,
                     width = 7, height = 5, dpi = .OUT_DPI) {
  stopifnot(grepl("^fig[0-9]+$", fig_id))          # figNN 命名を強制
  fig_dir <- here::here("outputs", "figures")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  base <- file.path(fig_dir, sprintf("%s_%s", fig_id, .fs_safe(desc)))
  for (ext in c("png", "pdf")) {                   # ★ 両形式を必ず出す
    path <- paste0(base, ".", ext)
    if (inherits(plot, "ggplot")) {
      ggplot2::ggsave(path, plot = plot + .FIG_THEME,
                      width = width, height = height, dpi = dpi)
    } else if (is.function(plot)) {
      if (ext == "png") {
        grDevices::png(path, width = width, height = height, units = "in", res = dpi)
      } else {
        grDevices::pdf(path, width = width, height = height)
      }
      # plot() が失敗してもデバイスを必ず 1 度だけ閉じる（null device 二重 close を避ける）
      tryCatch(plot(), finally = grDevices::dev.off())
    } else {
      stop("save_fig: plot は ggplot オブジェクトか描画関数を渡す")
    }
  }
  append_caption(fig_id, desc, script)
  invisible(base)
}

# --- 表の保存: outputs/tables/<tbl_id>_<desc>.tsv --------------------------------
save_table <- function(x, tbl_id, desc, script) {
  tbl_dir <- here::here("outputs", "tables")
  dir.create(tbl_dir, showWarnings = FALSE, recursive = TRUE)
  path <- file.path(tbl_dir, sprintf("%s_%s.tsv", tbl_id, .fs_safe(desc)))
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
  append_caption(tbl_id, desc, script)
  invisible(path)
}

# --- captions.tsv への 1 行 upsert（id/description/source_script）-----------------
# 同じ id は上書きし、図表番号とその由来スクリプトの索引を単一ファイルに集約する。
append_caption <- function(id, desc, script) {
  path <- here::here("captions.tsv")
  cols <- c("id", "description", "source_script")
  tab <- if (file.exists(path)) {
    utils::read.table(path, sep = "\t", header = TRUE,
                      stringsAsFactors = FALSE, quote = "", comment.char = "")
  } else {
    stats::setNames(
      data.frame(character(), character(), character(), stringsAsFactors = FALSE),
      cols)
  }
  tab <- tab[tab$id != id, , drop = FALSE]
  tab <- rbind(tab,
               data.frame(id = id, description = desc, source_script = script,
                          stringsAsFactors = FALSE))
  tab <- tab[order(tab$id), , drop = FALSE]
  utils::write.table(tab, path, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(path)
}
