# 出力ハーネス — 図表を outputs/ へ統一体裁で保存し captions.tsv に来歴を索引する。
# 生の ggsave/write.csv を直呼びせず、必ずこの save_fig / save_table を経由する。
# そうすることで (1) 図は PNG+PDF 両形式・固定 DPI・統一テーマ、(2) どのスクリプトが
# 生成したか（provenance）が captions.tsv に自動記録され、図表と由来の乖離を防ぐ。
suppressPackageStartupMessages({ library(here); library(ggplot2) })

SAVE_FIG_DPI <- 300                       # 固定 DPI（スクリプト間でばらつかせない）
theme_scrna  <- function() theme_bw(base_size = 11)  # 統一テーマ

.out_figures  <- function() here::here("outputs", "figures")
.out_tables   <- function() here::here("outputs", "tables")
.captions_tsv <- function() here::here("outputs", "captions.tsv")

# 図表番号 → 説明 → 由来スクリプト を captions.tsv に追記（無ければヘッダ付きで新規作成）。
.append_caption <- function(id, desc, script) {
  path <- .captions_tsv()
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  if (!file.exists(path)) cat("id\tdescription\tscript\n", file = path)
  cat(sprintf("%s\t%s\t%s\n", id, desc, script), file = path, append = TRUE)
  invisible(path)
}

# 図を PNG（閲覧用ラスタ）+ PDF（出版用ベクタ）両形式で固定 DPI・統一テーマ保存する。
# id="fig01" 形式・desc="umap_clusters" 等 → outputs/figures/fig01_umap_clusters.{png,pdf}
save_fig <- function(plot, id, desc, script, width = 7, height = 6) {
  dir.create(.out_figures(), showWarnings = FALSE, recursive = TRUE)
  base <- file.path(.out_figures(), sprintf("%s_%s", id, desc))  # figNN_<desc>
  p <- plot & theme_scrna()   # patchwork/ggplot 双方に効く（& は各パネルへ適用）
  ggplot2::ggsave(paste0(base, ".png"), p, width = width, height = height, dpi = SAVE_FIG_DPI)
  ggplot2::ggsave(paste0(base, ".pdf"), p, width = width, height = height)
  .append_caption(id, desc, script)
  invisible(base)
}

# 表を outputs/tables へ規定形式（CSV）で保存し captions.tsv に記録する。
# id="tbl01" 形式・desc="markers_all" 等 → outputs/tables/tbl01_markers_all.csv
save_table <- function(x, id, desc, script) {
  dir.create(.out_tables(), showWarnings = FALSE, recursive = TRUE)
  path <- file.path(.out_tables(), sprintf("%s_%s.csv", id, desc))
  utils::write.csv(x, path, row.names = FALSE)
  .append_caption(id, desc, script)
  invisible(path)
}
