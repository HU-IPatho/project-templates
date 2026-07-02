# 06: 図。主役=isoform 使用比プロット（DTU 上位 n_top 遺伝子・群別）、併走=DTE volcano（対比別）。
#   重い計算はしない（04/05 が保存した se_tx.rds / dtu.rds を読むだけ）。全て出力ハーネス経由。
# 実行（ルートから）:  Rscript analysis/06_figures.R
suppressPackageStartupMessages({
  library(here); library(yaml); library(ggplot2); library(SummarizedExperiment)
})
source(here::here("R", "helpers.R"))   # save_fig / %||%
source(here::here("R", "plots.R"))     # plot_isoform_proportions / plot_dte_volcano

cfg   <- yaml::read_yaml(here::here("config.yaml"))
se    <- readRDS(here::here("data", "processed", "se_tx.rds"))
res   <- readRDS(here::here("data", "processed", "dtu.rds"))
gcol  <- cfg$design$group_col
fdr   <- cfg$dtu$fdr %||% 0.05
n_top   <- cfg$dtu$n_top %||% 12
methods <- cfg$dtu$methods %||% c("drimseq")

# --- 主役図: isoform 使用比（primary エンジン=methods[1] の有意 DTU 遺伝子・上位 n_top）---
# プロット本体は se_tx の観測比率（engine 非依存）だが、載せる遺伝子の選抜は primary エンジンの
# gene-level 有意性に従う（primary が drimseq でなく dexseq/swish でも動くよう名前で参照）。
# gene スコア = drimseq/dexseq は $table$gene_padj、swish は tx-level qvalue(padx) の gene 内最小。
dtu_gene_scores <- function(engine_res) {
  do.call(rbind, lapply(engine_res, function(x) {
    if (is.data.frame(x)) {                        # swish: tx-level df（gene_id, padx=qvalue）
      data.frame(gene_id = x$gene_id, score = x$padx, stringsAsFactors = FALSE)
    } else {                                        # drimseq/dexseq: $table（gene_id, gene_padj）
      data.frame(gene_id = x$table$gene_id, score = x$table$gene_padj,
                 stringsAsFactors = FALSE)
    }
  }))
}

prim_dtu <- methods[1]
if (!is.null(res$dtu[[prim_dtu]])) {
  agg <- dtu_gene_scores(res$dtu[[prim_dtu]])
  agg <- agg[!is.na(agg$score), ]
  # 遺伝子ごとに（全対比中の）最小スコアで順位付け
  gene_rank <- stats::aggregate(score ~ gene_id, data = agg, FUN = min)
  gene_rank <- gene_rank[order(gene_rank$score), ]
  sig       <- gene_rank$gene_id[gene_rank$score < fdr]
  # 有意遺伝子があればそれを、無ければ順位上位を（図が空にならないよう）上位 n_top まで
  top_genes <- utils::head(if (length(sig) > 0) sig else gene_rank$gene_id, n_top)
  p_iso <- plot_isoform_proportions(
    se, top_genes, gcol,
    title = sprintf("Isoform usage: top %d DTU genes (%s)", length(top_genes), prim_dtu))
  save_fig(p_iso, num = 1, desc = sprintf("isoform_proportions_top%d", n_top),
           script = "analysis/06_figures.R", width = 9, height = 7)
}

# --- 併走: DTE volcano（選択した各 DTE エンジン・対比別・fig1NN バンドで衝突回避）---
# dte.method: both なら DESeq2/edgeR 両方の volcano を engine 修飾 desc で保存。
for (eng in names(res$dte_tidy)) {
  tid <- res$dte_tidy[[eng]]
  for (k in seq_along(tid)) {
    nm    <- names(tid)[k]
    p_vol <- plot_dte_volcano(tid[[nm]], fdr = cfg$dte$fdr %||% 0.05,
               title = sprintf("DTE volcano: %s (%s)", nm, eng))
    save_fig(p_vol, num = 10 + k, desc = sprintf("dte_volcano_%s_%s", tolower(eng), nm),
             script = "analysis/06_figures.R", width = 6, height = 5)
  }
}

cat("FIGURES_DONE\n")
