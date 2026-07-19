# 01: データ読込 → merge → doublet 取扱い → QC。設定=config.yaml、ローダ=R/load_data.R。
# 実行: プロジェクトルートから  Rscript analysis/01_load_qc.R  （数分以内なら前景で可）。
# パスは here::here() でルート相対解決するので、どこから叩いても同じに解決する。
#
# ★ambient 補正の前提（標準・filtered data 既知リスク）: SoupX / CellBender は raw/unfiltered droplet 行列を
#   必須とする。filtered（cell-called）行列のみで生 UMI が無い場合、ambient 補正は原理的に実行不能。
#   その場合は残留 ambient（高発現 housekeeping / 血液系遺伝子の低レベル発現）を marker 解釈から割り引く
#   運用にする（03/04 の marker 解釈で留意）。raw 行列の再取得は本雛形では要求しない。
suppressPackageStartupMessages({ library(here); library(Seurat); library(dplyr) })
source(here::here("analysis", "00_config.R"))
source(here::here("R", "load_data.R"))

cfg <- CONFIG
cfg$data_dir <- here::here(cfg$data_dir)   # ルート相対 → 絶対（実行位置非依存）
objs <- load_samples(cfg)
if (length(objs) == 1) {
  merged <- objs[[1]]
} else {
  merged <- merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs))
}
merged <- JoinLayers(merged)   # Seurat v5: sample ごとの layer を 1 つに結合
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = CONFIG$mito_pattern)
cat(sprintf("pre-QC: %d cells\n", ncol(merged)))

# --- doublet の取扱い（標準: 「有無・手法を記録」する。除去可否/手法/閾値の科学判断は grill ゲート）---
# method=none は除去せず rationale を記録するだけ。scDblFinder は per-sample にスコアし予測 doublet を除去する。
# どちらでも取扱いを merged@misc$doublet に機械可読で残す（下流 04 の confidence/監査が読む）。
dbl <- CONFIG$doublet
if (identical(dbl$method, "none")) {
  cat(sprintf("DOUBLET: method=none（除去なし）| rationale=%s\n", dbl$rationale))
} else if (identical(dbl$method, "scDblFinder")) {
  if (!requireNamespace("scDblFinder", quietly = TRUE))
    stop("doublet.method=scDblFinder には scDblFinder パッケージが必要です（renv::install('scDblFinder') 等で追加）。")
  n_before <- ncol(merged)
  sce <- Seurat::as.SingleCellExperiment(merged)
  # per-sample に走らせる（doublet は lane/sample 内で生じる）。batch_key 列を samples に渡す。
  sce <- scDblFinder::scDblFinder(sce, samples = as.character(merged@meta.data[[CONFIG$batch_key]]))
  merged$scDblFinder.class <- as.character(sce$scDblFinder.class)
  merged$scDblFinder.score <- as.numeric(sce$scDblFinder.score)
  merged <- subset(merged, subset = scDblFinder.class == "singlet")
  cat(sprintf("DOUBLET: method=scDblFinder | 除去 %d/%d cells（predicted doublet）\n",
              n_before - ncol(merged), n_before))
} else {
  stop("未知の doublet.method: ", dbl$method, "（\"none\" か \"scDblFinder\"）")
}
merged@misc$doublet <- list(method = dbl$method, rationale = dbl$rationale)

q <- CONFIG$qc
merged <- subset(merged, subset = nFeature_RNA > q$nFeature_min &
                                  nFeature_RNA < q$nFeature_max &
                                  percent.mt   < q$percent_mt_max)
cat(sprintf("post-QC: %d cells, %d %s\n",
            ncol(merged), length(unique(merged@meta.data[[CONFIG$batch_key]])), CONFIG$batch_key))

# QC 後の merged は data/raw から再生成可能な中間物 → data/interim（正準 object ではない）。
dir.create(here::here("data", "interim"), showWarnings = FALSE, recursive = TRUE)
saveRDS(merged, here::here("data", "interim", "seurat_qc.rds"))
cat("QC_DONE\n")
