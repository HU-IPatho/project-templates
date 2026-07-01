# ============================================================================
# 正準 object の構築 — counts 行列 + colData(config) → SummarizedExperiment
# ----------------------------------------------------------------------------
# ②固有パイプラインの正準 object は SummarizedExperiment（counts + colData を 1 object）。
# analysis/01_build_se.R がこれを呼び data/processed/se.rds に保存する（＝共有昇格の単位）。
# ============================================================================
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(S4Vectors)
})

# counts_file: gene×sample の raw counts 行列（TSV・1 列目=gene_id・ヘッダ=sample_id）。
# coldata: data.frame（行=sample・config の samples から構築。rownames=sample_id）。
# 返り値: assay "counts"（integer 行列）と colData を持つ SummarizedExperiment。
build_se <- function(counts_file, coldata, gene_id_col = 1) {
  if (!file.exists(counts_file)) {
    stop("build_se: counts_file が無い: ", counts_file,
         "（fetch_data.sh / python/fetch_counts.py で data/ に用意する）")
  }
  m <- utils::read.delim(counts_file, header = TRUE, row.names = gene_id_col,
                         check.names = FALSE)
  m <- round(as.matrix(m))            # counts は非負整数（tximport 等の推定値は丸める）
  storage.mode(m) <- "integer"

  # colData を counts の列順（sample_id）に厳密整合させる（順序ズレは解析を壊す）
  missing <- setdiff(colnames(m), rownames(coldata))
  if (length(missing) > 0) {
    stop("build_se: counts の列に対応する sample_id が config に無い: ",
         paste(missing, collapse = ", "))
  }
  coldata <- coldata[colnames(m), , drop = FALSE]

  SummarizedExperiment::SummarizedExperiment(
    assays  = list(counts = m),
    colData = S4Vectors::DataFrame(coldata)
  )
}
