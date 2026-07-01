# データ読込ローダ（CONFIG$loader で切替）。
# scRNA の生データは形式が複数あるため、形式ごとのローダをここに集約する。
# 返り値: named list（sample 名 → Seurat オブジェクト）。01_load_qc.R が merge する。
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

# --- 10x Cell Ranger 出力（最も一般的）---
# data_dir 直下に sample ごとのサブディレクトリがあり、各々に
# matrix.mtx.gz / features.tsv.gz / barcodes.tsv.gz を含む想定。
load_10x <- function(data_dir, min_cells, min_features) {
  sample_dirs <- list.dirs(data_dir, recursive = FALSE)
  if (length(sample_dirs) == 0)
    stop("10x: ", data_dir, " に sample サブディレクトリが無い（各 sample に matrix/features/barcodes.tsv.gz を置く）")
  objs <- list()
  for (d in sample_dirs) {
    sample <- basename(d)
    counts <- Read10X(d)   # features×cells の sparse matrix
    o <- CreateSeuratObject(counts = counts, project = sample,
                            min.cells = min_cells, min.features = min_features)
    o$sample <- sample
    objs[[sample]] <- o
    cat(sprintf("  [10x] %s: %d cells x %d genes\n", sample, ncol(o), nrow(o)))
  }
  objs
}

# --- dense な per-sample 行列（GEO processed matrix 等）---
# 各ファイル: tab 区切り gz。1 行目=バーコード、2 行目以降=遺伝子名 + 各細胞の値。
# ★ dogfood で判明した罠（この 3 点を外すと壊れる。自前で読み込みを書くときも同じ）:
#   1) ragged: header 行の列数 = データ行より 1 少ない → fread(header=TRUE) は列ズレで
#      バーコードが重複する → header を別途読んで割当て、本体は header=FALSE, skip=1 で読む。
#   2) 重複遺伝子名 → make.unique。ただし Seurat は feature 名の '_' を '-' に置換するので、
#      先に gsub('_','-') してから make.unique（後だと置換で重複が再発する）。
#   3) fread の gz 直読は R.utils 依存 → cmd="zcat ..." で回避する。
load_dense <- function(data_dir, pattern, sample_regex, min_cells, min_features) {
  if (!requireNamespace("data.table", quietly = TRUE)) renv::install("data.table")
  files <- sort(list.files(data_dir, pattern = pattern, full.names = TRUE))
  if (length(files) == 0)
    stop("dense: ", data_dir, " に pattern=", pattern, " のファイルが無い")
  objs <- list()
  for (f in files) {
    sample <- sub(sample_regex, "\\1", basename(f))
    con <- gzfile(f); barcodes <- strsplit(readLines(con, n = 1), "\t")[[1]]; close(con)  # 罠1
    dt <- data.table::fread(cmd = paste("zcat", shQuote(f)), sep = "\t",
                            header = FALSE, skip = 1)                                       # 罠1,3
    genes <- dt[[1]]; dt[[1]] <- NULL
    m <- as.matrix(dt)
    rownames(m) <- make.unique(gsub("_", "-", as.character(genes)))                        # 罠2
    colnames(m) <- make.unique(as.character(barcodes))
    m <- as(m, "CsparseMatrix")
    o <- CreateSeuratObject(counts = m, project = sample,
                            min.cells = min_cells, min.features = min_features)
    o$sample <- sample
    objs[[sample]] <- o
    cat(sprintf("  [dense] %s: %d cells x %d genes\n", sample, ncol(o), nrow(o)))
    rm(dt, m); gc(verbose = FALSE)
  }
  objs
}

load_samples <- function(cfg) {
  if (cfg$loader == "tenx") {
    load_10x(cfg$data_dir, cfg$min_cells, cfg$min_features)
  } else if (cfg$loader == "dense") {
    load_dense(cfg$data_dir, cfg$dense_pattern, cfg$dense_sample_regex,
               cfg$min_cells, cfg$min_features)
  } else {
    stop("未知の loader: ", cfg$loader, "（\"tenx\" か \"dense\"）")
  }
}
