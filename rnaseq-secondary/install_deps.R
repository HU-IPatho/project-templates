# bulk RNA-seq 2 次解析環境（DESeq2 + edgeR + SummarizedExperiment・Bioconductor）。
# renv で版を固定。renv.lock があれば「検証済みの版」に restore、無ければ / 不完全なら
# fresh install + snapshot で renv.lock を完成させる。
#
# 使い方（最初に一度だけ・プロジェクトルートから）:
#   Rscript install_deps.R
#
# 版ピン 3 系統（spec）: R+Bioconductor=renv.lock（Bioc リリース明示）/ Python=pyproject(uv.lock)
#   / システムツール=omics-dev イメージのタグ。本スクリプトは 1 系統目（renv.lock）を担う。
# 同梱の renv.lock は Bioconductor 3.20 系に pin した seed。fresh install はこの Bioc
# リリースへ揃えたうえで依存閉包込みの完全な renv.lock を snapshot し直す。
options(Ncpus = 4L)
ppm <- "https://packagemanager.posit.co/cran/__linux__/noble/latest"   # Linux バイナリで高速
options(repos = c(CRAN = ppm))
Sys.setenv(RENV_CONFIG_PPM_ENABLED = "TRUE")
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

# 本テンプレの直接依存（Bioc + CRAN）。
bioc_pkgs <- c("SummarizedExperiment", "DESeq2", "edgeR", "limma", "apeglm")
cran_pkgs <- c("here", "yaml", "ggplot2", "pheatmap", "matrixStats", "data.table",
               "knitr", "rmarkdown")   # reports/report.qmd（Quarto knitr エンジン）用

# renv.lock に明示された Bioconductor リリースを読む（無ければ NULL）。
lock_bioc_version <- function() {
  if (!file.exists("renv.lock")) return(NULL)
  tryCatch(renv::lockfile_read("renv.lock")$Bioconductor$Version,
           error = function(e) NULL)
}

fresh_install <- function() {
  message("[install_deps] fresh install + snapshot（renv.lock を完成させる）")
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  bioc_ver <- lock_bioc_version()
  if (!is.null(bioc_ver)) {
    BiocManager::install(version = bioc_ver, update = FALSE, ask = FALSE)
  }
  if (!file.exists("renv/activate.R")) renv::init(bare = TRUE, restart = FALSE)
  options(repos = c(CRAN = ppm))
  renv::install(c(cran_pkgs, paste0("bioc::", bioc_pkgs)))
  renv::snapshot(type = "all", prompt = FALSE)
}

# 直接依存がすべて load できるか（restore が閉包を満たしたかの検証）。
deps_loadable <- function() {
  all(vapply(c(bioc_pkgs, cran_pkgs),
             function(p) requireNamespace(p, quietly = TRUE), logical(1)))
}

if (file.exists("renv.lock")) {
  message("[install_deps] renv.lock を検出 → restore（版を固定して再現）")
  ok <- tryCatch({ renv::restore(prompt = FALSE); TRUE },
                 error = function(e) { message("restore 失敗: ", conditionMessage(e)); FALSE })
  if (!ok || !deps_loadable()) {
    message("[install_deps] restore が不完全 → fresh install で補完")
    fresh_install()
  }
} else {
  fresh_install()
}

suppressPackageStartupMessages({
  library(DESeq2); library(edgeR); library(SummarizedExperiment)
})
cat("DESeq2               ", as.character(packageVersion("DESeq2")), "\n")
cat("edgeR                ", as.character(packageVersion("edgeR")), "\n")
cat("SummarizedExperiment ", as.character(packageVersion("SummarizedExperiment")), "\n")
cat("INSTALL_DONE\n")
