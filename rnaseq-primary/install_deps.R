# rnaseq-primary 解析環境（R + Bioconductor）。renv で版を固定・PPM バイナリで高速化。
#
# 使い方: 最初に一度だけ実行する（omics-dev イメージ内）。
#   Rscript install_deps.R
# renv.lock があれば「検証済みの版」に固定して復元する（＝再現性）。
# 無ければ Bioconductor リリースを固定して新規 install → renv.lock を生成する。
options(Ncpus = 4L)
ppm <- "https://packagemanager.posit.co/cran/__linux__/noble/latest"  # Linux バイナリで高速
options(repos = c(CRAN = ppm))
Sys.setenv(RENV_CONFIG_PPM_ENABLED = "TRUE")
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

BIOC_VERSION <- "3.23"   # ★Bioconductor リリースを明示（renv.lock / settings.json と一致）

# CRAN + Bioconductor 双方から入れる依存（bioc:: 接頭辞で Bioc を指定）
pkgs <- c(
  "here", "yaml", "jsonlite", "ggplot2", "matrixStats",   # CRAN
  "bioc::tximeta", "bioc::tximport", "bioc::SummarizedExperiment"  # Bioconductor
)

fresh_install <- function() {
  message("[install_deps] 新規 install + snapshot（renv.lock を生成）")
  if (!file.exists("renv/activate.R")) renv::init(bare = TRUE, restart = FALSE)
  options(repos = c(CRAN = ppm))
  renv::install(pkgs)
  renv::snapshot(type = "all", prompt = FALSE)
}

# renv に Bioconductor リリースを認識させる（版ピンの要）
options(renv.bioconductor.version = BIOC_VERSION)

if (file.exists("renv.lock")) {
  message("[install_deps] renv.lock を検出 → restore（版を固定して再現）")
  ok <- tryCatch({ renv::restore(prompt = FALSE); TRUE },
                 error = function(e) { message("restore 失敗: ", conditionMessage(e)); FALSE })
  if (!ok) {
    message("[install_deps] restore に失敗 → 新規 install にフォールバック")
    fresh_install()
  }
} else {
  fresh_install()
}

suppressPackageStartupMessages({ library(tximeta); library(SummarizedExperiment) })
cat("tximeta",             as.character(packageVersion("tximeta")), "\n")
cat("SummarizedExperiment", as.character(packageVersion("SummarizedExperiment")), "\n")
cat("Bioconductor",         BIOC_VERSION, "\n")
cat("INSTALL_DONE\n")
