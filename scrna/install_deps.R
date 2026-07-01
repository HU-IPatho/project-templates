# scRNA-seq 解析環境（Seurat v5 + Harmony）。renv で版を固定・PPM バイナリで高速化。
#
# 使い方: 最初に一度だけ実行する。
#   Rscript install_deps.R
# renv.lock があれば「dogfood で検証済みの版」に固定して復元する（＝再現性）。
# 無ければ最新版を入れて renv.lock を新規生成する。
options(Ncpus = 4L)
ppm <- "https://packagemanager.posit.co/cran/__linux__/noble/latest"  # Linux バイナリで高速
options(repos = c(CRAN = ppm))
Sys.setenv(RENV_CONFIG_PPM_ENABLED = "TRUE")
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

# here=ルート相対パス解決 / yaml=config.yaml 読込（① 標準で必須）。
pkgs <- c("Seurat", "harmony", "dplyr", "patchwork", "ggplot2", "Matrix", "data.table",
          "here", "yaml")

fresh_install <- function() {
  message("[install_deps] 新規 install + snapshot（renv.lock を生成）")
  if (!file.exists("renv/activate.R")) renv::init(bare = TRUE, restart = FALSE)
  options(repos = c(CRAN = ppm))
  renv::install(pkgs)
  renv::snapshot(type = "all", prompt = FALSE)
}

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

suppressPackageStartupMessages({ library(Seurat); library(harmony) })
cat("Seurat",  as.character(packageVersion("Seurat")),  "\n")
cat("harmony", as.character(packageVersion("harmony")), "\n")
cat("INSTALL_DONE\n")
