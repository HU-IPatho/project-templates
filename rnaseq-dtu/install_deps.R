# rnaseq-dtu 解析環境（R + Bioconductor）。renv で版を固定・PPM バイナリで高速化。
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

# CRAN 依存（出力ハーネス・可視化・レポート）
cran_pkgs <- c("here", "yaml", "ggplot2", "pheatmap", "data.table", "knitr", "rmarkdown")
# Bioconductor 依存（DTE/DTU エンジン一式）:
#   import 基盤 = SummarizedExperiment / tximport / tximeta
#   DTU        = DRIMSeq / DEXSeq / stageR（二段 OFDR）
#   不確実性    = fishpond（swish・bootstrap 必須）
#   DTE        = DESeq2 / edgeR / limma / apeglm
bioc_pkgs <- c("SummarizedExperiment", "tximport", "tximeta",
               "DRIMSeq", "DEXSeq", "stageR", "fishpond",
               "DESeq2", "edgeR", "limma", "apeglm")
# 注: IsoformSwitchAnalyzeR は重い依存を引くため同梱しない。ISA を使う場合のみ各自で
#     renv::install("IsoformSwitchAnalyzeR") する（AGENTS.md「任意手法」参照）。
pkgs <- c(cran_pkgs, paste0("bioc::", bioc_pkgs))

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

suppressPackageStartupMessages({ library(tximport); library(DRIMSeq); library(stageR) })
cat("tximport", as.character(packageVersion("tximport")), "\n")
cat("DRIMSeq",  as.character(packageVersion("DRIMSeq")),  "\n")
cat("stageR",   as.character(packageVersion("stageR")),   "\n")
cat("Bioconductor", BIOC_VERSION, "\n")
cat("INSTALL_DONE\n")
