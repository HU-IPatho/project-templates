# 00: 設定ローダ。config.yaml（事実の SSoT）を読み CONFIG list を組む薄いローダ。
# ★編集するのは config.yaml だけ★（このファイルとパイプライン本体は触らなくてよい）。
# 01/02/03 は冒頭で source(here::here("analysis/00_config.R")) して CONFIG を得る。
suppressPackageStartupMessages({ library(here); library(yaml) })

CONFIG <- yaml::read_yaml(here::here("config.yaml"))

# 型の担保（YAML は数値/文字列を素直に返すが、下流が前提とする型を明示する）。
CONFIG$qc <- lapply(CONFIG$qc, as.numeric)
for (k in c("min_cells", "min_features", "n_variable", "n_pcs", "dims", "resolution"))
  CONFIG[[k]] <- as.numeric(CONFIG[[k]])

# marker_hint は named character（遺伝子 → 細胞種ラベル）に正規化する。
CONFIG$marker_hint <- unlist(CONFIG$marker_hint)
