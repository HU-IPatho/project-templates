# 00: 設定ローダ。config.yaml（事実の SSoT）を読み CONFIG list を組む薄いローダ。
# ★編集するのは config.yaml だけ★（このファイルとパイプライン本体は触らなくてよい）。
# 01/02/03/04 は冒頭で source(here::here("analysis/00_config.R")) して CONFIG を得る。
suppressPackageStartupMessages({ library(here); library(yaml) })

CONFIG <- yaml::read_yaml(here::here("config.yaml"))

# 型の担保（YAML は数値/文字列を素直に返すが、下流が前提とする型を明示する）。
CONFIG$qc <- lapply(CONFIG$qc, as.numeric)
for (k in c("min_cells", "min_features", "n_variable", "n_pcs", "dims", "resolution"))
  CONFIG[[k]] <- as.numeric(CONFIG[[k]])

# marker_hint は named character（遺伝子 → 細胞種ラベル）に正規化する。
CONFIG$marker_hint <- unlist(CONFIG$marker_hint)

# --- 標準 conform 由来キーの既定補完（古い config でも動く後方互換・欠落時は安全な既定へ）---
# 欠けたキーは「標準の既定 baseline（追加処理 off・記録のみ）」に落とす。既存プロジェクトの config を
# 壊さずに新しい 02/04 が走るようにするための defensive defaults。
`%||%` <- function(a, b) if (is.null(a)) b else a

CONFIG$sample_class <- CONFIG$sample_class %||% "mixed"
CONFIG$assay        <- CONFIG$assay        %||% "droplet_3prime"

# resolutions（走査解像度群）: 未指定なら primary resolution 単独に落とす。
CONFIG$resolutions <- as.numeric(unlist(CONFIG$resolutions %||% CONFIG$resolution))
CONFIG$resolutions <- sort(unique(c(CONFIG$resolution, CONFIG$resolutions)))

# doublet: method / rationale。未指定は none（除去なし）+ 既定 rationale。
CONFIG$doublet <- CONFIG$doublet %||% list()
CONFIG$doublet$method    <- CONFIG$doublet$method    %||% "none"
CONFIG$doublet$rationale <- CONFIG$doublet$rationale %||% "doublet 取扱い未記録（既定 baseline: 除去なし）。"

# cell_cycle: score / regress（論理値）。
CONFIG$cell_cycle <- CONFIG$cell_cycle %||% list()
CONFIG$cell_cycle$score  <- isTRUE(CONFIG$cell_cycle$score  %||% FALSE)
CONFIG$cell_cycle$regress <- isTRUE(CONFIG$cell_cycle$regress %||% FALSE)

# reference（two-track クロスチェック）。既定 off。
CONFIG$reference <- CONFIG$reference %||% list()
CONFIG$reference$enabled <- isTRUE(CONFIG$reference$enabled %||% FALSE)
CONFIG$reference$method  <- CONFIG$reference$method  %||% "none"

# cnv_gate（悪性 CNV 検証ゲート・carcinoma 限定）。既定 off。
CONFIG$cnv_gate <- CONFIG$cnv_gate %||% list()
CONFIG$cnv_gate$enabled <- isTRUE(CONFIG$cnv_gate$enabled %||% FALSE)
CONFIG$cnv_gate$tool    <- CONFIG$cnv_gate$tool    %||% "none"
CONFIG$cnv_gate$normal_reference <- CONFIG$cnv_gate$normal_reference %||% ""

# 純細胞株では malignancy が a priori 既知ゆえ CNV ゲートは意味を成さない（標準 annex）。
# enabled のまま pure_cell_line が来たら fail-safe で off にし、警告する（silent に走らせない）。
if (identical(CONFIG$sample_class, "pure_cell_line") && isTRUE(CONFIG$cnv_gate$enabled)) {
  warning("sample_class=pure_cell_line では CNV ゲートを省略します（標準の純細胞株 annex）。cnv_gate.enabled を無効化しました。")
  CONFIG$cnv_gate$enabled <- FALSE
}
