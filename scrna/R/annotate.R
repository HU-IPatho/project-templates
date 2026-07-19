# アノテーション・ロジック（04_annotate.R が source する再利用関数）。
# 研究室標準 scrna-annotation-standard（v0・方法論）に conform する scaffold。
# ここは「方法論の構造」（直交レイヤ・two-track・CNV ゲート・confidence/review routing）を用意し、
# 生物学固有の中身（marker panel・細胞種呼称・numeric 閾値・ツール pin）は config.yaml + /grill-me に委ねる。
#
# ★重要（標準）: 手動クラスタマーカー解釈が最終権威。本 scaffold は自動で identity を確定しない——
#   クラスタごとの「証拠」を組み立て、confidence と review-priority を付し、研究員の手動確定を助ける。
suppressPackageStartupMessages({ library(Seurat); library(dplyr) })

# --- 系統の当たり付け（marker_hint の従属補強・identity を単独確定しない）---
# クラスタの有意 positive marker のうち marker_hint に載る遺伝子を数え、hint ラベル別に集計する。
# 返り: data.frame(cluster, suggested_lineage, n_hint_hits, lineage_margin)。
#   suggested_lineage : hit 数最多の hint ラベル（0 hit は "unknown"）。あくまで示唆であり決定ではない。
#   lineage_margin    : 最多 hit と 2 位の差（0=曖昧・tie → review トリガー）。
suggest_lineage <- function(markers, marker_hint) {
  hint_genes <- names(marker_hint)
  clusters <- sort(unique(markers$cluster))
  out <- lapply(clusters, function(cl) {
    genes_cl <- markers$gene[markers$cluster == cl]
    hits <- marker_hint[intersect(hint_genes, genes_cl)]
    if (length(hits) == 0)
      return(data.frame(cluster = cl, suggested_lineage = "unknown",
                        n_hint_hits = 0L, lineage_margin = 0L, stringsAsFactors = FALSE))
    tab <- sort(table(as.character(hits)), decreasing = TRUE)
    margin <- if (length(tab) >= 2) as.integer(tab[1] - tab[2]) else as.integer(tab[1])
    data.frame(cluster = cl, suggested_lineage = names(tab)[1],
               n_hint_hits = as.integer(sum(tab)), lineage_margin = margin,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

# --- marker の発現特異性（effect specificity・p 値に依存しない）---
# クラスタごとに max(pct.1 - pct.2) を返す（1 に近いほど特異的マーカーを持つ）。標準: p 値を confidence に使わない。
marker_specificity <- function(markers) {
  markers %>%
    group_by(cluster) %>%
    summarise(top_marker_specificity = suppressWarnings(max(pct.1 - pct.2, na.rm = TRUE)), .groups = "drop") %>%
    mutate(top_marker_specificity = ifelse(is.finite(top_marker_specificity), top_marker_specificity, NA_real_))
}

# --- 組成の偏り（composition confound の review トリガー・固定閾値でなくヒューリスティック）---
# クラスタごとに「単一 sample 由来の細胞割合」を返す。高いほど単一サンプル駆動の疑い（design-specific・
# 標準: 固定閾値化は grill ゲート委譲。ここでは割合を記録し、既定トリガー閾値は緩く置く）。
composition_flags <- function(o, batch_key) {
  md <- o@meta.data
  cl <- as.character(md$seurat_clusters); smp <- as.character(md[[batch_key]])
  clusters <- sort(unique(cl))
  data.frame(
    cluster = clusters,
    dominant_sample_frac = vapply(clusters, function(c) {
      t <- table(smp[cl == c]); if (length(t) == 0) NA_real_ else max(t) / sum(t)
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
}

# --- cell-cycle の優勢 Phase（直交軸の記録・identity ではない）---
cellcycle_flags <- function(o) {
  if (!"Phase" %in% colnames(o@meta.data)) return(NULL)
  md <- o@meta.data; cl <- as.character(md$seurat_clusters); ph <- as.character(md$Phase)
  clusters <- sort(unique(cl))
  data.frame(cluster = clusters,
             dominant_phase = vapply(clusters, function(c) {
               t <- table(ph[cl == c]); if (length(t) == 0) NA_character_ else names(which.max(t))
             }, character(1)),
             stringsAsFactors = FALSE)
}

# --- 悪性判定の直交 CNV 検証ゲート（carcinoma 限定・tool-agnostic・非対称）---
# 標準: 上皮/EMT を marker 単独で「悪性」と呼ばず inferred CNV を検証ゲートに課す。ゲートは非対称
# （CNV 陰性→非悪性とは限らず「保留」）。ツールは pin しない。ここは構造（malignant_status 列）を用意する
# フック。既定 disabled では全クラスタ "not_assessed"。enabled 時は tool 実装を研究員が /grill-me で埋める。
cnv_gate_status <- function(o, cfg) {
  clusters <- sort(unique(as.character(o@meta.data$seurat_clusters)))
  na_out <- data.frame(cluster = clusters, malignant_status = "not_assessed",
                       cnv_note = "CNV ゲート未実施（cnv_gate.enabled=false or 非適用）", stringsAsFactors = FALSE)
  if (!isTRUE(cfg$cnv_gate$enabled)) return(na_out)
  tool <- cfg$cnv_gate$tool
  if (identical(tool, "none") || is.null(tool) || tool == "") {
    warning("cnv_gate.enabled=true ですが tool 未指定です。malignant_status=not_assessed のままにします。")
    return(na_out)
  }
  # tool-agnostic フック: 実際の CNV 推定（inferCNV / CopyKAT / SCEVAN 等）は非悪性 reference の指定と
  # gene-position 準備が要り dataset-specific ゆえ、研究員が /grill-me で本関数を完成させる。scaffold は
  # 未完成のまま「悪性」を捏造せず、明示的に停止する（silent に偽陽性を出さない・標準の非対称ゲート）。
  stop(sprintf(paste0(
    "cnv_gate.tool=%s の CNV 推定は未実装フックです。研究員が R/annotate.R の cnv_gate_status() を ",
    "完成させてください（非悪性 reference=%s を diploid baseline に、直交 CNV 証拠でクラスタ単位に判定。",
    "ゲートは非対称: CNV 陰性→非悪性でなく『保留(hold)』。near-diploid 腫瘍種は例外条項で非適用）。",
    "手順は /grill-me（grill ゲート）で確定する。"),
    tool, ifelse(nzchar(cfg$cnv_gate$normal_reference), cfg$cnv_gate$normal_reference, "未指定")))
}

# --- 参照ベース クロスチェック（two-track・適合アトラス存在時は必須 complementary）---
# 標準: 参照は無条件 primary/必須にはしない two-track。不一致クラスタを review へ routing する。
# ここは構造（reference_label / reference_agreement 列）を用意するフック。既定 disabled では NA。
reference_crosscheck <- function(o, cfg) {
  clusters <- sort(unique(as.character(o@meta.data$seurat_clusters)))
  na_out <- data.frame(cluster = clusters, reference_label = NA_character_,
                       reference_agreement = NA, stringsAsFactors = FALSE)
  if (!isTRUE(cfg$reference$enabled)) return(na_out)
  method <- cfg$reference$method
  if (identical(method, "none") || is.null(method) || method == "") {
    warning("reference.enabled=true ですが method 未指定です。reference_label=NA のままにします。")
    return(na_out)
  }
  # SingleR / Azimuth 等の適合アトラス指定は dataset-specific ゆえ研究員が /grill-me で本関数を完成させる。
  # scaffold は参照ラベルを捏造せず明示停止する（手動注釈との一致を後で計算する two-track の受け皿だけ用意）。
  stop(sprintf(paste0(
    "reference.method=%s の参照クロスチェックは未実装フックです。研究員が R/annotate.R の reference_crosscheck() を ",
    "完成させてください（適合アトラスを指定し per-cluster に参照ラベルを付与、手動注釈との一致を reference_agreement に。",
    "参照は必須 complementary の two-track であり無条件 primary にしない）。適合判定・アトラスは /grill-me で確定する。"),
    method))
}

# --- confidence + review-priority routing（標準の注釈監査層）---
# 各クラスタに操作的 confidence と review_flag を付す。標準の必須:「marker 遺伝子 p 値を confidence に使わない」。
# ここでは marker 特異性・系統 margin・組成偏り・（有効時）参照一致/CNV を統合したヒューリスティック。
# 閾値は「review へ回すトリガー」であり確定閾値ではない（skew_trigger は config review_triggers から・tune 可能）。
# sample_class / n_batches で composition トリガーの適用可否を切り替える（下記）。
cluster_confidence <- function(ev, sample_class = "mixed", n_batches = 2L, skew_trigger = 0.90) {
  # NA 安全化: marker を 1 つも持たないクラスタ（FindAllMarkers が閾値を通す positive marker を返さない）は
  # outer merge で n_hint_hits/lineage_margin が NA になる。NA を安全側（示唆なし＝review 対象）に倒す。
  hits   <- ifelse(is.na(ev$n_hint_hits),    0, ev$n_hint_hits)
  margin <- ifelse(is.na(ev$lineage_margin), 0, ev$lineage_margin)
  f_no_hint   <- hits == 0                                                   # 系統示唆なし（unknown）
  f_ambiguous <- margin == 0                                                 # 系統が曖昧（tie）
  # composition skew は「multi-replicate 比較設計」依存の design-specific（標準 R7/R8）。単一 sample/batch 設計と
  # 純細胞株（annex: composition フラグ無効化）では自動で無効化する（過検知回避・標準の非適用条件に一致）。
  skew_applicable <- (n_batches >= 2) && !identical(sample_class, "pure_cell_line")
  f_skew <- skew_applicable & !is.na(ev$dominant_sample_frac) & ev$dominant_sample_frac >= skew_trigger
  f_ref_mis <- !is.na(ev$reference_agreement) & ev$reference_agreement == FALSE  # 参照不一致
  ev$review_flag <- f_no_hint | f_ambiguous | f_skew | f_ref_mis
  ev$review_reason <- vapply(seq_len(nrow(ev)), function(i) {
    r <- c(if (f_no_hint[i])   "no_lineage_hint"    else NULL,
           if (f_ambiguous[i]) "ambiguous_lineage"  else NULL,
           if (f_skew[i])      "single_sample_skew" else NULL,
           if (f_ref_mis[i])   "reference_mismatch" else NULL)
    if (length(r) == 0) "ok" else paste(r, collapse = ";")
  }, character(1))
  ev
}
