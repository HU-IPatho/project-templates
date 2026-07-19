# 04: アノテーション（直交レイヤ・cluster-then-marker two-track・confidence/review routing）。
# 研究室標準 scrna-annotation-standard（v0・方法論）に conform する。ロジック本体は R/annotate.R。
# 実行: プロジェクトルートから  Rscript analysis/04_annotate.R
# 入力: data/processed/seurat_integrated.rds（02）+ outputs/tables/tbl01_markers_all.csv（03）
# 出力: outputs/tables/tbl03_cluster_evidence.csv・tbl04_annotation.csv + data/processed/seurat_annotated.rds
#
# ★方法論（標準・要点）:
#  - 手動クラスタマーカー解釈が最終権威。本段は「証拠」を組み立て confidence/review を付し手動確定を助ける
#    （自動で identity を確定しない）。最終ラベルは annotations.yaml（研究員が埋める）で与える。
#  - 直交レイヤ: 系統(lineage) / 悪性(malignant_status・CNV ゲート) / 状態(state) / cell-cycle（直交・identity でない）。
#  - two-track: 適合アトラスが在れば参照ベースを必須 complementary クロスチェックに（config reference.enabled）。
#  - 悪性は marker 単独で呼ばず直交 CNV 検証ゲート（carcinoma 限定・非対称: CNV 陰性→保留）で確定する。
suppressPackageStartupMessages({ library(here); library(Seurat); library(dplyr); library(ggplot2) })
source(here::here("analysis", "00_config.R"))
source(here::here("R", "helpers.R"))
source(here::here("R", "annotate.R"))

o  <- readRDS(here::here("data", "processed", "seurat_integrated.rds"))
o  <- JoinLayers(o)
Idents(o) <- "seurat_clusters"
mk <- read.csv(here::here("outputs", "tables", "tbl01_markers_all.csv"), stringsAsFactors = FALSE)

# --- クラスタ単位の証拠を組み立てる（直交レイヤ + 監査層の材料）---
tabs <- list(
  suggest_lineage(mk, CONFIG$marker_hint),   # 系統の当たり付け（従属補強）
  marker_specificity(mk),                    # marker 特異性（p 値に依存しない）
  composition_flags(o, CONFIG$batch_key),    # 組成の偏り（design-specific トリガー）
  cellcycle_flags(o),                        # cell-cycle 優勢 Phase（直交・NULL 可）
  cnv_gate_status(o, CONFIG),                # 悪性 CNV 検証ゲート（既定 not_assessed）
  reference_crosscheck(o, CONFIG)            # 参照 two-track（既定 NA）
)
tabs <- Filter(Negate(is.null), tabs)
tabs <- lapply(tabs, function(d) { d$cluster <- as.character(d$cluster); d })
ev <- Reduce(function(a, b) merge(a, b, by = "cluster", all = TRUE), tabs)
ev <- cluster_confidence(ev)                 # confidence + review_flag + review_reason
ev <- ev[order(suppressWarnings(as.numeric(ev$cluster)), ev$cluster), ]
save_table(ev, "tbl03", "cluster_evidence", "04_annotate.R")

# --- 手動確定ラベル（annotations.yaml・研究員が evidence を見て埋める）---
# 無ければ suggested_lineage を暫定に据え、確定は researcher に委ねる（標準: 手動解釈が最終権威）。
ann_path <- here::here("annotations.yaml")
final <- data.frame(cluster = ev$cluster,
                    lineage_final = ev$suggested_lineage,
                    malignant_final = ev$malignant_status,
                    state_final = NA_character_,
                    annotation_source = "suggested(auto)", stringsAsFactors = FALSE)
if (file.exists(ann_path)) {
  man <- yaml::read_yaml(ann_path)$clusters
  for (cl in names(man)) {
    i <- which(final$cluster == cl); if (length(i) != 1) next
    if (!is.null(man[[cl]]$lineage))          final$lineage_final[i]   <- man[[cl]]$lineage
    if (!is.null(man[[cl]]$malignant_status)) final$malignant_final[i] <- man[[cl]]$malignant_status
    if (!is.null(man[[cl]]$state))            final$state_final[i]     <- man[[cl]]$state
    final$annotation_source[i] <- "manual(annotations.yaml)"
  }
  cat(sprintf("annotations.yaml から %d クラスタの手動確定を反映\n", length(man)))
} else {
  cat("annotations.yaml 未作成 → suggested_lineage を暫定ラベルに使用。evidence(tbl03) を見て手動確定を埋めること。\n")
  cat("  例: annotations.yaml に  clusters:\\n    \"0\": { lineage: T, malignant_status: non_malignant }  の形で記録。\n")
}

annot <- merge(ev, final, by = "cluster", all = TRUE)
annot <- annot[order(suppressWarnings(as.numeric(annot$cluster)), annot$cluster), ]
save_table(annot, "tbl04", "annotation", "04_annotate.R")

# --- 確定ラベルを cell へ写像し注釈済み object を保存 ---
lut_lin <- setNames(final$lineage_final, final$cluster)
lut_mal <- setNames(final$malignant_final, final$cluster)
lut_rev <- setNames(ev$review_flag, ev$cluster)
cl_chr <- as.character(o$seurat_clusters)
o$annotation      <- unname(lut_lin[cl_chr])
o$malignant_status <- unname(lut_mal[cl_chr])
o$review_flag     <- unname(lut_rev[cl_chr])
saveRDS(o, here::here("data", "processed", "seurat_annotated.rds"))

# --- 図: 注釈 UMAP と review 対象の可視化 ---
save_fig(DimPlot(o, group.by = "annotation", label = TRUE) + ggtitle("annotation (lineage)"),
         "fig05", "umap_annotation", "04_annotate.R", width = 8, height = 6)
save_fig(DimPlot(o, group.by = "review_flag") + ggtitle("clusters flagged for review"),
         "fig06", "umap_review_flag", "04_annotate.R", width = 7, height = 6)

n_review <- sum(ev$review_flag, na.rm = TRUE)
cat(sprintf("ANNOTATE_DONE: %d clusters | review 要 %d（tbl04 の review_reason 参照）| malignant ゲート=%s\n",
            nrow(ev), n_review, ifelse(isTRUE(CONFIG$cnv_gate$enabled), CONFIG$cnv_gate$tool, "off")))
