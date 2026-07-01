# 04: salmon quant -> tximeta（来歴付インポート）-> gene-level SE -> data/processed/se.rds
#     さらに QC（library size / mapping rate / PCA）を計算し出力ハーネスで図表保存する。
#     data/processed/se.rds は下流 ②（rnaseq-secondary）の入力になる。
# 実行（ルートから）:  Rscript analysis/04_tximeta_se.R
suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(tximeta)
  library(tximport)
  library(SummarizedExperiment)
})
source(here::here("R", "helpers.R"))
source(here::here("R", "qc_se.R"))

cfg     <- yaml::read_yaml(here::here("config.yaml"))
ids     <- vapply(cfg$samples, function(s) s$id, character(1))
sf_dirs <- setNames(here::here("data", "interim", "salmon", ids), ids)
# files はサンプル名付きで持つ（file.path は名前属性を落とすため setNames で ids を付す）。
# こうすると tximport フォールバック時に txi$counts の列名がサンプル名になり QC の sample キーが一致する。
files   <- setNames(file.path(sf_dirs, "quant.sf"), ids)
if (!all(file.exists(files)))
  stop("salmon quant が未完（quant.sf 欠損）: ",
       paste(ids[!file.exists(files)], collapse = ", "), " — 先に 03_salmon_quant.sh を回す")

coldata <- data.frame(names = ids, files = files, stringsAsFactors = FALSE)

# --- 来歴担保: linkedTxome を登録（index を作った参照 FASTA/GTF と一致させる）---
# これで tximeta が transcriptome の由来（source/organism/release/genome）を SE に刻む。
tx <- cfg$reference$txome
makeLinkedTxome(
  indexDir = cfg$reference$salmon_index,
  source   = tx$source, organism = tx$organism, release = as.character(tx$release),
  genome   = tx$genome, fasta = tx$fasta, gtf = tx$gtf, write = FALSE
)

# --- transcript-level import（provenance 付）-> gene-level へ集約 ---
se <- tryCatch({
  se_tx <- tximeta(coldata)                 # tx レベル SE（来歴メタ付）
  summarizeToGene(se_tx)                     # gene レベルへ集約（linkedTxome の GTF 由来 tx2gene）
}, error = function(e) {
  # フォールバック: tximeta が解決できない環境では config.reference.tx2gene で tximport 集約
  message("tximeta 失敗 -> tximport(tx2gene) にフォールバック: ", conditionMessage(e))
  t2g <- utils::read.table(cfg$reference$tx2gene, header = FALSE,
                           col.names = c("tx", "gene"), stringsAsFactors = FALSE)
  txi <- tximport(files, type = "salmon", tx2gene = t2g)
  # 各 assay に列名(サンプル名)を明示する。file.path が名前属性を落とす場合でも
  # colnames(se) を primary 経路(tximeta が coldata$names を colnames に刻む)と揃え、
  # 直後の compute_se_qc() の sample キー(merge by="sample")を一致させる。
  colnames(txi$counts)    <- ids
  colnames(txi$abundance) <- ids
  rownames(coldata)       <- ids
  SummarizedExperiment(assays = list(counts = txi$counts, abundance = txi$abundance),
                       colData = DataFrame(coldata))
})

dir.create(here::here("data", "processed"), showWarnings = FALSE, recursive = TRUE)
saveRDS(se, here::here("data", "processed", "se.rds"))
cat(sprintf("SE: %d genes x %d samples -> data/processed/se.rds\n", nrow(se), ncol(se)))

# --- QC 計算（重い計算は script 側・保存済み object をレポートが読む）---
mapping <- read_salmon_mapping(sf_dirs)
qc      <- compute_se_qc(se, mapping_df = mapping)
saveRDS(qc, here::here("data", "processed", "qc.rds"))

# --- QC 図表を出力ハーネスで保存（PNG+PDF・captions.tsv 記録）---
save_fig(plot_lib_size(qc),     num = 1, desc = "library_size", script = "04_tximeta_se.R")
save_fig(plot_mapping_rate(qc), num = 2, desc = "mapping_rate", script = "04_tximeta_se.R")
save_fig(plot_pca(qc),          num = 3, desc = "pca",          script = "04_tximeta_se.R")
save_table(qc$summary,          num = 1, desc = "qc_summary",   script = "04_tximeta_se.R")

cat("TXIMETA_SE_DONE\n")
