# 04: salmon quant -> tximport(txOut=TRUE) を 2 通りの countsFromAbundance で実行
#     -> transcript-level SummarizedExperiment を data/processed/se_tx.rds に保存（正準 object）。
#
# なぜ txOut=TRUE / 2 通りの count 表現か（rnaseqDTU: Love/Soneson/Patro F1000）:
#   - txOut=TRUE         : gene 集約せず transcript 単位のまま保持する（DTE/DTU とも tx 単位で検定）。
#   - assay "counts"     = dtuScaledTPM   : TPM を「遺伝子内 isoform 間で比較可能な尺度」に整えた
#                          擬似 count。DRIMSeq/DEXSeq が期待する DTU（isoform 使用比）向け。gene 集約なし。
#   - assay "counts_dte" = lengthScaledTPM: 各転写産物を自身の長さで scale した count。DTE（発現
#                          「量」差）の妥当な count 表現。dtuScaledTPM は isoform 比率専用で DTE には使わない。
#   tx2gene（config.reference.tx2gene）で transcript_id -> gene_id を rowData に刻む（DRIMSeq の
#   gene グルーピングに使う）。② の gene-level se.rds とは別 object（名前 se_tx.rds で区別）。
# 実行（ルートから）:  Rscript analysis/04_import_tx.R
suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(tximport)
  library(SummarizedExperiment)
  library(S4Vectors)
})
source(here::here("R", "helpers.R"))     # `%||%`

cfg     <- yaml::read_yaml(here::here("config.yaml"))
ids     <- vapply(cfg$samples, function(s) s$id, character(1))
sf_dirs <- setNames(here::here("data", "interim", "salmon", ids), ids)
# files はサンプル名付きで持つ（file.path は名前属性を落とすため setNames で ids を付す）。
# こうすると txi$counts の列名がサンプル名になり colData の行順と整合する。
files   <- setNames(file.path(sf_dirs, "quant.sf"), ids)
if (!all(file.exists(files)))
  stop("salmon quant が未完（quant.sf 欠損）: ",
       paste(ids[!file.exists(files)], collapse = ", "), " — 先に 03_salmon_quant.sh を回す")

# tx2gene: transcript_id -> gene_id 対応表（DRIMSeq の gene グルーピングに使う）
t2g <- utils::read.table(cfg$reference$tx2gene, header = FALSE,
                         col.names = c("tx", "gene"), stringsAsFactors = FALSE)

# --- transcript-level import（gene 集約せず）---
# DTU 向け dtuScaledTPM: 「遺伝子内 isoform 間の median transcript length」で scale するため
# tx2gene が必須（gene グルーピングを渡す）。txOut=TRUE で出力は transcript 単位のまま。
txi <- tximport(files, type = "salmon", txOut = TRUE,
                countsFromAbundance = "dtuScaledTPM", tx2gene = t2g)
# DTE 向け lengthScaledTPM: 各転写産物を自身の長さで scale した count（DTE の妥当な count 表現）。
# 同じ tximport を countsFromAbundance だけ替えて呼び直す（tx2gene / txOut=TRUE は共通）。
txi_dte <- tximport(files, type = "salmon", txOut = TRUE,
                    countsFromAbundance = "lengthScaledTPM", tx2gene = t2g)

# --- colData: config.samples の group / covariate をそのまま刻む（fq1/fq2 は除く）---
coldata <- do.call(rbind, lapply(cfg$samples, function(s) {
  as.data.frame(s[setdiff(names(s), c("fq1", "fq2"))], stringsAsFactors = FALSE)
}))
rownames(coldata) <- coldata$id
coldata <- coldata[ids, , drop = FALSE]        # counts の列順（ids）に厳密整合

# --- rowData: transcript_id と対応 gene_id（tx2gene を counts の行順に合わせる）---
tx_ids  <- rownames(txi$counts)
gene_of <- t2g$gene[match(tx_ids, t2g$tx)]
if (any(is.na(gene_of)))
  warning("tx2gene に無い transcript が ", sum(is.na(gene_of)),
          " 件（gene_id=NA → DRIMSeq/DEXSeq では除外）。salmon index と tx2gene の整合を確認する。")
rowdata <- S4Vectors::DataFrame(tx_id = tx_ids, gene_id = gene_of)

# assay "counts"=dtuScaledTPM(DTU 用) / "counts_dte"=lengthScaledTPM(DTE 用) / abundance=TPM を共通保持。
# rowData(tx_id/gene_id)・colData は 3 assay 共通。counts_dte は txi の行/列順に厳密整合させる
# （同一 files/tx2gene/txOut ゆえ通常は同順だが防御的に再索引する）。
se <- SummarizedExperiment(
  assays  = list(counts     = txi$counts,
                 counts_dte = txi_dte$counts[rownames(txi$counts), colnames(txi$counts)],
                 abundance  = txi$abundance),
  colData = S4Vectors::DataFrame(coldata),
  rowData = rowdata)
colnames(se) <- ids

dir.create(here::here("data", "processed"), showWarnings = FALSE, recursive = TRUE)
saveRDS(se, here::here("data", "processed", "se_tx.rds"))
cat(sprintf("SE(tx): %d transcripts x %d samples -> data/processed/se_tx.rds\n",
            nrow(se), ncol(se)))
cat("IMPORT_TX_DONE\n")
