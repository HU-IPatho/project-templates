# ============================================================================
# DTU（differential transcript usage）— 遺伝子内 isoform「使用比」の群間変化
# ----------------------------------------------------------------------------
# 既定エンジン = DRIMSeq + stageR 二段確認（rnaseqDTU 準拠）。
#   Love, Soneson, Patro (2018) F1000Research "Swimming downstream" に忠実:
#     DRIMSeq: dmDSdata -> dmFilter -> dmPrecision -> dmFit -> dmTest
#     stageR : 二段の OFDR（overall false discovery rate）制御
#       screening    = 遺伝子（DRIMSeq gene-level p 値）
#       confirmation = transcript（DRIMSeq feature-level p 値）
#       method="dtu", alpha = config.dtu.fdr（screening を通った遺伝子内でのみ tx を確認）
#   count ベースなので bootstrap 不要（bootstrap が要るのは任意の swish のみ）。
#
# dmFilter の発現閾値は config.dtu（min_gene_expr / min_feature_expr / min_feature_prop）。
#   標本数依存の min_samps_*（何標本で閾値を満たすか）は実行時に n / n.small から決める
#   （rnaseqDTU の推奨式: n=総標本数, n.small=小さい群の標本数）。
#
# 各 contrast は [num, den] の 2 群に subset し、~condition（den 基準）の 2 群比較に落とす。
# 任意エンジン: swish（fishpond・inferential replicates 必須）/ dexseq（DEXSeq+stageR）。
# ============================================================================
suppressPackageStartupMessages({
  library(SummarizedExperiment)
})

# DRIMSeq/DEXSeq が testable でない gene/feature に返す NA p 値は 1 に置換（rnaseqDTU 準拠）。
.no_na <- function(x) ifelse(is.na(x), 1, x)

# --- contrast の 2 群を subset し dmFilter まで済ませた dmDSdata を返す（共通前処理）----
# DRIMSeq と DEXSeq で同一のフィルタ規則を共有し、drift を防ぐ。
.filtered_drim <- function(se, cfg, num, den) {
  gcol    <- cfg$design$group_col
  coldata <- as.data.frame(SummarizedExperiment::colData(se))
  rd      <- as.data.frame(SummarizedExperiment::rowData(se))
  sel     <- rownames(coldata)[coldata[[gcol]] %in% c(num, den)]
  cts     <- as.matrix(SummarizedExperiment::assay(se, "counts"))[, sel, drop = FALSE]

  # gene_id 不明（tx2gene に無い）transcript は DTU 対象外 → 除外
  keep_feat <- !is.na(rd$gene_id)
  counts_df <- data.frame(gene_id    = rd$gene_id[keep_feat],
                          feature_id = rd$tx_id[keep_feat],
                          cts[keep_feat, , drop = FALSE],
                          check.names = FALSE, stringsAsFactors = FALSE)
  # den を基準（reference）に置く → model.matrix の 2 列目が num の効果
  cond       <- factor(coldata[sel, gcol], levels = c(den, num))
  samples_df <- data.frame(sample_id = sel, condition = cond, stringsAsFactors = FALSE)

  d <- DRIMSeq::dmDSdata(counts = counts_df, samples = samples_df)
  n <- length(sel); n.small <- min(table(cond))     # 標本数依存の min_samps_* を決める
  d <- DRIMSeq::dmFilter(
    d,
    min_samps_gene_expr    = n,
    min_gene_expr          = cfg$dtu$min_gene_expr    %||% 10,
    min_samps_feature_expr = n.small,
    min_feature_expr       = cfg$dtu$min_feature_expr %||% 10,
    min_samps_feature_prop = n.small,
    min_feature_prop       = cfg$dtu$min_feature_prop %||% 0.1)
  list(d = d, cond = cond, sel = sel)
}

# --- stageR 二段確認の共通ヘルパー（DRIMSeq / DEXSeq 双方から呼ぶ）-----------------
# pScreen: 遺伝子名付き screening p（DRIMSeq=生 p→adjusted=FALSE / DEXSeq=perGeneQValue→TRUE）。
# pConf  : transcript 名を rownames に持つ confirmation p 値の 1 列行列。
.stager_dtu <- function(pScreen, pConf, tx2gene_df, fdr, screen_adjusted) {
  obj <- stageR::stageRTx(pScreen = pScreen, pConfirmation = pConf,
                          pScreenAdjusted = screen_adjusted, tx2gene = tx2gene_df)
  obj <- stageR::stageWiseAdjustment(obj, method = "dtu", alpha = fdr)
  # 列: geneID, txID, gene(screening OFDR), transcript(confirmation OFDR)
  stageR::getAdjustedPValues(obj, order = FALSE, onlySignificantGenes = FALSE)
}

# --- 既定エンジン: DRIMSeq + stageR（全 contrast）--------------------------------
run_dtu <- function(se, cfg) {
  suppressPackageStartupMessages({ library(DRIMSeq); library(stageR) })
  fdr <- cfg$dtu$fdr %||% 0.05

  out <- lapply(names(cfg$contrasts), function(nm) {
    spec <- cfg$contrasts[[nm]]; num <- spec[[1]]; den <- spec[[2]]
    f <- .filtered_drim(se, cfg, num, den)
    d <- f$d
    if (length(f$sel) < 4)
      warning("DTU ", nm, ": 標本 ", length(f$sel), " 件（各群 n>=2 推奨）")

    # den を基準に固定した design（samples(d) が condition を character 化しても正しい参照に）
    sdf <- DRIMSeq::samples(d)
    sdf$condition <- factor(sdf$condition, levels = c(den, num))
    design_full <- stats::model.matrix(~ condition, data = sdf)
    coef <- colnames(design_full)[2]                # "condition<num>"（num vs den の効果）

    set.seed(1)                                     # dmPrecision の初期化を再現可能に
    d <- DRIMSeq::dmPrecision(d, design = design_full)
    d <- DRIMSeq::dmFit(d, design = design_full)
    d <- DRIMSeq::dmTest(d, coef = coef)

    res_gene <- DRIMSeq::results(d)                    # screening（gene-level p）
    res_tx   <- DRIMSeq::results(d, level = "feature") # confirmation（transcript-level p）
    res_gene$pvalue <- .no_na(res_gene$pvalue)
    res_tx$pvalue   <- .no_na(res_tx$pvalue)

    # stageR 二段: screening=gene の DRIMSeq 生 p（adjusted=FALSE）→ confirmation=tx
    pScreen    <- stats::setNames(res_gene$pvalue, res_gene$gene_id)
    pConf      <- matrix(res_tx$pvalue, ncol = 1,
                         dimnames = list(res_tx$feature_id, NULL))
    tx2gene_df <- res_tx[, c("feature_id", "gene_id")]
    padj <- .stager_dtu(pScreen, pConf, tx2gene_df, fdr, screen_adjusted = FALSE)

    # DRIMSeq 生 p と stageR OFDR を transcript 行で結合した結果表
    tab <- merge(
      data.frame(gene_id    = res_tx$gene_id,
                 feature_id = res_tx$feature_id,
                 p_gene = res_gene$pvalue[match(res_tx$gene_id, res_gene$gene_id)],
                 p_tx   = res_tx$pvalue,
                 stringsAsFactors = FALSE),
      data.frame(gene_id    = padj$geneID,
                 feature_id = padj$txID,
                 gene_padj  = padj$gene,        # screening OFDR
                 tx_padj    = padj$transcript,  # confirmation OFDR
                 stringsAsFactors = FALSE),
      by = c("gene_id", "feature_id"), all.x = TRUE)
    tab$significant <- !is.na(tab$tx_padj) & tab$tx_padj < fdr
    tab$contrast    <- nm
    tab <- tab[order(tab$gene_padj, tab$tx_padj), ]

    # 主役図（06）用: screening OFDR で有意な DTU 遺伝子（gene_padj 昇順）
    sig_genes <- unique(tab$gene_id[!is.na(tab$gene_padj) & tab$gene_padj < fdr])
    list(table = tab, sig_genes = sig_genes)
  })
  names(out) <- names(cfg$contrasts)
  out
}

# --- 任意エンジン: swish（fishpond・inferential replicates が必須）-----------------
# salmon の bootstrap（config.salmon.num_bootstraps>0）由来の inferential replicates を使い、
# 不確実性を織り込んだ DTU（isoform 使用比）検定を行う。fishpond の DTU レシピに忠実:
#   gene_id を mcols に付与 → scaleInfReps → labelKeep → keep でフィルタ →
#   isoformProportions（gene 内 isoform 比率へ変換）→ swish（x=condition で比率の群差を検定）。
# swish の検定対象を「発現量」でなく「isoform 比率」にするのが DTU の要（isoformProportions が核）。
# bootstrap 未取得なら fail-fast（誤用防止）。
run_swish <- function(cfg) {
  if ((cfg$salmon$num_bootstraps %||% 0) < 1) {
    stop("swish は inferential replicates（salmon bootstrap）が必須。",
         "config.salmon.num_bootstraps を >0 にして analysis/03_salmon_quant.sh を再実行する。")
  }
  suppressPackageStartupMessages({
    library(tximeta); library(fishpond); library(SummarizedExperiment); library(here)
  })
  gcol <- cfg$design$group_col
  ids  <- vapply(cfg$samples, function(s) s$id, character(1))
  files <- stats::setNames(
    file.path(here::here("data", "interim", "salmon", ids), "quant.sf"), ids)
  coldata <- do.call(rbind, lapply(cfg$samples, function(s) {
    as.data.frame(s[setdiff(names(s), c("fq1", "fq2"))], stringsAsFactors = FALSE)
  }))
  rownames(coldata) <- coldata$id
  coldata$names <- coldata$id
  coldata$files <- files[coldata$id]
  # skipMeta=TRUE で linkedTxome 不要のまま inferential replicates 付きで tx-level import
  y0 <- tximeta::tximeta(coldata, type = "salmon", txOut = TRUE, skipMeta = TRUE)

  # isoformProportions に必須の gene_id を mcols へ付与（tx2gene を rownames に match）。
  # gene_id 欠損 tx（tx2gene に無い）は DTU 対象外 → 除外（04 の rowData 割当と同じ流儀）。
  t2g <- utils::read.table(cfg$reference$tx2gene, header = FALSE,
                           col.names = c("tx", "gene"), stringsAsFactors = FALSE)
  SummarizedExperiment::mcols(y0)$gene_id <- t2g$gene[match(rownames(y0), t2g$tx)]
  y0 <- y0[!is.na(SummarizedExperiment::mcols(y0)$gene_id), ]

  out <- lapply(names(cfg$contrasts), function(nm) {
    spec <- cfg$contrasts[[nm]]; num <- spec[[1]]; den <- spec[[2]]
    sel <- colnames(y0)[SummarizedExperiment::colData(y0)[[gcol]] %in% c(num, den)]
    y <- y0[, sel]
    y$condition <- factor(SummarizedExperiment::colData(y)[[gcol]], levels = c(den, num))
    # fishpond 標準 DTU 順序: scale してから label→keep フィルタ（scale が infRep を整える）。
    y   <- fishpond::scaleInfReps(y)
    y   <- fishpond::labelKeep(y)
    y   <- y[SummarizedExperiment::mcols(y)$keep, ]
    iso <- fishpond::isoformProportions(y)         # gene 内 isoform 比率へ変換（DTU の核）
    set.seed(1)                                    # swish の permutation を再現可能に
    iso <- fishpond::swish(iso, x = "condition")   # 比率の群差を検定（x=design の group 列）
    mc  <- as.data.frame(SummarizedExperiment::mcols(iso))
    # 他エンジンと揃う列（feature_id/gene_id/stat/pvalue/padx=qvalue）に整形
    data.frame(feature_id = rownames(iso),
               gene_id    = mc$gene_id,
               stat       = mc$stat,
               log2FC     = mc$log2FC,
               pvalue     = mc$pvalue,
               padx       = mc$qvalue,
               contrast   = nm, stringsAsFactors = FALSE)
  })
  names(out) <- names(cfg$contrasts)
  out
}

# --- 任意エンジン: DEXSeq + stageR（DRIMSeq と同じフィルタを共有）------------------
run_dexseq <- function(se, cfg) {
  suppressPackageStartupMessages({ library(DEXSeq); library(DRIMSeq); library(stageR) })
  fdr <- cfg$dtu$fdr %||% 0.05

  out <- lapply(names(cfg$contrasts), function(nm) {
    spec <- cfg$contrasts[[nm]]; num <- spec[[1]]; den <- spec[[2]]
    f  <- .filtered_drim(se, cfg, num, den)
    cd <- DRIMSeq::counts(f$d)                       # gene_id, feature_id, <sample 列...>
    sample.data <- DRIMSeq::samples(f$d)
    sample.data$condition <- factor(sample.data$condition, levels = c(den, num))
    count.data  <- round(as.matrix(cd[, -c(1, 2)]))

    dxd <- DEXSeq::DEXSeqDataSet(
      countData = count.data, sampleData = sample.data,
      design    = ~ sample + exon + condition:exon,
      featureID = cd$feature_id, groupID = cd$gene_id)
    dxd <- DEXSeq::estimateSizeFactors(dxd)
    dxd <- DEXSeq::estimateDispersions(dxd, quiet = TRUE)
    dxd <- DEXSeq::testForDEU(dxd, reducedModel = ~ sample + exon)
    dxr <- DEXSeq::DEXSeqResults(dxd, independentFiltering = FALSE)

    # stageR: screening は perGeneQValue（既に補正済み → pScreenAdjusted=TRUE）
    pConf      <- matrix(.no_na(dxr$pvalue), ncol = 1,
                         dimnames = list(dxr$featureID, NULL))
    pScreen    <- DEXSeq::perGeneQValue(dxr)
    tx2gene_df <- data.frame(feature_id = dxr$featureID, gene_id = dxr$groupID,
                             stringsAsFactors = FALSE)
    padj <- .stager_dtu(pScreen, pConf, tx2gene_df, fdr, screen_adjusted = TRUE)

    tab <- data.frame(gene_id = padj$geneID, feature_id = padj$txID,
                      gene_padj = padj$gene, tx_padj = padj$transcript,
                      contrast = nm, stringsAsFactors = FALSE)
    tab$significant <- !is.na(tab$tx_padj) & tab$tx_padj < fdr
    tab <- tab[order(tab$gene_padj, tab$tx_padj), ]
    sig_genes <- unique(tab$gene_id[!is.na(tab$gene_padj) & tab$gene_padj < fdr])
    list(table = tab, sig_genes = sig_genes)
  })
  names(out) <- names(cfg$contrasts)
  out
}
