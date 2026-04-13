# RNA-seq 解析プロジェクト

## 概要

このプロジェクトは RNA-seq データの解析を行う。

## 解析パイプライン

1. QC: FastQC → MultiQC
2. アラインメント: STAR or HISAT2
3. 定量: featureCounts or Salmon
4. 差分発現解析: DESeq2
5. 可視化: ggplot2, pheatmap, EnhancedVolcano

## ディレクトリ規約

- `data/raw/` — 生データ（FASTQ）。git 管理外。`work-fetch` で SSD にコピーして使用
- `data/reference/` — リファレンスゲノム・アノテーション。git 管理外
- `data/processed/` — 前処理済みデータ（BAM, counts）。git 管理外
- `src/` — 解析スクリプト。番号付き（01_, 02_, ...）で実行順序を明示
- `results/` — 解析結果（テーブル、統計量）
- `reports/` — レポート・図表

## データ操作コマンド

```bash
work-fetch /data/shared/datasets/<data>/   # HDD → SSD コピー
work-save results/                         # SSD → HDD 保存
work-share results/se.rds --message "説明" # 共有リクエスト
work-status                                # SSD 使用量確認
```

## 言語・ツール

- R (renv で依存管理)
- Bioconductor: DESeq2, SummarizedExperiment, GenomicFeatures
- CLI: samtools, STAR, featureCounts
