# scRNA-seq 解析プロジェクト

## 概要

このプロジェクトは single-cell RNA-seq データの解析を行う。

## 解析パイプライン

1. Cell Ranger: FASTQ → count matrix
2. QC/Filter: nFeature, nCount, percent.mt によるフィルタリング
3. 正規化: SCTransform or LogNormalize → PCA → UMAP
4. クラスタリング: FindNeighbors → FindClusters
5. マーカー遺伝子: FindAllMarkers
6. 細胞型アノテーション: SingleR or 手動アノテーション

## ディレクトリ規約

- `data/raw/` — Cell Ranger 出力（filtered_feature_bc_matrix 等）。git 管理外
- `data/reference/` — リファレンスデータ（refdata-gex-GRCh38 等）。git 管理外
- `data/processed/` — 前処理済みデータ（Seurat/SCE オブジェクト）。git 管理外
- `src/` — 解析スクリプト。番号付きで実行順序を明示
- `results/` — 解析結果（マーカー遺伝子テーブル、クラスター情報）
- `reports/` — レポート・UMAP プロット等

## データ操作コマンド

```bash
work-fetch /data/shared/datasets/<data>/   # HDD → SSD コピー
work-save results/                         # SSD → HDD 保存
work-share results/seurat.rds --message "説明" # 共有リクエスト
work-status                                # SSD 使用量確認
```

## 言語・ツール

- R (renv で依存管理)
- Seurat v5, SingleCellExperiment
- Cell Ranger (10x Genomics)
- scran, scater, SingleR
