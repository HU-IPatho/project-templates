#!/usr/bin/env python3
"""counts 行列の取得・整形（最小 Python モジュール・analysis/ から呼ぶ or 単体実行）。

②パイプラインの正準 object（SummarizedExperiment）は R 側（analysis/01_build_se.R）が
gene×sample の counts 行列 TSV から構築する。本モジュールは「入手した多様な形式の
count データを、その TSV 形式（1 列目=gene_id・ヘッダ=sample_id）に正規化する」取得変換層。

想定入力の例:
  - featureCounts の出力（先頭にコメント行 + Geneid/Chr/Start/... 列 + sample 列）
  - salmon/tximport 由来の gene-level counts（別途 R で集計した行列）
  - GEO の supplementary count matrix（gene×sample の TSV/CSV）

出力: data/interim/counts.tsv（config.yaml の counts_file と一致させる）。
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def load_featurecounts(path: Path) -> pd.DataFrame:
    """featureCounts の出力を gene×sample の counts に整形する。

    先頭の '# Program:...' コメント行を飛ばし、Geneid を index に、
    メタ列（Chr/Start/End/Strand/Length）を落として sample 列だけ残す。
    """
    df = pd.read_csv(path, sep="\t", comment="#")
    meta_cols = ["Chr", "Start", "End", "Strand", "Length"]
    df = df.set_index("Geneid").drop(columns=[c for c in meta_cols if c in df.columns])
    # featureCounts の列名はしばしば bam パス → basename の stem に短縮
    df.columns = [Path(c).name.replace(".bam", "") for c in df.columns]
    return df


def load_matrix(path: Path, sep: str) -> pd.DataFrame:
    """既に gene×sample になっている count 行列（1 列目=gene_id）を読む。"""
    df = pd.read_csv(path, sep=sep, index_col=0)
    return df


def write_counts(df: pd.DataFrame, out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    df.index.name = "gene_id"
    df.to_csv(out, sep="\t")
    print(f"COUNTS_DONE: {df.shape[0]} genes x {df.shape[1]} samples -> {out}")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="count データを gene×sample TSV に正規化")
    p.add_argument("input", type=Path, help="入力ファイル")
    p.add_argument("-o", "--output", type=Path, default=Path("data/interim/counts.tsv"))
    p.add_argument("--format", choices=["featurecounts", "tsv", "csv"], default="tsv")
    args = p.parse_args(argv)

    if args.format == "featurecounts":
        df = load_featurecounts(args.input)
    else:
        sep = "," if args.format == "csv" else "\t"
        df = load_matrix(args.input, sep)

    write_counts(df, args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
