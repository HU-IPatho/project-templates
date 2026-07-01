#!/usr/bin/env bash
# データ取得テンプレート。counts 行列を用意し data/ に置く（.gitignore 対象・git 非追跡）。
# 生の入力は data/raw/（不変層）へ、正規化した gene×sample counts は data/interim/ へ。
# 最終的に config.yaml の counts_file（既定 data/interim/counts.tsv）を指すこと。
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p data/raw data/interim

# --- パターンA: 共有データ（HDD）から取り込む（推奨・最速）---
# 共有データセットは /data/shared/datasets/ にある。registry.tsv で探す。
#   work-fetch /data/shared/datasets/<dataset>/
#   cp -r /work/<dataset>/matrix/* data/raw/

# --- パターンB: 公開データを直接ダウンロード（例: GEO の count matrix）---
#   cd data/raw
#   wget -q -O counts_raw.tsv.gz "<URL>"
#   gunzip counts_raw.tsv.gz
#   cd -

# --- 正規化: 入手形式を gene×sample TSV へ（python/fetch_counts.py）---
# featureCounts 出力の例:
#   python python/fetch_counts.py data/raw/featurecounts.txt --format featurecounts \
#     -o data/interim/counts.tsv
# 既に gene×sample の TSV/CSV なら:
#   python python/fetch_counts.py data/raw/counts_raw.tsv --format tsv -o data/interim/counts.tsv

echo "data/interim の中身:"; ls -1 data/interim 2>/dev/null | head
echo "FETCH_DONE"
