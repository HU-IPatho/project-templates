#!/usr/bin/env bash
# データ取得テンプレート。自分のデータに合わせて下のいずれかを有効化する。
# 生データは data/raw/ に置く（.gitignore 対象・git には入らない）。
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p data/raw

# --- パターンA: 共有データ（HDD）から取り込む（推奨・最速）---
# 共有データセットは /data/shared/datasets/ にある。まず registry.tsv で探す。
#   work-fetch /data/shared/datasets/<dataset>/
#   cp -r /work/<dataset>/* data/raw/

# --- パターンB: 公開データを直接ダウンロード（例: GEO）---
# 例）GSE134520（dogfood 参照実装と同じ胃がん scRNA・dense 形式）:
#   cd data/raw
#   wget -q -O GSE134520_RAW.tar "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE134nnn/GSE134520/suppl/GSE134520_RAW.tar"
#   tar xf GSE134520_RAW.tar
#   cd -

echo "data/raw の中身:"; ls -1 data/raw 2>/dev/null | head
echo "FETCH_DONE"
