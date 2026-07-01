#!/usr/bin/env bash
# 02: fastp でアダプタトリム/品質フィルタ（data/raw -> data/interim）→ FastQC。
# system ツール（fastp/fastqc）は omics-dev イメージが版固定担体（config でなくイメージタグで pin）。
# 実行（どこから叩いてもルートに解決）:  bash analysis/02_qc_fastp.sh
set -euo pipefail

# --- ルートに解決（spec: ルート実行・実行位置非依存）---
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

CFG() { python3 python/config.py "$@"; }
THREADS="$(CFG get fastp.threads)"
EXTRA="$(CFG get fastp.extra_flags)"      # 例: --detect_adapter_for_pe
PAIRED="$(CFG get project.paired_end)"
# samples は先に変数へ捕捉する（プロセス置換 < <(CFG samples) は set -e の対象外で
# config.py の列挙失敗を握り潰すため使わない）。command substitution 代入なので
# config.py が非零終了すれば set -e が発火して fail-closed に即停止する。
SAMPLES="$(CFG samples)"

TRIM_DIR="data/interim/trimmed"
FASTP_DIR="data/interim/fastp"            # fastp の json/html（MultiQC が集約）
FASTQC_DIR="data/interim/fastqc"          # FastQC 出力（MultiQC が集約）
mkdir -p "$TRIM_DIR" "$FASTP_DIR" "$FASTQC_DIR"

# config.samples を id<TAB>fq1<TAB>fq2 で受け取り、1 サンプルずつ処理する
while IFS=$'\t' read -r ID FQ1 FQ2; do
  [ -z "${ID:-}" ] && continue
  echo "[fastp] $ID"
  if [ "$PAIRED" = "true" ]; then
    fastp --thread "$THREADS" $EXTRA \
      -i "$FQ1" -I "$FQ2" \
      -o "$TRIM_DIR/${ID}_R1.trimmed.fastq.gz" \
      -O "$TRIM_DIR/${ID}_R2.trimmed.fastq.gz" \
      --json "$FASTP_DIR/${ID}.fastp.json" \
      --html "$FASTP_DIR/${ID}.fastp.html"
    fastqc -t "$THREADS" -o "$FASTQC_DIR" \
      "$TRIM_DIR/${ID}_R1.trimmed.fastq.gz" "$TRIM_DIR/${ID}_R2.trimmed.fastq.gz"
  else
    fastp --thread "$THREADS" $EXTRA \
      -i "$FQ1" \
      -o "$TRIM_DIR/${ID}.trimmed.fastq.gz" \
      --json "$FASTP_DIR/${ID}.fastp.json" \
      --html "$FASTP_DIR/${ID}.fastp.html"
    fastqc -t "$THREADS" -o "$FASTQC_DIR" "$TRIM_DIR/${ID}.trimmed.fastq.gz"
  fi
done <<< "$SAMPLES"

echo "QC_FASTP_DONE"
