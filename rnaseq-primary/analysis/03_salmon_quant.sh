#!/usr/bin/env bash
# 03: salmon で定量（selective alignment + decoy-aware index + --gcBias）→ MultiQC 集約。
#   - 参照 index（decoy-aware）は config.reference.salmon_index（管理者が /work/shared に事前構築・RO）。
#   - 定量は 1 サンプルが長時間 → job-run で detach 実行し job-wait で待つ（前景直実行は不可）。
#   - 最後に MultiQC で fastp/FastQC/salmon のログを集約し reports/ に HTML を出す。
# 実行（どこから叩いてもルートに解決）:  bash analysis/03_salmon_quant.sh
set -euo pipefail

# --- ルートに解決（spec: ルート実行・実行位置非依存）---
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

CFG() { python3 python/config.py "$@"; }
INDEX="$(CFG get reference.salmon_index)"
LIBTYPE="$(CFG get salmon.libtype)"
THREADS="$(CFG get salmon.threads)"
EXTRA="$(CFG get salmon.extra_flags)"     # --gcBias --seqBias --validateMappings（契約: --gcBias 必須）
PAIRED="$(CFG get project.paired_end)"
# samples は先に変数へ捕捉する（プロセス置換 < <(CFG samples) は set -e の対象外で
# config.py の列挙失敗を握り潰すため使わない）。command substitution 代入なので
# config.py が非零終了すれば set -e が発火して fail-closed に即停止する。
SAMPLES="$(CFG samples)"

TRIM_DIR="data/interim/trimmed"
QUANT_DIR="data/interim/salmon"
mkdir -p "$QUANT_DIR"

[ -d "$INDEX" ] || echo "[warn] salmon index が見つからない: $INDEX（管理者が /work/shared に構築済みか確認）" >&2

# サンプルごとに salmon quant を job-run で直列実行（同一リソースを --exclusive で直列化）
while IFS=$'\t' read -r ID FQ1 FQ2; do
  [ -z "${ID:-}" ] && continue
  OUT="$QUANT_DIR/$ID"
  if [ "$PAIRED" = "true" ]; then
    READS=(-1 "$TRIM_DIR/${ID}_R1.trimmed.fastq.gz" -2 "$TRIM_DIR/${ID}_R2.trimmed.fastq.gz")
  else
    READS=(-r "$TRIM_DIR/${ID}.trimmed.fastq.gz")
  fi
  echo "[salmon] $ID -> $OUT"
  # selective alignment: -i decoy-aware index / -l A(自動判定) / $EXTRA に --gcBias を含む
  JID="$(job-run --exclusive salmon-quant --label "salmon-$ID" -- \
    salmon quant -i "$INDEX" -l "$LIBTYPE" "${READS[@]}" \
      -p "$THREADS" $EXTRA -o "$OUT")"
  # job-wait: 124=未完→再待機 / 0=成功 / それ以外=失敗
  # 終了コードは job-wait 直後に捕捉する（`if ...; then break; fi` の後で $? を拾うと
  # else 無し if 複合の終了ステータス 0 を拾ってしまい 124 再待機分岐が死ぬため）。
  while :; do
    rc=0
    job-wait "$JID" --timeout 570 || rc=$?
    [ "$rc" -eq 0 ]   && break
    [ "$rc" -eq 124 ] && continue
    echo "[salmon] $ID 失敗 (rc=$rc)"; job-logs "$JID" -n 100 || true; exit 1
  done
done <<< "$SAMPLES"

# MultiQC: fastp/FastQC/salmon のログを 1 つの HTML に集約 → reports/（契約: 集約 HTML は reports 配下）
mkdir -p reports/multiqc
multiqc -f -o reports/multiqc data/interim
echo "[multiqc] reports/multiqc/multiqc_report.html"

echo "SALMON_QUANT_DONE"
