"""公開 RNA-seq FASTQ の取得ロジック（SRA/ENA accession → data/raw）。

方針: 取り込んだ生データは data/raw（不変入力層）に置く（spec の 3 層データレイアウト）。
既に data/raw に FASTQ があるなら fetch は不要。共有領域（/data/shared）由来の場合は
work-fetch で data/raw に取り込む（AGENTS.md「共有連携」参照）— 本モジュールは公開 DB 用。

実体の取得は omics-dev イメージ同梱の CLI に委譲する:
  - ffq        : accession → FASTQ の DL URL 解決
  - pysradb    : SRA メタデータ/実行 accession 解決
  - prefetch / fasterq-dump (sra-tools) : SRR の FASTQ 取得
ここでは「取得計画の解決」と「data/raw への配置」を薄くラップし、ネットワーク実行は
呼出し側（長時間なら job-run）に委ねる（テンプレは環境非依存の骨格に留める）。
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from python.config import ROOT, load_config

RAW_DIR = ROOT / "data" / "raw"


def _have(tool: str) -> bool:
    return shutil.which(tool) is not None


def resolve_commands(accessions: list[str]) -> list[list[str]]:
    """accession ごとの取得コマンド列を組み立てて返す（実行はしない＝dry-run 可能）。

    sra-tools が使えれば prefetch→fasterq-dump、無ければ ffq で URL を解決する。
    """
    cmds: list[list[str]] = []
    for acc in accessions:
        if _have("fasterq-dump"):
            cmds.append(["fasterq-dump", "--split-files", "--outdir", str(RAW_DIR), acc])
        elif _have("ffq"):
            cmds.append(["ffq", "--ftp", acc])
        else:
            raise RuntimeError(
                "fasterq-dump も ffq も見つからない。omics-dev イメージ内で実行しているか確認する。"
            )
    return cmds


def fetch(accessions: list[str], dry_run: bool = False) -> None:
    """accession を data/raw に取得する。dry_run=True なら計画を表示するだけ。"""
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    if not accessions:
        print("[fetch] config.fetch.accessions が空 — data/raw の既存 FASTQ を使う")
        return
    for cmd in resolve_commands(accessions):
        print("[fetch]", " ".join(cmd))
        if not dry_run:
            subprocess.run(cmd, check=True)


def main(dry_run: bool = False) -> None:
    cfg = load_config()
    accessions = list(cfg.get("fetch", {}).get("accessions", []) or [])
    fetch(accessions, dry_run=dry_run)
    print("FETCH_DONE")
