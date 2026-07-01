#!/usr/bin/env python3
"""01: 公開 FASTQ を data/raw に取得する（config.fetch.accessions を使う）。

実行（プロジェクトルートから）:
    python3 analysis/01_fetch.py            # 取得実行
    python3 analysis/01_fetch.py --dry-run  # 取得計画だけ表示

生データは data/raw（不変入力層）に着地する。共有領域由来のデータは代わりに
work-fetch で data/raw に取り込む（AGENTS.md「共有連携」）。取得が長時間なら job-run で回す:
    JID=$(job-run --label fetch -- python3 analysis/01_fetch.py); job-wait "$JID" --timeout 570
"""
import sys
from pathlib import Path

# プロジェクトルートを import パスに追加（ルート実行前提・実行位置非依存）
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from python.fetch import main  # noqa: E402

if __name__ == "__main__":
    main(dry_run="--dry-run" in sys.argv)
