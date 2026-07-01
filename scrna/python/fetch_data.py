"""データ取得モジュール（取得系は python/ に集約する＝① polyglot レイアウト）。

GEO series の supplementary（RAW.tar 等）を data/raw へ取得・展開する。標準ライブラリのみ
（urllib / tarfile）で依存を増やさない。analysis/ から import しても CLI から叩いてもよい。

CLI:
    python python/fetch_data.py GSE134520          # data/raw/ へ RAW.tar を取得・展開
    python python/fetch_data.py --url <URL> --out data/raw/foo.tar

共有領域（/data/shared）由来を使う場合は、ダウンロードせず work-fetch で data/raw に配置する
（AGENTS-base 参照）。本モジュールは公開 GEO を直接取得する経路を担う。
"""
from __future__ import annotations

import argparse
import sys
import tarfile
import urllib.request
from pathlib import Path

# プロジェクトルート（python/ の 1 つ上）をルート相対解決の起点にする（実行位置非依存）。
ROOT = Path(__file__).resolve().parent.parent


def geo_supplement_url(geo_id: str) -> str:
    """GEO series ID から supplementary RAW.tar の FTP URL を組む。

    例) GSE134520 → https://ftp.ncbi.nlm.nih.gov/geo/series/GSE134nnn/GSE134520/suppl/GSE134520_RAW.tar
    """
    if not geo_id.startswith("GSE"):
        raise ValueError(f"GEO series ID は GSE で始まる想定: {geo_id!r}")
    num = geo_id[3:]
    stub = (num[:-3] + "nnn") if len(num) > 3 else "nnn"
    return (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE{stub}/{geo_id}"
        f"/suppl/{geo_id}_RAW.tar"
    )


def download(url: str, dest: Path) -> Path:
    """url を dest へダウンロードする（親ディレクトリは自動作成）。"""
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"[fetch] {url} -> {dest}")
    urllib.request.urlretrieve(url, dest)  # noqa: S310 (公開 https/ftp のみ想定)
    return dest


def extract_tar(tar_path: Path, out_dir: Path) -> None:
    """tar を out_dir へ展開する（GEO の RAW.tar 想定）。"""
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[fetch] extract {tar_path} -> {out_dir}")
    with tarfile.open(tar_path) as tf:
        # filter="data" はパストラバーサル対策（PEP 706・Python 3.12+ / 3.11.4+）。
        # 古い 3.11 パッチ版（3.11.0–3.11.3）は filter 引数を持たず TypeError になるため
        # フォールバックする（requires-python は >=3.11.4 で下限を切っている）。
        try:
            tf.extractall(out_dir, filter="data")
        except TypeError:
            tf.extractall(out_dir)


def fetch_geo(geo_id: str, raw_dir: Path | None = None) -> Path:
    """GEO series の RAW.tar を data/raw へ取得・展開して展開先を返す。"""
    raw_dir = raw_dir or (ROOT / "data" / "raw")
    tar_path = raw_dir / f"{geo_id}_RAW.tar"
    download(geo_supplement_url(geo_id), tar_path)
    extract_tar(tar_path, raw_dir)
    return raw_dir


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="GEO/URL から生データを data/raw へ取得する")
    ap.add_argument("geo_id", nargs="?", help="GEO series ID（例 GSE134520）")
    ap.add_argument("--url", help="任意 URL を直接ダウンロード（geo_id の代わり）")
    ap.add_argument("--out", default="data/raw", help="出力ディレクトリ（既定 data/raw）")
    args = ap.parse_args(argv)

    out_dir = ROOT / args.out if not Path(args.out).is_absolute() else Path(args.out)
    if args.url:
        dest = download(args.url, out_dir / Path(args.url).name)
        if dest.suffix in (".tar", ".tgz") or dest.name.endswith(".tar.gz"):
            extract_tar(dest, out_dir)
    elif args.geo_id:
        fetch_geo(args.geo_id, out_dir)
    else:
        ap.error("geo_id か --url のいずれかを指定する")
    print("FETCH_DONE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
