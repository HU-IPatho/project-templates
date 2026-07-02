"""config.yaml リーダ（層 3 の事実データへの単一アクセス点）。

Python からは ``from python.config import load_config`` で dict を得る。
shell スクリプト（02/03）からは CLI で値を取り出す:

    python3 python/config.py get salmon.threads         # スカラを 1 行出力
    python3 python/config.py get salmon.num_bootstraps   # bootstrap 数（swish 用・既定 0）
    python3 python/config.py samples                     # id<TAB>fq1<TAB>fq2 を 1 サンプル 1 行

パスはプロジェクトルート相対で解決する（実行位置非依存 = spec「ルート実行」要件）。
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

# このファイル（python/config.py）の親の親 = プロジェクトルート
ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH = ROOT / "config.yaml"


def load_config(path: Path | str | None = None) -> dict[str, Any]:
    """config.yaml を読んで dict を返す。"""
    import yaml  # 遅延 import（PyYAML は pyproject の依存）

    p = Path(path) if path else CONFIG_PATH
    if not p.exists():
        raise FileNotFoundError(f"config.yaml が見つからない: {p}")
    with p.open("r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh) or {}
    return cfg


def _dig(cfg: dict[str, Any], dotted: str) -> Any:
    """'salmon.threads' のようなドット区切りキーで値を取り出す。"""
    cur: Any = cfg
    for part in dotted.split("."):
        if not isinstance(cur, dict) or part not in cur:
            raise KeyError(f"config.yaml にキーが無い: {dotted}")
        cur = cur[part]
    return cur


def _emit_scalar(val: Any) -> str:
    if isinstance(val, bool):
        return "true" if val else "false"
    if isinstance(val, (list, tuple)):
        return " ".join(str(v) for v in val)
    return str(val)


def _main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(__doc__, file=sys.stderr)
        return 2
    cmd = argv[1]
    cfg = load_config()
    if cmd == "get":
        if len(argv) < 3:
            print("usage: config.py get <dotted.key>", file=sys.stderr)
            return 2
        print(_emit_scalar(_dig(cfg, argv[2])))
        return 0
    if cmd == "samples":
        paired = bool(cfg.get("project", {}).get("paired_end", True))
        for s in cfg.get("samples", []):
            fq2 = s.get("fq2", "") if paired else ""
            print(f"{s['id']}\t{s.get('fq1', '')}\t{fq2}")
        return 0
    print(f"unknown command: {cmd}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(_main(sys.argv))
