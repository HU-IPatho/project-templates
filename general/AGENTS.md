# データ解析プロジェクト

## 概要

このプロジェクトはデータ解析を行う汎用テンプレートから作成された。

## ディレクトリ規約

- `data/raw/` — 生データ。git 管理外。`work-fetch` で SSD にコピーして使用
- `data/reference/` — リファレンスデータ。git 管理外
- `data/processed/` — 前処理済みデータ。git 管理外
- `src/` — 解析スクリプト。番号付き（01_, 02_, ...）で実行順序を明示
- `results/` — 解析結果（テーブル、統計量）
- `reports/` — レポート・図表

## データ操作コマンド

```bash
work-fetch /data/shared/datasets/<data>/   # HDD → SSD コピー
work-save results/                         # SSD → HDD 保存
work-share results/output.csv --message "説明" # 共有リクエスト
work-status                                # SSD 使用量確認
```

## 言語・ツール

プロジェクトに応じて R または Python を使用。
- R: renv で依存管理
- Python: uv で依存管理
