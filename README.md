# HU-IPatho Project Templates

北海道大学 IPatho 研究室の解析プロジェクト用テンプレートリポジトリ。

## 概要

`rnaseq`, `scrna`, `general` の3種類のプロジェクトテンプレートを提供する。ipatho1 の Coder 環境（omics-dev workspace）上で `project-init` コマンドから使用される。

## テンプレート

| テンプレート | 用途 |
|---|---|
| `rnaseq/` | バルク RNA-seq 解析（DESeq2/edgeR ベース） |
| `scrna/` | 単一細胞 RNA-seq 解析（Seurat/Scanpy ベース） |
| `general/` | 汎用解析（言語・パイプライン自由） |

## 使い方

### GitHub Template Repository として

本リポジトリは [Template Repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-template-repository) として設定されており、以下のコマンドで新規プロジェクトを作成できる。

```bash
gh repo create HU-IPatho/<project-name> --template HU-IPatho/project-templates --private
```

ただし推奨は後述の `project-init` 経由。

### omics-dev workspace の `project-init` コマンド

ipatho1 の omics-dev workspace で利用可能。

```bash
project-init <project-name> --template rnaseq
```

これにより以下が実行される:
1. 選択したテンプレート（`rnaseq`/`scrna`/`general`）を展開
2. `AGENTS.md` / `README.md.template` の変数置換（`{{PROJECT_NAME}}`, `{{DATE}}`, `{{USER}}`）
3. `git init` + 初回コミット
4. `gh repo create HU-IPatho/<project-name> --private`

## 各テンプレートの構成

全テンプレート共通:
- `AGENTS.md` — AI コーディングエージェント（Claude Code / GitHub Copilot CLI）向け指示書
- `CLAUDE.md` — `AGENTS.md` への後方互換 symlink
- `README.md.template` — プロジェクト README 雛形
- `.gitignore` — 大容量データ除外（`data/raw/` 等）
- `src/`, `results/`, `reports/` ディレクトリ（`.gitkeep` 付き）

テンプレート固有のディレクトリ構造は各 `AGENTS.md` を参照。

## メンテナンス

- このリポジトリは `HU-IPatho/coder` プロジェクトから分離された独立リポジトリ
- テンプレートの追加・修正は PR ベース
- `project-init` コマンド側の参照は `HU-IPatho/coder` の `scripts/project-init` を参照

## 関連

- [HU-IPatho/coder](https://github.com/HU-IPatho/coder) — Coder サーバー + workspace image 管理
