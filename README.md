# HU-IPatho Project Templates

北海道大学 IPatho 研究室の解析プロジェクト用テンプレートリポジトリ。ipatho1 の Coder 環境（omics-dev workspace）上で `project-init` コマンドから使用される。

## テンプレート

| テンプレート | 用途 |
|---|---|
| `rnaseq-primary/` | バルク RNA-seq **1 次解析**（FASTQ → fastp/FastQC/MultiQC → Salmon → tximeta → SummarizedExperiment） |
| `rnaseq-secondary/` | バルク RNA-seq **2 次解析**（SE/counts → DESeq2 / edgeR で DEG・n=1 は housekeeping 遺伝子で分散推定） |
| `scrna/` | 単一細胞 RNA-seq 解析（Seurat ベース・Harmony 統合） |
| `general/` | 汎用解析（言語・パイプライン自由） |

`rnaseq-primary` / `rnaseq-secondary` / `scrna` は **解析テンプレ標準**（`HU-IPatho/coder` の `specs/analysis-template-standard/spec.md`・ADR-0023）に conform する — 3 層データ `data/{raw,interim,processed}`、`analysis/`+`R/`+`python/` の polyglot 配置、`outputs/{tables,figures}`（PNG+PDF）+`reports/`（Quarto）、出力ハーネス（`save_fig`/`save_table`/`captions.tsv`）、`renv.lock`/`uv.lock`/イメージタグの 3 系統版ピン。`general` は旧来の自由構成。

## 使い方

### omics-dev workspace の `project-init` コマンド（推奨）

```bash
project-init --list                # 利用可能なテンプレートを一覧
project-init <template> <name>     # 例: project-init rnaseq-secondary my-analysis
```

これにより以下が実行される:

1. 選択したテンプレートを `/work/projects/<name>` に展開
2. `*.template` / `AGENTS.md` の変数置換（`{{PROJECT_NAME}}`, `{{AUTHOR}}`, `{{DATE}}`）
3. `git init` + 初回コミット + `gh repo create HU-IPatho/<name> --private`（gh 認証時）

### GitHub Template Repository として

```bash
gh repo create HU-IPatho/<name> --template HU-IPatho/project-templates --private
```

## メンテナンス

- `HU-IPatho/coder` プロジェクトから分離された独立リポジトリ。テンプレの追加・修正は PR ベース。
- `project-init` の実体は `HU-IPatho/coder` の `templates/omics-dev/main.tf`（inline heredoc）。テンプレートは**動的に列挙**される（新テンプレはディレクトリを追加するだけで `--list` に載る・名前の hardcode は無い）。

## 関連

- [HU-IPatho/coder](https://github.com/HU-IPatho/coder) — Coder サーバー + workspace image 管理・解析テンプレ標準（`specs/analysis-template-standard`）
