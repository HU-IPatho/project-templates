# scRNA-seq 解析プロジェクト規約（Seurat v5 + Harmony）

single-cell RNA-seq を**再現性高く**解析する雛形。参照実装＝`HU-IPatho/gse134520-scrna-dogfood`
（胃がん scRNA を実走）で判明した「ハマりどころ」を封入してある。

この AGENTS.md は**指示 3 層**（① analysis-template-standard）の**層 3＝このプロジェクト固有**にあたり、
**scRNA に固有で非自明な規約だけ**を書く。workspace 共通（`work-*` / `job-run` / 共有データ / git 作法）は
**層 1＝`AGENTS-base.md`**（omics-dev workspace が与える）に委ね、ここでは重複させない。
機械可読な「事実」は**`config.yaml`**（層 3 の事実）にある。`CLAUDE.md` は本ファイルへの symlink。

> **上位層へのポインタ**: `work-fetch/work-save/work-share`・長時間ジョブの `job-run/job-wait`・
> 共有データ（`/data/shared`・`registry.tsv`・`README.yaml`）の探し方・git コミット規約・R=renv / Python=uv の
> ツール規約は、すべて **`AGENTS-base.md`** を参照する（本ファイルには再掲しない）。

## 編集するのは `config.yaml` だけ

パイプライン本体（`analysis/01,02,03`・`R/load_data.R`・`R/helpers.R`）は触らない。
**`config.yaml`**（層 3 の事実）の値を変えるだけで動く。初回は同梱 scaffold `/scaffold-config` が
`/grill-me` を scRNA 固有論点で走らせて対話的に埋める（`.claude/skills/scaffold-config/`）。

固有で**間違えやすい**設定（非自明ゆえここに明記する）:
- **`mito_pattern`**: ヒト=`^MT-` / マウス=`^mt-`。**種を間違えると percent.mt が全部 0** になる。
- **`batch_key`**: Harmony でバッチ効果を除く軸（通常 `sample`）。**癌種などの生物学的変数を batch key に
  しない**（本物の生物差まで消える）。
- **`loader`**: `tenx`（10x Cell Ranger 出力）/ `dense`（per-sample 密行列 gz・GEO processed 等）。
- **`qc`**: `nFeature`/`percent.mt` 閾値。まず violin plot で分布を見てから決める。

## パイプラインの順序の勘所（崩すと結果が壊れる・scRNA 固有）

1. **`analysis/01_load_qc`**: 読込 → `merge` → `JoinLayers` → `percent.mt` → QC → `data/interim/seurat_qc.rds`
2. **`analysis/02_integrate`**: `Normalize`→`FindVariableFeatures`→`Scale`→`RunPCA`→
   **`RunHarmony(group.by.vars=batch_key)`**→`FindNeighbors/Clusters/UMAP`（`reduction="harmony"`）
   → `data/processed/seurat_integrated.rds`（**正準 object**）＋ UMAP 図
3. **`analysis/03_markers`**: `JoinLayers`→`FindAllMarkers`→top10 → `outputs/tables/`

- **Harmony は PCA の後・Neighbors/UMAP の前**。統合後は必ず `reduction="harmony"`（`"pca"` のままだと
  バッチ効果が残る）。**バッチが 1 つなら Harmony はスキップ**し PCA を使う（02 が自動判定）。
- **`FindAllMarkers` の前に `JoinLayers`**（v5 は sample ごとに layer が分かれ、結合しないと marker 検定が
  正しく走らない）。

## ハマりどころ（`R/load_data.R` に封入済み・自前で読み込む時も同じ罠）

- **gz 直読は `data.table::fread(cmd="zcat ...")`**（`fread` の gz 直読は R.utils 依存を招く）。
- **GEO の dense 行列は ragged**（header 行=バーコードのみで 1 列少ない）。`fread(header=TRUE)` は列ズレで
  バーコードが重複する → **header 行を別途読んで割当て、本体は `header=FALSE, skip=1`**。
- **重複遺伝子名は `make.unique`。ただし Seurat が feature 名の `_`→`-` を置換するので、先に
  `gsub("_","-")` してから `make.unique`**（順序が逆だと置換で重複が再発する）。

## ディレクトリ規約（① analysis-template-standard 準拠）

| パス | 用途 |
|---|---|
| `config.yaml` | ★編集するのはここ（層 3 の事実・データ形式/QC/統合パラメータ） |
| `analysis/00_config.R` | config.yaml を読む薄いローダ（触らない） |
| `analysis/01,02,03_*.R` | パイプライン本体（ルートから番号順に実行・`here::here()` でパス解決） |
| `R/load_data.R` / `R/helpers.R` | ローダ（gotcha 封入）/ 出力ハーネス（触らない） |
| `python/fetch_data.py` | データ取得（GEO/URL → `data/raw`） |
| `data/raw/` | 生データ（**不変**入力層・git 非追跡）。`python/fetch_data.py` か `work-fetch` で置く |
| `data/interim/` | 再生成可能な中間物（`seurat_qc.rds` 等・git 非追跡） |
| `data/processed/` | **正準 object**（`seurat_integrated.rds`・共有昇格の単位・git 非追跡） |
| `outputs/tables/` `outputs/figures/` | 表 / 図（figNN・PNG+PDF・git 追跡） |
| `reports/` | Quarto レポート（`report.qmd`・保存済み object を読むだけ・git 追跡） |

- **出力は必ず出力ハーネス経由**: `save_fig`（PNG+PDF・固定 DPI・統一テーマ）/ `save_table` を使い、
  生の `ggsave`/`write.csv` を直呼びしない（`outputs/captions.tsv` に図表→由来スクリプトが記録される）。
- **`results/` は使わない**（① で廃止・`outputs/` に一本化）。
- **共有昇格**は `data/processed` の正準 object を `work-share` で `_inbox` へ（詳細は AGENTS-base）。

## 再現性（renv で版を固定）

- `renv.lock` が dogfood 検証済みの版（Seurat 5.5.1 / harmony 2.0.5 等）を固定する。`install_deps.R` が
  `renv::restore()` で復元する。**`renv::update()` を安易に実行しない**（版が動くと結果が変わりうる）。
- Python は `pyproject.toml`（uv）。既定は依存ゼロ（取得系は標準ライブラリのみ）。
