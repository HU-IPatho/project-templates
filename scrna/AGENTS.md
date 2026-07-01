# scRNA-seq 解析プロジェクト規約（Seurat v5 + Harmony）

このテンプレートは single-cell RNA-seq を **再現性高く** 解析するための雛形。
`HU-IPatho/gse134520-scrna-dogfood`（胃がん scRNA を実走した参照実装）と同じ設計で、
そこで判明した「ハマりどころ」を封じ込めてある。**まず本ファイルを最後まで読むこと。**

## クイックスタート（この順で実行する）

```bash
# 1) 解析環境を用意（初回だけ・renv.lock から版を固定して復元）
Rscript install_deps.R                       # → INSTALL_DONE / Seurat 5.x / harmony 2.x

# 2) データを data/raw/ に置く（fetch_data.sh を自分のデータ用に編集して実行）
bash fetch_data.sh

# 3) 自分のデータに合わせて 00_config.R を編集（★編集するのは基本ここだけ★）

# 4) パイプライン実行（02 は長時間 → 必ず job-run。下記「長時間ジョブ」参照）
Rscript 01_load_qc.R                          # → QC_DONE
JID=$(job-run --label integrate -- Rscript 02_integrate.R)
job-wait "$JID" --timeout 570                 # 完了(0)/失敗(1)/未完(124→もう一度)
Rscript 03_markers.R                          # → MARKERS_DONE

# 5) 結果を HDD に退避（/work は揮発的・セッション終了前に必須）
work-save results/ && work-save reports/
```

## 自分のデータに合わせる — 編集するのは `00_config.R` だけ

スクリプト本体（01/02/03・R/load_data.R）は触らない。`00_config.R` の値を変えるだけで動く。

- **`loader`**: データ形式で選ぶ。
  - `"tenx"` … 10x Cell Ranger 出力（各 sample ディレクトリに `matrix.mtx.gz` / `features.tsv.gz` / `barcodes.tsv.gz`）。**最も一般的。**
  - `"dense"` … sample ごとの密行列 gz（GEO の processed matrix 等）。参照実装 GSE134520 はこれ。
- **`batch_key`**: 統合（Harmony）でバッチ効果を除く軸。通常は `"sample"`。**癌種などの生物学的変数を batch key にしない**（本物の差まで消える）。
- **`mito_pattern`**: ヒト=`"^MT-"` / マウス=`"^mt-"`。**種を間違えると percent.mt が全部 0 になる。**
- **`qc`**: `nFeature`/`percent.mt` 閾値。まず分布（violin plot）を見てから決めるのが定石。テンプレ既定は汎用値。

## 長時間ジョブは job-run で回す（Harmony 統合など・必須）

数分を超える処理（`02_integrate.R` の Harmony・大きな `FindAllMarkers`）は **前景で直接実行しない**。
前景実行は (a) ツールがタイムアウトし (b) 通信を切るとジョブごと死ぬ。必ず `job-run` で detach 実行し `job-wait` で完了を待つ。

```bash
JID=$(job-run --label integrate -- Rscript 02_integrate.R)   # 即 JOB_ID を返す（tmux 内で継続）
job-wait "$JID" --timeout 570   # ★ツール側 timeout は最大(=600000ms)で呼ぶ。--timeout は必ずそれ未満(570)
#   exit 124 = 未完 → 間を空けず同じ job-wait をもう一度（イベント駆動で完了の瞬間に返る）
#   exit 0   = 成功 / exit 1 = 失敗
job-logs "$JID" -n 50                # 結果の末尾を確認（失敗時は -n 200）
job-status "$JID" --env              # exit_code / runtime
job-list                             # 全ジョブ俯瞰
```

- **Claude Code**: `job-wait` ループの代わりに Monitor tool に登録すると呼び出しを圧縮できる（任意）。
- **Codex**: 長時間コマンドで固まることがあるため `--timeout` は短め(60)にするか `job-status` を数分間隔で単発確認する。
- 同じリソースを使うジョブは `job-run --exclusive <名前> -- …` で直列化する。

## パイプラインの中身と順序（ベストプラクティス）

1. **01_load_qc**: 読込 → `merge` → `JoinLayers` → `percent.mt` 算出 → QC フィルタ → `data/processed/seurat_qc.rds`
2. **02_integrate**: `NormalizeData` → `FindVariableFeatures` → `ScaleData` → `RunPCA` → **`RunHarmony(group.by.vars=batch_key)`** → `FindNeighbors/FindClusters/RunUMAP`（いずれも `reduction="harmony"`）→ `seurat_integrated.rds` + UMAP
3. **03_markers**: `JoinLayers` → `FindAllMarkers` → top10 → `results/markers_*.csv`

順序の勘所（崩すと結果が壊れる）:
- **Harmony は PCA の後・Neighbors/UMAP の前**。統合後は必ず `reduction="harmony"` を使う（`"pca"` のままだとバッチ効果が残る）。
- **バッチが 1 つなら Harmony は不要**（02 は自動でスキップし PCA を使う）。
- **`FindAllMarkers` の前に `JoinLayers`**（v5 は sample ごとに layer が分かれており、結合しないと marker 検定が正しく走らない）。

## ハマりどころ（`R/load_data.R` に封入済み。自前で読み込む時も同じ罠）

- **gz 直読は `data.table::fread(cmd="zcat ...")`**（`fread` の gz 直読は R.utils 依存を招く）。
- **GEO の dense 行列は ragged**（header 行=バーコードのみで 1 列少ない）。`fread(header=TRUE)` は列ズレでバーコードが重複する → **header 行を別途読んで割当て、本体は `header=FALSE, skip=1`**。
- **重複遺伝子名は `make.unique`。ただし Seurat が feature 名の `_`→`-` を置換するので、先に `gsub("_","-")` してから `make.unique`**（順序が逆だと置換で重複が再発する）。

## 再現性（renv で版を固定）

- `renv.lock` が dogfood で検証済みの版（Seurat 5.5.1 / harmony 2.0.5 等）を固定する。`install_deps.R` がこれを `renv::restore()` で復元する。
- **`renv::update()` を安易に実行しない**（版が動くと結果が変わりうる）。更新したら必ず `renv::snapshot()` して `renv.lock` を commit する。
- パッケージ取得は PPM バイナリ（`packagemanager.posit.co/.../noble/latest`）で高速。

## ディレクトリ規約

| パス | 用途 |
|---|---|
| `data/raw/` | 生データ（git 管理外）。`fetch_data.sh` でここへ置く |
| `data/processed/` | 中間 Seurat オブジェクト（`*.rds`・git 管理外） |
| `00_config.R` / `R/load_data.R` | 設定 と ローダ（編集するのは基本 config だけ）|
| `01/02/03_*.R` | パイプライン本体（番号順に実行）|
| `results/` | マーカーテーブル等（git 追跡）|
| `reports/` | UMAP 等の図（git 追跡）|

## データ操作（ワークスペース共通ヘルパー）

```bash
work-fetch /data/shared/datasets/<dataset>/   # HDD → SSD（解析前）
work-save  results/ reports/                   # SSD → HDD（解析後・必須）
work-share results/markers_top10.csv -m "説明" # 共有 _inbox へ投稿
work-status                                     # SSD 使用量
```

## エージェント共通の作法

- 推測せず調べてから答える。2 回失敗したら場当たり修正せず原因を分解する。
- 生データ（`/data/shared`）は read-only。書けるのは `/work` と `_inbox` のみ。
- 変更後は `git status && git diff` を確認し、Conventional Commits（`feat:`/`fix:`/`docs:`）でコミットする。
