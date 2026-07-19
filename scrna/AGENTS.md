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

パイプライン本体（`analysis/01,02,03,04`・`R/load_data.R`・`R/helpers.R`・`R/annotate.R`）は触らない
（例外: `R/annotate.R` の CNV/参照フックは有効化時に研究員が完成させる＝下記アノテ方法論）。
**`config.yaml`**（層 3 の事実）の値を変え、最終ラベルは **`annotations.yaml`**（研究員が手動確定）に書く。
初回は同梱 scaffold `/scaffold-config` が `/grill-me` を scRNA 固有論点で走らせて対話的に埋める
（`.claude/skills/scaffold-config/`）。

固有で**間違えやすい**設定（非自明ゆえここに明記する）:
- **`mito_pattern`**: ヒト=`^MT-` / マウス=`^mt-`。**種を間違えると percent.mt が全部 0** になる。
- **`batch_key`**: Harmony でバッチ効果を除く軸（通常 `sample`）。**癌種などの生物学的変数を batch key に
  しない**（本物の生物差まで消える）。
- **`loader`**: `tenx`（10x Cell Ranger 出力）/ `dense`（per-sample 密行列 gz・GEO processed 等）。
- **`qc`**: `nFeature`/`percent.mt` 閾値。まず violin plot で分布を見てから決める。
- **`sample_class`**: `mixed`（混在組織）/ `pure_cell_line`（純細胞株）。アノテ方法論の適用条件が分岐する（下記）。

## パイプラインの順序の勘所（崩すと結果が壊れる・scRNA 固有）

1. **`analysis/01_load_qc`**: 読込 → `merge` → `JoinLayers` → `percent.mt` → **doublet 取扱い（記録/除去）** → QC
   → `data/interim/seurat_qc.rds`
2. **`analysis/02_integrate`**: `Normalize`→**（cell-cycle スコア・直交軸）**→`FindVariableFeatures`→`Scale`→`RunPCA`→
   **`RunHarmony(group.by.vars=batch_key)`**→`FindNeighbors`→**複数解像度 `FindClusters`**→`UMAP`（`reduction="harmony"`）
   → `data/processed/seurat_integrated.rds`（**正準 object**）＋ UMAP 図
3. **`analysis/03_markers`**: `JoinLayers`→`FindAllMarkers`→top10 → `outputs/tables/`（marker 段）
4. **`analysis/04_annotate`**: 直交レイヤ注釈（系統/悪性/状態/cell-cycle）+ confidence/review routing
   → `outputs/tables/tbl03,tbl04` + `data/processed/seurat_annotated.rds`

- **Harmony は PCA の後・Neighbors/UMAP の前**。統合後は必ず `reduction="harmony"`（`"pca"` のままだと
  バッチ効果が残る）。**バッチが 1 つなら Harmony はスキップ**し PCA を使う（02 が自動判定）。
- **`FindAllMarkers` の前に `JoinLayers`**（v5 は sample ごとに layer が分かれ、結合しないと marker 検定が
  正しく走らない）。

## アノテーションの方法論（研究室標準 scrna-annotation-standard v0 に conform）

このプロジェクトのアノテーションは研究室標準（方法論のみ・numeric/ツール/marker panel は自分のデータで
`/grill-me` = grill ゲートで確定）に従う。**config.yaml のノブ**でパイプラインが標準の構造を踏む。要点:

- **sample-class 軸を先に決める**（`config.yaml: sample_class`）: `mixed`（腫瘍+正常/免疫/間質の混在組織・v0 の
  確定操作 scope）/ `pure_cell_line`（malignancy a priori 既知 → CNV ゲート省略・EMT を実 program 扱い・
  反復単位=culture replicate・composition 無効）。方法論の適用条件が分岐する。
- **cluster-then-marker（クラスタ先行→クラスタ単位手動マーカー解釈）**: **手動解釈が最終権威**。04 は自動で
  identity を確定せず、クラスタごとの証拠（top marker・hint 一致・cell-cycle・組成・（有効時）CNV/参照）を
  `tbl03_cluster_evidence` に組み、**confidence と review-priority** を付す。最終ラベルは **`annotations.yaml`**
  に研究員が手で埋める（evidence を見て確定）。**クラスタ＝細胞型ではない**——クラスタ内異質性があり、
  composition が変わるとラベルは不安定化しうる。曖昧/境界クラスタは subset 再クラスタか cell-level クロスチェックで
  再検証する。
- **marker p 値を絶対視しない**（同一データで clustering と marker を決める double-dipping で過大化）。effect size
  (`avg_log2FC`)・特異性(`pct.1/pct.2`)・文献照合・標本横断の再現性で総合判断する。
- **複数解像度をスキャン**（`config.yaml: resolutions`）: 単一固定にせず over/under-clustering を評価（02 が
  各解像度でクラスタし比較図を出す。primary=`resolution` を `seurat_clusters` に採用）。曖昧/境界クラスタは
  **反復 subset 再クラスタ**で再検証してよい（optional）。使う場合は**3 ガードレール必須**: (1) 再クラスタ解像度を
  統計的区別可能性で正当化（恣意的 resolution で細分を状態と決めつけない）、(2) double-dipping 回避（同一データの
  subtype marker は独立検証か null 補正を要する）、(3) subset 依存の score/integration 落とし穴の回避（batch 過補正
  回避・rank-based score・global と refined 間でスコアを直接比較しない）。低信頼ラベルは強引に細分せず粗ラベル/
  unknown へ後退する（reject option）。
- **直交レイヤで記述**: 系統/細胞型（marker/参照）・悪性性（CNV ゲート）・状態（program）・cell-cycle（直交軸）。
  **cell-cycle は identity と直交に「保持」**し（`config.yaml: cell_cycle.score`）、**無条件の regression 除去はしない**
  （`cell_cycle.regress` は identity 歪みが実証された時だけ TRUE）。
- **module/score は marker 従属の補強証拠**（identity を単独で上書きしない）。腫瘍層で高間葉スコアを EMT 癌細胞と
  即断しない（CAF/内皮で最高値になり腫瘍純度を追跡する）。score を足すなら **composition 非依存の rank-based
  手法（AUCell/UCell 等）を優先**し、**bulk 由来 ssGSEA/GSVA を single-cell identity の単独判定に使わない**。
- **悪性は marker 単独で呼ばない**（`config.yaml: cnv_gate`）: aneuploidy を持ちうる固形上皮 carcinoma では
  inferred CNV を直交検証ゲートに課す。**ゲートは非対称**（CNV 陽性→悪性は強いが、CNV 陰性→非悪性とは限らず
  「保留(hold)」）。**near-diploid 腫瘍種（前立腺癌・甲状腺癌・ccRCC・小児/造血器腫瘍・sarcoma 等）と純細胞株は
  非適用**（例外条項）。ツールは pin しない（inferCNV/CopyKAT/SCEVAN 等は researcher 選択・`R/annotate.R` の
  `cnv_gate_status()` を完成させる）。
- **参照ベースは two-track**（`config.yaml: reference`）: 適合アトラスが在るとき SingleR/Azimuth 等を**必須の
  complementary クロスチェック**として併走し、不一致クラスタを review へ routing する（無条件 primary/必須には
  しない）。`R/annotate.R` の `reference_crosscheck()` を完成させる。
- **confidence + review-priority routing**（04 が付す）: marker p 値を confidence に使わない（禁則）。示唆なし/曖昧・
  単一 sample 偏り・参照不一致を **review トリガー**（固定閾値でなくヒューリスティック・tune 可）にし
  `tbl04_annotation` の `review_flag/review_reason` に出す。低信頼ラベルは強引に細分せず粗ラベル/unknown へ後退。
- **群間 DE は sample 単位**（design-specific・注釈確定「後」の別層）: 細胞レベル検定でなく biological replicate 間
  分散を勘定（既定 pseudobulk、希少集団は mixed model の two-track）。反復単位は sample-class 別（組織=patient、
  細胞株=culture replicate/line）。本雛形は注釈段まで。群間 DE を足すときはこの規律に従う。

## ambient / doublet（filtered data の既知リスク）

- **ambient 補正（SoupX/CellBender）は raw/unfiltered droplet 行列が必須**。GEO の filtered(cell-called) 行列のみで
  生 UMI が無い場合は原理的に実行不能 → **既知リスクとして記録**し、残留 ambient（高発現 housekeeping/血液系の
  低レベル発現）を marker 解釈から割り引く（01 冒頭コメント参照）。
- **doublet の取扱いは必ず記録**（`config.yaml: doublet.method` と `rationale`）。除去可否・手法・閾値の科学判断は
  `/grill-me`（grill ゲート）。既定 baseline は除去なし（filtered 行列では検出力が落ちる）。

## ハマりどころ（`R/load_data.R` に封入済み・自前で読み込む時も同じ罠）

- **gz 直読は `data.table::fread(cmd="zcat ...")`**（`fread` の gz 直読は R.utils 依存を招く）。
- **GEO の dense 行列は ragged**（header 行=バーコードのみで 1 列少ない）。`fread(header=TRUE)` は列ズレで
  バーコードが重複する → **header 行を別途読んで割当て、本体は `header=FALSE, skip=1`**。
- **重複遺伝子名は `make.unique`。ただし Seurat が feature 名の `_`→`-` を置換するので、先に
  `gsub("_","-")` してから `make.unique`**（順序が逆だと置換で重複が再発する）。

## ディレクトリ規約（① analysis-template-standard 準拠）

| パス | 用途 |
|---|---|
| `config.yaml` | ★編集するのはここ（層 3 の事実・データ形式/QC/統合/アノテ方法論のノブ） |
| `annotations.yaml` | ★手動確定ラベル（04 の evidence を見てクラスタ→系統/悪性/状態を埋める・任意） |
| `analysis/00_config.R` | config.yaml を読む薄いローダ（触らない） |
| `analysis/01,02,03,04_*.R` | パイプライン本体（ルートから番号順に実行・`here::here()` でパス解決） |
| `R/load_data.R` / `R/helpers.R` | ローダ（gotcha 封入）/ 出力ハーネス（触らない） |
| `R/annotate.R` | アノテ方法論ロジック（直交レイヤ/confidence/CNV・参照フック。有効化時のみ研究員が完成） |
| `python/fetch_data.py` | データ取得（GEO/URL → `data/raw`） |
| `data/raw/` | 生データ（**不変**入力層・git 非追跡）。`python/fetch_data.py` か `work-fetch` で置く |
| `data/interim/` | 再生成可能な中間物（`seurat_qc.rds` 等・git 非追跡） |
| `data/processed/` | **正準 object**（`seurat_integrated.rds` / `seurat_annotated.rds`・共有昇格の単位・git 非追跡） |
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
