# bulk RNA-seq 2 次解析プロジェクト規約（SE → DEG・DESeq2/edgeR）

このテンプレートは、既に得られた **counts 行列**から下流（正準 object 化 → DEG → 図表 →
レポート）を **再現性高く**回すための雛形。1 次解析（FASTQ→counts）は別テンプレ（③ bulk 1 次）。
**まず本ファイルを最後まで読むこと。** 環境・共通作法は上位層（`AGENTS-base.md`）を参照し、
ここには本テンプレ固有の非自明な規約のみを置く（重複記載しない＝drift 防止）。

## この解析の骨格（analysis-template-standard 準拠）

- **正準 object = SummarizedExperiment**（counts + colData を 1 object）。`analysis/01_build_se.R`
  が `data/processed/se.rds` に保存する。これが共有領域への昇格単位（`work-share`）。
- **パイプラインは番号順**: `01_build_se` → `02_de` → `03_figures`。各スクリプトはプロジェクト
  ルートから実行し、パスは `here::here()` でルート相対解決（実行位置非依存）。
- **図表は必ず出力ハーネス経由**（`R/helpers.R`）: `save_fig`（1 呼出しで PNG+PDF・固定 DPI・
  統一テーマ）/ `save_table`（`outputs/tables` へ TSV）/ `captions.tsv`（図表→説明→由来script）。
  生の `ggsave`/`write.csv` を直呼びしない（両形式・体裁・来歴索引が抜ける）。
- **レポートは再計算しない**: `reports/report.qmd` は `se.rds`/`de.rds` と `outputs/` を読むだけ。
  重い計算は `analysis/` が担い object に保存する（数値/図の乖離防止）。

## 編集するのは `config.yaml` だけ（事実層）

スクリプト本体（`analysis/`・`R/`）は触らない。`config.yaml` の値だけ変えて動かす。着手時に
**`/scaffold-rnaseq-secondary`** を実行すると `/grill-me` が固有論点を対話で引き出し、
`config.yaml` と本ファイル「このプロジェクト固有」節を埋める。主なキー:

- **`counts_file`**: gene×sample の raw counts TSV（1 列目=gene_id・ヘッダ=sample_id）。
- **`method`**: `edger`（n=1 はこれ）/ `deseq2`（複製ありの標準）/ `both`（一致度も出力）。
- **`design`**: `group_col`（主比較変数）・`covariates`（バッチ等の共変量・単一レベルは自動除外）・
  `reference_level`（対照）。**バッチはモデル共変量で補正**する（scRNA の Harmony=埋め込み補正とは別レイヤ）。
- **`contrasts`**: `名前: [numerator, denominator]`（group のレベル名で多対比を宣言）。
- **`samples`**: サンプルシート（行=sample・`sample_id` は counts のヘッダと一致）。

## ★n=1 DEG（各群 1 サンプル）— 本テンプレの主眼

複製が無いと分散を各遺伝子から推定できない。そこで **housekeeping 遺伝子**（群間で発現不変が
前提）の群間ばらつきを分散の代理に使い、edgeR の `common.dispersion` を推定する（`R/de_edger.R`）:

1. group を全 1 に潰した DGEList の HK 部分集合で
   `estimateDisp(HK, trend.method="none", tagwise=FALSE)` → `common.dispersion`
2. それを全遺伝子の DGEList へ移植（trended/tagwise は付けない）
3. `glmFit`/`glmLRT` で多対比（NC vs sh1, NC vs sh2 …）

- **HK リストは同梱**（`resources/housekeeping/hk_human.txt` / `hk_mouse.txt`・Eisenberg &
  Levanon 2013 + 古典 HK の curated 部分集合）。`config$organism` で自動選択、
  `config$edger.hk_gene_file` で任意リスト（全 Eisenberg 等）へ差替可。
- HK が発現行列に `min_hk_genes` 未満しか無い / 推定失敗 → **BCV 固定へフォールバック**
  （`bcv_fallback`・human 0.4 / 細胞株 0.1 / technical 0.01）。
- **n=1 の結果は探索的**。FDR/有意性の解釈は限定的（`de.rds$edgeR$disp_source` に由来を記録）。
  スクリプトは実行時に warning でこれを明示する。**複製がある実験を n=1 と偽って回さないこと。**

## 前処理・正規化（エンジン差）

- **edgeR**: `filterByExpr` → TMM 正規化（`calcNormFactors(method="TMM")`）。
- **DESeq2**: 低発現 prefilter（`prefilter.min_count_sum`）→ median-of-ratios（`DESeq` 内部）。
- **両走時**（`method: both`）は `compare_methods` が対比ごとに有意遺伝子の Jaccard・logFC 相関・
  符号一致率を `outputs/tables/method_agreement_*.tsv` に出す。DESeq2 は複製が要るので、
  n=1 では `edger` を使う（両走は複製ありのとき）。

## 再現性（版ピン 3 系統）

- **R+Bioconductor**: `renv.lock`（Bioconductor **3.20** 明示・DESeq2/edgeR/SummarizedExperiment/
  apeglm）。`install_deps.R` が `renv::restore()`、不完全なら fresh install + snapshot で完成させる。
- **Python**: `pyproject.toml`（`uv lock && uv sync` で `uv.lock` を生成・追跡）。取得変換のみの最小構成。
- **システムツール**: omics-dev イメージのタグ（`AGENTS-base.md` / ADR-0020）。
- `renv::update()` を安易に実行しない（版が動くと結果が変わる）。更新したら `renv::snapshot()`。

## データ操作（ワークスペース共通ヘルパー・詳細は AGENTS-base.md）

```bash
work-fetch /data/shared/datasets/<dataset>/   # HDD → data/raw（不変入力層）
work-save  data/processed outputs reports      # SSD → HDD（永続・必須）
work-share data/processed/se.rds -m "説明"     # 正準 object を共有 _inbox へ昇格
```

- `data/` は git 非追跡（`raw`=不変入力 / `interim`=再生成可能 / `processed`=正準 object）。
- `outputs/`・`reports/` は git 追跡（PR/diff でレビュー）。生データ（`/data/shared`）は read-only。

---

## このプロジェクト固有

<!-- 層 3（narrative）: 着手時に /scaffold-rnaseq-secondary（/grill-me 駆動）が
     非自明な固有規約のみを追記する。事実（design/contrast/samples 等）は config.yaml へ。
     共通・自明な内容はここに書かず上位層（AGENTS-base.md・本ファイル上部）を指す。 -->

_（未設定。`/scaffold-rnaseq-secondary` を実行して埋める。）_
