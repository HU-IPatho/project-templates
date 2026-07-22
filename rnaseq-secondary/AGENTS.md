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
- **`organism` / `library_prep`**: scope（human/mouse ＋ full-length polyA bulk）。範囲外は fail-closed。
- **`provenance`**: `source`（cell_line/patient/technical・BCV バンド整合）・`replication_unavailable_reason`
  （n=1 経路を使うなら必須）。
- **`method`**: 複製あり lane の engine。`edger`（glmQLFTest）/ `deseq2` / `both`（concordance sanity check）。
  ※ n=1（min(group)<2）対比は engine に依らず必ず edgeR screening 経路。
- **`design`**: `group_col`（主比較変数）・`covariates`（バッチ等・**加法補正は複製あり lane 限定**・
  単一レベルは自動除外）・`reference_level`（対照）。
- **`replicate_independence`**: 複製の独立性（`biological`/`technical`/`pseudo`）。technical のみは n=1 扱い。
- **`contrasts`**: `名前: [numerator, denominator]`（group のレベル名で多対比を宣言）。
- **`edger.test`**: `QL`（既定）/ `LRT`（逸脱時 `test_deviation_reason` 必須）。**`deseq2.shrink`**:
  `ashr`（既定・contrast 対応）/ `apeglm`（単一係数のみ）/ `none`。
- **`screening`**: n=1 経路のゲート（G1–G6）・HK 検証・hairpin_map（G5）・master_regulator_file（G2）等（数値は例）。
- **`bcv.band`**: provenance.source 別 BCV バンド（感度スイープ用・数値は例示）。
- **`samples`**: サンプルシート（行=sample・`sample_id` は counts のヘッダと一致）。

## 差次発現（DEG）は研究室標準 bulk-secondary-deg-standard に conform

DEG 段の方法論は研究室標準 **bulk-secondary-deg-standard**（v0・方法論）に conform する
（`HU-IPatho/coder` の `specs/bulk-secondary-deg-standard/spec.md`・ADR-0034）。標準の
domain/dataset 固有の **numeric**（HK membership・BCV バンド値・各種閾値・G2 リスト・対比の中身）は
**TBD**＝config.yaml の「default 例」であり、自分のデータで診断を見てから確定する（確定は `/grill-me`
＝grill ゲート）。**科学的標準そのものを勝手に書き換えない**（変更は grill 批准が要る・ADR-0028）。

### 複製構造による routing（対比単位 `min(group)`）

経路分岐は **対比ごとの最小群サイズ `min(table(group))`** で決める（`R/routing.R`）。`min<2` の対比は
**n=1 スクリーニング経路**へ、`min>=2` は**複製あり経路**へ振られる（部分複製デザインの silent 漏れ防止）。
複製の**独立性**を `config$replicate_independence`（`biological`/`technical`/`pseudo`）で必ず申告する。
**technical のみの群は biological n=1 扱い**（technical 分散は生物変動を計上せず偽陽性が膨張する）。

### 複製あり経路（較正済み FDR lane）

- **edgeR**: `filterByExpr` → TMM → `estimateDisp(design)` → **`glmQLFit`+`glmQLFTest` が既定**
  （`edger.test: QL`・very reliable FDR control）。`glmLRT`（旧経路）は複製ありでは非推奨で、使うなら
  `edger.test: LRT` + `test_deviation_reason` に理由必須（意図的逸脱）。
- **DESeq2**: median-of-ratios・NC を reference relevel。低発現 prefilter は**任意**（速度/メモリ）で、
  FDR 統制は `results()` の **independent filtering** が担う（別概念）。LFC 収縮は `deseq2.shrink`:
  **非参照間対比（sh1_vs_sh2 等）は `ashr`（contrast 対応）**、参照との単一係数対比のみ `apeglm` 可。
  収縮できず生 LFC に落ちたら出力に `shrink=none` と**明示ラベル**（silent 化しない）。
- **バッチは複製あり lane 限定で加法補正**（`~batch+group`・残差自由度＝複製を要する）。n=1 lane では
  batch を足さず素の `~0+group`。**バッチが group と交絡すると（各バッチに一部の群しか無い等）加法補正でも
  補正不能**（`model.matrix` が rank-deficient）。交絡はデザイン段で各バッチに全群を配置してバランス化する
  のが上流要件（補正でなく設計で防ぐ）。
- **複製数の指針**（SHALL・Love 2014 / Schurch 2016）: 最低 **3**（実務下限）・望ましく **6 以上**・
  DE 全体網羅なら **12 以上**。少数複製でも較正・再現性は劣化する（n=2 を妥当と過信しない）。

### ★n=1 KD スクリーニングスタンダード（screening-grade）

複製取得不能時の n=1 KD 経路は **screening standard**（documented degraded path）。出力は
**screening-grade（候補絞り・仮説生成）**で、**較正済み FDR を確定的主張として出さない**。ヒットは
fold-change / 収縮 LFC のランキングに限り、**記述解析（FC ランキング + MDS・有意性主張なし）を必ず併走**
（edgeR §2.13 option 1）。`config$provenance$replication_unavailable_reason` に**複製が取れない理由を必ず記録**。

分散は **検証済 HK 由来 `common.dispersion`（BCV バンドより優先）**、不足時は **BCV バンド代表値**へ
フォールバック（`R/de_edger.R` / `R/screening.R`）。**anti-conservative 方向性バイアス**（HK の低分散を
高分散遺伝子へ移植 → 偽陽性が生物学的に興味深い変動遺伝子に集中）に注意（`disp_source` と診断表に surface）。

**この経路の妥当性は biology-conditional**（「当該摂動が大域 RNA 組成シフトを起こさない」かつ「HK 群が
その摂動下で非DE」の precondition 下でのみ成立）。無条件に「統計機構ゆえ組織/遺伝子非依存」と述べない。

#### n=1 経路の必須ゲート G1–G6（`R/screening.R`・数値は config/TBD）

- **G1 大域シフト（hard）**: HK 群の群間 logFC 分布を診断。中央値/分散が前提破綻を示せば DEG 破棄 or 記述降格。
- **G2 既知グローバル制御因子（適用禁止＋リスト外 warning）**: 標的が master regulator 短リスト
  （`resources/master_regulators/`・既定は空）に載れば HK 経路 適用禁止（spike-in/直交検証要求）。
- **G3 データ内 HK 非応答検証（hard）**: 汎用 HK の**無検証使用禁止**。当該データから経験的 control
  （低 |群間 logFC|・高発現・低 CV）を導出し、供給 HK を検証。既知 CNV アーム上遺伝子は除外。
- **G4 shRNA seed off-target スクリーン（既定 ON 推奨）**: `screening.seed_offtarget.enabled`。未実装で
  enabled なら fail-closed 停止。
- **G5 cross-hairpin concordance（第一信頼フィルタ）**: `screening.hairpin_map`（target→hairpins）を宣言
  すると「sh 間で符号一致かつ両者候補閾値超」を第一信頼に。1 標的 2 本以上推奨・単一は低信頼ダウンランク。
- **G6 直交検証（既定 ON 推奨）**: `screening.orthogonal_validation.required`。確定扱い前に qPCR/独立ハーピン等。

### BCV 運用（バンド＋感度スイープ）

固定単一 BCV を確定値運用しない（最感度パラメータ）。`config$bcv$band`（provenance.source 別・数値は例示）で
**感度スイープ**し、ヒットの順位安定性を出力（`outputs/tables/bcv_sensitivity_*.tsv`）。`provenance.source` に
対応バンドが無ければ fail-closed。HK 由来推定が成功すれば固定 BCV より優先。

### cross-engine concordance（補助 sanity check）

`method: both`（複製あり）では `compare_methods` が **Spearman 順位相関**主指標・Jaccard・符号一致・
非対称 overlap を **effect-size bin と方向で層別**して出す（`outputs/tables/concordance_*.tsv`）。これは
**頑健性の保証でなく補助 sanity check**（両エンジンは failure mode が相関）。**n=1 には適用しない**
（DESeq2 は複製なしで null）。n=1 の robustness proxy は cross-hairpin（G5）と BCV スイープの順位安定性。

### 適用範囲（scope・fail-closed）

v0 scope は **human/mouse ＋ full-length polyA bulk**。`config$organism`・`config$library_prep` を検証し
範囲外は **fail-closed で停止**（`R/scope.R`）。HK リスト不在の種・3'-tag/UMI/FFPE 等は silent に適用しない
（「要 recalibration」）。

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
