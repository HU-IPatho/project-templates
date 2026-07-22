---
name: scaffold-rnaseq-secondary
description: >-
  bulk RNA-seq 2 次解析テンプレ（rnaseq-secondary）の着手時 scaffold。/grill-me を
  テンプレ固有の論点で走らせ、config.yaml の事実（design/contrast/method/samples/HK/
  organism）と AGENTS.md「このプロジェクト固有」節（非自明な narrative）を対話で埋める。
  入力が共有領域（/data/shared）由来なら README.yaml から config を下書きする。
  Use when: このテンプレで新規プロジェクトを始める / config を埋める / 解析設計を確定する。
---

# scaffold: rnaseq-secondary の層 3 を埋める

このテンプレ（analysis-template-standard 準拠）の**層 3＝プロジェクト固有**を確定する薄い
scaffold。層 1（`AGENTS-base.md`）・層 2（本テンプレ `AGENTS.md` 上部）は固定物。ここで埋めるのは
**事実（`config.yaml`）**と**非自明な固有規約（`AGENTS.md`「このプロジェクト固有」節）**のみ。

## 手順

1. **共有由来データの下書き（あれば）**: 入力が `/data/shared/datasets/<dataset>/` 由来なら、
   その `README.yaml`（provenance マニフェスト・`specs/shared-data-curation`）を読み、
   `organism` / `counts_file`（`matrix/` の位置）/ サンプル一覧（`metadata/samples.tsv`）を
   `config.yaml` の下書きに反映する（手入力を減らし来歴を継承する）。

2. **`/grill-me` を下記の固有論点で起動**し、ユーザーと対話して各分岐を 1 つずつ解消する。
   各質問の前に論点を平易に説明してから問う（grill-me の作法）。

   - **適用範囲 scope**: `organism`（human/mouse のみ）・`library_prep`（v0 は `polya_fulllength`=
     full-length polyA bulk のみ）。3'-tag/UMI/FFPE/他種は範囲外で fail-closed（要 recalibration・grill 論点）。
   - **provenance**: `provenance.source`（cell_line/patient/technical・BCV バンド整合に使う）。
   - **実験デザイン**: 群（`group_col` の水準）は何か。対照（`reference_level`）はどれか。
     **各群のサンプル数と複製の独立性**（`replicate_independence`: biological / technical / pseudo）——
     対比ごとの `min(実効群サイズ)` で lane が決まる（`min<2` → n=1 screening・`min>=2` → 複製あり）。
     **technical のみの群は biological n=1 扱い**（生物変動を計上しないため）。
   - **複製あり lane の検定**: `edger.test`（既定 `QL`＝glmQLFTest）。`LRT` にするなら理由必須（意図的逸脱）。
     `deseq2.shrink`（既定 `ashr`＝非参照間対比も可・単一係数対比のみ `apeglm`）。
   - **共変量**: バッチ・処理日・ドナー等の交絡はあるか（`covariates`）。無ければ `[]`。
     ※ **加法バッチ補正は複製あり lane 限定**（n=1 lane では素の ~0+group）。単一レベルは自動除外。
   - **対比（contrasts）**: どの群 vs どの群を見たいか（多対比を `名前: [num, den]` で列挙）。
   - **counts 入力**: `counts_file` の場所と形式（featureCounts / matrix TSV）。
     `python/fetch_counts.py` で `data/interim/counts.tsv` に正規化する導線を確認。
   - **n=1 screening ゲート（該当時）**: `provenance.replication_unavailable_reason`（複製が取れない理由・必須）。
     HK 候補（`screening.hk_gene_file`・要データ内検証 G3）。KD 標的の master regulator 該当（G2・
     `screening.master_regulator_file`）。cross-hairpin（G5・`screening.hairpin_map` に target→hairpins）。
     seed off-target（G4）・直交検証（G6）の要否。**出力は screening-grade で確定 FDR を主張しない**旨を共有。
   - **BCV バンド**: `bcv.band`（provenance.source 別・数値は例示。自分のデータで感度スイープを見て確定）。
   - **非自明な固有規約**: このデータ固有の落とし穴（正規化の例外・除外サンプル・特殊な
     アノテーション整合など）があれば `AGENTS.md`「このプロジェクト固有」節に narrative で残す。

3. **書き込み**: 合意した事実を `config.yaml` に、非自明な規約のみを `AGENTS.md`
   「このプロジェクト固有」節に反映する。共通・自明な内容は書かず上位層を指すに留める
   （重複＝drift 防止）。事実を `AGENTS.md` に、規約を `config.yaml` に混ぜないこと。

4. **確認**: `Rscript analysis/01_build_se.R` が通る最小状態（counts_file と samples の
   sample_id 整合）まで導く。
