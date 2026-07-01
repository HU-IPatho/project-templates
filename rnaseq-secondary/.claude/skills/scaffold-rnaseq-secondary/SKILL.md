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

   - **実験デザイン**: 群（`group_col` の水準）は何か。対照（`reference_level`）はどれか。
     **各群のサンプル数**（複製あり / 各群 n=1 か）——これで `method` が決まる:
     - 各群 n=1 → `method: edger`（HK 分散が必須。下の HK 論点へ）。
     - 複製あり → `deseq2` / `edger` / `both`（一致度を見るなら `both`）。
   - **共変量**: バッチ・処理日・ドナー等の交絡はあるか（`covariates`）。無ければ `[]`。
     ※ バッチはモデル共変量で補正する（埋め込み補正ではない）。単一レベルは自動除外。
   - **対比（contrasts）**: どの群 vs どの群を見たいか（多対比を `名前: [num, den]` で列挙）。
   - **counts 入力**: `counts_file` の場所と形式（featureCounts / matrix TSV）。
     `python/fetch_counts.py` で `data/interim/counts.tsv` に正規化する導線を確認。
   - **生物種と HK**: `organism`（human/mouse で同梱 HK を自動選択）。n=1 なら HK を
     自分の count 行列のアノテーションに合わせる必要（不一致は無視）。特殊なら
     `edger.hk_gene_file` に差替リストを指す。`bcv_fallback`（HK 不足時）の妥当性も確認。
   - **非自明な固有規約**: このデータ固有の落とし穴（正規化の例外・除外サンプル・特殊な
     アノテーション整合など）があれば `AGENTS.md`「このプロジェクト固有」節に narrative で残す。

3. **書き込み**: 合意した事実を `config.yaml` に、非自明な規約のみを `AGENTS.md`
   「このプロジェクト固有」節に反映する。共通・自明な内容は書かず上位層を指すに留める
   （重複＝drift 防止）。事実を `AGENTS.md` に、規約を `config.yaml` に混ぜないこと。

4. **確認**: `Rscript analysis/01_build_se.R` が通る最小状態（counts_file と samples の
   sample_id 整合）まで導く。
