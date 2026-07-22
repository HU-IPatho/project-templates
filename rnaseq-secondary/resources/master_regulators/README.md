# master regulator 短リスト（G2・既知グローバル制御因子）

n=1 KD screening lane のゲート **G2**（`specs/bulk-secondary-deg-standard/spec.md`）で使う、
既知の **master regulator（グローバル制御因子）** の curated 短リストの置き場。

## なぜ必要か

KD 標的が主要 TF・クロマチン制御・スプライシング・翻訳・グローバル増幅因子（例: MYC）のような
**大域的 RNA 組成シフトを起こす制御因子**の場合、HK-dispersion 経路の前提（「大多数の遺伝子は非DE」）が
破綻する。このとき HK-dispersion 経路は**適用禁止**とし、スパイクイン（ERCC 等）正規化 and/or 直交検証を
要求する。網羅は不可能ゆえ「短リスト適用禁止 ＋ リスト外 warning」の二段構え。

## 使い方

1. curated 短リスト（1 行 1 遺伝子シンボル・`#` 始まりと空行は無視）を作る。例: `master_regulators_human.txt`。
2. `config.yaml` の `screening.master_regulator_file` にそのルート相対パスを指す。
3. KD 標的（`config.yaml` の `screening.hairpin_map` のキー）がリストに載れば `02_de.R` が
   HK 経路を **適用禁止（停止）** し、リスト外なら warning に留める。

## 例（そのまま採用しない・データ/文脈依存・grill ゲートで確定）

```
# 例: グローバル制御因子（採用前に自分の実験文脈で curate すること）
MYC
MYCN
EZH2
BRD4
```

**既定は未整備（`master_regulator_file: ""`）** ＝ 常に「リスト未整備」warning。リスト整備後に機械強制化する
（リスト membership の確定は grill ゲート委譲＝科学判断）。
