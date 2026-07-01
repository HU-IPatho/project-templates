# ============================================================================
# scRNA-seq パイプライン設定 — ★自分のデータに合わせて編集するのはこのファイルだけ★
# 01/02/03 のスクリプトはここの CONFIG を読んで動く（スクリプト本体は触らなくてよい）。
# ============================================================================
CONFIG <- list(

  # --- データ形式（ローダの選択）---
  # "tenx"  : 10x Cell Ranger 出力（各 sample ディレクトリに
  #           matrix.mtx.gz / features.tsv.gz / barcodes.tsv.gz）。最も一般的。
  # "dense" : sample ごとの密行列 gz（GEO の processed matrix 等。tab 区切り、
  #           1 行目=バーコード、2 行目以降=遺伝子+値）。dogfood 参照実装と同形式。
  loader   = "tenx",
  data_dir = "data/raw",         # 生データの場所（fetch_data.sh でここへ置く）

  # --- dense ローダ専用（loader="dense" のときだけ使う）---
  dense_pattern      = "\\.txt\\.gz$",             # 生データファイル名にマッチする正規表現
  dense_sample_regex = "^(.*)\\.txt\\.gz$",        # ファイル名から sample 名を取り出す（\\1 が sample）

  # --- 細胞フィルタの基準（作成時）---
  min_cells = 3, min_features = 200,

  # --- QC 閾値（データを見て調整。まず violin plot で分布を見るのが定石）---
  qc = list(nFeature_min = 200, nFeature_max = 6000, percent_mt_max = 20),
  mito_pattern = "^MT-",         # ミトコンドリア遺伝子: ヒト="^MT-" / マウス="^mt-"

  # --- 統合（Harmony）と次元削減 ---
  batch_key   = "sample",        # バッチ効果を除くキー（複数 sample を統合する軸）
  n_variable  = 2000,            # 高変動遺伝子の数
  n_pcs       = 30,              # PCA 次元数
  dims        = 30,              # 近傍/UMAP に使う次元数
  resolution  = 0.5,             # クラスタ解像度（大きいほど細かく分かれる）

  # --- アノテーションの当たり付け（存在すればマーカー数を報告するだけの目安）---
  marker_hint = c(EPCAM = "epithelial", PTPRC = "immune", CD3D = "T", CD79A = "B",
                  LYZ = "myeloid", PECAM1 = "endothelial", ACTA2 = "fibro/smooth")
)
