# EXPRORER_MSMD 使用ガイド

## 1. プロジェクト概要

EXPRORER_MSMDは、共溶媒分子動力学（Mixed-Solvent Molecular Dynamics; MSMD）シミュレーションエンジンと解析ツールをまとめたリポジトリです。このプロジェクトは、GROMACSを用いたMSMDを自動的に行うためのシステムを提供しています。

### 主要な機能

- **MSMDシミュレーション**: GROMACSを使用した自動MSMDシミュレーション
- **タンパク質ホットスポット探索**: タンパク質表面のどの部分にプローブ分子が存在しやすいかを解析
- **プローブ分子周辺残基環境の取得**: プローブ分子の周囲にどのような残基が存在しやすいかを可視化
- **Quantitative Inverse MSMD**: 分子の特定の原子を置換した場合の結合親和性の変化を見積もる

### プロジェクトの全体構造

プロジェクトは以下のような構造になっています：

- `exprorer_msmd`: MSMDシミュレーションを実行するメインスクリプト
- `protein_hotspot`: タンパク質のホットスポット探索を行うスクリプト
- `probe_profile`: プローブ分子周辺残基環境を取得するスクリプト
- `script/`: 各種ユーティリティスクリプトを含むディレクトリ
  - `quantitative_inverse_msmd.py`: Quantitative Inverse MSMDを実行するスクリプト
  - `utilities/`: 各種ユーティリティモジュールを含むディレクトリ
    - `molecular/`: 分子操作関連のモジュール
    - `probe/`: プローブ分子関連のモジュール
- `example/`: サンプルファイルを含むディレクトリ
- `docs/`: ドキュメントを含むディレクトリ

## 2. 環境構築

### 必要なソフトウェアとライブラリ

EXPRORER_MSMDを使用するには、以下のソフトウェアとライブラリが必要です：

- **Python 3**
  - NumPy
  - RDKit
  - Matplotlib
  - PIL (Pillow)
  - PyYAML
  - その他の依存ライブラリ（`.devcontainer/Dockerfile`を参照）
- **AmberTools 20** (`$TLEAP`, `$CPPTRAJ`)
- **Gromacs 2021.5** (`$GMX`)
- **Packmol 18.169** (`$PACKMOL`)

### Dockerを使用した環境構築方法

Dockerを使用すると、必要な環境を簡単に構築できます。

1. Dockerをインストールします（[Docker公式サイト](https://www.docker.com/get-started)からダウンロード）

2. リポジトリをクローンします：
   ```bash
   git clone https://github.com/your-username/exprorer_msmd.git
   cd exprorer_msmd
   ```

3. Dockerイメージをビルドします：
   ```bash
   docker build -t exprorer_msmd -f .devcontainer/Dockerfile .
   ```

4. Dockerコンテナを実行します：
   ```bash
   docker run -it --rm -v $(pwd):/workdir exprorer_msmd
   ```

これにより、必要なすべてのソフトウェアとライブラリがインストールされた環境が構築されます。

### 手動での環境構築方法

手動で環境を構築する場合は、以下の手順に従ってください：

1. Python 3をインストールします（[Python公式サイト](https://www.python.org/downloads/)からダウンロード）

2. 必要なPythonライブラリをインストールします：
   ```bash
   pip install numpy rdkit-pypi matplotlib pillow pyyaml
   ```

3. AmberTools 20をインストールします（[AmberMD公式サイト](https://ambermd.org/GetAmber.php)を参照）

4. Gromacs 2021.5をインストールします（[Gromacs公式サイト](https://www.gromacs.org/Downloads)を参照）

5. Packmolをインストールします（[Packmol公式サイト](https://m3g.github.io/packmol/)を参照）

6. 環境変数を設定します：
   ```bash
   export TLEAP=/path/to/tleap
   export CPPTRAJ=/path/to/cpptraj
   export GMX=/path/to/gmx
   export PACKMOL=/path/to/packmol
   ```

7. 動作確認を行います：
   ```bash
   ./exprorer_msmd --help
   ```

## 3. 基本的な使用方法

### MSMDシミュレーションの実行方法

MSMDシミュレーションを実行するには、タンパク質、共溶媒（プローブ）分子、シミュレーションプロトコルを定義したYAMLファイルを用意し、`exprorer_msmd`コマンドを実行します。

#### YAMLファイルの作成方法

YAMLファイルには以下の情報を記述します：

```yaml
general:
  name: "project_name"
  workdir: "./output"
  iter_index: "1-3"  # 3回の独立試行
  multiprocessing: true

input:
  protein:
    pdb: "example/protein.pdb"
  probe:
    sdf: "example/probe.sdf"
    cid: "probe_id"

simulation:
  # シミュレーションパラメータ
  temperature: 300
  pressure: 1.0
  timestep: 0.002
  nsteps: 500000  # 1 ns
```

詳細な設定項目については、`example/example_protocol.yaml`を参照してください。

#### exprorer_msmdコマンドの使用方法

```bash
./exprorer_msmd example/example_protocol.yaml
```

#### オプション

- `--skip-preprocess`: 前処理をスキップします
- `--skip-simulation`: シミュレーションをスキップします
- `--skip-postprocess`: 後処理をスキップします

#### 出力ファイル

MSMDシミュレーションを実行すると、以下のファイルが生成されます：

- `output/[プロジェクト名]/system/`: シミュレーションシステムのディレクトリ
  - `[プロジェクト名]_woWAT_10ps.pdb`: 水分子を除いたトラジェクトリ（10ps間隔）
  - その他のシミュレーション関連ファイル

### タンパク質ホットスポット探索の方法

タンパク質表面のどの部分にプローブ分子が存在しやすいかを探索するには、`protein_hotspot`コマンドを使用します。

```bash
./protein_hotspot example/example_protocol.yaml
```

#### 出力ファイル

- `output/[プロジェクト名]/maxPMAP_[プロジェクト名]_nV.dx`: OpenDXフォーマットのボクセルファイル

#### PyMOLでの可視化方法

PyMOLで以下のコマンドを実行することで、ホットスポットを可視化できます：

```python
load example/protein.pdb
load output/[プロジェクト名]/maxPMAP_[プロジェクト名]_nV.dx
isomesh mesh, maxPMAP_[プロジェクト名]_nV, 0.5
```

### プローブ分子周辺残基環境の取得方法

プローブ分子の周囲にどのような残基が存在しやすいかを可視化するには、`probe_profile`コマンドを使用します。

```bash
./probe_profile example/example_protocol.yaml
```

#### 出力ファイル

- `output/[プロジェクト名]/[プロジェクト名]_[プローブ名]_mesh_anion.dx`: 陰イオン性残基の分布
- `output/[プロジェクト名]/[プロジェクト名]_[プローブ名]_mesh_cation.dx`: 陽イオン性残基の分布
- その他6種類のOpenDXフォーマットのボクセルファイル

#### 可視化方法

PyMOLで以下のコマンドを実行することで、残基環境を可視化できます：

```python
load output/[プロジェクト名]/alignedresenv_[プロジェクト名].pdb
load output/[プロジェクト名]/[プロジェクト名]_[プローブ名]_mesh_anion.dx
isomesh mesh_anion, [プロジェクト名]_[プローブ名]_mesh_anion, 0.5
```

## 4. Quantitative Inverse MSMDの使用方法

### 概要と目的

Quantitative Inverse MSMDは、分子の特定の原子を置換した場合の結合親和性の変化を見積もるための手法です。2つの極めて構造が類似した化合物間の結合自由エネルギーの差を求めることができます。

### モジュール構造

Quantitative Inverse MSMDモジュールは、以下のディレクトリ構造に整理されています：

```
script/
├── quantitative_inverse_msmd.py  # メインスクリプト
└── utilities/
    ├── analysis/                 # 解析関連のモジュール
    │   ├── __init__.py
    │   ├── profile.py            # プロファイル関連のクラスと関数
    │   ├── matching_score.py     # マッチングスコア計算関連の関数
    │   └── compare_results.py    # 結果比較関連の関数
    ├── molecular/                # 分子操作関連のモジュール
    │   ├── mcs_extractor.py      # 非共通部分抽出関連の関数
    │   ├── molecule_superimposer.py # 分子重ね合わせ関連の関数
    │   ├── simple_atom_replacer.py  # 原子置換関連の関数
    │   └── extract_substructure.py  # 部分構造抽出関連の関数
    └── probe/                    # プローブ関連のモジュール
        └── probe_designer.py     # プローブ設計関連の関数
```

### 実行手順の詳細

#### 入力ファイルの準備

1. 元の分子のSDFファイルを用意します
2. 置換後の分子のSDFファイルを用意します（または`script/utilities/molecular/simple_atom_replacer.py`を使用して作成）
3. タンパク質のPDBファイルを用意します
4. プローブライブラリを用意します

#### simple_atom_replacer.pyを使用した原子置換

```bash
python script/utilities/molecular/simple_atom_replacer.py --mol example/ZINC000000330081.sdf --atom-idx 3 --new-atom N --output-dir output
```

このコマンドにより、以下のファイルが生成されます：
- `output/original_molecule.sdf`: 元の分子
- `output/replaced_molecule.sdf`: 4番目のCをNに置換した分子

#### quantitative_inverse_msmd.pyスクリプトの使用方法

```bash
python script/quantitative_inverse_msmd.py \
  --compound1 output/original_molecule.sdf \
  --compound2 output/replaced_molecule.sdf \
  --protein example/protein.pdb \
  --probe-library probe_library \
  --output-dir output/qimsmd \
  --only-non-common-rings \
  --separate-rings
```

#### オプション

- `--compound1`: 1つ目の化合物のSDFファイルパス（必須）
- `--compound2`: 2つ目の化合物のSDFファイルパス（必須）
- `--protein`: タンパク質構造のPDBファイルパス（必須）
- `--probe-library`: プローブライブラリのディレクトリパス（必須）
- `--output-dir`: 出力ディレクトリ（デフォルト: "./output"）
- `--only-non-common-rings`: 複合環のうち非共通な環のみをプローブ化します
- `--separate-rings`: 非共通な環を別々のプローブとして抽出します

### 結果の解釈方法

Quantitative Inverse MSMDを実行すると、以下のファイルが生成されます：

- `output/qimsmd/quantitative_inverse_msmd_result.txt`: 結果ファイル
  - 化合物1と化合物2の合致度スコア
  - 結合強度差
  - 結論（どちらの化合物が強く結合するか）

結合強度差が正の場合は化合物1の方が強く結合し、負の場合は化合物2の方が強く結合すると予測されます。

### 処理フロー

Quantitative Inverse MSMDの処理フローは以下の通りです：

1. **非共通部分の抽出**：2つの化合物間の非共通部分を抽出します
2. **プローブの設計**：非共通部分に基づいてプローブを設計します
3. **分子の重ね合わせ**：元の化合物とプローブを重ね合わせます
4. **MSMDシミュレーション**：プローブとタンパク質を使用してMSMDシミュレーションを実行します
5. **プロファイル作成**：MSMDシミュレーションの結果からプロファイルを作成します
6. **合致度計算**：プロファイルとタンパク質構造の合致度を計算します
7. **結合強度の推定**：合致度スコアの差から結合強度の差を推定します

## 5. 各ユーティリティツールの使用方法と入出力ファイル形式

### 分子の重ね合わせ自動化ツール（molecule_superimposer.py）

#### 機能の説明

このツールは、化合物とプローブの原子対応関係を自動生成し、Kabsch algorithmを使用して最適な重ね合わせを計算します。

#### コマンドライン引数の詳細

```bash
python script/utilities/molecular/molecule_superimposer.py \
  --compound example/compound.sdf \
  --probe example/probe.sdf \
  --output output/transformed_probe.sdf
```

#### 入力ファイル形式

- `compound.sdf`: 化合物のSDFファイル（3D座標を含む）
- `probe.sdf`: プローブのSDFファイル（3D座標を含む）

SDFファイルは、分子の原子座標や結合情報を含む標準的な化学情報ファイル形式です。

#### 出力ファイル形式と内容

- `transformed_probe.sdf`: 重ね合わせ後のプローブのSDFファイル
- `output/compound_mapping.png`: 化合物の原子対応関係を示す画像
- `output/probe_mapping.png`: プローブの原子対応関係を示す画像
- `output/molecules_comparison.png`: 両方の分子を並べて表示した画像
- `output/superimposed.png`: 重ね合わせ結果の可視化画像

### MCS計算と非共通部分抽出ツール（mcs_extractor.py）

#### 機能の説明

このツールは、2つの化合物間のMaximum Common Substructure（MCS）を計算し、非共通部分を抽出します。

#### コマンドライン引数の詳細

```bash
python script/utilities/molecular/mcs_extractor.py \
  --mol1 example/compound1.sdf \
  --mol2 example/compound2.sdf \
  --output-dir output/mcs \
  --only-non-common-rings \
  --separate-rings
```

#### 入力ファイル形式

- `compound1.sdf`: 1つ目の化合物のSDFファイル
- `compound2.sdf`: 2つ目の化合物のSDFファイル

#### 出力ファイル形式と内容

- `output/mcs/non_common_mol1.sdf`: 1つ目の化合物の非共通部分
- `output/mcs/non_common_mol2.sdf`: 2つ目の化合物の非共通部分

`--separate-rings`オプションを使用した場合：
- `output/mcs/non_common_mol1_ring1.sdf`, `non_common_mol1_ring2.sdf`, ...
- `output/mcs/non_common_mol2_ring1.sdf`, `non_common_mol2_ring2.sdf`, ...

### プローブ設計ツール（probe_designer.py）

#### 機能の説明

このツールは、非共通部分に合うプローブ分子を設計します。

#### コマンドライン引数の詳細

```bash
python probe_designer.py \
  --mol1 example/compound1.sdf \
  --mol2 example/compound2.sdf \
  --probe-library probe_library \
  --output-dir output/probes
```

#### 入力ファイル形式

- `compound1.sdf`, `compound2.sdf`: 化合物のSDFファイル
- `probe_library/`: プローブライブラリのディレクトリ（SDFファイルを含む）

#### 出力ファイル形式と内容

- `output/probes/probe1.sdf`: 1つ目の化合物の非共通部分に対応するプローブ
- `output/probes/probe2.sdf`: 2つ目の化合物の非共通部分に対応するプローブ

### 部分構造抽出ツール（extract_substructure_with_drawing.py）

#### 機能の説明

このツールは、分子から指定した原子の周りの部分構造を抽出し、PNG画像として描画します。

#### コマンドライン引数の詳細

```bash
python extract_substructure_with_drawing.py \
  --mol example/molecule.sdf \
  --atom-idx 3 \
  --radius 2 \
  --output-dir output/substructure
```

#### 入力ファイル形式

- `molecule.sdf`: 分子のSDFファイル

#### 出力ファイル形式と内容

- `output/substructure/molecule_substructure.sdf`: 抽出された部分構造のSDFファイル
- `output/substructure/molecule_substructure.png`: 部分構造の画像
- `output/substructure/molecule_highlighted.png`: 元の分子の中で部分構造をハイライトした画像
- `output/substructure/molecule_mapping.txt`: 元の分子のインデックスと部分構造のインデックスのマッピング

### SMILES表記からSDFファイルを生成するツール（smiles_to_sdf.py）

#### 機能の説明

このツールは、SMILES表記から3D座標を持つSDFファイルを生成します。

#### コマンドライン引数の詳細

```bash
python script/utilities/molecular/smiles_to_sdf.py \
  --smiles "CC(=O)Oc1ccccc1C(=O)O" \
  --output output/molecule.sdf
```

#### 入力SMILES形式の説明

SMILES（Simplified Molecular Input Line Entry System）は、分子構造を表現するための文字列表記法です。例えば、アスピリンのSMILESは`CC(=O)Oc1ccccc1C(=O)O`です。

#### 出力SDFファイル形式の説明

- `output/molecule.sdf`: 生成されたSDFファイル（3D座標を含む）

#### 3D座標生成と最適化の詳細

このツールは、RDKitを使用してSMILES表記から分子を生成し、3D座標を計算します。その後、UFFまたはMMFFを使用してエネルギー最小化を行います。

## 6. 実行例とサンプルコマンド

### 一般的なワークフロー例

#### 基本的なMSMDシミュレーションの実行例

```bash
# MSMDシミュレーションの実行
./exprorer_msmd example/example_protocol.yaml

# シミュレーション結果の確認
ls output/example_project/system/
```

#### タンパク質ホットスポット探索の実行例

```bash
# ホットスポット探索の実行
./protein_hotspot example/example_protocol.yaml

# 結果の確認
ls output/example_project/
```

#### Quantitative Inverse MSMDの実行例

```bash
# 原子置換
python script/utilities/molecular/simple_atom_replacer.py \
  --mol example/ZINC000000330081.sdf \
  --atom-idx 3 \
  --new-atom N \
  --output-dir output

# Quantitative Inverse MSMDの実行
python script/quantitative_inverse_msmd.py \
  --compound1 output/original_molecule.sdf \
  --compound2 output/replaced_molecule.sdf \
  --protein example/protein.pdb \
  --probe-library probe_library \
  --output-dir output/qimsmd \
  --only-non-common-rings

# 結果の確認
cat output/qimsmd/quantitative_inverse_msmd_result.txt
```

### コマンドラインオプションの詳細説明

#### exprorer_msmdのオプション

- `--skip-preprocess`: 前処理をスキップします
- `--skip-simulation`: シミュレーションをスキップします
- `--skip-postprocess`: 後処理をスキップします
- `--help`: ヘルプメッセージを表示します

#### quantitative_inverse_msmd.pyのオプション

- `--compound1`: 1つ目の化合物のSDFファイルパス（必須）
- `--compound2`: 2つ目の化合物のSDFファイルパス（必須）
- `--protein`: タンパク質構造のPDBファイルパス（必須）
- `--probe-library`: プローブライブラリのディレクトリパス（必須）
- `--output-dir`: 出力ディレクトリ（デフォルト: "./output"）
- `--only-non-common-rings`: 複合環のうち非共通な環のみをプローブ化する
- `--separate-rings`: 非共通な環を別々のプローブとして抽出する
- `--probe-library`: プローブライブラリのディレクトリパス
- `--output-dir`: 出力ディレクトリ
- `--only-non-common-rings`: 複合環のうち非共通な環のみをプローブ化する
- `--separate-rings`: 非共通な環を別々のプローブとして抽出する

### 入出力ファイルの詳細

#### YAMLファイル形式

YAMLファイルは、階層構造を持つ設定ファイルです。以下は基本的な構造です：

```yaml
general:
  # 一般的な設定
input:
  # 入力ファイルの設定
simulation:
  # シミュレーションの設定
```

#### SDFファイル形式

SDFファイルは、分子の原子座標や結合情報を含む標準的な化学情報ファイル形式です。以下は基本的な構造です：

```
分子名
  コメント行
  次のコメント行
  原子数 結合数 0 0 0 0 0 0 0 0 0 V2000
    x1    y1    z1 原子種1 0  0  0  0  0  0  0  0  0  0  0  0
    x2    y2    z2 原子種2 0  0  0  0  0  0  0  0  0  0  0  0
    ...
  原子1 原子2 結合種 0  0  0  0
  原子1 原子3 結合種 0  0  0  0
  ...
M  END
$$$$
```

#### PDBファイル形式

PDBファイルは、タンパク質や核酸の3D構造を記述するための標準的なファイル形式です。以下は基本的な構造です：

```
HEADER    タイトル
ATOM      1  N   ALA A   1      11.804  18.255  17.872  1.00  0.00           N
ATOM      2  CA  ALA A   1      11.804  18.255  17.872  1.00  0.00           C
...
END
```

#### OpenDXファイル形式

OpenDXファイルは、3D格子データを記述するためのファイル形式です。PMAPやホットスポット情報を格子データとして保存します。

## 7. トラブルシューティング

### 一般的な問題と解決方法

#### 環境構築時の問題

- **問題**: Dockerイメージのビルドに失敗する
  - **解決策**: Dockerのインストールを確認し、十分なディスク容量があることを確認してください

- **問題**: 必要なライブラリのインストールに失敗する
  - **解決策**: Pythonのバージョンを確認し、必要に応じて仮想環境を使用してください

#### 実行時のエラー

- **問題**: `exprorer_msmd`の実行に失敗する
  - **解決策**: 環境変数（`$TLEAP`, `$CPPTRAJ`, `$GMX`, `$PACKMOL`）が正しく設定されているか確認してください

- **問題**: シミュレーションが途中で停止する
  - **解決策**: メモリ不足の可能性があります。より小さなシステムで試すか、より多くのメモリを割り当ててください

#### 結果の解釈に関する問題

- **問題**: ホットスポットが表示されない
  - **解決策**: 閾値を調整してみてください。例えば、`isomesh mesh, maxPMAP_[プロジェクト名]_nV, 0.3`のように閾値を下げてみてください

- **問題**: Quantitative Inverse MSMDの結果が不安定
  - **解決策**: 複数回実行して結果の平均を取るか、`--only-non-common-rings`オプションを試してみてください

### エラーメッセージの解釈

- **エラー**: `ModuleNotFoundError: No module named 'rdkit'`
  - **原因**: RDKitがインストールされていない
  - **対処法**: `pip install rdkit-pypi`を実行してRDKitをインストールしてください

- **エラー**: `Command 'gmx' not found`
  - **原因**: Gromacsがインストールされていないか、パスが通っていない
  - **対処法**: Gromacsをインストールし、環境変数`$GMX`を設定してください

- **エラー**: `ValueError: 分子ファイル ... を読み込めませんでした`
  - **原因**: SDFファイルの形式が不正か、ファイルが存在しない
  - **対処法**: ファイルの存在と形式を確認してください

## 8. 参考資料

### 関連論文と文献

- **Keisuke Yanagisawa**, Yoshitaka Moriwaki, Tohru Terada, Kentaro Shimizu. "EXPRORER: Rational Cosolvent Set Construction Method for Cosolvent Molecular Dynamics Using Large-Scale Computation", *Journal of Chemical Information and Modeling*, **61**: 2744-2753, 2021/06. DOI: [10.1021/acs.jcim.1c00134](https://doi.org/10.1021/acs.jcim.1c00134)

- **Keisuke Yanagisawa**, Ryunosuke Yoshino, Genki Kudo, Takatsugu Hirokawa. "Inverse Mixed-Solvent Molecular Dynamics for Visualization of the Residue Interaction Profile of Molecular Probes", *International Journal of Molecular Sciences*, **23**: 4749, 2022/04. DOI: [10.3390/ijms23094749](https://doi.org/10.3390/ijms23094749)

### 詳細なドキュメントへのリンク

- [基本的な使い方](docs/ja/user_guide/basic.md)
- [高度な使用方法](docs/ja/user_guide/advanced.md)
- [プローブ分子の準備](docs/ja/user_guide/probe_preparation.md)
- [MSMDの実装](docs/ja/impl/msmd.md)
- [PMAPの詳細](docs/ja/impl/pmap.md)