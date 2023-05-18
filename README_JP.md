# EXPRORER_MSMD

共溶媒分子動力学 (mixed-solvent molecular dynamics; MSMD) シミュレーションエンジンと、解析ツールをまとめたリポジトリ

## システムの説明

GROMACSを用いたMSMDを自動的に行うためのシステム。

## 簡単な使い方

### 環境構築
このgithubにはDockerファイルが含まれており、dockerを用いることで実行可能な環境を容易に作ることができる。

中で利用されている特筆すべきアプリは以下の通り。

- python 3 (`$PYTHON`)
  - 利用しているライブラリは `.devcontainer/Dockerfile` を参照のこと。
- AmberTools 20 (`$TLEAP`, `$CPPTRAJ`)
- Gromacs 2021.5 (`$GMX`)
- packmol 18.169 (`$PACKMOL`)

### MSMDシミュレーションの実施 ( `exprorer_msmd` )

タンパク質、共溶媒（プローブ）分子、シミュレーションプロトコルを定義したyamlファイルを用意し（例えば `example/example_protocol.yaml` ）、
以下のコマンドを実行することで系の構築からシミュレーションの実施、確率密度マップ (spatial probability distribution map; PMAP) の作成を自動的に行う。

```
./exprorer_msmd example/example_protocol.yaml
```

なお、このシミュレーションの手順は 

- MSMDの系を作成する **preprocess**
- GROMACSを用いたシミュレーションを実施する **simulation**
- シミュレーション結果からPMAPファイル等を作成する **postprocess**

の3つに分かれており、それぞれは `--skip-preprocess` 、 `--skip-simulation` 、 `--skip-postprocess` で飛ばすことができる。

### MSMDシミュレーション結果の解析

#### シミュレーションのトラジェクトリを確認する

`./exprorer_msmd` を実行すると、出力ディレクトリの中に独立試行ごとに `system` フォルダが作成される。この中に、`[プロジェクト名]_woWAT_10ps.pdb` など、いくつかのPDBファイルが作成されるが、これがMSMDのトラジェクトリである。PyMOLで当該ファイルを開くことで、シミュレーション中のタンパク質とプローブの動きを見ることができる。（水分子は除外されている）

#### タンパク質のホットスポットを探索する ( `protein_hotspot` )

タンパク質表面のどの部分にプローブ分子が存在しやすいか？という
ホットスポット探索を行う。

```
./protein_hotspot example/example_protocol.yaml
```

これによって、出力ディレクトリ直下に `maxPMAP_[プロジェクト名]_nV.dx` などのOpenDXフォーマットのボクセルが作成される。このボクセルは、入力タンパク質構造に重なるように作成されている。これはPyMOLで読み込むことができ、
PyMOLの `isomesh` コマンドを用いることで描画することができる。

実行結果は以下のようなものになる（これは複数のプローブの計算結果を複数重ね合わせて表示している）。
![ホットスポット探索](https://i.imgur.com/bzxz0K6.png)

EXPRORER[^1] はこの結果を用いて、ボクセル間の類似度を計算、プローブ同士の類似度を計算している。

[^1]:**Keisuke Yanagisawa**, Yoshitaka Moriwaki, Tohru Terada, Kentaro Shimizu. "EXPRORER: Rational Cosolvent Set Construction Method for Cosolvent Molecular Dynamics Using Large-Scale Computation", *Journal of Chemical Information and Modeling*, **61**: 2744-2753, 2021/06. DOI: [10.1021/acs.jcim.1c00134](https://doi.org/10.1021/acs.jcim.1c00134)

#### プローブ分子周辺残基環境[^2]を取得する ( `probe_profile` )

前述のホットスポット解析とは反対に、プローブ分子の周囲にどのような残基が存在しやすいか？を描画する。

```
./protein_hotspot example/example_protocol.yaml
```

これによって、出力ディレクトリ直下に `[プロジェクト名]_[プローブ名]_mesh_anion.dx` など8種類のOpenDXフォーマットのボクセルが作成される。これは `alignedresenv_[プロジェクト名].pdb` に重なるように作成されており、これもPyMOLの `isomesh` コマンドで描画できる。

実行結果は以下のようなものになる。
![プローブ分子周辺残基環境](https://i.imgur.com/4QIZxhW.png)

[^2]:**Keisuke Yanagisawa**, Ryunosuke Yoshino, Genki Kudo, Takatsugu Hirokawa. "Inverse Mixed-Solvent Molecular Dynamics for Visualization of the Residue Interaction Profile of Molecular Probes", *International Journal of Molecular Sciences*, **23**: 4749, 2022/04. DOI: [10.3390/ijms23094749](https://doi.org/10.3390/ijms23094749)

## yamlファイルの書き換え

すべての設定は yaml ファイルに記述されている。
（工事中）