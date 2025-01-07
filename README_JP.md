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

## 詳細な動作説明

### `exprorer_msmd`

前述のように、以下のようなコマンドを実行することで、MSMDのシミュレーションを自動的に行うことができる。

```
./exprorer_msmd example/example_protocol.yaml
```

このコマンドを実行した場合に、どのような動作をするのか、その詳細を以下に示す。ただし、以下のことを前提とする。

- yamlファイル内の設定を参照する場合は `yaml[general][workdir]` のように記述する。
- 全ての相対パスは、yamlファイルが存在するディレクトリを基準とする。
- 「一時ファイル」と明記されている場合は、 `/tmp` ディレクトリ等に作成され、処理後に自動で破棄される。
- `example_protocol.yaml` は、 `yaml[general][iter_index]` が `0,1` となっており、以下の処理が 実行ID が 0 の場合と 1 の場合の計2回行われるが、以下の説明ではそのうちの1回分の処理に注目して記述する。

![exprorer_msmdの流れ図](https://i.imgur.com/EVEO8PV.png)

#### 1. 前処理 (preprocess)

まず、MSMDの系を作成する前処理 (preprocess) が実行される。
- **事前準備**
  - `yaml[general][workdir]` を参照し、各種ファイルを出力するための `../output/example_protocol` ディレクトリを作成する。
- **プローブの準備**
  - `yaml[input][probe][cid]` を参照し、プローブの構造を示す mol2 ファイルのパス `A11.mol2` を決定する。
  - `parmchk2` コマンドを利用して .frcmod 一時ファイルを作成する。
    - この際GAFF (general Amber force field) を用いるか GAFF2 を用いるかは `yaml[input][probe][atomtype]` で指定する。
- **タンパク質の準備**
  - `yaml[input][protein][pdb]` を参照し、タンパク質構造を記録している pdb ファイルのパス `1F47B.pdb` を決定する。
  - `1F47B.pdb` から、C末端の酸素原子 OXT を削除し、更新後の pdb ファイルを一時ファイルとして保存する。
  - タンパク質と水分子からなる系を `tleap` を用いて構築し、作成された系の直方体のサイズを記録する。
    - この際、タンパク質から 10 Å 離れた位置まで水分子で充填する。
  - 先述の直方体のサイズのうち、最も長い辺の長さを求め、1辺がその長さであるような立方体領域をMSMDの系のサイズとする。
    - タンパク質が長細い形状である場合など、境界をまたいでタンパク質同士が干渉することを防ぐために、立方体領域に拡張している。
- **MSMDの系の作成**
  - `yaml[input][probe][molar]` と、先述した立方体のサイズから、プローブの分子数を決定する。
    - このケースでは、約50個のプローブ分子が含まれる。
  - `packmol` コマンドを用いて、タンパク質と複数のプローブ分子のみからなる系を作成し、一時ファイルとして保存する。
    - この際、乱数のシード値を異なる値に設定することで、異なるプローブの初期配置を得る。
  - 作成された系に対して、 `tleap` を用いて水分子やイオン（Na+およびCl-）を充填し、水分子を含んだ系を作成する。
    - 既に系のサイズは決定されているため、水分子の充填はそのサイズに合わせて行う。
    - タンパク質およびプローブ分子が持つ電荷を中和する必要最低限のイオンを追加する。
  - 上記の系を、一旦 .parm7 および .rst7 一時ファイルとして保存する。
  - .parm7 および .rst7 ファイルを、 .top および .gro 一時ファイルに変換する。
  - プローブ分子間にLJポテンシャルによる仮想反発項を与える為に、プローブの重心座標に仮想原子を配置するように .gro ファイルを更新する。
    - プローブの重心座標を計算する際にはプローブが持つ水素原子も考慮する。
  - .top ファイルに、LJポテンシャルによる仮想反発項を与えるための設定を追加する。
    - 仮想原子の位置は GROMACS のトポロジーファイルの `virtual_sitesn` を使って固定する。
    - 仮想反発項のパラメータは $\epsilon_{LJ} = 10^{-6} \mathrm{kcal/mol}, R_{min} = 20.0 \mathrm{\AA}$ とする。
  - .top ファイルに、平衡化のための位置拘束設定を追加する。
    - 位置拘束はタンパク質の全ての重原子に設定する。
    - 位置拘束は $1000, 500, 200, 100, 50, 20, 10, 0 [\mathrm{kcal/mol/\AA^2}]$ の8種類を作成する。
- **前処理後のファイルの保存**
  - 作成された .top および .gro ファイルを、 `../output/example_protocol/system*/prep` ディレクトリに保存する。
    - `*` 部分には、前述した実行ID が入る。
    - ファイル名は `yaml[general][name]` に基づいて、 `TEST_PROJECT.top` および `TEST_PROJECT.gro` となる。
  - 後処理 postprocess の為に、 `TEST_PROJECT.gro` に対応する `TEST_PROJECT.pdb` ファイルを作成する。

#### 2. MSMD シミュレーションの実施 (simulation)

続いて、GROMACSを用いたMSMDのシミュレーションを実施する。
- `yaml[general][workdir]` を参照し、シミュレーション関係のファイルを出力するための `../output/example_protocol/system*/simulation` ディレクトリを作成する。
- `gmx make_ndx` コマンドを利用して、系内の原子群をグループ化するための `index.ndx` ファイルを作成する。
- `yaml[exprorer_msmd][sequence]` を参照して、 .mdp ファイルを作成する。
  - .mdp ファイルの名前は、 `yaml[exprorer_msmd][sequence]` の各ステップの `name` に基づいて、 `min1.mdp`, `heat.mdp`, `equil1.mdp`, `pr.mdp` などとなる。
  - 各ステップは 4種類のいずれかの `type` である。
    - minimization: エネルギー最小化を行う
    - heating: 系の温度を上昇させる
    - equilibration: 系の平衡化を行う
    - production: シミュレーションを実施する
- 各ステップを順に実行するためのシェルスクリプト `mdrun.sh` を作成し、実行することでシミュレーションを実施する。
  - 各ステップは直前の .gro ファイルを入力として、トラジェクトリファイル .xtc や 各ステップ終了時の全原子の位置および速度が記録された .gro ファイルなどを出力する。
    - 各ステップを通して、シミュレーションのステップ幅は `yaml[exprorer_msmd][general][dt]` で指定された値を用いる。デフォルトでは 0.002 ps = 2 fs となっている。
    - トラジェクトリの出力頻度は `yaml[exprorer_msmd][sequence][nstxtcout]` で指定された値を用いる。この例では、production runにおける `nstxtcout` は 5000 となっており、これは 2 fs * 5000 steps = 10 ps ごとにフレームを出力する事に相当する。
  - 最後のステップが出力した .xtc ファイル（この例では `pr.xtc`）のシンボリックリンク `TEST_PROJECT.xtc` を作成する。

#### 3. 後処理 (postprocess)

最後に、MSMDシミュレーション結果から、プローブ分子の空間分布を示す確率密度マップ (spatial probability distribution map, PMAP) や grid free energy (GFE) などを計算する。

- `yaml[map][snapshot]` を参照して、PMAPを作成するためのスナップショットを決定する。
  - この例では、 `TEST_PROJECT.xtc` の2001フレーム目から4000フレーム目を利用する。前述の通り、production runでは 10 ps ごとにフレームを書き出しているので、 20--40 ns のシミュレーション結果を利用することに相当する。
- `cpptraj` の `grid` コマンドを用いて、各原子の存在頻度を計算し、存在頻度マップとして保存する。
  - この時点では頻度であり、確率化していないため、 PMAP ではない事に注意せよ。
  - タンパク質構造は、初期構造の `../output/example_protocol/system*/prep/TEST_PROJECT.pdb` を参照構造として構造重ね合わせを行う。
  - 存在頻度を計算する領域のサイズは、 `yaml[map][map_size]` で指定された値を用いる。この例では、 1辺が$80\mathrm{\AA}$ の立方体領域としている。
  - 存在頻度を計算する原子群の種類は `yaml[map][maps]` に設定することが出来る。この例では、プローブ分子の重原子のみを考慮した　`nVH`、プローブ分子の全原子を考慮した `nV` の2種類を作成する。
    - `@VIS` はプローブの仮想原子を指す。仮想原子はプローブの重心であるから、 `@VIS` の存在頻度をカウントすることは、プローブの重心がタンパク質表面上の存在頻度を求めることに相当する。
    - `suffix` として `nVH` が指定された場合は、 `../output/example_protocol/system*/TEST_PROJECT_nVH.dx` というファイル名で原子の存在頻度が保存される。
- 存在頻度マップをもとに、確率化を行い、PMAPを作成する。
  - `yaml[map][normalization]` の設定によって、3種類の異なるマップを作成する。
    - `total` の場合、全ての原子の存在頻度の合計が1になるように規格化し、 `PMAP_TEST_PROJECT_nVH.dx` などのファイルとして出力する。
    - `snapshot` の場合、存在頻度をスナップショットの数で割ることで、各座標におけるプローブの存在確率を求めて、 `PMAP_TEST_PROJECT_nVH.dx` などのファイルとして出力する。
    - `GFE` の場合、PMAP ではなく grid free energy (GFE) を計算し、 `GFE_PMAP_TEST_PROJECT_nVH.dx` などのファイルとして出力する。
