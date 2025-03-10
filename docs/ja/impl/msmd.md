# 実装内容

Mixed-Solvent Molecular Dynamics (MSMD) は、タンパク質と複数の溶媒分子（プローブ）の相互作用を解析するための分子動力学シミュレーション手法です。

## 目次
1. [MSMDの基本](#msmdの基本)
   - [プローブ分子](#プローブ分子)
   - [プローブ分子間の仮想的な反発力](#プローブ分子間の仮想的な反発力)
2. [MSMDシミュレーションプロトコル](#msmdシミュレーションプロトコル)
   - [エネルギー最小化](#エネルギー最小化)
   - [系の加熱](#系の加熱)
   - [位置拘束の段階的緩和](#位置拘束の段階的緩和)
   - [プロダクションラン](#プロダクションラン)
3. [実装の詳細](#実装の詳細)
   - [前処理](#1-前処理-preprocess)
   - [シミュレーション実行](#2-シミュレーション実行-execute_single_simulation)
   - [後処理](#3-後処理-postprocess)

関連ドキュメント：
- [PMAPの詳細](pmap.md)
- [基本的な使い方](../user_guide/basic.md)
- [プローブ分子の準備](../user_guide/probe_preparation.md)

## MSMDの基本

MSMDでは、以下の特徴を持つ分子動力学シミュレーションを実行します。

1. 系の中に複数のプローブ分子（溶媒分子）が存在する
2. プローブ分子間に仮想的な反発力が働いている
3. 比較的短いシミュレーションを複数回実行し、その結果を統合する

### プローブ分子

プローブ分子とは、MSMDにおける解析用分子であり、ベンゼン、イソプロパノール等のような、極めて小さな有機化合物が採用されます。このプローブ分子がタンパク質と相互作用することで、タンパク質の構造変化の促進（cryptic結合部位探索）や、タンパク質表面上の分子の結合の解析（ホットスポット解析、結合親和性予測）などを行うことが出来ます。

### プローブ分子間の仮想的な反発力

しばしば、疎水的なプローブ分子を利用することがあり、そのような場合にはプローブ分子が凝集してしまいます。
このような問題を解決するため、プローブ分子間に仮想的な反発力を導入します。

プローブ分子間の仮想反発力は、Lennard-Jonesポテンシャルを用いて表現されます。

1. 各プローブ分子の重心に仮想原子 VIS を配置する
2. 仮想原子 VIS 同士の間にLennard-Jonesポテンシャルを適用する

$U_{LJ}(r) = 4\epsilon_{LJ} \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6} \right]$

- $\epsilon_{LJ} = 10^{-6} \mathrm{kcal/mol}$
- $R_{min} = 20.0 \mathrm{\AA}$ ($\sigma = R_{min}/2^{1/6}$)

この設定は、Guvenchらの研究[^1] で示されていたプローブ分子の炭素原子間距離をよく再現するようなパラメータを手動で探索し、定めたものになっています。

[^1]: Olgun Guvench and Alexander D. MacKerell Jr. "Computational Fragment-Based Binding Site Identification by Ligand Competitive Saturation", *PLoS Computational Biology*, **5**: e1000435, 2009. DOI: [10.1371/journal.pcbi.1000435](https://doi.org/10.1371/journal.pcbi.1000435)


## MSMDシミュレーションプロトコル

MSMDでは、以下の手順のシミュレーションを複数回、系の初期構造を変化させながら行います。

### エネルギー最小化
系の中に存在している歪みを多少解消するために、最急降下法によるエネルギー最小化を行います。
この際、まずは位置拘束を加えた状態で最小化を行い、次に拘束を解除して最小化を行います。

### 系の加熱
系を目的の温度まで加熱します。

- NVTアンサンブル
- タンパク質に位置拘束
- 温度を徐々に上昇

### 位置拘束の段階的緩和
ここから先は、平衡化ステップとして、位置拘束を段階的に緩和していきます。

1. 初期の強い拘束から開始（1000 kcal/mol/Å²）
2. 段階的に拘束を弱める
3. 最終的に拘束を完全に解放

拘束の強さは以下に示すような段階で変化します。

1. 1000 kcal/mol/Å²
2. 500 kcal/mol/Å²
3. 200 kcal/mol/Å²
4. 100 kcal/mol/Å²
5. 50 kcal/mol/Å²
6. 20 kcal/mol/Å²
7. 10 kcal/mol/Å²
8. 0 kcal/mol/Å²

これらは圧力一定（NPTアンサンブル）で行われます。

### プロダクションラン
最終的なシミュレーションを実行します。

- NPTアンサンブル
- 位置拘束なし
- 一定時間のシミュレーション（標準：40 ns）


## 実装の詳細

`./exprorer_msmd` は、3つの主要な処理段階で構成されています。

1. 前処理　`preprocess()`
2. シミュレーション実行 `execute_single_simulation()`
3. 後処理 `postprocess()`

なお、今後の説明では、シミュレーションの反復IDが `42` である場合を例に説明します。

### 1. 前処理 `preprocess()`

前処理では、MSMDシミュレーションの実行に必要な系の構築を行います。最終的に出力されるデータは GROMACS で扱える形式に変換されます。
実行結果は、`PATH/TO/OUTPUT/DIR/system42/prep` ディレクトリに保存されます。

#### 1-1. プローブの準備

- プローブのmol2ファイルを参照し、`parmchk2` コマンドを利用して不足している力場パラメータを補うfrcmodファイルを生成します。

#### 1-2. タンパク質の準備

- 入力PDBファイルを読み込み、ANISOUレコード（異方性温度因子の情報）とC末端の酸素原子（OXT）を削除します。
- タンパク質立体構造から、系のサイズを決定します。
   - `tleap` を使ってタンパク質と水分子だけからなる系を作成し、その系の最長辺の長さを記録します。
   - この長さを基準とした立方体を、MSMDシミュレーションの系のサイズとします。
     - タンパク質が長細い形状である場合など、境界をまたいでタンパク質同士が干渉することを防ぐために立方体領域としました。

#### 1-3. MSMDシステムの構築

- タンパク質とプローブ分子のみからなる系を作成します。
  - 系内部に置けるプローブ分子の数を、指定されたプローブのモル濃度と、1-2で定めた系のサイズから計算します。
  - `packmol` を用いて、タンパク質とプローブ分子を、系内部に配置します。この際、プローブ分子の初期配置はシミュレーションのIDをシード値とした乱数に基づいて各シミュレーション試行で異なるようにします。
- 水分子とイオンの追加
  - `tleap` を用いて系内部に水分子を充填し、系の電荷を中和するイオンを追加します。
    - Na+イオンとCl-イオンは必要最低限追加することとし、結果としてNa+もCl-も追加されない場合もあります。
- トポロジーファイルの変換
  - Amber形式のトポロジーファイル（.parm7/.rst7）をGROMACS形式（.top/.gro）に変換します。
- 仮想原子の追加
   - プローブ分子の重心位置に仮想原子 VIS を配置します。
   - VIS の座標を gro ファイルに追加します
   - VIS の位置を固定するために、 top ファイルに `[virtual_sitesn]` を設定します。
   - VIS 間のLennard-Jonesポテンシャルを適用するため、 top ファイルにVIS-VIS間の非共有結合相互作用を追記します。
- 位置拘束の設定
   - タンパク質の重原子に対して、8段階の位置高速を設定できるように top ファイルに追記を行います。

### 2. シミュレーション実行 `execute_single_simulation()`

前処理で作成された系を用いて、MSMDシミュレーションを実行します。シミュレーションエンジンとしては GROMACS を使用します。
実行結果は、`PATH/TO/OUTPUT/DIR/system42/simulation` ディレクトリに保存されます。

- インデックスファイルの生成
   - `gmx make_ndx` によって原子グループを定義します。exprorer_msmd では、 `Water` と `Non-Water` の2つのグループを利用します。
- mdp ファイルの作成
   - YAMLの設定に基づいて、各ステップ（最小化、加熱、平衡化、プロダクション）に対応する mdp ファイルを自動生成します。
      - 各ステップのテンプレートは　`./script/template/production.mdp` などに存在しています。
- シミュレーションの実行
   - 各シミュレーションステップを直列に実行する `mdrun.sh` を自動で作成し、実行します。
   - 最終ステップ（通常は pr (production run)）で得られたトラジェクトリファイルは、 `[プロジェクト名].xtc` という名前のシンボリックリンクが作成されます。

### 3. 後処理 `postprocess()`

このステップは [PMAPの詳細](pmap.md) に詳細が記載されていますので、そちらを参照してください。