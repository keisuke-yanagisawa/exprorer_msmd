# PMAP（確率密度マップ）

PMAP (spatial Probability distribution MAP) は、プローブ分子の空間分布を表現する3次元の確率密度マップです。

## 目次
1. [PMAPの基本](#pmapの基本)
2. [計算方法](#計算方法)
   - [基本的な手順](#1-基本的な手順)
   - [グリッドの設定](#2-グリッドの設定)
   - [原子の選択](#3-原子の選択)
   - [数値の変換](#4-数値の変換)
3. [実装の詳細](#実装の詳細)
   - [単一シミュレーション結果の処理](#1-単一シミュレーション結果に対する処理gen_pmap-関数)
   - [複数シミュレーション結果の統合](#2-複数シミュレーション結果の統合)

関連ドキュメント：
- [MSMDの実装](msmd.md)
- [基本的な使い方](../user_guide/basic.md)
- [高度な使用方法](../user_guide/advanced.md)

## PMAPの基本

PMAPは、共溶媒分子動力学 (mixed-solvent molecular dynamics; MSMD) シミュレーションにおいて観測されたプローブ分子の空間的存在確率を表したものです。

この情報は、定性的および定量的に活用することが可能です。

- プローブが存在しやすい領域は化合物が結合しやすいだろう（定性的、ホットスポット解析）
- プローブの存在確率が高い領域はより強く結合しているだろう（定量的、結合親和性予測）

## 計算方法

### 1. 基本的な手順

1. シミュレーションボックスを3次元グリッドに分割
2. 各グリッドポイントでのプローブ原子の存在頻度をカウント
3. 頻度を確率に変換

### 2. グリッドの設定

グリッドの作成は Amber Tools の cpptraj を使用して行っています。

- グリッド中心：タンパク質の重心（この座標を原点とする）
- デフォルトのボックスサイズ：80 Å × 80 Å × 80 Å
- グリッド間隔：1 Å

これにより、1 Å の間隔で 80 × 80 × 80 のグリッドが作成され、
各グリッドポイントにおけるプローブ原子の存在 **頻度** をカウントします。

### 3. 原子の選択

グリッドを作成する際に用いる原子は、yamlファイルで指定されたセレクタによって選択されます。
以下に3つの例を示します。

1. 重原子のみ
   - セレクタ：`(!@VIS)&(!@H*)`
   - 水素原子と仮想原子を除外
2. 全原子
   - セレクタ：`(!@VIS)`
   - 仮想原子のみを除外
3. プローブ重心
   - セレクタ：`@VIS`
   - プローブの仮想原子のみを選択
     - プローブの仮想原子はそのプローブの重心に配置されているので、仮想原子のみを選択することでプローブの重心を選択出来ます。

### 4. 数値の変換

前述の通り、ここまでに作成したグリッドはプローブ原子の存在 **頻度** をカウントしたものです。
当然、シミュレーション長（正確にはスナップショット数）によってカウントされる数値が変わってしまうため、
これを **確率** に変換する必要があります。この確率化は、以下の2種類が提供されています。
デフォルトは **total 正規化** です（が、snapshot 正規化の方が適切な場合も多いです）。

1. total 正規化
   $P(r) = N(r) / \sum N(r)$

   - $N(r)$：位置 $r \in \mathbb{R}^3$ での出現頻度
   - 全グリッドポイントの確率の総和が1

2. snapshot 正規化
   $P(r) = N(r) / N_{\mathrm{frames}}$

   - $N_{\mathrm{frames}}$：使用したスナップショット数
   - 各位置での存在確率を直接表現

また、この確率を使って、タンパク質とプローブ分子との相互作用エネルギー (grid free energy; GFE) を見積もる手法も提供されています。

$$\mathrm{GFE}(r) = -RT \ln(P(r)/P_{\mathrm{bulk}})$$

- $R$：気体定数
- $T$：温度
- $P_{\mathrm{bulk}}$：バルクの確率（全体の平均確率）

## 実装の詳細

### 1. 単一シミュレーション結果に対する処理（`gen_pmap()` 関数）

#### cpptrajによる頻度グリッド計算

cpptrajを使用したグリッド計算は以下の手順で実行されます。
以下の実行を起こすテンプレートは　`script/utilities/executables/template/cpptraj_pmap.in` に存在しますので、合わせて参照してください。

1. トラジェクトリの前処理
   トラジェクトリを読み込み、主鎖炭素原子を使ってシミュレーション前の初期構造との重ね合わせを実施します。
   ```
   trajin [プロジェクト名].xtc [開始フレーム] [終了フレーム]
   reference [初期構造].pdb
   rms reference @CA
   ```

2. 頻度グリッドの作成
   cpptrajのgridコマンドを用いて、指定された原子セレクタに基づいてグリッドを作成します。

   ```
   grid [保存するマップファイルパス].dx 
     [X方向サイズ(Å)] [X方向間隔(Å)] [Y方向サイズ(Å)] [Y方向間隔(Å)] [Z方向サイズ(Å)] [Z方向間隔(Å)] 
     gridcenter 0.0 0.0 0.0 [原子セレクタ]
   ```
   - デフォルトサイズ：80×80×80 Å
   - デフォルト間隔：1×1×1 Å

#### PMAPファイルの生成

cpptrajで作成した頻度グリッドを読み込み、PMAPを生成します。

1. 確率への変換
   - total正規化の場合：
     ```python
     total_count = np.sum(frequency_map)
     probability_map = frequency_map / total_count
     ```
   - snapshot正規化の場合：
     ```python
     probability_map = frequency_map / n_frames
     ```

2. GFEの計算（オプション） `convert_to_gfe()`
   ```python
   R = 0.001987  # kcal/mol/K
   T = 300       # K
   bulk_prob = np.mean(probability_map)
   gfe_map = -R * T * np.log(probability_map / bulk_prob)
   ```

   このとき、 `gfe_map` は無限大の値を持つ可能性がある（probability_mapの値が0であるような座標はそうなる）ので、
   `gfe_map` の各値が 3.0 kcal/mol 以上であれば、3.0 kcal/mol にクリップします。

### 2. 複数シミュレーション結果の統合

under construction