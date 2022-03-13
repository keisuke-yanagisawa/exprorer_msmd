# EXPRORER_MSMD

共溶媒分子動力学 (mixed-solvent molecular dynamics; MSMD) シミュレーションエンジンと、解析ツールをまとめたリポジトリ

## システムの説明

GROMACSを用いたMSMDを自動的に行うためのシステム。

## 簡単な使い方

### 環境構築
このgithubにはDockerファイルが含まれており、dockerを用いることで実行可能な環境を容易に作ることができる。

### MSMDシミュレーションの実施 ( `exprorer_msmd` )

タンパク質、共溶媒（プローブ）分子、シミュレーションプロトコルを定義したyamlファイルを用意し（例えば `example/msmd_protocol.yaml` ）、
以下のコマンドを実行することで系の構築からシミュレーションの実施、確率密度マップ (spatial probability distribution map; PMAP) の作成を自動的に行う。

```
./exprorer_msmd example/msmd_protocol.yaml
```

### MSMDシミュレーション結果の解析

現在、2つの解析手法が準備されている。

#### タンパク質のホットスポットを探索する ( `protein_hotspot` )

タンパク質表面のどの部分にプローブ分子が存在しやすいか？という
ホットスポット探索を行う。

```
./protein_hotspot example/msmd_protocol.yaml
```

実行結果は以下のようなものになる。
![ホットスポット探索](https://i.imgur.com/bzxz0K6.png)

EXPRORER[^1] はこの結果を用いて手法を提案している。

[^1]:**Keisuke Yanagisawa**, Yoshitaka Moriwaki, Tohru Terada, Kentaro Shimizu. "EXPRORER: Rational Cosolvent Set Construction Method for Cosolvent Molecular Dynamics Using Large-Scale Computation", *Journal of Chemical Information and Modeling*, **61**: 2744-2753, 2021/06. DOI: [10.1021/acs.jcim.1c00134](https://doi.org/10.1021/acs.jcim.1c00134)

#### プローブ分子周辺残基環境を取得する ( `probe_profile` )

前述のホットスポット解析とは反対に、プローブ分子の周囲にどのような残基が存在しやすいか？を描画する。

```
./protein_hotspot example/msmd_protocol.yaml
```

実行結果は以下のようなものになる。
![プローブ分子周辺残基環境](https://i.imgur.com/4QIZxhW.png)

## yamlファイルの書き換え

すべての設定は yaml ファイルに記述されている。