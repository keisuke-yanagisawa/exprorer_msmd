# 基本的な使い方

## はじめに

EXPRORER_MSMDは、共溶媒分子動力学（MSMD）シミュレーションを実行し、解析するためのツールです。

## 環境構築

このリポジトリにはDockerファイルが含まれており、以下のツールを含む実行環境を簡単に構築できます。

- Python 3
- AmberTools 20
- Gromacs 2021.5
- Packmol 18.169

## 基本的な使用手順

### 1. プロトコルファイルの準備

シミュレーションの設定は、YAMLファイルで定義します。

```yaml
general:
  name: TEST_PROJECT         # プロジェクト名
  workdir: ./PATH/TO/WORKDIR # 出力ディレクトリ

input:
  protein:
    pdb: protein.pdb       # タンパク質構造
  probe:
    cid: PROBE            # プローブ分子名
    molar: 0.25          # 濃度（mol/L）
```

YAMLファイルの記述については、 `example/example_protocol.yaml` も参照してください。

### 2. シミュレーションの実行

```bash
./exprorer_msmd protocol.yaml
```

これを行うことで、系の構築、シミュレーション、PMAPの作成が自動的に行われます。

### 3. 結果の解析

#### トラジェクトリの確認

シミュレーション結果は自動的にPDBファイルに変換され、
`./PATH/TO/WORKDIR/system*/[プロジェクト名]_woWAT_10ps.pdb` として保存されています。
なお、 `*` の部分は、実行時のシミュレーションIDであり、複数の独立試行がある場合にはそれぞれのIDが付与されます。

このトラジェクトリは、PyMOLなどで開くことで、シミュレーション中のタンパク質とプローブの動きを確認できます。

#### ホットスポット解析
```bash
./protein_hotspot protocol.yaml
```

Under construction

#### プローブ環境解析
```bash
./probe_profile protocol.yaml
```

Under construction
