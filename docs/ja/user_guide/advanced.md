# 高度な使用方法

## シミュレーション設定のカスタマイズ

### シミュレーションステップの制御

特定のフェーズをスキップするオプションが用意されています。

```bash
# 前処理をスキップ
./exprorer_msmd protocol.yaml --skip-preprocess

# シミュレーションをスキップ
./exprorer_msmd protocol.yaml --skip-simulation

# 後処理をスキップ
./exprorer_msmd protocol.yaml --skip-postprocess
```

### 独立試行の設定

```yaml
general:
  iter_index: 0,1,2  # 3回の独立試行を実行
  # 表記法:
  # "1-3" => 1,2,3
  # "5-9:2" => 5,7,9
  # "1-3,5-9:2" => 1,2,3,5,7,9
```

### シミュレーション条件の調整

```yaml
exprorer_msmd:
  general:
    dt: 0.002          # 時間刻み幅（ps）
    temperature: 300   # 温度（K）
    pressure: 1.0      # 圧力（bar）
    pbc: xyz          # 周期境界条件

  sequence:
    - name: pr        # プロダクションラン
      type: production
      nsteps: 20000000  # ステップ数（40 ns）
      nstxtcout: 5000   # 出力頻度（10 ps）
```

## 解析設定のカスタマイズ

### PMAPの計算設定

```yaml
map:
  type: pmap
  snapshot: 2001-4001:1  # 使用するスナップショット範囲
  maps:
    - suffix: nVH        # 重原子のみ
      selector: (!@VIS)&(!@H*)
    - suffix: nV         # 全原子
      selector: (!@VIS)
  map_size: 80           # マップサイズ（Å）
  normalization: total   # 正規化方法, snapshotおよびGFEが指定可能
```

### Inverse MSMD 関連の設定

```yaml
probe_profile:
  resenv:
    map: nVH            # 使用するマップ
    threshold: 0.001    # 確率閾値
  profile:
    types:
      - name: anion # マップの名称
        atoms:
          - ["ASP", " CB "]
          - ["GLU", " CB "]
```