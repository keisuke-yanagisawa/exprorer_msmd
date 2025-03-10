# プローブ分子の準備

プローブ分子の準備は、以下の手順に沿って行います。

1. 初期立体構造生成
2. プローブ分子の構造最適化
3. 電子密度計算
4. RESP法による部分電荷の割り当て

なお、プローブ分子は、**内部自由度を持たない** 事が望ましいです。特に、タンパク質表面におけるプローブの存在確率を計算する場合、プローブ分子自身の配座が大きな影響を及ぼします。プローブ分子に内部自由度を持たせないことで、プローブ分子の配座にかかるサンプリングの効率低下を防ぐことが出来ます。

もし内部自由度を持つプローブ分子を使用する場合には、MSMDシミュレーションの総量を増加させつつ、結果が十分に収束しているかを確認するようにしましょう。

例として登録されている `./example/A11.mol2` および `./example/A11.pdb` は、内部自由度を持たない isopropanol プローブになっています。

## 1. 初期立体構造生成

RDKitやOpenBabelなどのツールを使用して、プローブ分子の初期立体構造を生成します。
この初期構造を起点に量子化学計算を行うため、生成された初期立体構造が想定通りになっていることを良く確認してください。（例：キラル中心を持つプローブ分子、内部自由度を持つプローブ分子）

```sh
# OpenBabelを使用して、isopropanolの初期立体構造を生成する例
#（あらかじめisopropanolのSMILES表記を用意しておく）
obabel -ismi isopropanol.smi -omol2 -Oisopropanol.mol2 --gen3D
```

## 2. プローブ分子の構造最適化

続いて、量子化学計算を用いてプローブ分子の構造最適化を行います。この構造最適化が不十分だと、次のステップで行う部分電荷割り当てに影響を及ぼす可能性があります。

量子化学計算を行うには、Gaussian や GAMESS などの計算ソフトウェアを存在しますが、Gaussian 16 を使用した場合の例を以下に示します。

先ず、Gaussian 16 の入力ファイルを準備します。exprorer_msmdでは、計算レベルは B3LYP/6-31G(d) を利用しています。

`./example/A11_opt.gjf`
```
%chk=opt.chk
%mem=10GB
# p opt=(maxcycle=999, maxstep=5, tight) b3lyp/6-31g(d) scf=(qc,tight)
opt

0 1
0  1
C           1.06720        -0.04599        -0.01121 # ここから先は各原子の座標情報
C           2.58785        -0.07121        -0.02775
...
```

この入力ファイルを Gaussian 16 に入力します。

```sh
g16 < A11_opt.gjf > A11_opt.log
```

このようにすると、実行ログが `A11_opt.log` に出力されます。この中で、`Normal termination` という文字列が存在していれば、計算が正常に終了しています。

## 3. 電子密度計算

次に、最適化された構造について、電子密度を計算します。電子密度の計算には、Gaussian 16 の `pop=mk` オプションを使用します。この際には、B3LYP/6-31g(d) レベル **ではなく** 、 HF/6-31g(d) レベルでの計算を行います。


`./example/A11_elec.gjf`
```
%chk=hf.chk
%mem=10GB
# p HF/6-31g(d) iop(6/33=2,6/41=10,6/42=17) pop=(mk,readradii) scf=tight
hf

0 1
0  1
C           1.06720        -0.04599        -0.01121 # ここから先は各原子の座標情報
C           2.58785        -0.07121        -0.02775 # 構造最適化後の座標を記述する（この例は最適化前の座標なので注意）
...
```

```sh
g16 < A11_elec.gjf > A11_elec.log
```

これで、電子密度が計算され、その結果が log ファイルに記載されます。

## 4. RESP法による部分電荷の割り当て

最後に、 RESP 法を用いて部分電荷を割り当てます。この際には、`antechamber` を使用します。

```sh
# mol2 と pdb ファイルをそれぞれ出力するために、 antechamber を2回実行します。
# -nc は分子全体の電荷を指定します。Iso-propanol は中性分子なので 0 です。
# -rn は残基名を指定します。他の分子と名前が衝突しないようにしましょう。
antechamber \
  -fi gout -i A11_elec.log \
  -fo mol2 -o A11_resp.mol2 \
  -c resp \
  -nc 0 \
  -at gaff2 \
  -pf y \
  -rn A11

antechamber \
  -fi gout -i A11_elec.log \
  -fo pdb -o A11_resp.pdb \
  -c resp \
  -nc 0 \
  -at gaff2 \
  -pf y \
  -rn A11
```

なお、`antechamber` で mol2 ファイルを出力させる際に、部分電荷の合計値が0にならない事があります。MSMDシミュレーション実行時にエラーが発生しますので、その場合には合計値が0になるように手動で修正してください。