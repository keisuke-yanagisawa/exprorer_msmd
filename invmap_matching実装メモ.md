1. probe_profile の計算をCa原子だけに限定する

```py
create_residue_interaction_profile(pdb, [" CA "], resname)
```
※そもそもalignする前にCaだけが存在する状態にトラジェクトリを変換しておくべき

2. probe_profile の計算を全probe moleculesに対して行うようにする 
  距離を極めて大きくすればOKな話

```py
resenv(...., threshold=0) # 閾値を無くす
```

3. 全空間に均等に残基Caがあるとしたときの存在確率を計算

タンパク質の全Ca原子を取るように修正
```py
resenv(..., env_distance=10000) # 十分に遠い距離を指定する
```

続いて、全空間の体積を取得
```py
probe_volume = uPDB.estimate_exclute_volume(probe_molecules_obj)
volume = cpptraj_obj.last_volume - probe_volume
probability = cnt_residues(protein, "Arg") / volume * voxel_size # probability in A^3
```

4. probe_profile の出力値をバルク中の存在確率 probability で割り算することでオッズ比に変換

```py
ret = resenv(...)
ret /= probability # odds ratioに変換
```

5. ある中心・向きを指定されたときのprobe interaction profileとタンパク質構造との類似度評価を計算
  範囲内のタンパク質残基が存在する場所のオッズ比を全て持ってきて、積を取る


