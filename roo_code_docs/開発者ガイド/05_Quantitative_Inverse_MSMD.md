# Quantitative Inverse MSMD

Quantitative Inverse MSMDモジュールは、分子の特定の原子を置換した場合の結合親和性の変化を見積もるための機能を提供します。このモジュールは、創薬研究において重要な構造活性相関（SAR）の理解や、リード最適化に役立ちます。

## モジュール構造と責任

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

主要なモジュールの責任は以下の通りです：

- **quantitative_inverse_msmd.py**：Quantitative Inverse MSMD全体を管理するメインスクリプト
  - 非共通部分の抽出
  - プローブの設計
  - 分子の重ね合わせ
  - MSMDシミュレーションの設定
  - 結合親和性の変化の見積もり

- **utilities/analysis/profile.py**：プロファイル関連のクラスと関数
  - PreProfileクラス：各残基に対するbulk_probasumとprobasum_gridを持つ
  - Profileクラス：各残基に対する、bulkとの存在確率比を保持する
  - プロファイル生成関連の関数

- **utilities/analysis/matching_score.py**：マッチングスコア計算関連の関数
  - preprocessing：トラジェクトリの前処理
  - calculate_matching_score：マッチングスコアの計算

- **utilities/analysis/compare_results.py**：結果比較関連の関数
  - parse_msmd_result：MSMDシミュレーション結果の解析
  - compare_results：2つのMSMDシミュレーション結果の比較

- **utilities/molecular/mcs_extractor.py**：非共通部分抽出関連の関数
  - extract_non_common_parts：2つの分子の非共通部分を抽出

- **utilities/molecular/molecule_superimposer.py**：分子重ね合わせ関連の関数
  - generate_atom_mapping：原子対応関係の生成
  - superimpose_molecules：分子の重ね合わせ

- **utilities/molecular/simple_atom_replacer.py**：原子置換関連の関数
  - replace_atom_in_molecule：分子の特定の原子を別の原子に置換

- **utilities/molecular/extract_substructure.py**：部分構造抽出関連の関数
  - extract_substructure_around_atom：指定した原子の周りの部分構造を抽出

- **utilities/probe/probe_designer.py**：プローブ設計関連の関数
  - design_probe：非共通部分に基づいてプローブを設計

## 主要な関数と機能

```python
def semi_automatic_quantitative_inverse_msmd(
    compound1_file, compound2_file, protein_file, probe_library_dir, output_dir,
    only_non_common_rings=False, separate_rings=False
):
    """
    半自動Quantitative Inverse MSMDを実行する

    Parameters
    ----------
    compound1_file : str
        1つ目の化合物のSDFファイルパス
    compound2_file : str
        2つ目の化合物のSDFファイルパス
    protein_file : str
        タンパク質構造のPDBファイルパス
    probe_library_dir : str
        プローブライブラリのディレクトリパス
    output_dir : str
        出力ディレクトリ
    only_non_common_rings : bool, optional
        Trueの場合、複合環のうち非共通な環のみをプローブ化する, by default False
    separate_rings : bool, optional
        Trueの場合、非共通な環を別々のプローブとして抽出する, by default False

    Returns
    -------
    dict
        結果情報
    """
    # 出力ディレクトリを作成
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. 2つの化合物間の非共通部分を抽出
    print("ステップ1: 非共通部分の抽出")
    non_common_parts = extract_non_common_parts(
        compound1_file, compound2_file, only_non_common_rings, separate_rings
    )
    
    # 2. 各非共通部分に対してプローブを設計
    print("ステップ2: プローブの設計")
    probe1 = design_probe(non_common_parts['mol1'], probe_library_dir)
    probe2 = design_probe(non_common_parts['mol2'], probe_library_dir)
    
    # 3. プローブをファイルに保存
    probe1_file = os.path.join(output_dir, "probe1.sdf")
    probe2_file = os.path.join(output_dir, "probe2.sdf")
    
    writer1 = Chem.SDWriter(probe1_file)
    writer1.write(probe1)
    writer1.close()
    
    writer2 = Chem.SDWriter(probe2_file)
    writer2.write(probe2)
    writer2.close()
    
    # 4. ユーザーによる量子化学計算（手動ステップ）
    print("ステップ3: プローブの量子化学計算")
    print("  プローブ1とプローブ2の量子化学計算を実行してください")
    print("  計算が完了したら、任意のキーを押して続行してください")
    input()
    
    # 5. 化合物とプローブの重ね合わせ
    print("ステップ4: 分子の重ね合わせ")
    
    # 化合物1とプローブ1の原子対応関係を自動生成
    print("  化合物1とプローブ1の原子対応関係を自動生成中...")
    atom_mapping1 = generate_atom_mapping(compound1_file, probe1_file)
    atom_mapping1 = validate_and_refine_atom_mapping(compound1_file, probe1_file, atom_mapping1)
    
    # 化合物1とプローブ1の重ね合わせ
    transformed_probe1 = superimpose_molecules(compound1_file, probe1_file, atom_mapping1)
    transformed_probe1_file = os.path.join(output_dir, "transformed_probe1.sdf")
    
    writer1 = Chem.SDWriter(transformed_probe1_file)
    writer1.write(transformed_probe1)
    writer1.close()
    
    # 化合物2とプローブ2の原子対応関係を自動生成
    print("  化合物2とプローブ2の原子対応関係を自動生成中...")
    atom_mapping2 = generate_atom_mapping(compound2_file, probe2_file)
    atom_mapping2 = validate_and_refine_atom_mapping(compound2_file, probe2_file, atom_mapping2)
    
    # 化合物2とプローブ2の重ね合わせ
    transformed_probe2 = superimpose_molecules(compound2_file, probe2_file, atom_mapping2)
    transformed_probe2_file = os.path.join(output_dir, "transformed_probe2.sdf")
    
    writer2 = Chem.SDWriter(transformed_probe2_file)
    writer2.write(transformed_probe2)
    writer2.close()
    
    # 6. inverse MSMDシミュレーションの実行
    print("ステップ5: inverse MSMDシミュレーションの実行")
    
    # プローブ1のinverse MSMDシミュレーション
    print("  プローブ1のinverse MSMDシミュレーションを実行中...")
    probe1_settings = setup_msmd_simulation(transformed_probe1_file, protein_file, output_dir)
    probe1_simulation_dir = run_msmd_simulation(probe1_settings)
    
    # プローブ2のinverse MSMDシミュレーション
    print("  プローブ2のinverse MSMDシミュレーションを実行中...")
    probe2_settings = setup_msmd_simulation(transformed_probe2_file, protein_file, output_dir)
    probe2_simulation_dir = run_msmd_simulation(probe2_settings)
    
    # 7. プロファイル作成
    print("ステップ6: プロファイルの作成")
    
    # プローブ1のプロファイル作成
    probe1_residue_env = extract_residue_environment(probe1_simulation_dir)
    profile1 = create_odds_ratio_profile(probe1_residue_env)
    
    # プローブ2のプロファイル作成
    probe2_residue_env = extract_residue_environment(probe2_simulation_dir)
    profile2 = create_odds_ratio_profile(probe2_residue_env)
    
    # 8. 合致度計算
    print("ステップ7: 合致度計算と結合強度推定")
    
    # 各プロファイルの合致度を計算
    score1 = calculate_matching_score(protein_file, profile1)
    score2 = calculate_matching_score(protein_file, profile2)
    
    # 9. 結合強度の推定
    strength_diff = estimate_binding_strength(score1, score2)
    
    # 10. 結果の出力
    print("ステップ8: 結果出力")
    print(f"  化合物1の合致度スコア: {score1}")
    print(f"  化合物2の合致度スコア: {score2}")
    print(f"  結合強度差: {strength_diff}")
    
    if strength_diff > 0:
        print("  結論: 化合物1の方が強く結合すると予測されます")
    elif strength_diff < 0:
        print("  結論: 化合物2の方が強く結合すると予測されます")
    else:
        print("  結論: 両化合物の結合強度は同程度と予測されます")
    
    # 結果をファイルに保存
    result_file = os.path.join(output_dir, "quantitative_inverse_msmd_result.txt")
    with open(result_file, "w") as f:
        f.write(f"化合物1: {compound1_file}\n")
        f.write(f"化合物2: {compound2_file}\n")
        f.write(f"タンパク質: {protein_file}\n")
        f.write(f"プローブ1: {transformed_probe1_file}\n")
        f.write(f"プローブ2: {transformed_probe2_file}\n")
        f.write(f"化合物1の合致度スコア: {score1}\n")
        f.write(f"化合物2の合致度スコア: {score2}\n")
        f.write(f"結合強度差: {strength_diff}\n")
        
        if strength_diff > 0:
            f.write("結論: 化合物1の方が強く結合すると予測されます\n")
        elif strength_diff < 0:
            f.write("結論: 化合物2の方が強く結合すると予測されます\n")
        else:
            f.write("結論: 両化合物の結合強度は同程度と予測されます\n")
    
    print(f"結果を {result_file} に保存しました")
    
    return {
        'probe1': probe1,
        'probe2': probe2,
        'profile1': profile1,
        'profile2': profile2,
        'score1': score1,
        'score2': score2,
        'strength_diff': strength_diff
    }
```

## 他モジュールとの連携ポイント

Quantitative Inverse MSMDモジュールは、以下のモジュールと連携しています：

- **MCS計算と非共通部分抽出モジュール**：`script/utilities/molecular/mcs_extractor.py`で実装されています。
- **プローブ設計モジュール**：`script/utilities/probe/probe_designer.py`で実装されています。
- **分子の重ね合わせ自動化モジュール**：`script/utilities/molecular/molecule_superimposer.py`で実装されています。
- **MSMDシミュレーションエンジン**：`script/quantitative_inverse_msmd.py`内の関数で設定されています。
- **プロファイル解析モジュール**：`script/utilities/analysis/profile.py`で実装されています。
- **マッチングスコア計算モジュール**：`script/utilities/analysis/matching_score.py`で実装されています。
- **結果比較モジュール**：`script/utilities/analysis/compare_results.py`で実装されています。

## Quantitative Inverse MSMDの処理フロー

Quantitative Inverse MSMDの処理フローは以下の通りです：

1. **非共通部分の抽出**：
   - 2つの化合物のSDFファイルを読み込みます。
   - 最大共通部分構造（MCS）を計算します。
   - MCSに基づいて非共通部分を抽出します。

2. **プローブの設計**：
   - 非共通部分に基づいてプローブを設計します。
   - プローブの化学的性質（極性、電荷、サイズなど）を調整します。

3. **分子の重ね合わせ**：
   - 元の化合物とプローブの原子対応関係を生成します。
   - Kabsch algorithmを使用して分子を重ね合わせます。

4. **MSMDシミュレーションの設定**：
   - プローブとタンパク質を使用してMSMDシミュレーションを設定します。
   - シミュレーションパラメータを最適化します。

5. **残基環境の抽出**：
   - MSMDシミュレーションの結果からプローブ周辺の残基環境を抽出します。
   - 残基環境からプロファイルを生成します。

6. **オッズ比プロファイルの作成**：
   - 残基環境のオッズ比プロファイルを作成します。
   - オッズ比プロファイルを正規化します。

7. **合致度の計算**：
   - プロファイルとタンパク質構造の合致度を計算します。
   - 合致度スコアを正規化します。

8. **結合強度の推定**：
   - 合致度スコアの差から結合強度の差を推定します。
   - 推定結果を検証します。

## MCS計算と非共通部分抽出

MCS（Maximum Common Substructure）計算と非共通部分抽出は、Quantitative Inverse MSMDの重要なステップです。以下のようなアルゴリズムを使用しています：

1. **MCSの計算**：
   - RDKitのFMCSアルゴリズムを使用してMCSを計算します。
   - MCSの計算パラメータ（原子タイプの一致条件、結合タイプの一致条件など）を最適化します。

2. **非共通部分の抽出**：
   - MCSに基づいて、2つの分子の非共通部分を抽出します。
   - 非共通部分の化学的妥当性をチェックします。

3. **非共通部分の最適化**：
   - 非共通部分の構造を最適化して、化学的に安定な構造にします。
   - 非共通部分の電荷や水素原子の状態を調整します。

## プローブ設計

プローブ設計は、非共通部分に基づいてプローブ分子を設計するステップです。以下のような方法を使用しています：

1. **プローブテンプレートの選択**：
   - 非共通部分の性質に合わせてプローブテンプレートを選択します。
   - テンプレートには、メチルベンゼン、ジメチルベンゼン、トリメチルベンゼンなどがあります。

2. **プローブの修飾**：
   - プローブテンプレートに非共通部分を組み込みます。
   - プローブの化学的性質を調整します。

3. **プローブの最適化**：
   - プローブの構造を最適化して、化学的に安定な構造にします。
   - プローブの電荷や水素原子の状態を調整します。

## 分子の重ね合わせ

分子の重ね合わせは、元の化合物とプローブを重ね合わせるステップです。以下のような方法を使用しています：

1. **原子対応関係の生成**：
   - 元の化合物とプローブの原子対応関係を生成します。
   - 対応関係の生成には、原子タイプ、結合情報、環構造などを考慮します。

2. **Kabsch algorithmによる重ね合わせ**：
   - 原子対応関係に基づいて、Kabsch algorithmを使用して分子を重ね合わせます。
   - 重ね合わせの精度を評価して、必要に応じて調整します。

3. **重ね合わせ結果の最適化**：
   - 重ね合わせ結果を最適化して、より良い重ね合わせを得ます。
   - 最適化には、原子間距離の最小化、エネルギー最小化などを使用します。

## 結合親和性の変化の見積もり

結合親和性の変化の見積もりは、MSMDシミュレーションの結果から結合親和性の変化を推定するステップです。以下のような方法を使用しています：

1. **残基環境プロファイルの比較**：
   - 2つのプローブの残基環境プロファイルを比較します。
   - プロファイルの差異を定量化します。

2. **合致度スコアの計算**：
   - プロファイルとタンパク質構造の合致度を計算します。
   - 合致度スコアを正規化します。

3. **結合強度の差の推定**：
   - 合致度スコアの差から結合強度の差を推定します。
   - 推定結果を検証します。

## 実装上の注意点

### パフォーマンス最適化

Quantitative Inverse MSMDは計算コストが高いため、以下のような最適化が行われています：

1. **MCS計算の最適化**：MCS計算のパラメータを最適化して、計算を高速化しています。
2. **並列処理**：複数の分子ペアの解析を並列化して、処理を高速化しています。
3. **キャッシュ**：頻繁に使用するデータをキャッシュして、計算を高速化しています。

### 精度と感度のバランス

Quantitative Inverse MSMDでは、精度と感度のバランスが重要です。以下のような工夫が行われています：

1. **MCS計算パラメータの最適化**：MCS計算のパラメータを最適化して、適切なMCSを得ています。
2. **プローブ設計の最適化**：プローブ設計のパラメータを最適化して、適切なプローブを得ています。
3. **結合親和性推定モデルの最適化**：結合親和性推定モデルのパラメータを最適化して、精度の高い推定を行っています。

### 拡張性

Quantitative Inverse MSMDモジュールは、以下のような拡張性を持っています：

1. **新しいMCSアルゴリズムの追加**：新しいMCS計算アルゴリズムを簡単に追加できるように設計されています。
2. **新しいプローブテンプレートの追加**：新しいプローブテンプレートを簡単に追加できるように設計されています。
3. **新しい結合親和性推定モデルの追加**：新しい結合親和性推定モデルを簡単に追加できるように設計されています。

## 使用例

以下は、Quantitative Inverse MSMDモジュールの使用例です：

```python
# Quantitative Inverse MSMDの実行
from script.quantitative_inverse_msmd import semi_automatic_quantitative_inverse_msmd

# 実行
result = semi_automatic_quantitative_inverse_msmd(
    'path/to/compound1.sdf',
    'path/to/compound2.sdf',
    'path/to/protein.pdb',
    'path/to/probe_library',
    'output/qimsmd',
    only_non_common_rings=True,
    separate_rings=False
)

# 結果の表示
print("結合親和性の変化:")
print(f"  化合物1の合致度スコア: {result['score1']}")
print(f"  化合物2の合致度スコア: {result['score2']}")
print(f"  結合強度差: {result['strength_diff']}")

# 結果の可視化
# （省略）
```

## 応用例

Quantitative Inverse MSMDモジュールは、以下のような応用例があります：

1. **構造活性相関（SAR）の解析**：化合物の構造変化と活性の関係を解析します。
2. **リード最適化**：リード化合物の最適化に役立つ情報を提供します。
3. **ファーマコフォアモデルの構築**：活性に重要な構造的特徴を特定します。
4. **新規化合物の設計**：結合親和性の高い新規化合物を設計します。