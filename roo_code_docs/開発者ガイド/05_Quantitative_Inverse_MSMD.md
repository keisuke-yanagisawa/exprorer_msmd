# Quantitative Inverse MSMD

Quantitative Inverse MSMDモジュールは、分子の特定の原子を置換した場合の結合親和性の変化を見積もるための機能を提供します。このモジュールは、創薬研究において重要な構造活性相関（SAR）の理解や、リード最適化に役立ちます。

## クラス構造と責任

Quantitative Inverse MSMDモジュールは、以下の主要なクラスから構成されています：

- **QuantitativeInverseMSMD**：Quantitative Inverse MSMD全体を管理するクラス
  - 非共通部分の抽出
  - プローブの設計
  - 分子の重ね合わせ
  - MSMDシミュレーションの設定
  - 結合親和性の変化の見積もり

- **BindingAffinityEstimator**：結合親和性の見積もりを担当するクラス
  - プロファイルとタンパク質構造の合致度を計算
  - 合致度スコアの差から結合強度の差を推定

## 主要なメソッドと機能

```python
class QuantitativeInverseMSMD:
    def __init__(self, config):
        """
        Quantitative Inverse MSMDクラスを初期化する

        Parameters
        ----------
        config : dict
            設定情報
        """
        self.config = config
        self.mcs_extractor = MCSExtractor()
        self.probe_designer = ProbeDesigner()
        self.molecule_superimposer = MoleculeSuperimposer()
        self.binding_affinity_estimator = BindingAffinityEstimator()

    def run(self, compound1_file, compound2_file, protein_file):
        """
        Quantitative Inverse MSMDを実行する

        Parameters
        ----------
        compound1_file : str
            1つ目の化合物のSDFファイルパス
        compound2_file : str
            2つ目の化合物のSDFファイルパス
        protein_file : str
            タンパク質構造のPDBファイルパス

        Returns
        -------
        dict
            結果情報
        """
        # 非共通部分の抽出
        non_common_parts = self.mcs_extractor.extract_non_common_parts(
            compound1_file, compound2_file
        )
        
        # プローブの設計
        probe1 = self.probe_designer.design_probe(non_common_parts['mol1'])
        probe2 = self.probe_designer.design_probe(non_common_parts['mol2'])
        
        # 分子の重ね合わせ
        atom_mapping1 = self.molecule_superimposer.generate_atom_mapping(
            compound1_file, probe1
        )
        transformed_probe1 = self.molecule_superimposer.superimpose_molecules(
            compound1_file, probe1, atom_mapping1
        )
        
        atom_mapping2 = self.molecule_superimposer.generate_atom_mapping(
            compound2_file, probe2
        )
        transformed_probe2 = self.molecule_superimposer.superimpose_molecules(
            compound2_file, probe2, atom_mapping2
        )
        
        # MSMDシミュレーションの設定
        msmd_settings1 = self._setup_msmd_simulation(
            transformed_probe1, protein_file
        )
        msmd_settings2 = self._setup_msmd_simulation(
            transformed_probe2, protein_file
        )
        
        # 結合親和性の変化の見積もり
        binding_affinity_change = self.binding_affinity_estimator.estimate(
            msmd_settings1, msmd_settings2
        )
        
        # 結果の管理
        result = {
            'non_common_parts': non_common_parts,
            'probe1': probe1,
            'probe2': probe2,
            'transformed_probe1': transformed_probe1,
            'transformed_probe2': transformed_probe2,
            'msmd_settings1': msmd_settings1,
            'msmd_settings2': msmd_settings2,
            'binding_affinity_change': binding_affinity_change
        }
        
        return result
```

## 他モジュールとの連携ポイント

Quantitative Inverse MSMDモジュールは、以下のモジュールと連携しています：

- **MCS計算と非共通部分抽出モジュール**：非共通部分の抽出に使用します。
- **プローブ設計モジュール**：プローブの設計に使用します。
- **分子の重ね合わせ自動化モジュール**：分子の重ね合わせに使用します。
- **MSMDシミュレーションエンジン**：MSMDシミュレーションの設定に使用します。

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
# 設定情報
config = {
    'mcs': {
        'atom_compare': 'elements',  # 原子タイプの一致条件
        'bond_compare': 'any',  # 結合タイプの一致条件
        'ring_matches_ring_only': True,  # 環構造の一致条件
        'complete_rings_only': True  # 環構造の完全一致条件
    },
    'probe': {
        'template': 'methylbenzene',  # プローブテンプレート
        'optimize': True  # プローブの最適化
    },
    'superimpose': {
        'algorithm': 'kabsch',  # 重ね合わせアルゴリズム
        'optimize': True  # 重ね合わせの最適化
    },
    'binding_affinity': {
        'model': 'odds_ratio',  # 結合親和性推定モデル
        'normalize': True  # 結果の正規化
    }
}

# Quantitative Inverse MSMDクラスの初期化
qimsmd = QuantitativeInverseMSMD(config)

# Quantitative Inverse MSMDの実行
result = qimsmd.run(
    'path/to/compound1.sdf',
    'path/to/compound2.sdf',
    'path/to/protein.pdb'
)

# 結果の表示
print("非共通部分:")
print(f"  化合物1: {result['non_common_parts']['mol1']}")
print(f"  化合物2: {result['non_common_parts']['mol2']}")

print("プローブ:")
print(f"  プローブ1: {result['probe1']}")
print(f"  プローブ2: {result['probe2']}")

print("結合親和性の変化:")
print(f"  ΔΔG: {result['binding_affinity_change']['ddg']} kcal/mol")
print(f"  信頼度: {result['binding_affinity_change']['confidence']}")

# 結果の可視化
# （省略）
```

## 応用例

Quantitative Inverse MSMDモジュールは、以下のような応用例があります：

1. **構造活性相関（SAR）の解析**：化合物の構造変化と活性の関係を解析します。
2. **リード最適化**：リード化合物の最適化に役立つ情報を提供します。
3. **ファーマコフォアモデルの構築**：活性に重要な構造的特徴を特定します。
4. **新規化合物の設計**：結合親和性の高い新規化合物を設計します。