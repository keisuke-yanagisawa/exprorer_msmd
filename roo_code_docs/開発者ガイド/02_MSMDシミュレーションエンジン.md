# MSMDシミュレーションエンジン

MSMDシミュレーションエンジンは、EXPRORER_MSMDの中核となるコンポーネントで、GROMACSを使用してMixed-Solvent Molecular Dynamics（MSMD）シミュレーションを実行します。このモジュールは、前処理、シミュレーション、後処理の3つのサブモジュールから構成されています。

## 前処理（preprocess）モジュール

### クラス構造と責任

前処理モジュールは、MSMDシミュレーションの前処理を担当します。主なクラスと責任は以下の通りです：

- **PreprocessManager**：前処理全体を管理するクラス
  - YAMLファイルの解析
  - タンパク質構造の準備
  - プローブ分子の準備
  - シミュレーションシステムの構築

- **ProteinPreprocessor**：タンパク質構造の前処理を担当するクラス
  - PDBファイルの読み込み
  - 水素原子の追加
  - プロトン化状態の最適化
  - 力場パラメータの割り当て

- **ProbePreprocessor**：プローブ分子の前処理を担当するクラス
  - SDFファイルの読み込み
  - 力場パラメータの割り当て
  - プローブ分子の配置

- **SystemBuilder**：シミュレーションシステムの構築を担当するクラス
  - 溶媒の追加
  - イオンの追加
  - 周期境界条件の設定
  - エネルギー最小化の設定

### 主要なメソッドと機能

```python
class PreprocessManager:
    def __init__(self, yaml_file):
        """
        前処理マネージャーを初期化する

        Parameters
        ----------
        yaml_file : str
            YAMLファイルのパス
        """
        self.yaml_file = yaml_file
        self.config = self._parse_yaml()
        self.protein_preprocessor = ProteinPreprocessor(self.config)
        self.probe_preprocessor = ProbePreprocessor(self.config)
        self.system_builder = SystemBuilder(self.config)

    def _parse_yaml(self):
        """
        YAMLファイルを解析する

        Returns
        -------
        dict
            設定情報
        """
        with open(self.yaml_file, 'r') as f:
            config = yaml.safe_load(f)
        return config

    def run(self):
        """
        前処理を実行する
        """
        # タンパク質構造の準備
        protein = self.protein_preprocessor.preprocess()
        
        # プローブ分子の準備
        probe = self.probe_preprocessor.preprocess()
        
        # シミュレーションシステムの構築
        system = self.system_builder.build(protein, probe)
        
        return system
```

### 他モジュールとの連携ポイント

前処理モジュールは、以下のモジュールと連携しています：

- **シミュレーションモジュール**：前処理モジュールが生成したシミュレーションシステムを使用してシミュレーションを実行します。
- **外部ツール**：AmberToolsやPackmolなどの外部ツールを使用して、タンパク質構造の準備やシミュレーションシステムの構築を行います。

## シミュレーション（simulation）モジュール

### クラス構造と責任

シミュレーションモジュールは、GROMACSを使用してMSMDシミュレーションを実行します。主なクラスと責任は以下の通りです：

- **SimulationManager**：シミュレーション全体を管理するクラス
  - シミュレーションの設定
  - GROMACSコマンドの実行
  - シミュレーション結果の管理

- **GromacsRunner**：GROMACSコマンドの実行を担当するクラス
  - エネルギー最小化
  - 平衡化
  - 本番シミュレーション

- **SimulationMonitor**：シミュレーションの進行状況を監視するクラス
  - シミュレーションの進行状況の表示
  - エラーの検出
  - リソース使用状況の監視

### 主要なメソッドと機能

```python
class SimulationManager:
    def __init__(self, config):
        """
        シミュレーションマネージャーを初期化する

        Parameters
        ----------
        config : dict
            設定情報
        """
        self.config = config
        self.gromacs_runner = GromacsRunner(config)
        self.monitor = SimulationMonitor(config)

    def run(self, system):
        """
        シミュレーションを実行する

        Parameters
        ----------
        system : dict
            シミュレーションシステム

        Returns
        -------
        dict
            シミュレーション結果
        """
        # エネルギー最小化
        self.gromacs_runner.minimize(system)
        
        # 平衡化
        self.gromacs_runner.equilibrate(system)
        
        # 本番シミュレーション
        trajectory = self.gromacs_runner.production(system)
        
        # シミュレーション結果の管理
        result = {
            'trajectory': trajectory,
            'energy': self.gromacs_runner.get_energy(),
            'log': self.gromacs_runner.get_log()
        }
        
        return result
```

### 他モジュールとの連携ポイント

シミュレーションモジュールは、以下のモジュールと連携しています：

- **前処理モジュール**：前処理モジュールが生成したシミュレーションシステムを使用してシミュレーションを実行します。
- **後処理モジュール**：シミュレーション結果を後処理モジュールに渡して解析します。
- **外部ツール**：GROMACSを使用してシミュレーションを実行します。

## 後処理（postprocess）モジュール

### クラス構造と責任

後処理モジュールは、シミュレーション結果の解析を担当します。主なクラスと責任は以下の通りです：

- **PostprocessManager**：後処理全体を管理するクラス
  - シミュレーション結果の読み込み
  - トラジェクトリの解析
  - PMAPの生成

- **TrajectoryAnalyzer**：トラジェクトリの解析を担当するクラス
  - トラジェクトリの読み込み
  - プローブ分子の位置の抽出
  - 統計解析

- **PmapGenerator**：PMAPの生成を担当するクラス
  - プローブ分子の位置からPMAPを生成
  - PMAPの正規化
  - PMAPの保存

### 主要なメソッドと機能

```python
class PostprocessManager:
    def __init__(self, config):
        """
        後処理マネージャーを初期化する

        Parameters
        ----------
        config : dict
            設定情報
        """
        self.config = config
        self.trajectory_analyzer = TrajectoryAnalyzer(config)
        self.pmap_generator = PmapGenerator(config)

    def run(self, simulation_result):
        """
        後処理を実行する

        Parameters
        ----------
        simulation_result : dict
            シミュレーション結果

        Returns
        -------
        dict
            後処理結果
        """
        # トラジェクトリの解析
        probe_positions = self.trajectory_analyzer.analyze(simulation_result['trajectory'])
        
        # PMAPの生成
        pmap = self.pmap_generator.generate(probe_positions)
        
        # 後処理結果の管理
        result = {
            'pmap': pmap,
            'statistics': self.trajectory_analyzer.get_statistics()
        }
        
        return result
```

### 他モジュールとの連携ポイント

後処理モジュールは、以下のモジュールと連携しています：

- **シミュレーションモジュール**：シミュレーションモジュールが生成したシミュレーション結果を使用して解析を行います。
- **タンパク質ホットスポット探索モジュール**：後処理モジュールが生成したPMAPを使用してホットスポットを探索します。
- **プローブ分子周辺残基環境取得モジュール**：後処理モジュールが生成したPMAPを使用して残基環境を解析します。
- **外部ツール**：CPPTRAJなどの外部ツールを使用してトラジェクトリの解析を行います。

## 実装上の注意点

### パフォーマンス最適化

MSMDシミュレーションは計算コストが高いため、以下のような最適化が行われています：

1. **並列計算**：GROMACSの並列計算機能を活用して、シミュレーションを高速化しています。
2. **GPUアクセラレーション**：GPUを使用してシミュレーションを高速化しています。
3. **効率的なデータ構造**：大量のデータを効率的に処理するために、適切なデータ構造を選択しています。
4. **メモリ管理**：大きなトラジェクトリデータを扱うために、効率的なメモリ管理を行っています。

### エラーハンドリング

MSMDシミュレーションでは、様々なエラーが発生する可能性があります。以下のようなエラーハンドリングが実装されています：

1. **入力チェック**：入力ファイルの形式や内容をチェックして、エラーを早期に検出します。
2. **シミュレーションモニタリング**：シミュレーション中にエラーが発生した場合、適切に対応します。
3. **リカバリーメカニズム**：シミュレーションが途中で中断された場合、再開できるようにチェックポイントを保存します。
4. **ログ記録**：エラーや警告を詳細にログに記録して、デバッグを容易にします。

### 拡張性

MSMDシミュレーションエンジンは、以下のような拡張性を持っています：

1. **新しいプローブの追加**：新しいプローブ分子を簡単に追加できるように設計されています。
2. **シミュレーションパラメータのカスタマイズ**：YAMLファイルを通じて、シミュレーションパラメータをカスタマイズできます。
3. **解析機能の拡張**：新しい解析機能を追加するためのインターフェースが用意されています。
4. **外部ツールの統合**：新しい外部ツールを統合するための仕組みが用意されています。