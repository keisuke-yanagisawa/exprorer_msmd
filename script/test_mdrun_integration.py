import os
import shutil
import pytest
from pathlib import Path
import subprocess

class TestMDRunIntegration:
    @pytest.fixture(scope="class")
    def check_dependencies(self):
        """外部依存ソフトウェアのチェック"""
        dependencies = {
            "gmx": "gromacs",
            "packmol": "packmol",
            "tleap": "tleap",
            "cpptraj": "cpptraj"
        }
        missing = []
        for cmd, name in dependencies.items():
            if os.system(f"which {cmd} > /dev/null 2>&1") != 0:
                missing.append(name)
        if missing:
            pytest.skip(f"Missing dependencies: {', '.join(missing)}")

    @pytest.fixture
    def minimal_msmd_yaml(self, tmp_path):
        """最小構成のMSMD設定（MDシミュレーションのみ）"""
        yaml_content = """
general:
  iter_index: 0
  workdir: .
  name: test_run
  executables:
    python: python
    gromacs: gmx
    packmol: packmol
    tleap: tleap
    cpptraj: cpptraj

input:
  protein:
    pdb: tripeptide.pdb
  probe:
    cid: A11
    atomtype: gaff2
    molar: 0.25

exprorer_msmd:
  title: Test protocol
  general:
    seed: 42
    pbc: xyz
  sequence:
    - name: min
      type: minimization
      nsteps: 10
      nstlog: 1
    - name: heat
      type: heating
      nsteps: 10
      nstlog: 1
    - name: pr
      type: production
      nsteps: 10
      nstlog: 1
"""
        yaml_path = tmp_path / "test_msmd.yaml"
        yaml_path.write_text(yaml_content)

        # テストデータのコピー
        test_data_dir = Path("script/test_data")
        for file in ["tripeptide.pdb", "A11.mol2", "A11.pdb"]:
            shutil.copy2(test_data_dir / file, tmp_path)

        return yaml_path

    def test_msmd_minimal_run(self, check_dependencies, minimal_msmd_yaml, tmp_path):
        """最小構成でのMDシミュレーション実行テスト"""
        # コマンド実行
        cmd = f"timeout 2m ./exprorer_msmd {minimal_msmd_yaml} --skip-postprocess"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        # return_codeの確認
        assert result.returncode == 0, f"Command failed with:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        
        # シミュレーションディレクトリのパス
        sim_dir = tmp_path / "system0" / "simulation"
        
        # 各ステージの出力確認
        for stage in ["min", "heat", "pr"]:
            assert (sim_dir / f"{stage}.mdp").exists(), f"{stage}.mdp was not created"
            assert (sim_dir / f"{stage}.tpr").exists(), f"{stage}.tpr was not created"
            assert (sim_dir / f"{stage}.log").exists(), f"{stage}.log was not created"

    def test_msmd_missing_files(self, check_dependencies, tmp_path):
        """必要なファイルが存在しない場合のテスト"""
        yaml_content = """
general:
  iter_index: 0
  workdir: .
  name: test_run
  executables:
    python: python
    gromacs: gmx
    packmol: packmol
    tleap: tleap
    cpptraj: cpptraj

input:
  protein:
    pdb: nonexistent.pdb
  probe:
    cid: A11
    atomtype: gaff2
    molar: 0.25

exprorer_msmd:
  title: Test protocol
  general:
    seed: 42
    pbc: xyz
  sequence:
    - name: min
      type: minimization
      nsteps: 10
      nstlog: 1
"""
        yaml_path = tmp_path / "test_msmd.yaml"
        yaml_path.write_text(yaml_content)

        # コマンド実行
        cmd = f"timeout 2m ./exprorer_msmd {yaml_path} --skip-postprocess"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        # エラーが発生することを確認
        assert result.returncode != 0, "Command should fail with missing files"