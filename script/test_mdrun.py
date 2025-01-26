import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock

from script.mdrun import (
    gen_mdp,
    gen_mdrun_job,
    prepare_sequence,
    prepare_md_files
)

@pytest.fixture
def general_settings():
    """一般的なMD設定を提供するフィクスチャ"""
    return {
        "temperature": 300,
        "pressure": 1.0,
        "dt": 0.002,
        "pbc": "xyz"
    }

@pytest.fixture
def mock_jinja():
    """Jinja2環境のモックを提供するフィクスチャ"""
    with patch('jinja2.Environment') as mock_env:
        mock_template = MagicMock()
        mock_template.render.return_value = "rendered template content"
        mock_env.return_value.get_template.return_value = mock_template
        yield mock_env, mock_template

# gen_mdpのテスト
class TestGenMdp:
    @pytest.mark.parametrize("protocol", [
        # 代表的なプロトコルのみをテスト
        {
            "type": "minimization",
            "name": "min",
            "title": "Energy minimization",
            "define": "-DPOSRES1000",
            "nsteps": 200,
            "nstlog": 1,
            "dt": 0.002
        },
        {
            "type": "heating",
            "name": "heat",
            "title": "Heating system",
            "define": "-DPOSRES1000",
            "nsteps": 100000,
            "nstxtcout": 500,
            "nstlog": 500,
            "nstenergy": 500,
            "temperature": 300,
            "initial_temp": 100,
            "target_temp": 300,
            "pcoupl": "no",
            "dt": 0.002
        },
        {
            "type": "production",
            "name": "pr",
            "title": "Production Run",
            "define": "",
            "nsteps": 20000000,
            "nstxtcout": 5000,
            "nstenergy": 5000,
            "nstlog": 5000,
            "pcoupl": "Parrinello-Rahman",
            "dt": 0.002
        }
    ])
    def test_gen_mdp_protocols(self, tmp_path, mock_jinja, protocol):
        """代表的なプロトコルに基づくMDPファイル生成テスト"""
        mock_env, mock_template = mock_jinja
        
        gen_mdp(protocol, tmp_path)
        
        # テンプレートの選択を確認
        mock_env.return_value.get_template.assert_called_once_with(
            f"./template/{protocol['type']}.mdp"
        )
        
        # プロトコルパラメータが正しく渡されることを確認
        mock_template.render.assert_called_once_with(protocol)
        
        # heatingタイプの場合、追加パラメータを確認
        if protocol["type"] == "heating":
            expected_target_temp = protocol.get("target_temp", protocol["temperature"])
            expected_initial_temp = protocol.get("initial_temp", 0)
            expected_duration = protocol["nsteps"] * protocol["dt"]
            
            assert protocol["target_temp"] == expected_target_temp
            assert protocol["initial_temp"] == expected_initial_temp
            assert protocol["duration"] == pytest.approx(expected_duration)

    @pytest.mark.parametrize("invalid_type", ["invalid", None, ""])
    def test_gen_mdp_invalid_types(self, tmp_path, invalid_type):
        """無効なシミュレーションタイプでのエラー処理テスト"""
        protocol = {
            "type": invalid_type,
            "name": "test"
        }
        
        with pytest.raises(ValueError, match=f"Invalid simulation type: {invalid_type}"):
            gen_mdp(protocol, tmp_path)

# gen_mdrun_jobのテスト
class TestGenMdrunJob:
    @pytest.fixture
    def job_script(self, tmp_path):
        """ジョブスクリプトのパスを提供するフィクスチャ"""
        return tmp_path / "mdrun.sh"

    def test_gen_mdrun_job_sequence(self, job_script, mock_jinja):
        """標準的なシミュレーションシーケンスでのジョブスクリプト生成テスト"""
        mock_env, mock_template = mock_jinja
        
        sequence = ["min", "heat", "pr"]
        name = "TEST_PROJECT"
        top = Path("system.top")
        gro = Path("system.gro")
        out_traj = Path("traj.xtc")
        post_comm = "echo 'Simulation completed'"
        
        gen_mdrun_job(sequence, name, job_script, top, gro, out_traj, post_comm)
        
        render_data = mock_template.render.call_args[0][0]
        assert render_data["NAME"] == name
        assert render_data["TOP"] == top
        assert render_data["GRO"] == gro
        assert render_data["OUT_TRAJ"] == out_traj
        assert render_data["POST_COMMAND"] == post_comm
        assert render_data["STEP_NAMES"] == " ".join(sequence)

# prepare_sequenceのテスト
class TestPrepareSequence:
    def test_prepare_sequence_protocol(self, general_settings):
        """標準的なMDシーケンスの準備テスト"""
        sequence = [
            {
                "type": "minimization",
                "name": "min",
                "define": "-DPOSRES1000",
                "nsteps": 200
            },
            {
                "type": "production",
                "name": "pr",
                "define": "",
                "nsteps": 20000000,
                "temperature": 310  # 一般設定をオーバーライド
            }
        ]
        
        result = prepare_sequence(sequence, general_settings)
        
        assert len(result) == len(sequence)
        # 最初のステップで一般設定が適用されていることを確認
        assert result[0]["temperature"] == general_settings["temperature"]
        assert result[0]["pressure"] == general_settings["pressure"]
        # 2番目のステップで設定がオーバーライドされていることを確認
        assert result[1]["temperature"] == 310
        assert result[1]["pressure"] == general_settings["pressure"]

# prepare_md_filesのテスト
class TestPrepareMdFiles:
    def test_prepare_md_files_protocol(self, tmp_path, mock_jinja):
        """標準的なMDプロトコルでのファイル準備テスト"""
        sequence = [
            {
                "type": "minimization",
                "name": "min",
                "define": "-DPOSRES1000"
            },
            {
                "type": "production",
                "name": "pr",
                "define": ""
            }
        ]
        
        index = 42
        jobname = "TEST_PROJECT"
        top = Path("system.top")
        gro = Path("system.gro")
        out_traj = Path("traj.xtc")
        
        prepare_md_files(index, sequence, tmp_path, jobname, top, gro, out_traj)
        
        # シード値の設定を確認
        for step in sequence:
            assert step["seed"] == index
        
        # gen_mdpの呼び出し回数を確認
        mock_env, _ = mock_jinja
        assert mock_env.return_value.get_template.call_count == len(sequence) + 1  # +1 for mdrun.sh template
