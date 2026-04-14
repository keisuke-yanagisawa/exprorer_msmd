import subprocess

import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock

from script.mdrun import (
    detect_mdrun_thread_flag,
    _detect_mpi_type,
    gen_mdp,
    gen_mdrun_job,
    prepare_sequence,
    prepare_md_files,
    run_md_sequence,
)

@pytest.fixture
def general_settings():
    """Fixture providing general MD settings"""
    return {
        "temperature": 300,
        "pressure": 1.0,
        "dt": 0.002,
        "pbc": "xyz"
    }

@pytest.fixture
def mock_jinja():
    """Fixture providing mock for Jinja2 environment"""
    with patch('jinja2.Environment') as mock_env:
        mock_template = MagicMock()
        mock_template.render.return_value = "rendered template content"
        mock_env.return_value.get_template.return_value = mock_template
        yield mock_env, mock_template

# Test only representative protocols
class TestGenMdp:
    @pytest.mark.parametrize("protocol", [
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
        """Test MDP file generation based on representative protocols"""
        mock_env, mock_template = mock_jinja
        
        gen_mdp(protocol, tmp_path)
        
        # Verify template selection
        mock_env.return_value.get_template.assert_called_once_with(
            f"./template/{protocol['type']}.mdp"
        )
        
        # Verify protocol parameters are passed correctly
        mock_template.render.assert_called_once_with(protocol)
        
        # For heating type, verify additional parameters
        if protocol["type"] == "heating":
            expected_target_temp = protocol.get("target_temp", protocol["temperature"])
            expected_initial_temp = protocol.get("initial_temp", 0)
            expected_duration = protocol["nsteps"] * protocol["dt"]
            
            assert protocol["target_temp"] == expected_target_temp
            assert protocol["initial_temp"] == expected_initial_temp
            assert protocol["duration"] == pytest.approx(expected_duration)

    @pytest.mark.parametrize("invalid_type", ["invalid", None, ""])
    def test_gen_mdp_invalid_types(self, tmp_path, invalid_type):
        """Test error handling for invalid simulation types"""
        protocol = {
            "type": invalid_type,
            "name": "test"
        }
        
        with pytest.raises(ValueError, match=f"Invalid simulation type: {invalid_type}"):
            gen_mdp(protocol, tmp_path)

# Test gen_mdrun_job
class TestGenMdrunJob:
    @pytest.fixture
    def job_script(self, tmp_path):
        """Fixture providing job script path"""
        return tmp_path / "mdrun.sh"

    def test_gen_mdrun_job_sequence(self, job_script, mock_jinja):
        """Test job script generation for standard simulation sequence"""
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

# Test prepare_sequence
class TestPrepareSequence:
    def test_prepare_sequence_protocol(self, general_settings):
        """Test preparation of standard MD sequence"""
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
                "temperature": 310  # Override general settings
            }
        ]
        
        result = prepare_sequence(sequence, general_settings)
        
        assert len(result) == len(sequence)
        # Verify general settings are applied in first step
        assert result[0]["temperature"] == general_settings["temperature"]
        assert result[0]["pressure"] == general_settings["pressure"]
        # Verify settings are overridden in second step
        assert result[1]["temperature"] == 310
        assert result[1]["pressure"] == general_settings["pressure"]

# Test prepare_md_files
class TestPrepareMdFiles:
    def test_prepare_md_files_protocol(self, tmp_path, mock_jinja):
        """Test file preparation for standard MD protocol"""
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
        
        # Verify seed value settings
        for step in sequence:
            assert step["seed"] == index
        
        # Verify number of gen_mdp calls
        mock_env, _ = mock_jinja
        assert mock_env.return_value.get_template.call_count == len(sequence) + 1  # +1 for mdrun.sh template


class TestDetectMdrunThreadFlag:
    """Tests for GROMACS MPI build type detection"""

    def setup_method(self):
        _detect_mpi_type.cache_clear()

    def test_thread_mpi_build(self):
        """thread-MPI build should return -nt"""
        with patch('shutil.which', return_value='/usr/bin/gmx'), \
             patch('subprocess.run') as mock_run:
            mock_run.return_value = MagicMock(
                stdout="GROMACS version: 2024\nMPI library:      thread_mpi\n",
                stderr=""
            )
            assert detect_mdrun_thread_flag(Path("gmx")) == "-nt"

    def test_real_mpi_build(self):
        """real MPI build should return -ntomp"""
        with patch('shutil.which', return_value='/usr/bin/gmx_mpi'), \
             patch('subprocess.run') as mock_run:
            mock_run.return_value = MagicMock(
                stdout="GROMACS version: 2024\nMPI library:         MPI (GPU-aware: CUDA)\n",
                stderr=""
            )
            assert detect_mdrun_thread_flag(Path("gmx_mpi")) == "-ntomp"

    def test_executable_not_found(self):
        """Missing executable should fall back to -nt"""
        with patch('shutil.which', return_value=None):
            assert detect_mdrun_thread_flag(Path("nonexistent_gmx")) == "-nt"

    def test_version_timeout(self):
        """Timeout on --version should fall back to -nt"""
        with patch('shutil.which', return_value='/usr/bin/gmx'), \
             patch('subprocess.run', side_effect=subprocess.TimeoutExpired('gmx', 30)):
            assert detect_mdrun_thread_flag(Path("gmx_timeout")) == "-nt"

    def test_no_mpi_line_in_output(self):
        """Missing MPI library line should fall back to -nt"""
        with patch('shutil.which', return_value='/usr/bin/gmx'), \
             patch('subprocess.run') as mock_run:
            mock_run.return_value = MagicMock(
                stdout="GROMACS version: 2024\nSome other info\n",
                stderr=""
            )
            assert detect_mdrun_thread_flag(Path("gmx_no_mpi_line")) == "-nt"


class TestRunMdSequence:
    """Tests for run_md_sequence THREAD_FLAG passing"""

    def test_thread_flag_passed_to_shell(self):
        """Verify THREAD_FLAG environment variable is set in os.system command"""
        with patch('script.mdrun.detect_mdrun_thread_flag', return_value='-ntomp'), \
             patch('os.system') as mock_system:
            run_md_sequence(0, Path("/tmp/sim"), Path("gmx_mpi"), 4, "test")
            call_args = mock_system.call_args[0][0]
            assert 'THREAD_FLAG=-ntomp' in call_args
