import math
import sys
import subprocess
import pytest
from unittest.mock import patch, MagicMock

from script.utilities.GPUtil import get_gpuids, is_mps_control_running


@pytest.fixture
def mock_gpu_env():
    """Basic mock fixture for GPU testing"""
    with patch('GPUtil.getAvailable') as mock_available, \
         patch('os.getenv') as mock_getenv, \
         patch('script.utilities.GPUtil.logger') as mock_logger:
        yield {
            'available': mock_available,
            'getenv': mock_getenv,
            'logger': mock_logger
        }


def test_get_gpuids_no_gpu(mock_gpu_env):
    """
    Test for when no GPU exists
    Expected: Returns CPU only mode ([-1])
    """
    mock_gpu_env['available'].return_value = []
    mock_gpu_env['getenv'].return_value = None
    
    gpuids = get_gpuids()
    
    assert gpuids == [-1]
    mock_gpu_env['logger'].warn.assert_any_call("No GPU is allowed/existed to use")
    mock_gpu_env['logger'].warn.assert_any_call("Switch to CPU-only mode, it greatly decreases the simulation speed")


def test_get_gpuids_single_gpu(mock_gpu_env):
    """
    Test for when one GPU is available (close to real environment)
    Expected: Returns [0]
    """
    mock_gpu_env['available'].return_value = [0]
    mock_gpu_env['getenv'].return_value = None
    
    gpuids = get_gpuids()
    
    assert gpuids == [0]
    mock_gpu_env['logger'].info.assert_any_call("1 GPUs are detected")
    mock_gpu_env['logger'].info.assert_any_call(f"GPU IDs of {[0]} will be used")


def test_get_gpuids_multiple_gpus(mock_gpu_env):
    """
    Test for when multiple GPUs are available
    Expected: Returns a list of all available GPU IDs
    """
    mock_gpu_env['available'].return_value = [0, 1, 2]
    mock_gpu_env['getenv'].return_value = None  # CUDA_VISIBLE_DEVICES is undefined
    
    gpuids = get_gpuids()
    
    assert gpuids == [0, 1, 2]
    mock_gpu_env['logger'].info.assert_any_call("3 GPUs are detected")


def test_get_gpuids_with_cuda_visible_devices(mock_gpu_env):
    """
    Test for when CUDA_VISIBLE_DEVICES is set
    Expected: Returns only GPU IDs specified in CUDA_VISIBLE_DEVICES
    """
    mock_gpu_env['available'].return_value = [0, 1, 2, 3]
    mock_gpu_env['getenv'].side_effect = lambda x, default=None: "1,2" if x == "CUDA_VISIBLE_DEVICES" else default
    
    gpuids = get_gpuids()
    
    assert set(gpuids) == {1, 2}
    mock_gpu_env['logger'].info.assert_any_call("CUDA_VISIBLE_DEVICES detected")


def test_get_gpuids_ignore_cuda_visible_devices(mock_gpu_env):
    """
    Test for when ignoring CUDA_VISIBLE_DEVICES
    Expected: Ignores CUDA_VISIBLE_DEVICES setting and returns all available GPU IDs
    """
    mock_gpu_env['available'].return_value = [0, 1, 2]
    mock_gpu_env['getenv'].side_effect = lambda x, default=None: "1" if x == "CUDA_VISIBLE_DEVICES" else default
    
    gpuids = get_gpuids(ignore_cuda_visible_devices=True)
    
    assert gpuids == [0, 1, 2]
    # Verify that CUDA_VISIBLE_DEVICES log is not output
    assert not any("CUDA_VISIBLE_DEVICES detected" in str(call) 
                  for call in mock_gpu_env['logger'].info.call_args_list)


def test_get_gpuids_empty_cuda_visible_devices(mock_gpu_env):
    """
    Test for when CUDA_VISIBLE_DEVICES is an empty string
    Expected: Returns CPU only mode ([-1])
    """
    mock_gpu_env['available'].return_value = [0, 1, 2]
    mock_gpu_env['getenv'].side_effect = lambda x, default=None: "" if x == "CUDA_VISIBLE_DEVICES" else default
    
    gpuids = get_gpuids()
    
    assert gpuids == [-1]  # Do not use GPU
    mock_gpu_env['logger'].info.assert_any_call("CUDA_VISIBLE_DEVICES detected")
    mock_gpu_env['logger'].warn.assert_any_call("No GPU is allowed/existed to use")


def test_get_gpuids_maxload_maxmemory(mock_gpu_env):
    """
    Test for GPUtil.getAvailable call parameters
    Expected: GPUtil.getAvailable is called with correct parameters
    """
    mock_gpu_env['available'].return_value = [0, 1]
    mock_gpu_env['getenv'].return_value = None
    
    gpuids = get_gpuids()
    
    mock_gpu_env['available'].assert_called_once_with(
        maxLoad=math.inf,
        maxMemory=math.inf,
        limit=sys.maxsize
    )
    assert gpuids == [0, 1]


@pytest.mark.parametrize("output,expected,error", [
    ("100", True, None),  # When operating normally
    ("Cannot find MPS control daemon process", False, None),  # When daemon is not found
    ("Unexpected output", True, None),  # When output is unexpected
    (None, False, subprocess.CalledProcessError(1, 'cmd')),  # When error occurs
])
def test_is_mps_control_running(output, expected, error):
    """
    Test for MPS control state check
    Parameterized test for various cases
    """
    with patch('subprocess.check_output') as mock_check_output:
        if error:
            mock_check_output.side_effect = error
        else:
            mock_check_output.return_value = output
        
        result = is_mps_control_running()
        assert result == expected

        if not error:
            mock_check_output.assert_called_once_with(
                "echo get_default_active_thread_percentage | nvidia-cuda-mps-control",
                shell=True,
                text=True
            )