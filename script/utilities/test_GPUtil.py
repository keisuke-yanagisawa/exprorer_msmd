import math
import sys
import subprocess
import pytest
from unittest.mock import patch, MagicMock

from script.utilities.GPUtil import get_gpuids, is_mps_control_running


@pytest.fixture
def mock_gpu_env():
    """GPUテスト用の基本的なモックフィクスチャ"""
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
    GPUが存在しない場合のテスト
    期待値: CPU onlyモード（[-1]）を返す
    """
    mock_gpu_env['available'].return_value = []
    mock_gpu_env['getenv'].return_value = None
    
    gpuids = get_gpuids()
    
    assert gpuids == [-1]
    mock_gpu_env['logger'].warn.assert_any_call("No GPU is allowed/existed to use")
    mock_gpu_env['logger'].warn.assert_any_call("Switch to CPU-only mode, it greatly decreases the simulation speed")


def test_get_gpuids_single_gpu(mock_gpu_env):
    """
    1基のGPUが利用可能な場合のテスト（実環境に近い状況）
    期待値: [0]を返す
    """
    mock_gpu_env['available'].return_value = [0]
    mock_gpu_env['getenv'].return_value = None
    
    gpuids = get_gpuids()
    
    assert gpuids == [0]
    mock_gpu_env['logger'].info.assert_any_call("1 GPUs are detected")
    mock_gpu_env['logger'].info.assert_any_call(f"GPU IDs of {[0]} will be used")


def test_get_gpuids_multiple_gpus(mock_gpu_env):
    """
    複数のGPUが利用可能な場合のテスト
    期待値: 利用可能なすべてのGPU IDのリストを返す
    """
    mock_gpu_env['available'].return_value = [0, 1, 2]
    mock_gpu_env['getenv'].return_value = None  # CUDA_VISIBLE_DEVICESが未定義
    
    gpuids = get_gpuids()
    
    assert gpuids == [0, 1, 2]
    mock_gpu_env['logger'].info.assert_any_call("3 GPUs are detected")


def test_get_gpuids_with_cuda_visible_devices(mock_gpu_env):
    """
    CUDA_VISIBLE_DEVICESが設定されている場合のテスト
    期待値: CUDA_VISIBLE_DEVICESで指定されたGPU IDのみを返す
    """
    mock_gpu_env['available'].return_value = [0, 1, 2, 3]
    mock_gpu_env['getenv'].side_effect = lambda x, default=None: "1,2" if x == "CUDA_VISIBLE_DEVICES" else default
    
    gpuids = get_gpuids()
    
    assert set(gpuids) == {1, 2}
    mock_gpu_env['logger'].info.assert_any_call("CUDA_VISIBLE_DEVICES detected")


def test_get_gpuids_ignore_cuda_visible_devices(mock_gpu_env):
    """
    CUDA_VISIBLE_DEVICESを無視する場合のテスト
    期待値: CUDA_VISIBLE_DEVICESの設定を無視し、すべての利用可能なGPU IDを返す
    """
    mock_gpu_env['available'].return_value = [0, 1, 2]
    mock_gpu_env['getenv'].side_effect = lambda x, default=None: "1" if x == "CUDA_VISIBLE_DEVICES" else default
    
    gpuids = get_gpuids(ignore_cuda_visible_devices=True)
    
    assert gpuids == [0, 1, 2]
    # CUDA_VISIBLE_DEVICESのログが出力されていないことを確認
    assert not any("CUDA_VISIBLE_DEVICES detected" in str(call) 
                  for call in mock_gpu_env['logger'].info.call_args_list)


def test_get_gpuids_empty_cuda_visible_devices(mock_gpu_env):
    """
    CUDA_VISIBLE_DEVICESが空文字列の場合のテスト
    期待値: CPU onlyモード（[-1]）を返す
    """
    mock_gpu_env['available'].return_value = [0, 1, 2]
    mock_gpu_env['getenv'].side_effect = lambda x, default=None: "" if x == "CUDA_VISIBLE_DEVICES" else default
    
    gpuids = get_gpuids()
    
    assert gpuids == [-1]  # GPUを利用しない
    mock_gpu_env['logger'].info.assert_any_call("CUDA_VISIBLE_DEVICES detected")
    mock_gpu_env['logger'].warn.assert_any_call("No GPU is allowed/existed to use")


def test_get_gpuids_maxload_maxmemory(mock_gpu_env):
    """
    GPUtil.getAvailableの呼び出しパラメータのテスト
    期待値: 正しいパラメータでGPUtil.getAvailableが呼び出される
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
    ("100", True, None),  # 正常に動作している場合
    ("Cannot find MPS control daemon process", False, None),  # デーモンが見つからない場合
    ("Unexpected output", True, None),  # 予期しない出力の場合
    (None, False, subprocess.CalledProcessError(1, 'cmd')),  # エラーが発生した場合
])
def test_is_mps_control_running(output, expected, error):
    """
    MPSコントロールの状態チェックのテスト
    様々なケースをパラメータ化してテスト
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