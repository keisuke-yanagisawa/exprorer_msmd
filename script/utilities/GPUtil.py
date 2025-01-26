import math
import os
import subprocess
import sys

import GPUtil

from .logger import logger


def get_gpuids(ignore_cuda_visible_devices=False):
    """
    Retrieve available GPU IDs, respecting the CUDA_VISIBLE_DEVICES environment variable if applicable.
    This function checks for GPUs that are currently available using GPUtil, logging how many are found. 
    If ignore_cuda_visible_devices is False and the CUDA_VISIBLE_DEVICES environment variable is set:
        • An empty CUDA_VISIBLE_DEVICES string implies no GPU usage.
        • Otherwise, only the GPUs specified in CUDA_VISIBLE_DEVICES are used.
    If no GPUs remain after applying these filters, a warning is logged and the function falls back to CPU-only mode by returning [-1].
    Args:
            ignore_cuda_visible_devices (bool): 
                    If True, ignores the CUDA_VISIBLE_DEVICES environment variable 
                    and uses all available GPUs. Defaults to False.
    Returns:
            list of int:
                    A list of GPU IDs to be used, or [-1] if no GPUs are available.
    """
    gpuids = set(GPUtil.getAvailable(maxLoad=math.inf, maxMemory=math.inf, limit=sys.maxsize))
    logger.info(f"{len(gpuids)} GPUs are detected")
    
    if (ignore_cuda_visible_devices is False) and (os.getenv("CUDA_VISIBLE_DEVICES") is not None):
        logger.info("CUDA_VISIBLE_DEVICES detected")
        cuda_devices = os.getenv("CUDA_VISIBLE_DEVICES", default="")
        if cuda_devices == "":  # 空文字列の場合はGPUを利用しない
            gpuids = set()
        else:
            cvd = [int(s) for s in cuda_devices.split(",") if s.strip()]
            gpuids &= set(cvd)
    
    if len(gpuids) == 0:
        logger.warn("No GPU is allowed/existed to use")
        logger.warn("Switch to CPU-only mode, it greatly decreases the simulation speed")
    else:
        logger.info(f"GPU IDs of {list(gpuids)} will be used")

    gpuids = list(gpuids)
    if len(gpuids) == 0:  # there is only CPU
        gpuids = [-1]

    return gpuids


def is_mps_control_running() -> bool:
    # check if nvidia-cuda-mps-control process is running
    try:
        output = subprocess.check_output(
            "echo get_default_active_thread_percentage | nvidia-cuda-mps-control",
            shell=True,
            text=True
        )
        if "Cannot find MPS control daemon process" in output:
            return False
        else:
            return True
    except subprocess.CalledProcessError:
        return False
