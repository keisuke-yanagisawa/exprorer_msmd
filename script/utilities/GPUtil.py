import math
import os
import sys
import subprocess

import GPUtil


from .logger import logger


def get_gpuids(ignore_cuda_visible_devices=False):
    gpuids = set(GPUtil.getAvailable(maxLoad=math.inf, maxMemory=math.inf, limit=sys.maxsize))
    logger.info(f"{len(gpuids)} GPUs are detected")
    if (ignore_cuda_visible_devices is False) and (os.getenv("CUDA_VISIBLE_DEVICES") is not None):
        logger.info("CUDA_VISIBLE_DEVICES detected")
        cvd = [int(s) for s in os.getenv("CUDA_VISIBLE_DEVICES", default="").split(",")]
        gpuids &= set(cvd)
        gpuids = list(gpuids)

    if len(gpuids) == 0:
        logger.warn("No GPU is allowed/existed to use")
        logger.warn("Switch to CPU-only mode, it greatly decreases the simulation speed")
    else:
        logger.info(f"GPU IDs of {gpuids} will be used")

    ngpus = len(gpuids)
    if ngpus == 0:  # there is only CPU
        gpuids = [-1]
        ngpus = 1

    return list(gpuids)

def is_mps_control_running() -> bool:
    # check if nvidia-cuda-mps-control process is running
    try:
        output = subprocess.check_output("ps x | grep nvidia-cuda-mps-control", shell=True, text=True)
        lines = output.strip().split('\n')
        if len(lines) >= 3:
            return True
        else:
            return False
    except subprocess.CalledProcessError:
        return False
