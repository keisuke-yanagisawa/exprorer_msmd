import math
import os

import GPUtil


from .logger import logger


def get_gpuids(ignore_cuda_visible_devices=False):
    gpuids = set(GPUtil.getAvailable(limit=math.inf))
    logger.info(f"{len(gpuids)} GPUs are detected")
    if (ignore_cuda_visible_devices==False) and (os.getenv("CUDA_VISIBLE_DEVICES") is not None):
        logger.info(f"CUDA_VISIBLE_DEVICES detected")
        cvd = [int(s) for s in os.getenv("CUDA_VISIBLE_DEVICES", default="").split(",")]
        gpuids &= set(cvd)
        gpuids = list(gpuids)

    if len(gpuids) == 0:
        logger.warn(f"No GPU is allowed/existed to use")
        logger.warn(f"Switch to CPU-only mode, it greatly decreases the simulation speed")
    else:
        logger.info(f"GPU IDs of {gpuids} will be used")
    
    ngpus = len(gpuids)
    if ngpus == 0: # there is only CPU
        gpuids = [-1]
        ngpus = 1
    
    return list(gpuids)
