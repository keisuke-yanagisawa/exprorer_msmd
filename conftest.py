import glob
import json
import os
import sys


def get_ambertools_version() -> str:
    """Get the installed AmberTools version from conda metadata.

    Searches CONDA_PREFIX first, then falls back to sys.prefix and
    sys.base_prefix for environments where CONDA_PREFIX is unset
    (e.g., devcontainers activated via PATH rather than conda activate).
    """
    search_paths: list = []
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        search_paths.append(conda_prefix)
    search_paths.append(sys.prefix)
    if sys.base_prefix != sys.prefix:
        search_paths.append(sys.base_prefix)

    for prefix in search_paths:
        meta_files = glob.glob(f"{prefix}/conda-meta/ambertools-*.json")
        if meta_files:
            with open(meta_files[0]) as f:
                return json.load(f).get("version", "unknown")
    return "unknown"


AMBERTOOLS_VERSION = get_ambertools_version()
