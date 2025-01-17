import os
from pathlib import Path
from typing import Any, Dict, List, Optional, TypedDict, Union, cast

import yaml

from .utilities.util import expandpath, getabsolutepath, update_dict


class ProteinConfig(TypedDict):
    """Protein configuration settings."""
    pdb: Path
    ssbond: List[str]

class ProbeConfig(TypedDict):
    """Probe molecule configuration settings."""
    cid: str
    mol2: Optional[Path]
    pdb: Optional[Path]
    atomtype: str
    molar: float

class GeneralConfig(TypedDict):
    """General configuration settings."""
    workdir: Path
    multiprocessing: int
    num_process_per_gpu: int
    yaml: Path

class MapConfig(TypedDict):
    """Map generation configuration settings."""
    snapshot: str
    valid_dist: float
    map_size: int
    normalization: str
    aggregation: str
    maps: List[Dict[str, str]]

class MSMDConfig(TypedDict):
    """MSMD simulation configuration settings."""
    dt: float
    temperature: float
    pressure: float


class Configuration:
    """Configuration manager for MSMD simulations."""
    
    def __init__(self, yaml_path: Path) -> None:
        """Initialize configuration from YAML file.
        
        Args:
            yaml_path: Path to YAML configuration file
            
        Raises:
            FileNotFoundError: If YAML file does not exist
            IsADirectoryError: If path points to directory
            ValueError: If file extension is not .yaml
        """
        self._yaml_path = getabsolutepath(yaml_path)
        self._yaml_dir = self._yaml_path.parent
        self._validate_yaml_path()
        self._config = self._load_default_config()
        self._load_yaml()
        self._ensure_compatibility_v1_1()
        self._resolve_paths()


    def _validate_yaml_path(self) -> None:
        """Validate YAML file path."""
        if not self._yaml_path.exists():
            raise FileNotFoundError(f"YAML file not found: {self._yaml_path}")
        if self._yaml_path.is_dir():
            raise IsADirectoryError(f"Given YAML file path {self._yaml_path} is a directory")
        if not os.path.splitext(self._yaml_path)[1][1:] == "yaml":
            raise ValueError(f"YAML file must have .yaml extension: {self._yaml_path}")
            
    def _load_default_config(self) -> Dict[str, Any]:
        """Load default configuration settings."""
        return {
            "general": {
                "workdir": Path(""),
                "multiprocessing": -1,
                "num_process_per_gpu": 1,
            },
            "input": {
                "protein": {
                    "pdb": Path(""),
                },
                "probe": {
                    "cid": "",
                },
            },
            "exprorer_msmd": {
                "general": {
                    "dt": 0.002,
                    "temperature": 300,
                    "pressure": 1.0,
                },
            },
            "map": {
                "snapshot": "",
                "valid_dist": 5.0,
                "map_size": 80,
                "normalization": "total",
                "aggregation": "max",
                "maps": [
                    {
                        "suffix": "nVH",
                        "selector": "(!@VIS)&(!@H*)",
                    }
                ],
            },
            "probe_profile": {
                "threshold": 0.001,
                "resenv": {
                    "env_dist": 4.0,
                },
            },
        }
        
    def _load_yaml(self) -> None:
        """Load and merge YAML configuration."""
        with self._yaml_path.open() as fin:
            yaml_dict = yaml.safe_load(fin)
            if yaml_dict is not None:
                update_dict(self._config, yaml_dict)
                
    def _ensure_compatibility_v1_1(self) -> None:
        """Ensure compatibility with exprorer_msmd v1.1."""
        if "pmap" in self._config["exprorer_msmd"]:
            update_dict(self._config["map"], self._config["exprorer_msmd"]["pmap"])
        self._config["map"]["snapshot"] = self._config["map"]["snapshot"].split("|")[-1]


    def _resolve_paths(self) -> None:
        """Resolve relative paths in configuration."""
        # Set probe paths if not specified
        if "mol2" not in self._config["input"]["probe"] or self._config["input"]["probe"]["mol2"] is None:
            self._config["input"]["probe"]["mol2"] = self._config["input"]["probe"]["cid"] + ".mol2"
        if "pdb" not in self._config["input"]["probe"] or self._config["input"]["probe"]["pdb"] is None:
            self._config["input"]["probe"]["pdb"] = self._config["input"]["probe"]["cid"] + ".pdb"
            
        # Resolve workdir path
        workdir = str(expandpath(Path(self._config["general"]["workdir"])))
        if not (workdir.startswith("/") or workdir.startswith("$HOME") or workdir.startswith("~")):
            workdir = str(self._yaml_dir / workdir)
        self._config["general"]["workdir"] = Path(workdir)
        
        # Resolve protein PDB path
        pdb_path = str(expandpath(Path(self._config["input"]["protein"]["pdb"])))
        if not (pdb_path.startswith("/") or pdb_path.startswith("$HOME") or pdb_path.startswith("~")):
            pdb_path = str(self._yaml_dir / pdb_path)
        self._config["input"]["protein"]["pdb"] = Path(pdb_path)
        
        # Resolve probe paths
        for path_key in ["pdb", "mol2"]:
            probe_path = str(expandpath(Path(self._config["input"]["probe"][path_key])))
            if not (probe_path.startswith("/") or probe_path.startswith("$HOME") or probe_path.startswith("~")):
                probe_path = str(self._yaml_dir / probe_path)
            self._config["input"]["probe"][path_key] = Path(probe_path)
            
        # Set default ssbond if not specified
        if "ssbond" not in self._config["input"]["protein"] or self._config["input"]["protein"]["ssbond"] is None:
            self._config["input"]["protein"]["ssbond"] = []
            
        # Store YAML path
        self._config["general"]["yaml"] = self._yaml_path
        
    @property
    def config(self) -> Dict[str, Any]:
        """Get complete configuration dictionary."""
        return self._config
        
    @property
    def general(self) -> GeneralConfig:
        """Get general configuration settings."""
        return cast(GeneralConfig, self._config["general"])
        
    @property
    def protein(self) -> ProteinConfig:
        """Get protein configuration settings."""
        return cast(ProteinConfig, self._config["input"]["protein"])
        
    @property
    def probe(self) -> ProbeConfig:
        """Get probe configuration settings."""
        return cast(ProbeConfig, self._config["input"]["probe"])
        
    @property
    def map(self) -> MapConfig:
        """Get map configuration settings."""
        return cast(MapConfig, self._config["map"])
        
    @property
    def msmd(self) -> MSMDConfig:
        """Get MSMD simulation settings."""
        return cast(MSMDConfig, self._config["exprorer_msmd"]["general"])


def parse_yaml(yamlpath: Path) -> Dict[str, Any]:
    """Parse YAML configuration file (wrapper for backward compatibility)."""
    config = Configuration(yamlpath)
    return config.config
