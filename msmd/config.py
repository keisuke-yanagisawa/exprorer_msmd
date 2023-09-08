import os
import yaml
from typing import Final, Dict, Any, List, Tuple
from .variable import Command, Path, Name
from .parameter import IterIndices
from .protein.parameter import Protein
from .probe.parameter import Probe, AtomType
from .simulation import interface
from .simulation.minimization import MinimizationStep
from .simulation.heating import HeatingStep
from .simulation.equilibration import EquilibrationStep
from .simulation.production import ProductionStep
from .simulation.sequence import SimulationSequence
from .analysis.analysis_parameter import CpptrajSelector, PMAPNormalization
from .unit import Angstrom, Molar


class _ExecutableConfig:
    def __init__(self, executable_map: Dict[str, str]):
        self.PYTHON: Final[Command] = Command(executable_map["python"])
        self.GROMACS: Final[Command] = Command(executable_map["gromacs"])
        self.PACKMOL: Final[Command] = Command(executable_map["packmol"])
        self.TLEAP: Final[Command] = Command(executable_map["tleap"])
        self.CPPTRAJ: Final[Command] = Command(executable_map["cpptraj"])


class _GeneralConfig:
    def __init__(self, general_map: Dict[str, Any]):
        self.ITER_INDICES: Final[IterIndices] = IterIndices(general_map["iter_index"])
        self.WORKING_DIRECTORY: Final[Path] = Path(general_map["workdir"])
        self.PROJECT_NAME: Final[Name] = Name(general_map["name"])
        self.EXECUTABLE: Final[_ExecutableConfig] = _ExecutableConfig(general_map["executables"])


class InputConfig:

    @staticmethod
    def __fill_default(input_map: Dict[str, Any]) -> Dict[str, Any]:
        ret_map = input_map

        ret_map["protein"]["ssbond"] = ret_map["protein"].get("ssbond", None)
        ret_map["probe"]["mol2"] = ret_map["probe"].get("mol2", None)
        ret_map["probe"]["pdb"] = ret_map["probe"].get("pdb", None)

        if ret_map["protein"]["ssbond"] is None:
            ret_map["protein"]["ssbond"] = []

        return ret_map

    @staticmethod
    def __parse_protein(protein_map: Dict[str, Any]) -> Protein:
        __pdb: Path = Path(protein_map["pdb"], "pdb")
        __ssbond: List[Tuple[int, int]] = protein_map["ssbond"]
        return Protein(__pdb, __ssbond)

    @staticmethod
    def __parse_probe(probe_map: Dict[str, Any]) -> Probe:
        __cid = probe_map["cid"]
        __mol2 = probe_map["mol2"]
        __pdb = probe_map["pdb"]
        __atom_type = AtomType(probe_map["atomtype"])
        __concentration = Molar(probe_map["molar"])
        return Probe(__cid, __mol2, __pdb, __atom_type, __concentration)

    def __init__(self, input_map: Dict[str, Any]):
        tmp_map = self.__fill_default(input_map)
        self.PROTEIN: Final[Protein] = self.__parse_protein(tmp_map["protein"])
        self.PROBE: Final[Probe] = self.__parse_probe(tmp_map["probe"])


class _SimulationStepFactory:
    @staticmethod
    def __merge_config(step_config: Dict[str, Any], default_parameters: Dict[str, Any]) -> Dict[str, Any]:
        default_dt = default_parameters.get("dt", 0.002)
        default_temp = default_parameters.get("temperature", 300)
        default_pressure = default_parameters.get("pressure", 1.0)
        default_pbc = default_parameters.get("pbc", "xyz")

        ret_map = step_config
        ret_map["dt"] = ret_map.get("dt", default_dt)
        ret_map["temperature"] = ret_map.get("temperature", default_temp)
        ret_map["pressure"] = ret_map.get("pressure", default_pressure)
        ret_map["pbc"] = ret_map.get("pbc", default_pbc)
        return ret_map

    def __init__(self, default_parameters: Dict[str, Any]):
        self.__default_parameters: Final[Dict[str, Any]] = default_parameters

    def create(self, step_config: Dict[str, Any]) -> interface.SimulationInterface:
        config = self.__merge_config(step_config, self.__default_parameters)
        if config["type"] == "minimization":
            return MinimizationStep(config)
        if config["type"] == "heating":
            return HeatingStep(config)
        if config["type"] == "equilibration":
            return EquilibrationStep(config)
        if config["type"] == "production":
            return ProductionStep(config)
        raise ValueError(f"The step_type: {config['type']} is not supported")


class _SimulationConfig:
    @staticmethod
    def __parse_sequence(sequence_map: List[Dict[str, Any]], default_parameters: Dict[str, Any]) -> SimulationSequence:
        factory = _SimulationStepFactory(default_parameters)
        simulation_seq: List[SimulationInterface] = []
        for step_config in sequence_map:
            simulation_seq.append(factory.create(step_config))
        return SimulationSequence(simulation_seq)

    def __init__(self, simulation_map: Dict[str, Any]):
        self.TITLE: Final[str] = simulation_map["title"]
        self.SEQUENCE: Final[SimulationSequence] = self.__parse_sequence(simulation_map["sequence"], simulation_map["general"])


class _PMAPTypeConfig:
    def __init__(self, pmap_type_map: Dict[str, Any]):
        self.SUFFIX: Final[Name] = Name(pmap_type_map["suffix"])
        self.SELECTOR: Final[CpptrajSelector] = CpptrajSelector(pmap_type_map["selector"])
        self.NORMALIZATION: Final[PMAPNormalization] = PMAPNormalization(pmap_type_map["normalization"])
        self.MAP_SIZE: Final[Angstrom] = Angstrom(pmap_type_map["map_size"])


class _PMAPConfig:

    @staticmethod
    def __parse_pmap_types(pmap_map: Dict[str, Any]) -> List[_PMAPTypeConfig]:
        __default_normalization = pmap_map.get("normalization", "snapshot")
        __default_map_size = pmap_map.get("map_size", 80)  # TODO: 自動的に判断されるべき

        ret_list: List[_PMAPTypeConfig] = []
        for pmap_type in pmap_map["maps"]:
            pmap_type["normalization"] = pmap_map.get("normalization", __default_normalization)
            pmap_type["map_size"] = pmap_type.get("map_size", __default_map_size)
            ret_list.append(_PMAPTypeConfig(pmap_type))
        return ret_list

    def __init__(self, pmap_map: Dict[str, Any]):
        self.TYPE: Final[str] = pmap_map["type"]
        # self.SNAPSHOT: Final[Frame] = Frame(pmap_map["snapshot"]) # TODO: フォーマットがどうなっているか確認
        self.MAPS: Final[List[_PMAPTypeConfig]] = self.__parse_pmap_types(pmap_map)
        self.MAPSIZE: Final[Angstrom] = Angstrom(pmap_map["map_size"])  # TODO: map_sizeは本来PMAPTypeConfigに含まれるべき。
        self.NORMALIZATION: Final[PMAPNormalization] = PMAPNormalization(pmap_map["normalization"])  # TODO: normalizationは本来PMAPTypeConfigに含まれるべき。


class _InverseMSMDTypeConfig:
    def __init__(self, type_map: Dict[str, Any]):
        self.NAME: Final[Name] = Name(type_map["name"])
        self.ATOMS: Final[List[Tuple[str, str]]] = type_map["atoms"]


class _InverseMSMDConfig:
    @staticmethod
    def __parse_types(type_list: List[Dict[str, Any]]) -> List[_InverseMSMDTypeConfig]:
        ret_list: List[_InverseMSMDTypeConfig] = []
        for type_map in type_list:
            ret_list.append(_InverseMSMDTypeConfig(type_map))
        return ret_list

    def __init__(self, profile_map: Dict[str, Any]):
        self.TYPES: Final[List[_InverseMSMDTypeConfig]] = self.__parse_types(profile_map["profile"]["types"])


class __Config:
    def __init__(self, config_map: Dict[str, Any]):
        self.GENERAL: Final[_GeneralConfig] = _GeneralConfig(config_map["general"])
        self.INPUT: Final[InputConfig] = InputConfig(config_map["input"])
        self.SIMULATION: Final[_SimulationConfig] = _SimulationConfig(config_map["exprorer_msmd"])
        self.PMAP: Final[_PMAPConfig] = _PMAPConfig(config_map["map"])
        self.INVERSE_MSMD: Final[_InverseMSMDConfig] = _InverseMSMDConfig(config_map["probe_profile"])


def __load_yaml(path: str) -> dict:
    if os.path.exists(path) is False:
        raise FileNotFoundError(f"The file: {path} is not found")
    with open(path, "r") as f:
        return yaml.safe_load(f)


CONFIG: __Config
""" This global constant is initialized only by load_config()"""


def load_config(yaml: str) -> None:
    global CONFIG
    CONFIG = __Config(__load_yaml(yaml))
