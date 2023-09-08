from .system import System, PDBPreparator, FrcmodCreator
from typing import Final
from .variable import Path, PDBString
from .unit import Angstrom
from scipy import constants
from .protein.parameter import Protein
from .probe.parameter import Probe
from .config import InputConfig


class MSMDSystemCreator:

    factor: Final[float] = 0.8
    """
    The generated system will shrink because of NPT ensemble. 
    This factor is used to adjust the concentration change.
    """

    @staticmethod
    def __calculate_system_size(INPUT_CONFIG: InputConfig) -> Angstrom:
        raise NotImplementedError()

    @staticmethod
    def __calculate_protein_volume(protein: Protein) -> float:
        raise NotImplementedError()

    @classmethod
    def __execute_packmol(cls, protein, probe, system_length) -> PDBString:
        # calculate specific number of probe molecules
        protein_volume: Final[float] = cls.__calculate_protein_volume(protein)
        system_volume: Final[float] = system_length.get() ** 3 - protein_volume
        num_probe_molecules: Final[int] = int(constants.N_A * probe.concentration.get() * system_volume * (10**-27) * cls.factor)

        # packmolを実行する
        raise NotImplementedError()

    @staticmethod
    def __execute_tleap(packmol_box: PDBString, probe: Probe, frcmod: Path) -> System:
        raise NotImplementedError()

    @classmethod
    def create(cls, INPUT_CONFIG: InputConfig) -> System:
        # prepare protein and probe molecules
        protein: Final[Protein] = INPUT_CONFIG.PROTEIN
        original_protein: Final[PDBString] = PDBPreparator.read_pdb_file(protein.pdbpath)
        prepared_protein: Final[PDBString] = PDBPreparator.prepare(original_protein)
        probe: Final[Probe] = INPUT_CONFIG.PROBE

        # execute packmol
        system_length: Final[Angstrom] = cls.__calculate_system_size(INPUT_CONFIG)
        packmol_system_box: Final[PDBString] = cls.__execute_packmol(prepared_protein, probe, system_length)

        frcmod_path: Final[Path] = FrcmodCreator.create(probe)

        system: Final[System] = cls.__execute_tleap(packmol_system_box, probe, frcmod_path)

        return system
