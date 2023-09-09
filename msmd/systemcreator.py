import tempfile

from .executable.packmol import Packmol
from .system import System, PDBPreparator, SystemGenerator
from typing import Final
from .variable import Path, PDBString
from .unit import Angstrom
from scipy import constants
from .protein.parameter import Protein
from .probe.parameter import Probe
from . import config
from .config import InputConfig
from .executable.tleap import TLeap
import parmed as pmd
from .jinja2 import render_file


class MSMDSystemCreator:

    factor: Final[float] = 0.8
    """
    The generated system will shrink because of NPT ensemble.
    This factor is used to adjust the concentration change.
    """

    @staticmethod
    def __determine_system_size(parm7: Path, rst7: Path) -> Angstrom:
        """
        get longest box size from rst7 file
        """
        tmp = pmd.load_file(parm7.path, rst7.path)
        box_xyz = tmp.box_vectors[:3]  # in Angstrom
        size = max([box_xyz[0][0], box_xyz[1][1], box_xyz[2][2]]).value_in_unit(pmd.unit.angstrom)
        # boxは常にx, y, z座標系であると仮定している。斜めな軸を持つboxは考慮していない
        return Angstrom(size)

    # box_size_str = gop(f"tail -n 1 {rst7} | cut -c -36")
    # box_size = [float(s) for s in box_size_str.split()]
    # box_size = max(box_size)
    # return box_size

    @ classmethod
    def __calculate_system_size(cls, INPUT_CONFIG: InputConfig) -> Angstrom:
        tmp_parm7 = Path(tempfile.mkstemp(suffix=".parm7")[1])
        tmp_rst7 = Path(tempfile.mkstemp(suffix=".rst7")[1])
        tleap_source_file: Final[Path] = render_file("calculate_system_size.in", pdbfile=INPUT_CONFIG.PROTEIN.pdbpath.get(), tmp_parm7=tmp_parm7.get(), tmp_rst7=tmp_rst7.get())

        # TODO: EXECUTABLEの取得方法を考える
        TLeap(config.CONFIG.GENERAL.EXECUTABLE.TLEAP).execute(tleap_source_file)
        return cls.__determine_system_size(tmp_parm7, tmp_rst7)

    @ classmethod
    def __execute_packmol(cls, protein: Protein, probe: Probe, system_length: Angstrom) -> PDBString:
        # calculate specific number of probe molecules
        protein_volume: Final[float] = Protein.compute_volume(protein)
        system_volume: Final[float] = system_length.get() ** 3 - protein_volume
        num_probe_molecules: Final[int] = int(constants.N_A * probe.concentration.get() * system_volume * (10**-27) * cls.factor)

        packmol = Packmol(config.CONFIG.GENERAL.EXECUTABLE.PACKMOL)
        packmol.set(protein, probe, system_length, num_probe_molecules)
        box: Final[PDBString] = packmol.run()
        return box

    @ staticmethod
    def __execute_tleap(protein: Protein, probe: Probe, system_length: Angstrom) -> System:
        tmp_parm7 = Path(tempfile.mkstemp(suffix=".parm7")[1])
        tmp_rst7 = Path(tempfile.mkstemp(suffix=".rst7")[1])

        data = {
            "LIGAND_PARAM": f"leaprc.{probe.atomtype.get()}",
            "SS_BONDS": protein.ssbonds,
            "PROBE_ID": probe.cid,
            "PROBE_MOL2": probe.mol2path.get(),
            "OUT_PARM7": tmp_parm7.get(),
            "OUT_RST7": tmp_rst7.get(),
            "SYSTEM_PATH": protein.pdbpath.get(),
            "PROBE_FRCMOD": probe.frcmod.get(),
            "SIZE": system_length.get()
        }

        tleap_source_file: Final[Path] = render_file("leap.in", **data)
        tleap = TLeap(config.CONFIG.GENERAL.EXECUTABLE.TLEAP)
        tleap_output = tleap.execute(tleap_source_file)

        # System中のchargeがどうなっているか確認。
        final_charge_info = [s.strip() for s in tleap_output.split("\n")
                             if s.strip().startswith("Total unperturbed charge")][0]
        final_charge_value = float(final_charge_info.split()[-1])
        if final_charge_value != 0:
            raise ValueError(f"Final charge is {final_charge_info}. It should be 0.")

        out_system = SystemGenerator().from_amber(tmp_parm7, tmp_rst7)
        return out_system

    @ classmethod
    def create(cls, INPUT_CONFIG: InputConfig) -> System:

        system_length: Final[Angstrom] = cls.__calculate_system_size(INPUT_CONFIG)

        # prepare protein and probe molecules
        protein: Final[Protein] = INPUT_CONFIG.PROTEIN
        original_protein_str: Final[PDBString] = PDBPreparator.read_pdb_file(protein.pdbpath)
        prepared_protein_str: Final[PDBString] = PDBPreparator.prepare(original_protein_str)
        prepared_protein: Final[Protein] = Protein.from_string(prepared_protein_str, protein.ssbonds)
        probe: Final[Probe] = INPUT_CONFIG.PROBE

        # execute packmol
        probe.generate_frcmod(config.CONFIG.GENERAL.EXECUTABLE.PARMCHK)
        packmol_system_box: Final[Protein] = Protein.from_string(cls.__execute_packmol(prepared_protein, probe, system_length), protein.ssbonds)

        system: Final[System] = cls.__execute_tleap(packmol_system_box, probe, system_length)

        return system
