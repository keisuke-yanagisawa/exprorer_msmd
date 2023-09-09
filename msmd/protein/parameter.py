import tempfile
from ..variable import PDBString, VariableInterface, Path
from typing import Final, List, Tuple
from ..biopython.Bio.PDB import get_structure
from Bio.PDB import SASA


class ResidueID(VariableInterface):
    @staticmethod
    def _validation(resid: int) -> None:
        """residue idが適切な値かどうかのチェック"""
        if resid <= 0:
            raise ValueError(f"The residue id: {resid} is not positive")

    def __init__(self, resid: int = 1):
        self.__resid: Final[int] = resid
        self._validation(self.__resid)

    def get(self) -> int:
        return self.__resid

    @property
    def resid(self) -> int:
        return self.__resid

    def __str__(self) -> str:
        return f"ResidueID({self.__resid})"


class SSBond(VariableInterface):
    @staticmethod
    def _validation(ssbond: Tuple[ResidueID, ResidueID]) -> None:
        if ssbond[0].get() == ssbond[1].get():
            raise ValueError(f"The same residue id: {ssbond[0]} is specified as a pair of ssbond residues")
        if ssbond[0].get() < 0:
            raise ValueError(f"The residue id: {ssbond[0]} is not positive")
        if ssbond[1].get() < 0:
            raise ValueError(f"The residue id: {ssbond[1]} is not positive")

    def __init__(self, resid_from: int, resid_to: int):
        self.__ssbond: Final[Tuple[ResidueID, ResidueID]] = (ResidueID(resid_from), ResidueID(resid_to))
        self._validation(self.__ssbond)

    def get(self) -> Tuple[ResidueID, ResidueID]:
        return self.__ssbond

    @property
    def ssbond(self) -> Tuple[ResidueID, ResidueID]:
        return self.__ssbond

    def __str__(self):
        return f"SSBond({self.__ssbond[0]}, {self.__ssbond[1]})"


class Protein:
    @staticmethod
    def _validation_struct(pdbpath: Path) -> None:
        struct = get_structure(pdbpath.get())
        if len([res for res in struct.get_residues()]) == 0:
            raise ValueError("The protein structure is empty")
        if len([res for res in struct.get_residues()]) == 1:
            raise ValueError("The protein structure has only one residue; is it a probe?")

        for res in struct.get_residues():
            atoms = [atom for atom in res.get_atoms()]
            atom_names = [atom.get_name() for atom in atoms]
            if len(atoms) != len(set(atom_names)):
                raise ValueError(f"The residue: {res} on the protein {pdbpath} has alternative location")

    @staticmethod
    def _validation_ssbonds(pdbpath: Path, ssbonds: List[SSBond]) -> None:
        struct = get_structure(pdbpath.get())
        residues = [res for res in struct.get_residues()]
        res_nums = [res.get_id()[1] for res in residues]
        res_num_set = set(res_nums)
        for ssbond in ssbonds:
            if ssbond.ssbond[0].get() not in res_num_set:
                raise ValueError(f"The residue id specified as a ssbond residue: {ssbond.ssbond[0]} is not in the protein {pdbpath}")
            if ssbond.ssbond[1].get() not in res_num_set:
                raise ValueError(f"The residue id specified as a ssbond residue: {ssbond.ssbond[1]} is not in the protein {pdbpath}")
            if ssbond.ssbond[0].get() == ssbond.ssbond[1].get():
                raise ValueError(f"The same residue id: {ssbond.ssbond[0]} is specified as a pair of ssbond residues")
            # TODO: ssbondに指定された残基はCYXになっているか

    def compute_volume(self) -> float:
        struct = get_structure(self.pdbpath.get())
        SASA.ShrakeRupley().compute(struct, level="S")  # calculate structure SASA
        return struct.sasa  # type: ignore
        # sasa is set by above compute() method

    def __init__(self, pdbpath: Path, ssbonds: List[SSBond] = []):
        self.__pdbpath: Final[Path] = pdbpath
        self.__ssbonds: Final[List[SSBond]] = ssbonds
        self._validation_struct(self.__pdbpath)
        self._validation_ssbonds(self.__pdbpath, self.__ssbonds)

    @property
    def pdbpath(self) -> Path:
        return self.__pdbpath

    @property
    def ssbonds(self) -> List[SSBond]:
        return self.__ssbonds

    @staticmethod
    def from_string(pdbstring: PDBString, ssbonds: List[SSBond] = []) -> "Protein":
        pdbpath: Final[Path] = Path(tempfile.mkstemp(suffix=".pdb")[1])
        with open(pdbpath.get(), "w") as f:
            f.write(pdbstring.get())
        return Protein(pdbpath, ssbonds)
