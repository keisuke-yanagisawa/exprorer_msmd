from ..variable import VariableInterface, Path
from typing import Final, List, Tuple
from ..biopython.Bio.PDB import get_structure


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
    def _validation_ssbonds(pdbpath: Path, ssbonds: List[Tuple[ResidueID, ResidueID]]) -> None:
        struct = get_structure(pdbpath.get())
        residues = [res for res in struct.get_residues()]
        res_nums = [res.get_id()[1] for res in residues]
        res_num_set = set(res_nums)
        for ssbond in ssbonds:
            if ssbond[0].get() not in res_num_set:
                raise ValueError(f"The residue id specified as a ssbond residue: {ssbond[0]} is not in the protein {pdbpath}")
            if ssbond[1].get() not in res_num_set:
                raise ValueError(f"The residue id specified as a ssbond residue: {ssbond[1]} is not in the protein {pdbpath}")
            if ssbond[0].get() == ssbond[1].get():
                raise ValueError(f"The same residue id: {ssbond[0]} is specified as a pair of ssbond residues")
            # TODO: ssbondに指定された残基はCYXになっているか

    @staticmethod
    def __parse_ssbonds(ssbonds: List[Tuple[int, int]]) -> List[Tuple[ResidueID, ResidueID]]:
        ret_ssbonds: List[Tuple[ResidueID, ResidueID]] = []
        for ssbond in ssbonds:
            ret_ssbonds.append((ResidueID(ssbond[0]), ResidueID(ssbond[1])))
        return ret_ssbonds

    def __init__(self, pdbpath: Path, ssbonds: List[Tuple[int, int]] = []):
        self.__pdbpath: Final[Path] = pdbpath
        self.__ssbonds: Final[List[Tuple[ResidueID, ResidueID]]] = self.__parse_ssbonds(ssbonds)
        self._validation_struct(self.__pdbpath)
        self._validation_ssbonds(self.__pdbpath, self.__ssbonds)

    @property
    def pdbpath(self) -> Path:
        return self.__pdbpath
    # getterが必要
    # def get(self) -> Structure:
