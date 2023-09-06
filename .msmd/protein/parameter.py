from ..variable import VariableInterface
from typing import Final, List, Tuple
from ..biopython.Bio.PDB import get_structure
from Bio.PDB.Structure import Structure
import itertools


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
    def _validation_struct(struct: Structure) -> None:
        if True:
            """alternative locationが存在してはいけない"""
            raise NotImplementedError()
            raise ValueError("The protein structure has alternative location")

    @staticmethod
    def _validation_ssbonds(struct: Structure, ssbonds: List[Tuple[ResidueID, ResidueID]]) -> None:
        if True:
            """ssbondsのresidueが存在しているかどうかのチェック"""
            raise NotImplementedError()
            for resid in itertools.chain(*ssbonds):
                raise NotImplementedError()
                raise ValueError(f"The residue {resid.get()} is not found in the protein structure")

    @staticmethod
    def __parse_ssbonds(ssbonds: List[Tuple[int, int]]) -> List[Tuple[ResidueID, ResidueID]]:
        """ssbondsの形式を変換する"""
        raise NotImplementedError()
        return ssbonds

    def __init__(self, pdb: str, ssbonds: List[Tuple[int, int]] = []):
        self.__struct: Final[Structure] = get_structure(pdb)
        self.__ssbonds: Final[List[Tuple[ResidueID, ResidueID]]] = self.__parse_ssbonds(ssbonds)
        self._validation_struct(self.__struct)
        self._validation_ssbonds(self.__struct, self.__ssbonds)

    # getterが必要
    # def get(self) -> Structure:
