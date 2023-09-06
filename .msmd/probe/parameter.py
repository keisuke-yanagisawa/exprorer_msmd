from typing import Literal, Optional
from ..variable import VariableInterface
from ..parameter import Path


class Probe:
    @staticmethod
    def __generate_mol2_path(cid, mol2_path) -> str:
        if mol2_path is None:
            return f"{cid}.mol2"
        return mol2_path

    @staticmethod
    def __generate_pdb_path(cid, pdb_path) -> str:
        if pdb_path is None:
            return f"{cid}.pdb"
        return pdb_path

    def __init__(self, cid: str = "A11", mol2_path: Optional[str] = None, pdb_path: Optional[str] = None):
        self.__cid: str = cid
        self.__mol2_path: Path = Path(self.__generate_mol2_path(cid, mol2_path))
        self.__pdb_path: Path = Path(self.__generate_pdb_path(cid, pdb_path))


class AtomType(VariableInterface):
    @staticmethod
    def _validation(atomtype: Literal["gaff", "gaff2"]) -> None:
        if atomtype not in ["gaff", "gaff2"]:
            raise ValueError(f"The atomtype: {atomtype} is not supported")

    def __init__(self, atomtype: Literal["gaff", "gaff2"] = "gaff2"):
        self.__atomtype: str = atomtype
        self._validation(self.__atomtype)

    def get(self) -> str:
        return self.__atomtype


class Molar(VariableInterface):
    @staticmethod
    def _validation(molar: float) -> None:
        if molar < 0:
            raise ValueError(f"The molar: {molar} is negative")

    def __init__(self, molar: float = 0.25):
        self.__molar: float = molar
        self._validation(self.__molar)

    def get(self) -> float:
        return self.__molar
