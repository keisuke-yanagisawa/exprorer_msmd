from typing import Dict, Final, Literal, Optional

from ..executable.parmchk import Parmchk
from ..variable import VariableInterface, Path
from ..unit import Molar
from ..executable.command import Executable


class AtomType(VariableInterface):
    at_indices: Dict[str, int] = {"gaff": 1, "gaff2": 2}

    @staticmethod
    def _validation(atomtype: Literal["gaff", "gaff2"]) -> None:
        if atomtype not in ["gaff", "gaff2"]:
            raise ValueError(f"The atomtype: {atomtype} is not supported")

    def __init__(self, atomtype: Literal["gaff", "gaff2"] = "gaff2"):
        self.__atomtype: Final[str] = atomtype
        self._validation(self.__atomtype)

    def get(self) -> str:
        return self.__atomtype

    def __str__(self) -> str:
        return f"AtomType('{self.__atomtype}')"

    @property
    def atomtype_index(self) -> int:
        return self.at_indices[self.__atomtype]


class Probe:
    @staticmethod
    def _validation(cid: str) -> None:
        """pdbファイルの中身が、単一の残基から構成されていることをチェック"""
        """probe濃度が極端に濃ければwarningを出す"""

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

    def generate_frcmod(self, parmchk: Executable) -> None:
        self.__frcmod = Parmchk(parmchk).set(self).run()

    def __init__(self, cid: str,
                 mol2_path: str,
                 pdb_path: str,
                 atomtype: AtomType,
                 concentration: Molar):
        self.__cid: Final[str] = cid
        self.__mol2_path: Final[Path] = Path(self.__generate_mol2_path(cid, mol2_path), "mol2")
        self.__pdb_path: Final[Path] = Path(self.__generate_pdb_path(cid, pdb_path), "pdb")
        self.__atomtype: Final[AtomType] = atomtype
        self.__concentration: Final[Molar] = concentration
        self.__frcmod: Optional[Path] = None

    @ property
    def concentration(self) -> Molar:
        return self.__concentration

    @ property
    def cid(self) -> str:
        return self.__cid

    @ property
    def mol2path(self) -> Path:
        return self.__mol2_path

    @ property
    def pdbpath(self) -> Path:
        return self.__pdb_path

    @property
    def atomtype(self) -> AtomType:
        return self.__atomtype

    @property
    def frcmod(self) -> Path:
        if self.__frcmod is None:
            raise ValueError("frcmod is not generated")
        return self.__frcmod
