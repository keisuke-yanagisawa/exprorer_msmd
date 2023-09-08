from ..variable import VariableInterface, Path
from typing import Final, List, Tuple


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
        """
        チェックすべき項目
        ・複数残基から構成されている
        ・alternative locationが存在してはいけない
        """
        raise NotImplementedError()

    @staticmethod
    def _validation_ssbonds(pdbpath: Path, ssbonds: List[Tuple[ResidueID, ResidueID]]) -> None:
        """
        チェックすべき項目
        ・ssbondsのresidueがproteinに存在している
        """
        raise NotImplementedError()

    @staticmethod
    def __parse_ssbonds(ssbonds: List[Tuple[int, int]]) -> List[Tuple[ResidueID, ResidueID]]:
        """ssbondsの形式を変換する"""
        raise NotImplementedError()
        return ssbonds

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
