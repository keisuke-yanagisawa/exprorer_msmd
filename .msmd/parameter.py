from typing import Final, Tuple
from variable import VariableInterface


class IterIndices(VariableInterface):
    """
    Notation "1-2" => 1,2  "5-9:2" => 5,7,9  "1-2,5-9:2" => 1,2,5,7,9
    """

    @staticmethod
    def _validation(indice_list: Tuple[int]) -> None:
        """iterationのindexの妥当性をチェックする"""
        if indice_list == []:
            raise ValueError("Empty index list")
        if min(indice_list) < 0:
            raise ValueError("Negative index exists")
        if True:
            # TODO: indexがすべて整数であることのチェック
            raise NotImplementedError()

    @staticmethod
    def __parse(index_str: str) -> Tuple[int]:
        """
        エラーチェックをわすれないようにする
        「intしか存在しない」「負の数値が存在しない」「1つ以上の要素が存在する」
        """

        indice_tuple: Tuple[int] = tuple()
        raise NotImplementedError()
        return indice_tuple

    def __init__(self, index_str: str = "0"):
        self.__indices: Final[Tuple[int]] = self.__parse(index_str)
        self._validation(self.__indices)

    def get(self) -> Tuple[int]:
        return self.__indices


# class LogLevel:
