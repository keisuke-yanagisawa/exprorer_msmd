from typing import Final, List, Union
from .variable import VariableInterface


class IterIndices(VariableInterface):
    """
    Notation "1-2" => 1,2  "5-9:2" => 5,7,9  "1-2,5-9:2" => 1,2,5,7,9
    """

    @staticmethod
    # Notation 1-2 => 1,2  5-9:2 => 5,7,9
    def __expand_index(ind_info: str) -> List[int]:
        ret: List[int] = []
        elems = ind_info.split(",")
        for elem in elems:
            if elem.find("-") == -1:
                elem = int(elem)
                ret.append(elem)
            elif elem.find(":") == -1:
                st, ed = elem.split("-")
                st, ed = int(st), int(ed)
                ret.extend(list(range(st, ed + 1)))
            else:
                window, offset = elem.split(":")
                st, ed = window.split("-")
                st, ed, offset = int(st), int(ed), int(offset)
                ret.extend(list(range(st, ed + 1, offset)))
        return ret

    @staticmethod
    def _validation(indice_list: List[int]) -> None:
        """iterationのindexの妥当性をチェックする"""
        if indice_list == []:
            raise ValueError("Empty index list")
        if min(indice_list) < 0:
            raise ValueError("Negative index exists")
        for vartype in [type(index) for index in indice_list]:
            if vartype != int:
                raise ValueError("Non-integer index exists")

    @classmethod
    def __parse(cls, index_str: str) -> List[int]:
        indice_tuple: List[int] = cls.__expand_index(index_str)
        return indice_tuple

    def __init__(self, indices: Union[str, int] = "0"):
        indices_str = str(indices)
        self.__indices: Final[List[int]] = self.__parse(indices_str)
        self._validation(self.__indices)

    def get(self) -> List[int]:
        return self.__indices

    def __str__(self) -> str:
        index_str = ", ".join([str(index) for index in self.__indices])
        return f"IterIndices('[{index_str}]')"

# class LogLevel:
