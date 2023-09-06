from typing import Final, Tuple
from .variable import VariableInterface


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


class Path(VariableInterface):
    @staticmethod
    def _validation(path: str) -> None:
        """ファイルパスとして適切な文字列か否かのチェックが必要"""
        raise NotImplementedError()
        if True:
            raise ValueError(f"The path: {path} cannot be used as path")

    def __init__(self, path: str = "."):
        self.__path: Final[str] = path
        self._validation(self.__path)

    def get(self) -> str:
        return self.__path


class ProjectName(VariableInterface):
    @staticmethod
    def _validation(name: str) -> None:
        """ファイル名の一部として適切な文字列か否かのチェックが必要"""
        raise NotImplementedError()
        if True:
            raise ValueError(f"The project name: {name} cannot be used as path")

    def __init__(self, name: str = "TEST_PROJECT"):
        self.__name: Final[str] = name
        self._validation(self.__name)

    def get(self) -> str:
        return self.__name


class Command(VariableInterface):
    @staticmethod
    def _validation(command: str) -> None:
        """`which` コマンドで実行可能かどうかをチェック"""
        raise NotImplementedError()
        if True:
            raise ValueError(f"The command: {command} cannot be executed")

    def __init__(self, command: str = "python"):
        self.__command: Final[str] = command
        self._validation(self.__command)

    def get(self) -> str:
        return self.__command

# class LogLevel:
