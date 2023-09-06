import abc
from typing import Final


class VariableInterface(abc.ABC):
    @staticmethod
    @abc.abstractmethod
    def _validation(var):
        pass

    @abc.abstractmethod
    def get(self):
        pass


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


class Name(VariableInterface):
    @staticmethod
    def _validation(name: str) -> None:
        """ファイル名の一部として適切な文字列か否かのチェックが必要"""
        raise NotImplementedError()
        if True:
            raise ValueError(f"The name: {name} cannot be used as path")

    def __init__(self, name: str):
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
