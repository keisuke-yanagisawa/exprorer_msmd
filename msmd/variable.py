import abc
from typing import Final, Optional
import subprocess
from .standard_library.os import getabsolutepath


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
    def _validation(path: str, ext: Optional[str] = None) -> None:
        if (ext is not None) and (path.endswith(ext) is False):
            raise ValueError(f"The path: {path} does not end with {ext}")
        # TODO: ファイルパスとして適切な文字列か否かのチェックが必要

    def __init__(self, path: str = ".", ext: Optional[str] = None):
        self.__path: Final[str] = getabsolutepath(path)
        self._validation(self.__path, ext)

    def get(self) -> str:
        return self.__path


class Name(VariableInterface):
    @staticmethod
    def _validation(name: str) -> None:
        # TODO: ファイル名の一部として適切な文字列か否かのチェックが必要
        pass

    def __init__(self, name: str):
        self.__name: Final[str] = name
        self._validation(self.__name)

    def get(self) -> str:
        return self.__name


class PDBString(VariableInterface):
    """PDB形式の文字列を表すクラス"""
    @staticmethod
    def _validation(pdbstring: str) -> None:
        pass

    def __init__(self, pdbstring: str):
        self.__pdbstring: Final[str] = pdbstring
        self._validation(self.__pdbstring)

    def get(self) -> str:
        return self.__pdbstring


class Command(VariableInterface):
    @staticmethod
    def _validation(command: str) -> None:
        which_returncode = subprocess.run(["which", command], capture_output=True).returncode
        if which_returncode != 0:
            raise ValueError(f"The command: {command} cannot be executed")

    def __init__(self, command: str = "python"):
        self.__command: Final[str] = command
        self._validation(self.__command)

    def get(self) -> str:
        return self.__command
