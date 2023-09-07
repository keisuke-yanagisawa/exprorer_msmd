import abc
from typing import Optional
from .variable import Path


class SystemInterface(abc.ABC):
    @abc.abstractmethod
    def get_system(self) -> "System":
        pass


class System(SystemInterface):
    def get_system(self) -> "System":
        return self


class SystemString(SystemInterface):
    def get_system(self) -> System:
        raise NotImplementedError()


class Trajectory(SystemInterface):
    def get_system(self) -> System:
        raise NotImplementedError()

    def save(self, path: Optional[Path]) -> None:
        """トラジェクトリファイルを陽な場所に保存する。"""
        pass
