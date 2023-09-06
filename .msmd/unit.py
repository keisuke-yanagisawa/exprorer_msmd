from variable import VariableInterface
from typing import Final


class Molar(VariableInterface):
    @staticmethod
    def _validation(molar: float) -> None:
        if molar < 0:
            raise ValueError(f"The molar: {molar} is negative")

    def __init__(self, molar: float = 0.25):
        self.__molar: Final[float] = molar
        self._validation(self.__molar)

    def get(self) -> float:
        return self.__molar


class PicoSecond(VariableInterface):
    @staticmethod
    def _validation(ps: float) -> None:
        if ps <= 0:
            raise ValueError(f"The time length: {ps} ps is not positive")

    def __init__(self, ps: float = 0.002):
        self.__ps: Final[float] = ps
        self._validation(self.__ps)

    def get(self) -> float:
        return self.__ps


class Kelvin(VariableInterface):
    @staticmethod
    def _validation(kelvin: float) -> None:
        if kelvin < 0:
            raise ValueError(f"The temperature: {kelvin} K is negative")

    def __init__(self, kelvin: float = 300.0):
        self.__kelvin: Final[float] = kelvin
        self._validation(self.__kelvin)

    def get(self) -> float:
        return self.__kelvin


class Angstrom(VariableInterface):
    @staticmethod
    def _validation(distance: float) -> None:
        if distance <= 0:
            raise ValueError(f"The distance: {distance} angstrom is not positive")

    def __init__(self, distance: float):
        self.__distance: Final[float] = distance
        self._validation(self.__distance)

    def get(self) -> float:
        return self.__distance
