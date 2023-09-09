from typing import Final
from ..variable import VariableInterface


class PressureCoupling(VariableInterface):
    @staticmethod
    def _validation(pcoupl: str) -> None:
        if pcoupl not in ["no", "berendsen", "parrinello-rahman"]:
            raise ValueError(f"The pcoupl: {pcoupl} is not supported")

    def __init__(self, pcoupl: str):
        self.__pcoupl: Final[str] = pcoupl.lower()
        self._validation(self.__pcoupl)

    def get(self) -> str:
        return self.__pcoupl

    def __str__(self) -> str:
        return f"PressureCoupling('{self.__pcoupl}')"


class NumStep(VariableInterface):
    @staticmethod
    def _validation(steps: int) -> None:
        if steps <= 0:
            raise ValueError(f"The steps: {steps} is not positive")
        if isinstance(steps, int) is False:
            raise ValueError(f"The steps: {steps} is not integer")

    def __init__(self, steps: int = 1000000):
        self.__steps: Final[int] = steps
        self._validation(self.__steps)

    def get(self) -> int:
        return self.__steps

    def __str__(self):
        return f"NumStep('{self.__steps}')"
