from typing import Final
from ..variable import VariableInterface


class Steps(VariableInterface):
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
