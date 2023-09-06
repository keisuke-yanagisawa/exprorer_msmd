from ..variable import VariableInterface


class Temperature(VariableInterface):
    @staticmethod
    def _validation(temperature: int) -> None:
        if temperature < 0:
            raise ValueError(f"The temperature: {temperature} is negative")

    def __init__(self, temperature: int = 300):
        self.__temperature: int = temperature
        self._validation(self.__temperature)

    def get(self) -> int:
        return self.__temperature


class StepSize(VariableInterface):
    @staticmethod
    def _validation(stepsize: float) -> None:
        if stepsize <= 0:
            raise ValueError(f"The stepsize: {stepsize} is not positive")

    def __init__(self, stepsize: float = 0.002):
        self.__stepsize: float = stepsize
        self._validation(self.__stepsize)

    def get(self) -> float:
        return self.__stepsize


class Steps(VariableInterface):
    @staticmethod
    def _validation(steps: int) -> None:
        if steps <= 0:
            raise ValueError(f"The steps: {steps} is not positive")
        if isinstance(steps, int) is False:
            raise ValueError(f"The steps: {steps} is not integer")

    def __init__(self, steps: int = 1000000):
        self.__steps: int = steps
        self._validation(self.__steps)

    def get(self) -> int:
        return self.__steps
