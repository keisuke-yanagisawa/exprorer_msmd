import subprocess
from ..variable import VariableInterface
from ..standard_library.logging.logger import logger


class Executable(VariableInterface):
    @staticmethod
    def _validation(command: str) -> None:
        # subprocess.run raises CalledProcessError if the command is not found
        subprocess.run(["which", command], capture_output=True, check=True)

    def __init__(self, command: str = "python"):
        self.__command: Final[str] = command
        self._validation(self.__command)

    def __str__(self) -> str:
        return f"Executable('{self.__command}')"

    def get(self) -> str:
        return self.__command


class Command(VariableInterface):
    @staticmethod
    def _validation(comm: str) -> None:
        pass

    def __init__(self, comm: str):
        self.__comm = comm
        self._validation(self.__comm)

    def run(self) -> str:
        res = subprocess.run(self.__comm, shell=True, stdout=subprocess.PIPE, check=True)
        return res.stdout.decode("utf-8")

    def __str__(self) -> str:
        return f"Command('{self.__comm}')"

    def get(self) -> str:
        return self.__comm
