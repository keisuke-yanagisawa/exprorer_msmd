from typing import Final
from .command import Executable, Command
from ..standard_library.logging.logger import logger
from ..variable import Path


class TLeap(object):
    def __init__(self, exe: Executable, debug=False):
        self.__exe: Final[Executable] = exe
        self.__debug: Final[bool] = debug

    def execute(self, sourcefile: Path) -> str:
        """
        tleapを実行し、tleapの出力を文字列で返す
        """
        command = Command(f"{self.__exe.get()} -f {sourcefile.get()}")
        logger.debug(command)
        output = command.run()
        logger.info(output)
        return output
