from typing import List, Tuple

from msmd.system import System

_default_restraints = [
    ("POSRES1000", 1000),
    ("POSRES500", 500),
    ("POSRES200", 200),
    ("POSRES100", 100),
    ("POSRES50", 50),
    ("POSRES20", 20),
    ("POSRES10", 10),
    ("POSRES0", 0)
]


class PositionRestraintAdder:
    def __init__(self, restraints: List[Tuple[str, float]] = _default_restraints):
        pass

    def add(self, system: System) -> System:
        raise NotImplementedError()
