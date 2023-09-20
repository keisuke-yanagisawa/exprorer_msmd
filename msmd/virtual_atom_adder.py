from msmd.probe.parameter import Probe
from msmd.system import System


class VirtualAtomAdder:
    def __init__(self, sigma: float = 2, epsilon: float = 4.184e-6):
        # sigmaおよびepsilonは、VIS-VISのLJパラメータ。
        pass

    def add(self, system: System, target_probe: Probe) -> System:
        raise NotImplementedError()
