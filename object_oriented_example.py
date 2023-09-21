from msmd.positon_restraint_adder import PositionRestraintAdder
from msmd.systemcreator import MSMDSystemCreator
from msmd.system import System, Trajectory
from typing import Final
from msmd import config
from msmd.standard_library.logging.logger import logger
from msmd.variable import Name, Path
from msmd.virtual_atom_adder import VirtualAtomAdder


config.load_config("/workspaces/exprorer_msmd/example/debug_protocol.yaml")
print(config.CONFIG)
logger.setLevel("DEBUG")

for iter_index in config.CONFIG.GENERAL.ITER_INDICES:
    output_directory: Path = config.CONFIG.GENERAL.WORKING_DIRECTORY + Name(f"system{iter_index}")

    sys: System = MSMDSystemCreator.create(config.CONFIG.INPUT)
    sys = VirtualAtomAdder().add(sys, config.CONFIG.INPUT.PROBE)
    sys = PositionRestraintAdder().add(sys, config.CONFIG.INPUT.PROBE)
    sys.save(output_directory + Name("prep"), config.CONFIG.GENERAL.PROJECT_NAME)

    simulation_output_directory: Path = output_directory + Name("simulation")
    traj: Final[Trajectory] = config.CONFIG.SIMULATION.SEQUENCE.run(sys, outdir=simulation_output_directory)
    # うーん。CONFIGをそのままrunできるのは気持ち悪いな
