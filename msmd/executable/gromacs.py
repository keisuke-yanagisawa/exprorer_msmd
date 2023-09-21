import tempfile
from typing import Final

from msmd.system import Trajectory

from .command import Executable, Command
from ..variable import Name, Path


class Gromacs(object):
    def __init__(self, exe: Executable, debug=False):
        self.__exe: Final[Executable] = exe
        self.__debug: Final[bool] = debug

    def make_ndx(self, gro: Path) -> None:

        self.__ndx = Path(tempfile.mkstemp(suffix=".ndx")[1])

        comm = Command(f"""{self.__exe} make_ndx -f {gro.get()} -o {self.__ndx.get()}<< EOF
                           q
                           EOF""")
        comm.run()

    def grompp(self, name: Name, input_mdp: Path, gro: Path, top: Path, prev_cpt: Path):
        self.__mdp = Path(tempfile.mkstemp(suffix=".mdp")[1])
        comm = Command(f"""{self.__exe} grompp \
                           -f {input_mdp} -po {self.__mdp} \
                           -o {name.get()}.tpr -c {gro.get()} -p {top.get()} -r {gro.get()} \
                           -n {self.__ndx} {prev_cpt.get()}""")
        comm.run()

    def mdrun(self, name: Name) -> Trajectory:
        finished_step_list_file = Path("finished_step_list")
        self.__cpt = Path(tempfile.mkstemp(suffix=".cpt")[1])
        comm = Command(f"""{self.__exe} mdrun \
                           -reprod -v -s {name.get()} -deffnm {name.get()} \
                           -cpi {self.__cpt.get()} \
                           && echo {name.get()} >> {finished_step_list_file.get()}""")
        comm.run()
        raise NotImplementedError()


"""{exe} make_ndx -f {gro} << EOF
q
EOF
"""

"""{exe} grompp -f {mdp}. -po {output_mdp} -o {stepname}.tpr -c {previous}.gro -p {top} -r {previous}.gro -n index.ndx {previous_cpt}"""
""""""
