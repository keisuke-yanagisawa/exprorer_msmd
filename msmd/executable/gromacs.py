import os
import shutil
import tempfile
from typing import Final, Tuple
from msmd.system import System, Trajectory
from .command import Executable, Command
from ..variable import Name, Path


class Gromacs(object):
    def __init__(self, debug=False):
        from msmd.config import CONFIG
        self.__exe: Final[Executable] = CONFIG.GENERAL.EXECUTABLE.GROMACS
        self.__debug: Final[bool] = debug

    def run(self, name: Name, system: System, input_mdp: Path) -> Trajectory:
        self.__input_mdp = input_mdp
        self.__ndx = self.__make_ndx(system.gro)
        self.__tpr = self.__grompp(self.__input_mdp, system.gro, system.top, self.__ndx)
        cpt_path, edr_path, gro_path, log_path, xtc_path = self.__mdrun(name, self.__tpr)
        return Trajectory(top=system.top, gro=gro_path,
                          trj=xtc_path, edr=edr_path,
                          log=log_path, cpt=cpt_path)

    def save(self, basedirpath: Path, prefix: Name) -> None:
        path_prefix: Path = basedirpath + prefix
        os.makedirs(basedirpath.get(), exist_ok=True)
        shutil.copy(self.__ndx.get(), path_prefix.get() + ".ndx")
        shutil.copy(self.__tpr.get(), path_prefix.get() + ".tpr")
        shutil.copy(self.__input_mdp.get(), path_prefix.get() + ".mdp")

    def __make_ndx(self, gro: Path) -> Path:
        ndx_path = Path(tempfile.mkstemp(suffix=".ndx")[1])
        comm = Command(f"""{self.__exe.get()} make_ndx -f {gro.get()} -o {ndx_path.get()}<< EOF
                           q
                           EOF""")
        comm.run()
        return ndx_path

    def __grompp(self, input_mdp: Path, gro: Path, top: Path, ndx: Path) -> Path:
        """return: tpr_path"""
        tpr_path = Path(tempfile.mkstemp(suffix=".tpr")[1])
        comm = Command(f"""{self.__exe.get()} grompp \
                           -f {input_mdp.get()} -po /dev/null \
                           -o {tpr_path.get()} -c {gro.get()} -p {top.get()} -r {gro.get()} \
                           -n {ndx.get()}""")
        comm.run()
        return tpr_path

    def __mdrun(self, name: Name, tpr: Path) -> Tuple[Path, Path, Path, Path, Path]:
        """return (cpt_path, edr_path, gro_path, log_path, xtc_path)"""
        finished_step_list_file = Path("finished_step_list")  # これの処理どうする？Command自体がCalledProcessErrorを投げてくれるので、外で処理できる

        cpt_path = Path(tempfile.mkstemp(suffix=".cpt")[1])
        edr_path = Path(tempfile.mkstemp(suffix=".edr")[1])
        gro_path = Path(tempfile.mkstemp(suffix=".gro")[1])
        log_path = Path(tempfile.mkstemp(suffix=".log")[1])
        xtc_path = Path(tempfile.mkstemp(suffix=".xtc")[1])

        comm = Command(f"""{self.__exe.get()} mdrun -reprod -v -s {tpr.get()} \
                           -o /dev/null -cpo {cpt_path.get()} -x {xtc_path.get()} -c {gro_path.get()} \
                           -e {edr_path.get()} -g {log_path.get()} \
                           && echo {name.get()} >> {finished_step_list_file.get()}""")
        comm.run()
        return (cpt_path, edr_path, gro_path, log_path, xtc_path)

    # TODO: create_pdb は Gromacsの計算実施とは別軸の話なので、同じクラスに入っているのは違和感。
    def create_pdb(self, gro: Path) -> Path:
        pdb_path = Path(tempfile.mkstemp(suffix=".pdb")[1])
        comm = Command(f"""{self.__exe.get()} trjconv -s {gro.get()} \
                           -f {gro.get()} -o {pdb_path.get()}.pdb <<EOF
                            0
                            EOF""")
        comm.run()
        return pdb_path
