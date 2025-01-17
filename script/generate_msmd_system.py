#!/usr/bin/python3

import os
import tempfile
from pathlib import Path
from subprocess import getoutput as gop
from typing import Dict, List, Literal, Optional, Tuple, TypedDict, Union, cast

from script.utilities import const

class ProteinSettings(TypedDict):
    pdb: Union[str, Path]
    ssbond: List[str]

class ProbeSettings(TypedDict):
    mol2: Union[str, Path]
    pdb: Union[str, Path]
    cid: str
    atomtype: Literal["gaff", "gaff2"]
    molar: float

class InputSettings(TypedDict):
    protein: ProteinSettings
    probe: ProbeSettings
from script.utilities.executable import Packmol, Parmchk, TLeap
from script.utilities.logger import logger

VERSION = "2.0.0"

tmp_leap = """
source leaprc.protein.ff14SB
source leaprc.water.tip3p

prot = loadPDB {pdbfile}

addIons2 prot Na+ 0
addIons2 prot Cl- 0
solvateBox prot TIP3PBOX 10.0

center prot
charge prot

saveAmberParm prot {tmp_prefix}.parm7 {tmp_prefix}.rst7
quit
"""


def protein_pdb_preparation(pdbfile: Path) -> Path:
    """
    Remove OXT and ANISOU from pdbfile to avoid tleap error.
    
    Args:
        pdbfile: Path to PDB file
        
    Returns:
        Path to cleaned PDB file without OXT and ANISOU records
        
    Raises:
        FileNotFoundError: If input PDB file does not exist
    """
    tmp1 = Path(tempfile.mkstemp(prefix=const.TMP_PREFIX, suffix=const.EXT_PDB)[1])
    gop(f"grep -v OXT {pdbfile} | grep -v ANISOU > {tmp1}")
    return tmp1


def __calculate_boxsize(pdbfile: Path) -> float:
    tmpdir = tempfile.mkdtemp()
    tmp_prefix = f"{tmpdir}/{const.TMP_PREFIX}"
    with open(f"{tmp_prefix}.in", "w") as fout:
        fout.write(tmp_leap.format(pdbfile=str(pdbfile), tmp_prefix=tmp_prefix))
    logger.info(gop(f"tleap -f {tmp_prefix}.in | tee {tmp_prefix}.in.result"))

    try:
        size = calculate_boxsize(Path(f"{tmp_prefix}.rst7"))
    except ValueError as e:
        logger.error("====tleap input commands====")
        logger.error(gop(f"cat {tmp_prefix}.in"))
        logger.error("====tleap output====")
        logger.error(gop(f"cat {tmp_prefix}.in.result"))
        logger.error("failed to calculate box size: tleap error")
        raise e
    return size


def calculate_boxsize(rst7: Path) -> float:
    """
    get longest box size from rst7 file
    """

    box_size_str = gop(f"tail -n 1 {rst7} | cut -c -36")
    box_size = [float(s) for s in box_size_str.split()]
    box_size = max(box_size)
    return box_size


def _create_frcmod(
    mol2file: Path,
    atomtype: Literal["gaff", "gaff2"],
    debug: bool = False
) -> Path:
    """
    create frcmod file from mol2 file
    -----
    input
        mol2file: path to mol2 file
        atomtype: atom type (GAFF / GAFF2)
    output
        cfrcmod: path to frcmod file
    """
    cfrcmod = Path(tempfile.mkstemp(suffix=".frcmod")[1])
    Parmchk(debug=debug).set(mol2file, atomtype).run(frcmod=cfrcmod)
    return cfrcmod


def create_system(
    setting_protein: ProteinSettings,
    setting_probe: ProbeSettings,
    probe_frcmod: Path,
    debug: bool = False,
    seed: int = -1
) -> Tuple[Path, Path]:
    """
    create system from protein and probe
    """
    pdbpath = protein_pdb_preparation(Path(setting_protein["pdb"]))
    boxsize = __calculate_boxsize(pdbpath)
    ssbonds = setting_protein["ssbond"]
    cmol = Path(setting_probe["mol2"])
    cpdb = Path(setting_probe["pdb"])
    cid = setting_probe["cid"]
    atomtype = setting_probe["atomtype"]
    probemolar = float(setting_probe["molar"])

    box_pdb = Path(tempfile.mkstemp(suffix=".pdb")[1])
    Packmol(debug=debug).set(pdbpath, cpdb, boxsize, probemolar).run(box_pdb, seed=seed)

    tleap_obj = TLeap(debug=debug).set(cid, cmol, probe_frcmod, box_pdb, boxsize, ssbonds, atomtype)

    while True:
        _, fileprefix = tempfile.mkstemp(suffix="")
        os.system(f"rm {fileprefix}")
        tleap_obj.run(fileprefix)
        system_charge = tleap_obj._final_charge_value

        if system_charge == 0:
            break
        else:
            logger.warn("the system is not neutral. generate system again")

    return tleap_obj.parm7, tleap_obj.rst7


class MSMDSystemGenerator:
    """Generator class for Mixed-Solvent Molecular Dynamics systems."""
    
    def __init__(self, setting: Dict[str, InputSettings]) -> None:
        """Initialize the MSMD system generator.
        
        Args:
            setting: Dictionary containing system configuration
        """
        self.setting = setting
        self._protein_settings = cast(ProteinSettings, setting["input"]["protein"])
        self._probe_settings = cast(ProbeSettings, setting["input"]["probe"])
        
    def _prepare_protein_pdb(self, pdbfile: Path) -> Path:
        """Prepare protein PDB file by removing OXT and ANISOU records.
        
        Args:
            pdbfile: Path to input PDB file
            
        Returns:
            Path to cleaned PDB file
        """
        return protein_pdb_preparation(pdbfile)
        
    def _calculate_box_size(self, pdbfile: Path) -> float:
        """Calculate simulation box size from protein PDB.
        
        Args:
            pdbfile: Path to protein PDB file
            
        Returns:
            Box size in Angstroms
        """
        return __calculate_boxsize(pdbfile)
        
    def _create_frcmod(self, mol2file: Path, atomtype: Literal["gaff", "gaff2"], debug: bool = False) -> Path:
        """Create force field modification file.
        
        Args:
            mol2file: Path to mol2 file
            atomtype: Atom type (GAFF / GAFF2)
            debug: Enable debug output
            
        Returns:
            Path to generated frcmod file
        """
        return _create_frcmod(mol2file, atomtype, debug)
        
    def create_system(self, debug: bool = False, seed: int = -1) -> Tuple[Path, Path]:
        """Create MSMD system from current settings.
        
        Args:
            debug: Enable debug output
            seed: Random seed for system generation
            
        Returns:
            Tuple containing paths to (parm7, rst7) files
        """
        cfrcmod = self._create_frcmod(
            Path(self._probe_settings["mol2"]),
            self._probe_settings["atomtype"],
            debug=debug
        )
        return create_system(
            self._protein_settings,
            self._probe_settings,
            cfrcmod,
            debug=debug,
            seed=seed
        )

def generate_msmd_system(
    setting: Dict[str, InputSettings],
    debug: bool = False,
    seed: int = -1
) -> Tuple[Path, Path]:
    """Generate MSMD system from settings.
    
    This is a backward-compatible wrapper around MSMDSystemGenerator.
    For new code, prefer using MSMDSystemGenerator directly.
    
    Args:
        setting: Dictionary containing input settings
        debug: Enable debug output
        seed: Random seed for system generation
        
    Returns:
        Tuple containing paths to parm7 and rst7 files
    """
    generator = MSMDSystemGenerator(setting)
    return generator.create_system(debug=debug, seed=seed)
