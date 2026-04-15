import os
import tempfile
from typing import List

import numpy as np
import parmed as pmd
from scipy import constants

_NM_PER_ANGSTROM = 0.1
_ANGSTROM_PER_NM = 10.0


class GroAtom:
    """Atom data compatible with GRO format. Coordinates in nanometers."""

    def __init__(self):
        self.resi = -1
        self.resn = ""
        self.atomtype = ""
        self.atom_id = -1
        self.point = np.zeros((3,))
        self.comment = ""
        self.atomic_mass = 0.0

    @classmethod
    def from_parmed(cls, atom: pmd.Atom) -> "GroAtom":
        """Create GroAtom from a ParmEd Atom. Converts Angstrom to nm."""
        ga = cls()
        ga.resi = atom.residue.number
        ga.resn = atom.residue.name
        ga.atomtype = atom.name
        ga.atom_id = atom.idx + 1
        ga.point = np.array([atom.xx, atom.xy, atom.xz]) * _NM_PER_ANGSTROM
        ga.atomic_mass = atom.mass
        return ga


class Gro:
    """GRO file handler backed by ParmEd.

    Coordinates and box dimensions are stored in nanometers (GRO convention).
    ParmEd handles file I/O, converting between Angstrom and nm internally.
    """

    def __init__(self, path=""):
        self._structure: pmd.Structure = pmd.Structure()
        self._description: str = ""
        if path != "":
            self.parse(path)

    @property
    def description(self) -> str:
        return self._description

    @description.setter
    def description(self, value: str):
        self._description = value

    @property
    def natoms(self) -> int:
        return len(self._structure.atoms)

    @property
    def box_size(self) -> List[float]:
        if self._structure.box is not None:
            return [float(b * _NM_PER_ANGSTROM) for b in self._structure.box[:3]]
        return [0.0, 0.0, 0.0]

    @box_size.setter
    def box_size(self, value: List[float]):
        self._structure.box = np.array(
            [value[0] * _ANGSTROM_PER_NM, value[1] * _ANGSTROM_PER_NM, value[2] * _ANGSTROM_PER_NM, 90.0, 90.0, 90.0]
        )

    @property
    def atoms(self) -> List[GroAtom]:
        return [GroAtom.from_parmed(a) for a in self._structure.atoms]

    def parse(self, path):
        """Parse GRO file using ParmEd."""
        with open(path) as f:
            self._description = f.readline().rstrip()
            natoms = int(f.readline().strip())
            if natoms == 0:
                # ParmEd cannot parse empty GRO files; read box from last line
                box_line = f.readline().strip()
                box_vals = [float(s) for s in box_line.split()]
                self._structure = pmd.Structure()
                self.box_size = box_vals[:3]
                return
        self._structure = pmd.load_file(path)

    def get_atoms(self, resi=-1, resn="", atomtype="", atom_id=-1) -> List[GroAtom]:
        """Filter atoms by criteria. All coordinates in nm."""
        ret = []
        for atom in self._structure.atoms:
            if resi != -1 and atom.residue.number != resi:
                continue
            if resn != "" and atom.residue.name != resn:
                continue
            if atomtype != "" and atom.name != atomtype:
                continue
            if atom_id != -1 and (atom.idx + 1) != atom_id:
                continue
            ret.append(GroAtom.from_parmed(atom))
        return ret

    def add_atom(self, atom):
        """Add a GroAtom to the structure."""
        if not isinstance(atom, GroAtom):
            raise TypeError("the input is NON-GRO_ATOM")

        pmd_atom = pmd.Atom(name=atom.atomtype, mass=atom.atomic_mass)
        pmd_atom.xx = atom.point[0] * _ANGSTROM_PER_NM
        pmd_atom.xy = atom.point[1] * _ANGSTROM_PER_NM
        pmd_atom.xz = atom.point[2] * _ANGSTROM_PER_NM

        # Find existing residue or create new one
        target_res = None
        for res in self._structure.residues:
            if res.number == atom.resi and res.name == atom.resn:
                target_res = res
                break

        if target_res is not None:
            self._structure.add_atom_to_residue(pmd_atom, target_res)
        else:
            self._structure.add_atom(pmd_atom, atom.resn, atom.resi)

    def molar(self, resn):
        """Calculate molar concentration of a given residue type."""
        resis = {a.residue.number for a in self._structure.atoms if a.residue.name == resn}
        box = self.box_size
        volume = box[0] * box[1] * box[2]  # nm^3
        return (len(resis) / constants.N_A) / (volume * 1e-24)  # nm^3 -> cm^3

    def __repr__(self):
        """Serialize to GRO format string using ParmEd."""
        with tempfile.NamedTemporaryFile(suffix=".gro", delete=False) as f:
            tmppath = f.name
        try:
            self._structure.save(tmppath, overwrite=True, combine="all")
            with open(tmppath) as f:
                content = f.read()
        finally:
            os.unlink(tmppath)
        # Replace ParmEd's default title with our description
        lines = content.split("\n")
        if self._description:
            lines[0] = self._description
        return "\n".join(lines)
