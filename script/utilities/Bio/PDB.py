# coding: utf-8

"""
BioPytnonの Bio.PDB モジュールに関係する関数群。

version: 1.3.0
last update: 1 Jul, 2022
Authors: Keisuke Yanagisawa
"""
import collections
import gzip
import os
from typing import Any, Callable, List, Literal, Optional, Union
import warnings
import numpy as np
from Bio import PDB
from Bio.PDB import PDBExceptions
from collections.abc import Iterable
import tempfile
import io
from Bio.PDB.Atom import Atom
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
import numpy.typing as npt
from ..scipy.spatial_func import estimate_volume
from ..util import expandpath


class AtomSelector(PDB.Select):
    """
    構造から、原子を抜き出すためのクラス。
    毎回classを作るのが面倒なので、別途用意した。
    """

    def __init__(self, sele):
        self.sele = sele

    def accept_atom(self, atom):
        return self.sele(atom)


class MultiModelPDBReader(object):
    """
    多数のモデルが含まれるPDBファイルを
    省メモリで読むためのヘルパークラス。
    iteratorに対応し、1つずつ読んでくれる。
    一応get_modelもあるが、これの実装は雑なので注意。
    """

    def _init_fileobj(self):
        self.model_positions = []
        self.fileend = False
        self.file.seek(0)
        while True:
            line = self.file.readline()
            if line.startswith(self.header):
                self.model_positions.append(
                    self.file.tell() - len(line.encode())
                )
                break

    def __init__(self, file: str, header: str = "MODEL"):
        if os.path.splitext(file)[1] == ".gz":
            self.file = gzip.open(file, "rt")
        else:
            self.file = open(file)
        self.model_positions = []
        self.fileend = False
        self.header = header
        self._init_fileobj()

    def __del__(self):
        if hasattr(self, "file"):
            self.file.close()

    def get_model(self, idx: int):
        """
        get a model with 0-origin
        Note that it uses idx, not MODEL ID.
        It will do sequential search, thus the computation complexity is O(N), not O(1)
        """
        if idx < 0:
            raise IndexError(f"{self.header} index out of range")

        while (len(self.model_positions) <= idx + 1):
            if self.fileend:
                raise IndexError(f"{self.header} index out of range")
            self.model_positions.append(self._next())
        self.file.seek(self.model_positions[idx])

        with tempfile.NamedTemporaryFile("w") as f:
            n_bytes_to_be_read \
                = self.model_positions[idx + 1] - self.model_positions[idx]
            self.file.seek(self.model_positions[idx])
            f.write(self.file.read(n_bytes_to_be_read))
            f.flush()
            return get_structure(f.name)

    def _next(self) -> Union[int, None]:
        """
        get next STARTING point
        """

        if (self.fileend):
            return None

        self.file.seek(self.model_positions[-1])

        while True:
            line = self.file.readline()
            if line == "":
                self.fileend = True
                break
            elif line.startswith(self.header):
                cur = self.file.tell() - len(line.encode())
                if cur == self.model_positions[-1]:
                    continue
                self.file.seek(cur)
                break
        return self.file.tell()

    def __iter__(self):
        self._init_fileobj()
        return self

    def __next__(self):
        try:
            return self.get_model(len(self.model_positions) - 1)
        except IndexError:
            raise StopIteration


class PDBIOhelper():
    """
    多数のmodelを単一のPDBに登録する時のヘルパクラスです。
    PDBIOクラスは全てのモデルを一旦メモリ上に載せる設計ですが、
    これは1つずつモデルを保存していきます。
    """

    def __init__(self, path: str):
        self.path = expandpath(path)
        self.open()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def open(self):
        self.fo = open(self.path, "w")
        self.n_models = 0

    def close(self):
        self.fo.close()
        del self.fo

    def save(self, pdb_object):
        if "fo" not in dir(self):
            raise ValueError("File is not opened. Please call open() first.")

        model_lst = []
        if pdb_object.level == "S":
            for model in pdb_object:
                model_lst.append(model.copy())
        elif pdb_object.level == "M":
            model_lst.append(pdb_object.copy())
        else:
            exit(1)
            # TODO raise error

        out_structure = Structure("")
        for model in model_lst:
            self.n_models += 1
            model.serial_num = self.n_models
            model.id = model.serial_num - 1
            out_structure.add(model)

        pdbio = PDB.PDBIO(use_model_flag=True)
        pdbio.set_structure(out_structure)
        with io.StringIO() as f:
            pdbio.save(f, write_end=False)
            f.seek(0)
            self.fo.write(f.read())


def get_structure(filepath: str, structname="") -> Structure:
    """
    Read PDB file.

    Parameters
    ----------
    filepath : str
        filepath to a PDB file

    Returns
    -------
    Bio.PDB.Structure
        Structure object of the PDB file
    """
    filepath = expandpath(filepath)
    if filepath.endswith(".gz"):
        fileobj = gzip.open(filepath, "rt")
    else:
        fileobj = open(filepath)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBExceptions.PDBConstructionWarning)
        return PDB.PDBParser(QUIET=True).get_structure(structname, fileobj)


def get_atom_attr(atom: Atom,
                  attr: Literal["resid", "resname", "coord", "element", "fullname"]
                  ) -> Union[int, str, npt.NDArray[np.float_], tuple]:
    """
    Get attribute from Bio.PDB.Atom object.
    {"resid", "resname", "coord", "element", "fullname"}
    are only acceptable as ``attr`` so far.
    Other attributes raises NotImplementedError.

    Parameters
    ----------
    atom : Atom
        An atom object.
    attr : str
        An attribute name which will be obtained.

    Returns
    -------
    int or str or np.array or tuple
        An attribute of the atom.

    Raises
    ------
    NotImplementedError
        If the ``attr`` is not "resid", "resname", coord", "element", nor "fullname".
    """

    if attr == "resid":
        return get_resi(atom)
    elif attr == "resname":
        return get_resname(atom)
    elif attr == "coord":
        return atom.get_coord()
    elif attr == "element":
        return atom.element
    elif attr == "fullname":
        return atom.fullname
    else:
        raise NotImplementedError(f"Attribute {attr} is not supported yet.")


def get_attr(model: Union[Structure, Model],
             attr: Literal["resid", "resname", "coord", "element", "fullname"],
             sele: Optional[Callable[[Atom], bool]] = None
             ) -> npt.NDArray[Any]:
    """
    Get attribute from Bio.PDB.Model object.
    {"resid", "resname", "coord", "element", "fullname"} 
    are only acceptable as ``attr`` so far.
    Other attributes raises NotImplementedError.

    Parameters
    ----------
    model : Model
        A model object.
    attr : str
        An attribute name which will be obtained.
    sele : function, optional
        Atom selector function. all atoms will be selected if ``sele`` is not
        provided.

    Returns
    -------
    np.array
        An array of attributes of all atoms selected by ``sele`` function.

    Raises
    ------
    NotImplementedError
        If the ``attr`` is not "resid", "resname", coord", "element", nor "fullname".
    """

    data = []
    for atom in model.get_atoms():
        if sele is None or sele(atom):
            data.append(get_atom_attr(atom, attr))
    return np.array(data)


def get_resname(atom: Atom) -> str:
    """
    Get residue name from Bio.PDB.Atom.

    Parameters
    ----------
    atom : Bio.PDB.Atom
        An atom object.

    Returns
    -------
    str
        A residue name of the atom of interest.
    """
    return atom.get_parent().get_resname()  # type: ignore


def is_water(atom: Atom) -> bool:
    """
    Check whether the provided Bio.PDB.Atom instance is the one of water atoms.

    Parameters
    ----------
    atom : Bio.PDB.Atom
        An atom object.

    Returns
    -------
    bool
        Whether the atom of interest is the one of water atoms.
    """
    return get_resname(atom) == "WAT" or get_resname(atom) == "HOH"


def get_resi(atom: Atom) -> int:
    """
    Get residue sequence number (residue ID) from Bio.PDB.Atom.

    Parameters
    ----------
    atom : Bio.PDB.Atom
        An atom object.

    Returns
    -------
    int
        A residue ID of the atom of interest.
    """
    return atom.get_full_id()[3][1]  # type: ignore


def is_hetero(atom: Atom) -> bool:
    """
    Check whether the provided Bio.PDB.Atom instance is a hetero atom.

    Parameters
    ----------
    atom : Bio.PDB.Atom
        An atom object.

    Returns
    -------
    bool
        Whether the atom of interest is a hetero atom.
    """
    return atom.get_full_id()[3][0] != " "  # type: ignore


def is_hydrogen(atom: Atom) -> bool:
    """
    Check whether the provided Bio.PDB.Atom instance is hydrogen.

    Parameters
    ----------
    atom : Bio.PDB.Atom
        An atom object.

    Returns
    -------
    bool
        Whether the atom of interest is hydrogen.
    """
    return atom.get_fullname()[1] == "H"


def set_attr(model: Model, attr: str, lst: npt.NDArray, sele=None) -> None:
    """
    Set attribute to Bio.PDB.Model object.
    attr == "coord" is only acceptable so far.
    Other attributes raises NotImplementedError.

    Parameters
    ----------
    model : Bio.PDB.Model
    attr : str
    lst : array_like
    sele : function, optional

    Raises
    ------
    NotImplementedError
        If the ``attr`` is not "coord".
    """

    # TODO check the length of lst and the number of atoms.
    # if they are different, set_attr() must not assign new values.

    lst_idx = 0
    for atom in model.get_atoms():
        if sele is None or sele(atom):
            if attr == "coord":
                atom.set_coord(lst[lst_idx])
            else:
                raise NotImplementedError(f"set_attr(attr={attr}) is not implemented")
            lst_idx += 1


def save(structs, path) -> None:
    """
    Set attribute to Bio.PDB.Model object.
    attr == "coord" is only acceptable so far.
    Other attributes raises NotImplementedError.

    Parameters
    ----------
    structs : Bio.PDB.Struct or list of Bio.PDB.Struct
    path : str
    """
    path = expandpath(path)

    if not isinstance(structs, Iterable):
        structs = [structs]

    mod_structs = []
    for struct in structs:
        io = PDB.PDBIO()
        io.set_structure(struct)
        with tempfile.NamedTemporaryFile(suffix=".pdb") as fp:
            io.save(fp.name)
            mod_structs.append(get_structure(fp.name)[0])

    out_structure = Structure("")
    for struct in mod_structs:
        struct.id = len(out_structure)
        struct.serial_num = struct.id + 1
        out_structure.add(struct)

    io = PDB.PDBIO()
    io.set_structure(out_structure)
    io.save(path)


def concatenate_structures(structs: List[Structure]) -> Structure:
    """
    Concatenate structures.
    All structures are saved to a single structure with multiple models.
    """

    with tempfile.NamedTemporaryFile("w") as f:

        out_helper = PDBIOhelper(f.name)
        for struct in structs:
            out_helper.save(struct)
        out_helper.close()
        ret_structure = get_structure(f.name)
    return ret_structure


_ATOMIC_RADII = collections.defaultdict(lambda: 2.0)
_ATOMIC_RADII.update(
    {
        "H": 1.200,
        "HE": 1.400,
        "C": 1.700,
        "N": 1.550,
        "O": 1.520,
        "F": 1.470,
        "NA": 2.270,
        "MG": 1.730,
        "P": 1.800,
        "S": 1.800,
        "CL": 1.750,
        "K": 2.750,
        "CA": 2.310,
        "NI": 1.630,
        "CU": 1.400,
        "ZN": 1.390,
        "SE": 1.900,
        "BR": 1.850,
        "CD": 1.580,
        "I": 1.980,
        "HG": 1.550,
    }
)


def estimate_exclute_volume(prot: Union[Structure, Model]) -> float:
    """
    VdW半径に基づいてタンパク質の排除体積を計算
    タンパク質は原子種類に応じて処理し、
    溶媒は炭素原子1つ分の大きさであるとする
    """
    coords = []
    radii = []
    for atom in prot.get_atoms():
        if not atom.get_parent().resname in ["HOH", "WAT"]:
            coords.append(atom.get_coord())
            radii.append(_ATOMIC_RADII[atom.element])
    coords = np.array(coords)
    radii = np.array(radii)
    radii += _ATOMIC_RADII["C"]  # solvents' VdW radius: estimated by carbon radius.
    return estimate_volume(coords, radii)


class Selector(PDB.Select):
    def __init__(self, sele):
        self.sele = sele

    def accept_atom(self, atom):
        return self.sele(atom)


def extract_substructure(struct: Union[Structure, Model], sele: PDB.Select) -> Structure:
    pdbio = PDB.PDBIO()
    pdbio.set_structure(struct)
    with tempfile.NamedTemporaryFile(suffix=".pdb") as fp:
        pdbio.save(fp.name, select=sele)
        substruct = get_structure(fp.name)
    return substruct
