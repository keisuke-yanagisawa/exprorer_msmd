# coding: utf-8

"""
BioPytnonの Bio.PDB モジュールに関係する関数群。

version: 1.2.0
last update: 14 Sep, 2021
Authors: Keisuke Yanagisawa
"""
import gzip
import numpy as np
from Bio import PDB
from collections.abc import Iterable
import tempfile
import io

class PDBIOhelper():
    """
    多数のmodelを単一のPDBに登録する時のヘルパクラスです。
    PDBIOクラスは全てのモデルを一旦メモリ上に載せる設計ですが、
    これは1つずつモデルを保存していきます。
    """
    def __init__(self, path):
        self.path = path
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
        if not "fo" in dir(self):
            exit(1)
            # TODO raise error
            
        model_lst = []
        if pdb_object.level == "S":
            for model in pdb_object:
                model_lst.append(model.copy())
        elif pdb_object.level == "M":
            model_lst.append(pdb_object.copy())
        else:
            exit(1)
            # TODO raise error

        out_structure = PDB.Structure.Structure("")
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
        

def get_structure(filepath: str, structname="") -> PDB.Structure:
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
    if filepath.endswith(".gz"):
        fileobj = gzip.open(filepath, "rt")
    else:
        fileobj = open(filepath)

    return PDB.PDBParser(QUIET=True).get_structure(structname, fileobj)


def get_attr(model, attr, sele=None):
    """
    Get attribute from Bio.PDB.Model object.
    {"resid", "resname", "coord", "element", "fullname"} 
    are only acceptable as ``attr`` so far.
    Other attributes raises NotImplementedError.

    Parameters
    ----------
    model : Bio.PDB.Model
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
    method = None
    if attr == "resid":
        method = lambda a: get_resi(a)
    elif attr == "resname":
        method = lambda a: get_resname(a)
    elif attr == "coord":
        method = lambda a: a.get_coord()
    elif attr == "element":
        method = lambda a: a.element
    elif attr == "fullname":
        method = lambda a: a.fullname
    else:
        raise NotImplementedError(f"attr {attr} is not implemented")

    if isinstance(model, PDB.Atom.Atom):
        return method(model)
    else:
        data = []
        for atom in model.get_atoms():
            if sele is None or sele(atom):
                data.append(method(atom))
        return np.array(data)


def get_resname(atom: PDB.Atom) -> str:
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
    return atom.get_parent().get_resname()


def is_water(atom: PDB.Atom) -> bool:
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


def get_resi(atom: PDB.Atom) -> int:
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
    return atom.get_full_id()[3][1]


def is_hetero(atom: PDB.Atom) -> bool:
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
    return atom.get_full_id()[3][0] != " "


def is_hydrogen(atom: PDB.Atom) -> bool:
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


def set_attr(model: PDB.Model, attr: str, lst, sele=None):
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

    if not isinstance(structs, Iterable):
        structs = [structs]

    mod_structs = []
    for struct in structs:
        io = PDB.PDBIO()
        io.set_structure(struct)
        with tempfile.NamedTemporaryFile(suffix=".pdb") as fp:
            io.save(fp.name)
            mod_structs.append(get_structure(fp.name)[0])
        
    out_structure = PDB.Structure.Structure("")
    for struct in mod_structs:
        struct.id = len(out_structure)
        struct.serial_num = struct.id + 1
        out_structure.add(struct)
    
    io = PDB.PDBIO()
    io.set_structure(out_structure)
    io.save(path)
