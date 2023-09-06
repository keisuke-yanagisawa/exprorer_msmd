# coding: utf-8

"""
BioPytnonの Bio.PDB モジュールに関係する関数群。

version: 2.0.0
last update: 6 Sep, 2023
Authors: Keisuke Yanagisawa
"""
import gzip
import warnings
from Bio import PDB
from Bio.PDB import PDBExceptions
from Bio.PDB.Structure import Structure


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
    if filepath.endswith(".gz"):
        fileobj = gzip.open(filepath, "rt")
    else:
        fileobj = open(filepath)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PDBExceptions.PDBConstructionWarning)
        return PDB.PDBParser(QUIET=True).get_structure(structname, fileobj)
