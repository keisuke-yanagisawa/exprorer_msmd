import tempfile
from typing import List
from tqdm import tqdm
from Bio.PDB.Atom import Atom
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from script.utilities.Bio import PDB as uPDB
from script.utilities.Bio.sklearn_interface import SuperImposer

DESCRIPTION = """
superimpose structures in accordance with specific atoms
"""


def align_res_env(structs: List[Structure],
                  reference: Model,
                  resn: str,
                  focused: List[str] = [],
                  verbose: bool = False) -> Structure:
    """
    Align the residues of a list of structures to a reference model

    Parameters
    ----------
    structs : List[Structure]
        the list of structures to be aligned
        It is acceptable that the structures contain different number of models
    reference : Model
        the reference model to be aligned to
    resn : str
        the residue name of the residues to be aligned
    focused : List[str], optional
        the list of atom names to be aligned, by default []
        example: [" CB ", " CA ", " N  ", " C  "]
    verbose : bool, optional
        If True, display progress information, by default False

    Returns
    -------
    Structure
        the aligned structure contains the same number of models as the input
    """

    structs = [struct for struct in structs if len([a for a in struct.get_atoms()]) > 0]  # remove empty structures
    if len(structs) == 0:
        raise ValueError("No structure to align")

    def selector(a: Atom):
        cond1 = uPDB.get_atom_attr(a, "resname") == resn
        cond2 = len(focused) == 0 or uPDB.get_atom_attr(a, "fullname") in focused
        return cond1 and cond2

    ref_probe_c_coords = uPDB.get_attr(reference, "coord", sele=selector)
    if len(ref_probe_c_coords) == 0:
        raise ValueError("No reference atom to align")
    # print(ref_probe_c_coords)

    sup = SuperImposer()
    _, tmppdb = tempfile.mkstemp(suffix=".pdb")
    with uPDB.PDBIOhelper(tmppdb) as pdbio:
        for struct in structs:
            for model in tqdm(struct, desc="[align res. env.]", disable=not verbose):

                # print(struct, i)
                probe_coords = uPDB.get_attr(model, "coord", sele=selector)
                sup.fit(probe_coords, ref_probe_c_coords)
                all_coords = uPDB.get_attr(model, "coord")
                uPDB.set_attr(model, "coord", sup.transform(all_coords))

                pdbio.save(model)
                # print(len(pdbio))

    return uPDB.get_structure(tmppdb)
