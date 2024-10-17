from typing import Optional

import pytraj as pt


def align_by_single_residue(trajectory: pt.Trajectory, resi: int, ref_struct: Optional[pt.Trajectory] = None) -> None:
    """
    Align all snapshots in the trajectory so that the structures of the specified residue ID (`resi`) overlap.
    If a `ref_struct` is provided, use it as the reference structure.
    If not, use the first snapshot of the trajectory as the reference structure.

    The `trajectory` is modified in place.
    """

    if ref_struct is None:
        ref_struct = trajectory[f":{resi}"][0]

    # check residue id is the same as the probe
    ref_name: str = next(ref_struct.top.residues).name
    focused_name: str = next(trajectory[f":{resi}"].top.residues).name
    if ref_name != focused_name:
        raise ValueError(
            f"Focused residue name {focused_name} (resi {resi}) is different from reference structure name {ref_name}"
        )

    pt.fiximagedbonds(trajectory)
    pt.rmsd(trajectory, ref=ref_struct, mask=f":{resi}&!@VIS", ref_mask=f":*&!@/H")  # @VIS: virtual site
