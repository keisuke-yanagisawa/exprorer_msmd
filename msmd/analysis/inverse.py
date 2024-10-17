from typing import Optional

import numpy as np
import numpy.typing as npt
import pytraj as pt
from tqdm import tqdm

from msmd.pytraj import align_by_single_residue


def __enumerate_resis(traj: pt.Trajectory, resn: str) -> set[int]:
    """
    trajの中に含まれるresnのresidを列挙する
    この関数は、residue ID が1から始まることを前提としており、1-indexedで返す
    """
    tmp = traj[0:1][f":{resn}"]
    return set([res.index + 1 for res in tmp.top.residues])


def create_interaction_profile(
    trajectory: pt.Trajectory,
    probe_resn: str,
    aa_resns: list[str],
    sele_profiling_atoms: str = "@*",
    ref_probe_structure: Optional[pt.Trajectory] = None,
    verbose: bool = False,
    profile_size: int = 80,
    profile_pitch: int = 1,
) -> npt.NDArray[np.float64]:
    """create residue interaction profile from single trajectory

    Parameters
    ----------
    trajectory : pt.Trajectory
        trajectory data
    probe_resn : str
        probe residue name
    aa_resns : list[str]
        profiling residue names
    sele_profiling_atoms : str = "@*"
        selection string for profiling atoms
    ref_probe_structure : Optional[pt.Trajectory] = None
        reference structure for probe residue

        If none, the first frame of trajectory is used as reference
    verbose : bool = False
        verbose mode
    profile_size : int = 80
        profile size
    profile_pitch : int = 1
        profile pitch

    Returns
    -------
    npt.NDArray[np.float64]
        interaction profile

    Raises
    ------
    ValueError
        if probe residue is included in the profiling residues

        if profiling residues are not found in the trajectory

        if probe residue is not found in the trajectory

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.iterload(..., ...)
    >>> data = create_interaction_profile(traj, "FAA", ["HIS", "ARG", "LYS"])
    >>> print(data.shape)
    (80, 80, 80)
    """

    if probe_resn in aa_resns:
        raise ValueError("Probe residue should not be included in the profiling residues.")

    try:
        trajectory[f"(:{','.join(aa_resns)})"]
    except IndexError as e:  # 指定した残基がそのタンパク質に存在しなかった
        raise ValueError(f"Residue(s) {','.join(aa_resns)} is not existed in the trajectory.")

    probe_resis = __enumerate_resis(trajectory, probe_resn)
    if len(probe_resis) == 0:
        raise ValueError(f"Probe residue {probe_resn} is not found in the trajectory.")

    tmp = pt.strip(trajectory, ":WAT|:HOH|@/H")
    box_volume = np.prod(trajectory.unitcells[0][:3])
    bulk_ideal_cnt = 0
    map_cnt = np.zeros((profile_size, profile_size, profile_size))
    for probe_resi in tqdm(probe_resis, disable=not verbose):
        align_by_single_residue(tmp, probe_resi, ref_probe_structure)
        pt.center(
            tmp, mask=f":{probe_resi}", mass=True, center="origin"
        )  # TODO: centerを行うのではなく、 pt.gridの gridcenter x y z を指定するべき

        tmp = pt.image(tmp, "origin center byatom")
        # "byatom" によって、原子単位でimagingさせる。
        # これによってprobeからのprotein atomsの相対位置を、最も近い場所で定義する
        # （正確にはbox構造だから違うかもしれないけど、だいたいそうなっている）

        data = pt.grid(
            tmp,
            command=f"{profile_size} {profile_pitch} {profile_size} {profile_pitch} {profile_size} {profile_pitch} ({sele_profiling_atoms})&(:{','.join(aa_resns)})",
            dtype="dataset",
        ).to_ndarray()
        bulk_ideal_cnt += (
            np.sum(data) / box_volume
        )  # TODO: sum(data) ではなく、実際に何件の原子が存在したか？をカウントすべき。gridの外にあふれた原子がカウントされなくなってしまっている
        map_cnt += data
    return map_cnt / bulk_ideal_cnt
