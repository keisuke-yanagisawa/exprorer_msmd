#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
マッチングスコア計算関連の関数を提供するモジュール
"""

import numpy as np
from gridData import Grid


def preprocessing(TRAJ, ref_probe, probe_resi: int):
    """
    protein atomsがprobeの周りにくるようにして、probeの向きを整えたtrajectoryを作成する

    Parameters
    ----------
    TRAJ : pytraj.Trajectory
        トラジェクトリ
    ref_probe : pytraj.Trajectory
        参照プローブ
    probe_resi : int
        プローブの残基ID

    Returns
    -------
    pytraj.Trajectory
        前処理されたトラジェクトリ
    """
    import pytraj as pt

    # 水分子と水素原子を除去
    tmp = pt.strip(TRAJ, ":WAT|:HOH|@/H")
    
    # 周期境界条件の修正
    pt.fiximagedbonds(tmp)
    
    # プローブの向きを参照プローブに合わせる
    pt.rmsd(tmp, mask=f":{probe_resi}&!@VIS", ref=ref_probe, ref_mask=":*&!@/H")
    
    # プローブを中心に配置
    tmp = pt.center(tmp, mask=f":{probe_resi}", mass=True, center="origin")
    
    # 原子単位でイメージング
    tmp = pt.image(
        tmp, "origin center byatom"
    )  # 原子単位でimagingさせる。これによってprobeからのprotein atomsの相対位置を、最も近い場所で定義する
    
    return tmp


def calculate_matching_score(atoms_of_interest, profile_residues, probe_center, profiles, GAMMA=0.00):
    """
    マッチングスコアを計算する

    Parameters
    ----------
    atoms_of_interest : list
        対象となる原子のリスト
    profile_residues : list
        プロファイル残基のリスト
    probe_center : numpy.ndarray
        プローブの中心座標
    profiles : dict
        残基名をキー、プロファイルを値とする辞書
    GAMMA : float, optional
        重み付け係数, by default 0.00

    Returns
    -------
    float
        マッチングスコア
    """
    log_score = 0
    
    for atom in atoms_of_interest:
        resname = atom.get_parent().get_resname()
        
        if resname not in profile_residues:
            continue  # GLYなど、対象外の残基はスキップ
        
        coord = atom.get_coord()
        distance = np.linalg.norm(probe_center - coord)
        weight = np.exp(-GAMMA * (distance**2))
        
        log_score += (
            np.log(
                max(profiles[resname].interpolated([coord[0]], [coord[1]], [coord[2]])[0], profiles[resname].grid.min())
            )
            * weight
        )
    
    return log_score


def superimpose_protein_to_probe(protein_file, probe_file, atom_pairs_file, output_file):
    """
    タンパク質をプローブに重ね合わせる

    Parameters
    ----------
    protein_file : str
        タンパク質のPDBファイルパス
    probe_file : str
        プローブのPDBファイルパス
    atom_pairs_file : str
        原子対応関係のファイルパス
    output_file : str
        出力ファイルパス

    Returns
    -------
    None
    """
    from pathlib import Path

    import numpy as np

    from script.utilities.Bio import PDB
    from script.utilities.Bio.sklearn_interface import SuperImposer

    # タンパク質とプローブの読み込み
    protein = PDB.get_structure(Path(protein_file))
    probe = PDB.get_structure(Path(probe_file))
    
    # 座標の取得
    protein_coords = PDB.get_attr(protein, "coord")
    probe_coords = PDB.get_attr(probe, "coord")
    
    # 原子対応関係の読み込み
    atom_pairs = np.loadtxt(atom_pairs_file, int)
    
    # 対応する座標の取得
    probe_coords_target = probe_coords[atom_pairs[0]]
    protein_coords_target = protein_coords[atom_pairs[1]]
    
    # 重ね合わせ
    si = SuperImposer()
    si.fit(protein_coords_target, probe_coords_target)
    
    # タンパク質の座標を変換
    # モデルを取得して座標を設定
    model = protein.get_list()[0]  # 最初のモデルを取得
    PDB.set_attr(model, "coord", si.transform(protein_coords))
    
    # 結果を保存
    PDB.save(protein, output_file)


def calculate_matching_scores_for_probes(protein_file, probe_ids, profile_dir, output_file=None):
    """
    複数のプローブに対してマッチングスコアを計算する

    Parameters
    ----------
    protein_file : str
        タンパク質のPDBファイルパス
    probe_ids : list
        プローブIDのリスト
    profile_dir : str
        プロファイルディレクトリ
    output_file : str, optional
        出力ファイルパス, by default None

    Returns
    -------
    dict
        プローブIDをキー、マッチングスコアを値とする辞書
    """
    from pathlib import Path

    from gridData import Grid

    from script.utilities.Bio import PDB

    # 結果を格納する辞書
    results = {}
    
    # プロファイル残基のリスト
    profile_residues = "ALA CYS GLU ASP PHE HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split(" ")
    
    for probe_id in probe_ids:
        # タンパク質の読み込み
        protein = PDB.get_structure(Path(protein_file), f"protein_{probe_id}")
        
        # プローブの読み込み
        probe_file = f"probe_{probe_id}.pdb"
        probe = PDB.get_structure(Path(probe_file), f"probe_{probe_id}")
        
        # プローブの中心座標を計算
        probe_center = PDB.get_attr(probe, "coord").mean(axis=0)
        
        # 対象となる原子を取得（CB原子）
        atoms_of_interest = [a for a in protein.get_atoms() if a.get_name() == "CB"]
        
        # プロファイルの読み込み
        profiles = {
            res: Grid(f"{profile_dir}/{probe_id}_{res}_profile.dx") for res in profile_residues
        }
        
        # マッチングスコアの計算
        score = calculate_matching_score(atoms_of_interest, profile_residues, probe_center, profiles)
        
        # 結果を格納
        results[probe_id] = score
        
        print(f"{probe_id}: {score:.2f}")
    
    # 結果をファイルに保存
    if output_file:
        with open(output_file, "w") as f:
            for probe_id, score in results.items():
                f.write(f"{probe_id}\t{score:.2f}\n")
    
    return results