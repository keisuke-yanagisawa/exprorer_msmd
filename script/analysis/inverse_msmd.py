#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
逆MSMD解析モジュール

このモジュールは、タンパク質とマップの適合度計算、MCS計算と非共通部分抽出、
および化合物立体構造の重ね合わせに基づくdxファイル変換の機能を提供します。
"""

import os
import warnings
from pathlib import Path
from typing import Any, Iterable, Union

import gridData
import numpy as np
from Bio.PDB.Structure import Structure
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem.rdmolops import GetSymmSSSR
from scipy.ndimage import map_coordinates

from script.utilities.Bio.PDB import get_atom_attr
from script.utilities.Bio.sklearn_interface import SuperImposer
from script.utilities.molecule import classify_rings, create_substructure


# タンパク質とマップの適合度計算関連関数
def calculate_matching_score(atoms_of_interest, profiles):
    """
    マッチングスコアを計算する

    Parameters
    ----------
    atoms_of_interest : list
        対象となる原子のリスト
    profiles : dict
        残基名をキー、プロファイルを値とする辞書

    Returns
    -------
    float
        マッチングスコア

    Examples
    --------
    >>> atoms_of_interest = [a for a in protein.get_atoms() if a.get_name() == "CB"]
    >>> profile_files: dict[str, Union[str, Path]] = {
    >>>    res: Grid(f"{res}_profile.dx")
    >>>    for res in "ALA CYS GLU ASP PHE HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split(" ")
    >>> }
    >>> score = calculate_matching_score(atoms_of_interest, profile_files)
    """
    log_score = 0

    # 対象の残基だけ抽出
    atoms_of_interest = [a for a in atoms_of_interest if get_atom_attr(a, "resname") in profiles.keys()]

    for atom in atoms_of_interest:

        resname = get_atom_attr(atom, "resname")
        coord = get_atom_attr(atom, "coord")

        # TODO: interpolated はスプライン補間であり、 -1 との境界での値が大きく変化する可能性がある
        #       取り急ぎ 0 としたが、抜本的な対策が必要
        profiles[resname].grid[profiles[resname].grid == -1] = 0
        warnings.warn(
            f"グリッド {resname} には -1 の値があります。これは補間に予期しない動作を引き起こす可能性があります。",
            UserWarning,
        )
        log_score += np.log(
            max(profiles[resname].interpolated([coord[0]], [coord[1]], [coord[2]])[0], profiles[resname].grid.min())
        )

    return log_score


def calculate_protein_map_matching(
    protein: Structure,
    profiles: dict[str, gridData.Grid],
    target_atoms: list[str] = ["CB"],
) -> float:
    """
    タンパク質とマップの適合度を計算する

    Parameters
    ----------
    protein : Bio.PDB.Structure.Structure
        タンパク質構造オブジェクト
    profiles : dict
        残基名をキー、グリッドオブジェクトを値とする辞書
    target_atoms : list, optional
        対象となる原子のリスト, by default ["CB"]

    Returns
    -------
    float
        適合度スコア
    """
    # 適合度計算対象の原子を取得
    atoms_of_interest = [a for a in protein.get_atoms() if a.get_name() in target_atoms]

    # 適合度の計算
    score = calculate_matching_score(atoms_of_interest, profiles)

    return score


# MCS計算と非共通部分抽出関連関数
def get_connection_atoms(mol: Chem.Mol, atom_indices: Iterable[int], mcs_indices: Iterable[int]) -> set[int]:
    """
    非共有部分に隣接するMCSに含まれる原子のインデックスを取得する

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        分子
    atom_indices : Iterable[int]
        非共通部分の原子インデックスリスト
    mcs_indices : Iterable[int]
        MCSに対応する原子インデックスリスト

    Returns
    -------
    set
        接続点の原子インデックスリスト
    """
    connection_atoms = set()

    for atom_idx in atom_indices:
        atom = mol.GetAtomWithIdx(atom_idx)

        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()

            if neighbor_idx in mcs_indices:
                connection_atoms.add(neighbor_idx)

    return connection_atoms


def extract_ring_atoms(ring: Iterable[int], atom_indices: Iterable[int]) -> set[int]:
    """
    環の原子のうち、非共通部分に含まれる原子のみを抽出する

    Parameters
    ----------
    ring : Iterable[int]
        環の原子インデックスリスト
    atom_indices : Iterable[int]
        非共通部分の原子インデックスリスト

    Returns
    -------
    set
        環の原子のうち、非共通部分に含まれる原子のインデックスセット
    """
    ring_atom_indices = set()
    for atom_idx in ring:
        if atom_idx in atom_indices:
            ring_atom_indices.add(atom_idx)
    return ring_atom_indices


def process_separate_rings(mol, non_common_rings, atom_indices, mcs_indices):
    """
    非共通環を別々のサブ構造として処理する

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        元の分子
    non_common_rings : list
        非共通環のリスト
    atom_indices : list
        非共通部分の原子インデックスリスト
    mcs_indices : list
        MCSに対応する原子インデックスリスト

    Returns
    -------
    list
        サブ構造のリスト
    """
    substructures = []

    for ring in non_common_rings:
        # 環の原子のうち、非共通部分に含まれる原子のみを追加
        ring_atom_indices = extract_ring_atoms(ring, atom_indices)

        # 接続点を計算（環の原子に隣接する共通部分の原子のみを含める）
        ring_connection_atoms = get_connection_atoms(mol, ring_atom_indices, mcs_indices)

        # 環の原子と接続点を抽出
        if ring_atom_indices:
            ring_extraction_indices = list(ring_atom_indices) + ring_connection_atoms

            # サブ構造を作成
            substructure = create_substructure(mol, ring_extraction_indices)
            substructures.append({"mol": substructure, "atom_indices": list(ring_atom_indices)})

    return substructures


def process_non_common_rings(mol, non_common_rings, atom_indices, mcs_indices):
    """
    非共通環を単一のサブ構造として処理する

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        元の分子
    non_common_rings : list
        非共通環のリスト
    atom_indices : list
        非共通部分の原子インデックスリスト
    mcs_indices : list
        MCSに対応する原子インデックスリスト

    Returns
    -------
    list
        抽出する原子インデックスリスト
    """
    # 非共通部分の環に含まれる原子を集める
    ring_atom_indices = set()
    for ring in non_common_rings:
        # 環の原子のうち、非共通部分に含まれる原子のみを追加
        ring_atom_indices.update(extract_ring_atoms(ring, atom_indices))

    # 非共通部分の環に含まれる原子と接続点を抽出
    if ring_atom_indices:
        # 接続点を再計算（非共通環の原子に隣接する共通部分の原子のみを含める）
        new_connection_atoms = get_connection_atoms(mol, ring_atom_indices, mcs_indices)

        return list(ring_atom_indices) + new_connection_atoms

    return None


def get_substructure(
    mol: Chem.Mol,
    atom_indices: Iterable[int],
    mcs_indices: Iterable[int],
    only_non_common_rings: bool = False,
    separate_rings: bool = False,
) -> list[dict[str, Any]]:
    """
    指定された原子インデックスに基づいてサブ構造を抽出
    MCSとの接続点も含める
    TODO: MCSとの接続点を含めることで、環の中の原子だと問題になる可能性がある

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        元の分子
    atom_indices : Iterable[int]
        抽出する原子のインデックスリスト
    mcs_indices : Iterable[int]
        MCSに対応する原子のインデックスリスト
    only_non_common_rings : bool, optional
        Trueの場合、複合環のうち非共通な環のみを抽出する, by default False
    separate_rings : bool, optional
        Trueの場合、非共通な環を別々のプローブとして抽出する, by default False

    Returns
    -------
    list[dict[str, Any]]
        抽出されたサブ構造のリスト
        [{'mol': rdkit.Chem.Mol, 'atom_indices': list[int]}, ...]
    """
    # 環構造の情報を取得
    ring_atoms = GetSymmSSSR(mol)

    # 非共通部分の原子に加えて、MCSとの接続点も含める
    connection_atoms = get_connection_atoms(mol, atom_indices, mcs_indices)

    # 抽出する原子インデックスを結合
    extraction_indices = list(atom_indices) + list(connection_atoms)

    # 環の処理が必要な場合
    if only_non_common_rings:
        # 環を共通部分と非共通部分に分類
        non_common_rings, _ = classify_rings(ring_atoms, atom_indices, mcs_indices)

        if non_common_rings:
            # 環ごとに分割するかどうか
            if separate_rings:
                # 各環ごとに別々のサブ構造を作成
                substructures = process_separate_rings(mol, non_common_rings, atom_indices, mcs_indices)

                if substructures:
                    return substructures
            else:
                # 非共通環を単一のサブ構造として処理
                new_extraction_indices = process_non_common_rings(mol, non_common_rings, atom_indices, mcs_indices)

                if new_extraction_indices:
                    extraction_indices = new_extraction_indices

    # サブ構造を作成して返す
    substructure = create_substructure(mol, extraction_indices)
    return [{"mol": substructure, "atom_indices": atom_indices}]


def extract_non_common_parts(
    mol1: Chem.Mol,
    mol2: Chem.Mol,
    only_non_common_rings: bool = False,
    separate_rings: bool = False,
) -> dict:
    """
    2つの化合物間のMCSを計算し、非共通部分を抽出する

    Parameters
    ----------
    mol1 : Chem.Mol
        1つ目の化合物のRDKitの分子オブジェクト
    mol2 : Chem.Mol
        2つ目の化合物のRDKitの分子オブジェクト
    only_non_common_rings : bool, optional
        Trueの場合、複合環のうち非共通な環のみをプローブ化する, by default False
    separate_rings : bool, optional
        Trueの場合、非共通な環を別々のプローブとして抽出する, by default False

    Returns
    -------
    dict
        各化合物の非共通部分の情報
        {
            'mol1': [{'mol': rdkit.Chem.rdchem.Mol, 'atom_indices': list}, ...],
            'mol2': [{'mol': rdkit.Chem.rdchem.Mol, 'atom_indices': list}, ...]
        }
    """

    # MCSの計算
    mcs_result = rdFMCS.FindMCS([mol1, mol2], completeRingsOnly=True, ringMatchesRingOnly=True, matchValences=True)

    # MCSのSMARTSパターンを取得
    mcs_smarts = mcs_result.smartsString
    mcs_mol = Chem.MolFromSmarts(mcs_smarts)

    # 各分子でMCSに対応する原子インデックスを取得
    match1 = mol1.GetSubstructMatch(mcs_mol)
    match2 = mol2.GetSubstructMatch(mcs_mol)

    # 非共通部分の原子インデックスを取得
    all_atoms1 = set(range(mol1.GetNumAtoms()))
    all_atoms2 = set(range(mol2.GetNumAtoms()))

    non_common_atoms1 = list(all_atoms1 - set(match1))
    non_common_atoms2 = list(all_atoms2 - set(match2))

    # 非共通部分を含むサブ構造を抽出
    non_common_mol1 = get_substructure(mol1, non_common_atoms1, match1, only_non_common_rings, separate_rings)
    non_common_mol2 = get_substructure(mol2, non_common_atoms2, match2, only_non_common_rings, separate_rings)

    # separate_ringsがTrueの場合、複数のサブ構造を返す
    if separate_rings and only_non_common_rings:
        # 既に環ごとに分割されたサブ構造のリストが返されている
        return {"mol1": non_common_mol1, "mol2": non_common_mol2}  # リスト  # リスト
    else:
        # 単一のサブ構造が返されている
        return {"mol1": non_common_mol1, "mol2": non_common_mol2}


# 化合物立体構造の重ね合わせとdxファイル変換関連関数
def transform_grid(grid, sup, inverse=True):
    """
    SuperImposerオブジェクトを使用してグリッドデータを変換する

    Parameters:
    -----------
    grid : gridData.Grid
        変換するグリッドデータ
    sup : SuperImposer
        変換情報を持つSuperImposerオブジェクト
    inverse : bool, optional
        Trueの場合、逆変換を適用（新しいグリッドの各点から元のグリッドの対応する点を求める）

    Returns:
    --------
    gridData.Grid
        変換後のグリッドデータ
    """
    # グリッドの形状とサイズを取得
    grid_shape = grid.grid.shape

    # 新しいグリッドを作成
    new_grid = gridData.Grid()
    new_grid.grid = np.full(grid_shape, -1.0)  # 初期値を-1に設定
    new_grid.origin = grid.origin.copy()
    # TODO: SuperImposerは平行移動もさせるので、originは変化するはず。
    new_grid.delta = grid.delta.copy()

    # グリッド座標を生成
    x = np.arange(grid_shape[0])
    y = np.arange(grid_shape[1])
    z = np.arange(grid_shape[2])
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")

    # グリッド座標を実際の座標に変換
    real_coords = np.zeros((X.size, 3))
    for i in range(3):
        real_coords[:, i] = (np.array([X.ravel(), Y.ravel(), Z.ravel()])[i] * grid.delta[i]) + grid.origin[i]

    # 逆変換を適用して元のグリッドの対応する点を求める
    # TODO: どちらかが正しくて、どちらかは誤り。要調査。
    if inverse:
        transformed_coords = sup.inverse_transform(real_coords)
    else:
        transformed_coords = sup.transform(real_coords)

    # 元のグリッド座標系に戻す
    grid_coords = np.zeros((transformed_coords.shape[0], 3))
    for i in range(3):
        grid_coords[:, i] = (transformed_coords[:, i] - grid.origin[i]) / grid.delta[i]

    # 補間によって新しいグリッドデータを生成
    # 境界外の点は-1として処理
    values = map_coordinates(
        grid.grid, [grid_coords[:, 0], grid_coords[:, 1], grid_coords[:, 2]], order=1, mode="constant", cval=-1.0
    )

    new_grid.grid = values.reshape(grid_shape)

    return new_grid


def generate_atom_mapping(compound: Chem.Mol, probe: Chem.Mol) -> list:
    """
    化合物とプローブの原子対応関係を自動生成する

    Parameters
    ----------
    compound : Chem.Mol
        化合物のRDKitの分子オブジェクト
    probe : Chem.Mol
        プローブのRDKitの分子オブジェクト

    Returns
    -------
    list of tuple
        (化合物原子インデックス, プローブ原子インデックス)のリスト
    """
    warnings.warn(
        "generate_atom_mapping() 関数は不備があるので修正が必要です。プローブ切り出しの際にISO でタグをつけ、その情報だけを使って重ね合わせを行うように修正する予定です。"
    )

    # 最大共通部分構造（MCS）の計算
    mcs_result = rdFMCS.FindMCS([compound, probe], completeRingsOnly=True, ringMatchesRingOnly=True, matchValences=True)

    # MCSのSMARTSパターンを取得
    mcs_smarts = mcs_result.smartsString
    mcs_mol = Chem.MolFromSmarts(mcs_smarts)

    # 各分子でMCSに対応する原子インデックスを取得
    compound_match = compound.GetSubstructMatch(mcs_mol)
    probe_match = probe.GetSubstructMatch(mcs_mol)

    # 原子対応関係を生成
    atom_mapping = []
    for i in range(min(len(compound_match), len(probe_match))):
        atom_mapping.append((compound_match[i], probe_match[i]))

    # 結果を表示して確認
    print("自動生成された原子対応関係:")
    for compound_atom_idx, probe_atom_idx in atom_mapping:
        compound_atom = compound.GetAtomWithIdx(compound_atom_idx)
        probe_atom = probe.GetAtomWithIdx(probe_atom_idx)
        print(
            f"  化合物原子 {compound_atom_idx} ({compound_atom.GetSymbol()}) <-> プローブ原子 {probe_atom_idx} ({probe_atom.GetSymbol()})"
        )

    return atom_mapping


def transform_dx_based_on_molecule_alignment(
    mol1: Chem.Mol,
    mol2: Chem.Mol,
    grid: gridData.Grid,
    output_dx: Union[str, Path] = None,
) -> str:
    """
    化合物立体構造の重ね合わせに基づいてdxファイルの各値を移動させる

    Parameters
    ----------
    mol1 : Chem.Mol
        参照化合物のRDKitの分子オブジェクト
    mol2 : Chem.Mol
        ターゲット化合物のRDKitの分子オブジェクト
    grid : gridData.Grid
        変換するグリッドオブジェクト
    output_dx : Union[str, Path], optional
        出力dxファイルパス, by default None

    Returns
    -------
    str
        出力dxファイルパス
    """
    # 出力ファイルパスの設定
    if output_dx is None:
        output_dx = "transformed_grid.dx"

    print(f"出力dxファイル: {output_dx}")

    # 1. 化合物の重ね合わせ
    print("\n1. 化合物の重ね合わせを行います...")

    # 原子対応関係の自動生成
    print("  原子対応関係を自動生成しています...")
    atom_mapping = generate_atom_mapping(mol1, mol2)

    # 2. SuperImposerオブジェクトの作成
    print("\n2. 変換情報を取得しています...")

    # 対応する原子の座標を抽出
    coords1 = []
    coords2 = []

    for mol1_atom_idx, mol2_atom_idx in atom_mapping:
        coords1.append(mol1.GetConformer().GetAtomPosition(mol1_atom_idx))
        coords2.append(mol2.GetConformer().GetAtomPosition(mol2_atom_idx))

    # NumPy配列に変換
    coords1 = np.array([(p.x, p.y, p.z) for p in coords1])
    coords2 = np.array([(p.x, p.y, p.z) for p in coords2])

    # SuperImposerオブジェクトの作成
    sup = SuperImposer()
    sup.fit(coords2, coords1)  # mol2をmol1に合わせる

    print(f"  回転行列:\n{sup.rot_}")
    print(f"  平行移動ベクトル: {sup.tran_}")

    # 3. dxファイルの変換
    print("\n3. dxファイルを変換しています...")

    # グリッドデータの変換
    print("  グリッドデータを変換しています...")
    transformed_grid = transform_grid(grid, sup)

    # 変換後のグリッドデータの保存
    print(f"  変換後のグリッドデータを {output_dx} に保存しています...")
    transformed_grid.export(output_dx, type="double")

    print(f"\n変換後のdxファイルを保存しました: {output_dx}")

    return output_dx
