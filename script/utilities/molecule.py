#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
分子構造操作ユーティリティモジュール
"""

from typing import Iterable

from rdkit import Chem


def create_substructure(mol, extraction_indices):
    """
    指定された原子インデックスに基づいてサブ構造を作成

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        元の分子
    extraction_indices : list
        抽出する原子のインデックスリスト

    Returns
    -------
    rdkit.Chem.rdchem.Mol
        抽出されたサブ構造
    """
    # 新しい分子を作成
    edit_mol = Chem.EditableMol(Chem.Mol())

    # 抽出する原子を追加
    new_idx_map = {}
    for old_idx in extraction_indices:
        atom = mol.GetAtomWithIdx(old_idx)
        new_idx = edit_mol.AddAtom(atom)
        new_idx_map[old_idx] = new_idx

    # 結合を追加
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()

        if begin_idx in extraction_indices and end_idx in extraction_indices:
            edit_mol.AddBond(new_idx_map[begin_idx], new_idx_map[end_idx], bond.GetBondType())

    # 編集可能な分子から通常の分子に変換
    substructure = edit_mol.GetMol()

    # 3D座標を設定
    conf = Chem.Conformer(len(extraction_indices))
    for old_idx, new_idx in new_idx_map.items():
        old_pos = mol.GetConformer().GetAtomPosition(old_idx)
        conf.SetAtomPosition(new_idx, old_pos)

    substructure.AddConformer(conf)

    return substructure


def classify_rings(ring_atoms, atom_indices, mcs_indices):
    """
    環を共通部分と非共通部分に分類する

    Parameters
    ----------
    ring_atoms : list
        環構造の原子インデックスのリスト
    atom_indices : list
        非共通部分の原子インデックスリスト
    mcs_indices : list
        MCSに対応する原子インデックスリスト

    Returns
    -------
    tuple
        (非共通環のリスト, 共通環のリスト)
    """
    non_common_rings = []
    common_rings = []

    for ring in ring_atoms:
        # 環の原子が非共通部分に含まれるかチェック
        if any(atom_idx in atom_indices for atom_idx in ring):
            # 環全体が共通部分に含まれる場合は共通環とする
            if all(atom_idx in mcs_indices for atom_idx in ring):
                common_rings.append(ring)
            # それ以外は非共通環とする
            else:
                non_common_rings.append(ring)

    return non_common_rings, common_rings
