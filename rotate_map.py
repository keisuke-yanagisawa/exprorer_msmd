#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
DXファイル形式の3次元密度マップに対して、2つのPDB構造間のアラインメントに基づく
アフィン変換（回転と並進）を適用し、変換後のマップを新しいDXファイルとして出力するスクリプト
"""

import argparse
from pathlib import Path

import gridData
import numpy as np
from scipy.ndimage import map_coordinates

from script.utilities.Bio import PDB as uPDB
from script.utilities.Bio.sklearn_interface import SuperImposer


def select_atoms(model, resn, atom_names=None):
    """
    指定された残基名と原子名に基づいて原子を選択する
    
    Parameters:
    -----------
    model : Bio.PDB.Model
        原子を選択するモデル
    resn : str
        選択する残基名
    atom_names : list of str, optional
        選択する原子名のリスト。Noneの場合はすべての原子を選択
        
    Returns:
    --------
    numpy.ndarray
        選択された原子の座標
    """
    def selector(a):
        cond1 = uPDB.get_atom_attr(a, "resname") == resn
        cond2 = atom_names is None or uPDB.get_atom_attr(a, "fullname") in atom_names
        return cond1 and cond2
    
    return uPDB.get_attr(model, "coord", sele=selector)


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
    new_grid.delta = grid.delta.copy()
    
    # グリッド座標を生成
    x = np.arange(grid_shape[0])
    y = np.arange(grid_shape[1])
    z = np.arange(grid_shape[2])
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    # グリッド座標を実際の座標に変換
    real_coords = np.zeros((X.size, 3))
    for i in range(3):
        real_coords[:, i] = (np.array([X.ravel(), Y.ravel(), Z.ravel()])[i] * grid.delta[i]) + grid.origin[i]
    
    # 逆変換を適用して元のグリッドの対応する点を求める
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
    values = map_coordinates(grid.grid, [grid_coords[:, 0], grid_coords[:, 1], grid_coords[:, 2]], 
                            order=1, mode='constant', cval=-1.0)
    
    new_grid.grid = values.reshape(grid_shape)
    
    return new_grid


def main():
    parser = argparse.ArgumentParser(description='DXファイルに対してPDB構造間のアラインメントに基づく変換を適用する')
    parser.add_argument('input_dx', type=str, help='入力DXファイルのパス')
    parser.add_argument('ref_pdb', type=str, help='参照PDBファイルのパス')
    parser.add_argument('target_pdb', type=str, help='ターゲットPDBファイルのパス')
    parser.add_argument('--resn', type=str, default='LIG', help='アラインメントに使用する残基名')
    parser.add_argument('--atom_names', type=str, nargs='+', help='アラインメントに使用する原子名のリスト')
    parser.add_argument('--output_dx', type=str, help='出力DXファイルのパス（指定しない場合は入力ファイル名に_transformedを追加）')
    
    args = parser.parse_args()
    
    # 入力ファイルのパスを設定
    input_dx_path = Path(args.input_dx)
    ref_pdb_path = Path(args.ref_pdb)
    target_pdb_path = Path(args.target_pdb)
    
    # 出力ファイルのパスを設定
    if args.output_dx:
        output_dx_path = Path(args.output_dx)
    else:
        output_dx_path = input_dx_path.parent / f"{input_dx_path.stem}_transformed{input_dx_path.suffix}"
    
    # PDBファイルを読み込む
    ref_structure = uPDB.get_structure(ref_pdb_path)
    target_structure = uPDB.get_structure(target_pdb_path)
    
    # 最初のモデルを取得
    ref_model = ref_structure[0]
    target_model = target_structure[0]
    
    # アラインメントのための原子を選択
    ref_coords = select_atoms(ref_model, args.resn, args.atom_names)
    target_coords = select_atoms(target_model, args.resn, args.atom_names)
    
    if len(ref_coords) == 0 or len(target_coords) == 0:
        raise ValueError(f"選択された原子がありません。残基名: {args.resn}, 原子名: {args.atom_names}")
    
    if len(ref_coords) != len(target_coords):
        raise ValueError(f"参照構造とターゲット構造で選択された原子数が異なります: {len(ref_coords)} vs {len(target_coords)}")
    
    # SuperImposerを使用して変換を計算
    sup = SuperImposer()
    sup.fit(target_coords, ref_coords)
    
    # DXファイルを読み込む
    grid = gridData.Grid(str(input_dx_path))
    
    # グリッドデータを変換
    transformed_grid = transform_grid(grid, sup)
    
    # 変換後のグリッドデータを保存
    transformed_grid.export(str(output_dx_path), type="double")
    
    print(f"変換後のマップを保存しました: {output_dx_path}")


if __name__ == "__main__":
    main()