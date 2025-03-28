#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
分子から部分構造を抽出し、PNG画像として描画するスクリプト
"""

import os
import sys
from io import BytesIO

import matplotlib

matplotlib.use('Agg')  # グラフィカルディスプレイを必要としないバックエンドを設定
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D


def extract_substructure_around_atom(mol_file, atom_idx, radius=2, output_dir="output"):
    """
    指定した原子の周りの部分構造を抽出する
    
    Parameters
    ----------
    mol_file : str
        分子のSDFファイルパス
    atom_idx : int
        中心となる原子のインデックス
    radius : int, optional
        取り込む周辺原子の範囲（結合距離）, by default 2
    output_dir : str, optional
        出力ディレクトリ, by default "output"
    
    Returns
    -------
    str
        抽出された部分構造のSDFファイルパス
    """
    # 出力ディレクトリを作成
    os.makedirs(output_dir, exist_ok=True)
    
    # 分子の読み込み
    mol_supplier = Chem.SDMolSupplier(mol_file)
    mol = mol_supplier[0]
    
    if mol is None:
        raise ValueError(f"分子ファイル {mol_file} を読み込めませんでした")
    
    # 中心原子から指定した半径内の原子を取得
    atoms = set()
    atoms.add(atom_idx)
    
    # 指定した半径まで原子を追加
    current_atoms = {atom_idx}
    for r in range(radius):
        next_atoms = set()
        for atom_id in current_atoms:
            atom = mol.GetAtomWithIdx(atom_id)
            for neighbor in atom.GetNeighbors():
                next_atoms.add(neighbor.GetIdx())
        # 新しい原子を追加
        atoms.update(next_atoms)
        current_atoms = next_atoms
    
    # 環構造を保持するための処理
    # 中心原子が環構造に含まれる場合、その環全体を含める
    ring_info = mol.GetRingInfo()
    for ring_atoms in ring_info.AtomRings():
        if atom_idx in ring_atoms:
            atoms.update(ring_atoms)
    
    # 部分構造を抽出
    atoms_list = list(atoms)
    
    # 部分構造を抽出するためにRWMolを使用
    rwmol = Chem.RWMol()
    
    # 原子マッピングを作成
    atom_mapping = {}
    
    # 選択した原子を追加
    for old_idx in atoms_list:
        atom = mol.GetAtomWithIdx(old_idx)
        # 原子番号を取得
        atomic_num = atom.GetAtomicNum()
        # 新しい原子を追加
        new_idx = rwmol.AddAtom(Chem.Atom(atomic_num))
        atom_mapping[old_idx] = new_idx
        
        # 原子の特性をコピー
        new_atom = rwmol.GetAtomWithIdx(new_idx)
        new_atom.SetFormalCharge(atom.GetFormalCharge())
        new_atom.SetIsAromatic(atom.GetIsAromatic())
        new_atom.SetChiralTag(atom.GetChiralTag())
    
    # 結合を追加
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        
        if begin_idx in atoms_list and end_idx in atoms_list:
            rwmol.AddBond(
                atom_mapping[begin_idx],
                atom_mapping[end_idx],
                bond.GetBondType()
            )
    
    # RWMolから通常の分子に変換
    substructure = rwmol.GetMol()
    
    # 3D座標を設定
    if mol.GetNumConformers() > 0:
        conf = Chem.Conformer(substructure.GetNumAtoms())
        for old_idx, new_idx in atom_mapping.items():
            old_pos = mol.GetConformer().GetAtomPosition(old_idx)
            conf.SetAtomPosition(new_idx, old_pos)
        
        substructure.AddConformer(conf)
    
    # 分子を正規化（sanitize）
    try:
        Chem.SanitizeMol(substructure)
    except Exception as e:
        print(f"警告: 分子の正規化に失敗しました: {e}")
    
    # 部分構造を保存
    output_file = os.path.join(output_dir, os.path.basename(mol_file).replace(".sdf", "_substructure.sdf"))
    writer = Chem.SDWriter(output_file)
    writer.write(substructure)
    writer.close()
    
    print(f"部分構造を {output_file} に保存しました")
    
    # 部分構造を描画して保存（元の分子のインデックスを保持）
    draw_molecule_with_indices(substructure, output_file.replace(".sdf", ".png"), atom_mapping)
    
    # 元の分子の中で部分構造をハイライトして描画
    highlight_substructure(mol, atoms_list, os.path.join(output_dir, os.path.basename(mol_file).replace(".sdf", "_highlighted.png")))
    
    # 元の分子のインデックスと部分構造のインデックスのマッピングを保存
    mapping_file = os.path.join(output_dir, os.path.basename(mol_file).replace(".sdf", "_mapping.txt"))
    with open(mapping_file, "w") as f:
        f.write("元の分子のインデックス -> 部分構造のインデックス\n")
        for old_idx, new_idx in atom_mapping.items():
            f.write(f"{old_idx} -> {new_idx}\n")
    
    print(f"インデックスマッピングを {mapping_file} に保存しました")
    
    return output_file

def draw_molecule_with_indices(mol, output_file, atom_mapping=None):
    """
    分子構造を描画して保存する（原子インデックス付き）
    
    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        描画する分子
    output_file : str
        出力ファイルパス
    atom_mapping : dict, optional
        元の分子のインデックスと部分構造のインデックスのマッピング, by default None
    """
    # 2D座標を生成
    mol_2d = Chem.Mol(mol)
    AllChem.Compute2DCoords(mol_2d)
    
    # 原子にインデックスラベルを付ける
    atom_labels = {}
    
    if atom_mapping is None:
        # 通常のインデックスを表示（元素記号+インデックス）
        for atom in mol_2d.GetAtoms():
            idx = atom.GetIdx()
            symbol = atom.GetSymbol()
            atom_labels[idx] = f"{symbol}{idx}"
    else:
        # 元の分子のインデックスを表示（元素記号+インデックス）
        reverse_mapping = {v: k for k, v in atom_mapping.items()}
        for atom in mol_2d.GetAtoms():
            idx = atom.GetIdx()
            symbol = atom.GetSymbol()
            if idx in reverse_mapping:
                orig_idx = reverse_mapping[idx]
                atom_labels[idx] = f"{symbol}{orig_idx}"
            else:
                atom_labels[idx] = f"{symbol}{idx}"
    
    # rdMolDraw2Dを使用して分子を描画
    drawer = rdMolDraw2D.MolDraw2DCairo(800, 600)
    drawer.SetFontSize(14)  # フォントサイズを大きく
    
    # 分子を描画
    drawer.DrawMolecule(mol_2d)
    drawer.FinishDrawing()
    
    # 画像データを取得
    png_data = drawer.GetDrawingText()
    
    # PILを使用して画像を読み込む
    img = Image.open(BytesIO(png_data))
    
    # Matplotlibを使用して画像を表示し、原子ラベルを追加
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.imshow(np.array(img))
    
    # 原子の座標を取得し、ラベルを追加
    for atom_idx, label in atom_labels.items():
        atom_pos = mol_2d.GetConformer().GetAtomPosition(atom_idx)
        # RDKitの座標系からMatplotlibの座標系に変換
        # 注意: 座標変換は実際の画像サイズに合わせて調整が必要
        x = atom_pos.x * 30 + 400  # 画像の中心に合わせて調整
        y = -atom_pos.y * 30 + 300  # Y軸は反転
        ax.annotate(label, (x, y), fontsize=14, ha='center', va='center', color='black',
                   bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
    
    # 画像を保存
    plt.axis('off')  # 軸を非表示
    plt.tight_layout()
    plt.savefig(output_file, dpi=100, bbox_inches='tight')
    plt.close(fig)
    
    print(f"分子構造の画像を {output_file} に保存しました")


def highlight_substructure(mol, atom_indices, output_file):
    """
    分子の中で指定した原子をハイライトして描画する
    
    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        元の分子
    atom_indices : list
        ハイライトする原子のインデックスリスト
    output_file : str
        出力ファイルパス
    """
    # 2D座標を生成
    mol_2d = Chem.Mol(mol)
    AllChem.Compute2DCoords(mol_2d)
    
    # 原子にインデックスラベルを付ける（元素記号+インデックス）
    atom_labels = {}
    for idx in range(mol_2d.GetNumAtoms()):
        atom = mol_2d.GetAtomWithIdx(idx)
        symbol = atom.GetSymbol()
        atom_labels[idx] = f"{symbol}{idx}"
    
    # rdMolDraw2Dを使用して分子を描画（ハイライト付き）
    drawer = rdMolDraw2D.MolDraw2DCairo(800, 600)
    drawer.SetFontSize(14)  # フォントサイズを大きく
    
    # ハイライトする原子と色を設定
    highlight_atoms = atom_indices
    highlight_bonds = []
    atom_colors = {atom_idx: (1, 0, 0) for atom_idx in highlight_atoms}
    bond_colors = {}
    
    # 分子を描画
    drawer.DrawMolecule(mol_2d, highlightAtoms=highlight_atoms,
                        highlightBonds=highlight_bonds,
                        highlightAtomColors=atom_colors,
                        highlightBondColors=bond_colors)
    drawer.FinishDrawing()
    
    # 画像データを取得
    png_data = drawer.GetDrawingText()
    
    # PILを使用して画像を読み込む
    img = Image.open(BytesIO(png_data))
    
    # Matplotlibを使用して画像を表示し、原子ラベルを追加
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.imshow(np.array(img))
    
    # 原子の座標を取得し、ラベルを追加
    for atom_idx, label in atom_labels.items():
        atom_pos = mol_2d.GetConformer().GetAtomPosition(atom_idx)
        # RDKitの座標系からMatplotlibの座標系に変換
        # 注意: 座標変換は実際の画像サイズに合わせて調整が必要
        x = atom_pos.x * 30 + 400  # 画像の中心に合わせて調整
        y = -atom_pos.y * 30 + 300  # Y軸は反転
        color = 'red' if atom_idx in atom_indices else 'black'
        ax.annotate(label, (x, y), fontsize=14, ha='center', va='center', color=color,
                   bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
    
    # 画像を保存
    plt.axis('off')  # 軸を非表示
    plt.tight_layout()
    plt.savefig(output_file, dpi=100, bbox_inches='tight')
    plt.close(fig)
    
    print(f"ハイライト付き分子構造の画像を {output_file} に保存しました")


def main():
    """
    メイン関数
    """
    import argparse
    
    parser = argparse.ArgumentParser(description='分子から部分構造を抽出する')
    parser.add_argument('--mol', required=True, help='分子のSDFファイルパス')
    parser.add_argument('--atom-idx', type=int, required=True, help='中心となる原子のインデックス')
    parser.add_argument('--radius', type=int, default=2, help='取り込む周辺原子の範囲（結合距離）')
    parser.add_argument('--output-dir', default="output", help='出力ディレクトリ')
    args = parser.parse_args()
    
    # 部分構造の抽出
    output_file = extract_substructure_around_atom(
        args.mol, args.atom_idx, args.radius, args.output_dir
    )
    
    print("処理が完了しました")
    print(f"部分構造: {output_file}")


if __name__ == "__main__":
    main()