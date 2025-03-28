#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
分子の重ね合わせ自動化機能を提供するモジュール
"""

import os
import select
import sys
import time

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdDepictor, rdFMCS


def generate_atom_mapping(compound_file, probe_file):
    """
    化合物とプローブの原子対応関係を自動生成する

    Parameters
    ----------
    compound_file : str
        化合物のSDFファイルパス
    probe_file : str
        プローブのSDFファイルパス

    Returns
    -------
    list of tuple
        (化合物原子インデックス, プローブ原子インデックス)のリスト
    """
    # 分子の読み込み
    compound = Chem.SDMolSupplier(compound_file)[0]
    probe = Chem.SDMolSupplier(probe_file)[0]
    
    if compound is None:
        raise ValueError(f"化合物ファイル {compound_file} を読み込めませんでした")
    if probe is None:
        raise ValueError(f"プローブファイル {probe_file} を読み込めませんでした")
    
    # 最大共通部分構造（MCS）の計算
    mcs_result = rdFMCS.FindMCS(
        [compound, probe],
        completeRingsOnly=True,
        ringMatchesRingOnly=True,
        matchValences=True
    )
    
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
        print(f"  化合物原子 {compound_atom_idx} ({compound_atom.GetSymbol()}) <-> プローブ原子 {probe_atom_idx} ({probe_atom.GetSymbol()})")
    
    # 分子の可視化
    visualize_atom_mapping(compound, probe, atom_mapping)
    
    return atom_mapping


def visualize_atom_mapping(compound, probe, atom_mapping):
    """
    原子対応関係を可視化する

    Parameters
    ----------
    compound : rdkit.Chem.rdchem.Mol
        化合物分子
    probe : rdkit.Chem.rdchem.Mol
        プローブ分子
    atom_mapping : list of tuple
        (化合物原子インデックス, プローブ原子インデックス)のリスト
    """
    try:
        # 2D座標を生成
        compound_2d = Chem.Mol(compound)
        probe_2d = Chem.Mol(probe)
        
        # 2D座標を生成
        rdDepictor.Compute2DCoords(compound_2d)
        rdDepictor.Compute2DCoords(probe_2d)
        
        # 対応する原子をハイライト
        compound_atoms = [idx for idx, _ in atom_mapping]
        probe_atoms = [idx for _, idx in atom_mapping]
        
        # アスキーアートで表示
        print("\n化合物の構造（対応する原子に*付き）:")
        compound_ascii = create_ascii_depiction(compound_2d, compound_atoms)
        print(compound_ascii)
        
        print("\nプローブの構造（対応する原子に*付き）:")
        probe_ascii = create_ascii_depiction(probe_2d, probe_atoms)
        print(probe_ascii)
        
        # 対応関係を表示
        print("\n原子対応関係:")
        for i, (compound_idx, probe_idx) in enumerate(atom_mapping):
            compound_atom = compound.GetAtomWithIdx(compound_idx)
            probe_atom = probe.GetAtomWithIdx(probe_idx)
            print(f"  対応 {i+1}: 化合物原子 {compound_idx} ({compound_atom.GetSymbol()}) <-> プローブ原子 {probe_idx} ({probe_atom.GetSymbol()})")
        
        # 画像ファイルに保存
        try:
            # 対応する原子に番号を付ける
            compound_labels = {idx: str(i) for i, (idx, _) in enumerate(atom_mapping)}
            probe_labels = {idx: str(i) for i, (_, idx) in enumerate(atom_mapping)}
            
            # 分子を描画
            compound_img = Draw.MolToImage(compound_2d, highlightAtoms=compound_atoms, highlightColor=(1, 0, 0), 
                                          highlightBonds=[], atomLabels=compound_labels)
            probe_img = Draw.MolToImage(probe_2d, highlightAtoms=probe_atoms, highlightColor=(1, 0, 0), 
                                       highlightBonds=[], atomLabels=probe_labels)
            
            # 画像を保存
            os.makedirs("./output", exist_ok=True)
            compound_img.save("./output/compound_mapping.png")
            probe_img.save("./output/probe_mapping.png")
            
            print("\n分子の可視化画像を保存しました:")
            print("  化合物: ./output/compound_mapping.png")
            print("  プローブ: ./output/probe_mapping.png")
            
            # 両方の分子を並べて表示
            img = Draw.MolsToGridImage([compound_2d, probe_2d], 
                                      molsPerRow=2, 
                                      subImgSize=(300, 300), 
                                      legends=["化合物", "プローブ"],
                                      highlightAtomLists=[compound_atoms, probe_atoms],
                                      highlightBondLists=[[], []],
                                      highlightColor=(1, 0, 0))
            img.save("./output/molecules_comparison.png")
            print("  比較画像: ./output/molecules_comparison.png")
        except Exception as e:
            print(f"画像の保存中にエラーが発生しました: {e}")
    except Exception as e:
        print(f"可視化中にエラーが発生しました: {e}")


def create_ascii_depiction(mol, highlight_atoms=None):
    """
    分子のアスキーアート表現を生成する

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        分子
    highlight_atoms : list, optional
        ハイライトする原子のリスト, by default None

    Returns
    -------
    str
        アスキーアート表現
    """
    # 原子のシンボルを取得
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
    
    # 結合情報を取得
    bonds = []
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        bond_type = bond.GetBondType()
        bonds.append((begin_idx, end_idx, bond_type))
    
    # アスキーアート表現を生成
    ascii_art = []
    ascii_art.append("  Atoms:")
    for i, symbol in enumerate(symbols):
        highlight = "*" if highlight_atoms and i in highlight_atoms else " "
        ascii_art.append(f"    {i}: {symbol}{highlight}")
    
    ascii_art.append("  Bonds:")
    for begin_idx, end_idx, bond_type in bonds:
        ascii_art.append(f"    {begin_idx} -- {end_idx} ({bond_type})")
    
    return "\n".join(ascii_art)


def validate_and_refine_atom_mapping(compound_file, probe_file, atom_mapping):
    """
    自動生成された原子対応関係を検証し、必要に応じて修正する

    Parameters
    ----------
    compound_file : str
        化合物のSDFファイルパス
    probe_file : str
        プローブのSDFファイルパス
    atom_mapping : list of tuple
        (化合物原子インデックス, プローブ原子インデックス)のリスト

    Returns
    -------
    list of tuple
        検証・修正後の(化合物原子インデックス, プローブ原子インデックス)のリスト
    """
    # 分子の読み込み
    compound = Chem.SDMolSupplier(compound_file)[0]
    probe = Chem.SDMolSupplier(probe_file)[0]
    
    # 自動生成された原子対応関係を検証
    valid_mapping = []
    
    for compound_atom_idx, probe_atom_idx in atom_mapping:
        compound_atom = compound.GetAtomWithIdx(compound_atom_idx)
        probe_atom = probe.GetAtomWithIdx(probe_atom_idx)
        
        # 原子タイプが一致するか確認
        if compound_atom.GetSymbol() == probe_atom.GetSymbol():
            # 結合パターンが類似しているか確認
            compound_bonds = len(compound_atom.GetBonds())
            probe_bonds = len(probe_atom.GetBonds())
            
            if abs(compound_bonds - probe_bonds) <= 1:
                valid_mapping.append((compound_atom_idx, probe_atom_idx))
    
    # 検証結果を表示
    print(f"検証前の対応関係数: {len(atom_mapping)}")
    print(f"検証後の有効な対応関係数: {len(valid_mapping)}")
    
    # 分子の可視化
    if valid_mapping:
        visualize_atom_mapping(compound, probe, valid_mapping)
    
    # ユーザーに確認（10秒のタイムアウトを設定）
    print("自動生成された原子対応関係を修正しますか？ (y/n)")
    print("（10秒以内に入力がない場合は、自動的にnが選択されます）")
    
    # 標準入力からの入力を10秒間待機
    response = "n"  # デフォルト値
    
    # select.selectを使用して、タイムアウト付きの入力待機を実装
    rlist, _, _ = select.select([sys.stdin], [], [], 10)
    if rlist:
        response = sys.stdin.readline().strip().lower()
    else:
        print("タイムアウトしました。自動的にnが選択されました。")
    
    if response == "y":
        # ユーザーによる修正
        print("修正する対応関係を入力してください（空行で終了）:")
        print("形式: 化合物原子インデックス プローブ原子インデックス")
        
        while True:
            # 標準入力からの入力を10秒間待機
            rlist, _, _ = select.select([sys.stdin], [], [], 10)
            if not rlist:
                print("タイムアウトしました。修正を終了します。")
                break
            
            user_input = sys.stdin.readline().strip()
            if user_input == "":
                break
            
            try:
                compound_atom_idx, probe_atom_idx = user_input.split()
                valid_mapping.append((int(compound_atom_idx), int(probe_atom_idx)))
            except:
                print("無効な入力です。形式: 化合物原子インデックス プローブ原子インデックス")
        
        # 修正後の対応関係を可視化
        if valid_mapping:
            visualize_atom_mapping(compound, probe, valid_mapping)
    
    return valid_mapping


def calculate_optimal_transform(coords1, coords2):
    """
    Kabsch algorithmを使用して最適な回転行列と平行移動ベクトルを計算

    Parameters
    ----------
    coords1 : numpy.ndarray
        変換元の座標セット
    coords2 : numpy.ndarray
        変換先の座標セット

    Returns
    -------
    tuple
        (回転行列, 平行移動ベクトル)
    """
    # 重心を計算
    centroid1 = np.mean(coords1, axis=0)
    centroid2 = np.mean(coords2, axis=0)
    
    # 重心を原点に移動
    coords1_centered = coords1 - centroid1
    coords2_centered = coords2 - centroid2
    
    # 共分散行列を計算
    covariance_matrix = np.dot(coords1_centered.T, coords2_centered)
    
    # 特異値分解
    U, S, Vt = np.linalg.svd(covariance_matrix)
    
    # 回転行列を計算
    rotation_matrix = np.dot(Vt.T, U.T)
    
    # 行列式が1になるように調整（反転を防ぐ）
    if np.linalg.det(rotation_matrix) < 0:
        Vt[-1, :] *= -1
        rotation_matrix = np.dot(Vt.T, U.T)
    
    # 平行移動ベクトルを計算
    translation_vector = centroid2 - np.dot(rotation_matrix, centroid1)
    
    return rotation_matrix, translation_vector


def align_2d_coords(template_mol, target_mol, atom_mapping):
    """
    テンプレート分子に合わせて対象分子の2D座標を整列させる

    Parameters
    ----------
    template_mol : rdkit.Chem.rdchem.Mol
        テンプレート分子
    target_mol : rdkit.Chem.rdchem.Mol
        対象分子
    atom_mapping : list of tuple
        (テンプレート原子インデックス, 対象原子インデックス)のリスト

    Returns
    -------
    rdkit.Chem.rdchem.Mol
        整列後の対象分子
    """
    # テンプレート分子の2D座標を生成
    rdDepictor.Compute2DCoords(template_mol)
    
    # 対象分子の2D座標を生成（テンプレートに合わせる）
    aligned_mol = Chem.Mol(target_mol)
    
    # 対応する原子のマッピングを作成
    coordMap = {}
    for template_idx, target_idx in atom_mapping:
        # テンプレート分子の2D座標を取得
        template_pos = template_mol.GetConformer().GetAtomPosition(template_idx)
        # 対象分子の対応する原子に座標を設定
        coordMap[target_idx] = Chem.rdGeometry.Point2D(template_pos.x, template_pos.y)
    
    # 対象分子の2D座標を生成（テンプレートに合わせる）
    rdDepictor.Compute2DCoords(aligned_mol, coordMap=coordMap)
    
    return aligned_mol


def create_vertical_comparison_image(compound, probe, atom_mapping, output_path):
    """
    化合物とプローブの分子を縦方向に並べて比較画像を作成する

    Parameters
    ----------
    compound : rdkit.Chem.rdchem.Mol
        化合物分子
    probe : rdkit.Chem.rdchem.Mol
        プローブ分子
    atom_mapping : list of tuple
        (化合物原子インデックス, プローブ原子インデックス)のリスト
    output_path : str
        出力画像ファイルパス
    """
    from PIL import Image, ImageDraw

    # 2D座標を生成
    compound_2d = Chem.Mol(compound)
    rdDepictor.Compute2DCoords(compound_2d)
    
    # プローブの2D座標を生成（化合物に合わせる）
    probe_2d = align_2d_coords(
        compound_2d,
        probe,
        [(compound_idx, probe_idx) for compound_idx, probe_idx in atom_mapping]
    )
    
    # 対応する原子をハイライト
    compound_atoms = [idx for idx, _ in atom_mapping]
    probe_atoms = [idx for _, idx in atom_mapping]
    
    # 分子を描画
    compound_img = Draw.MolToImage(compound_2d, size=(300, 300), highlightAtoms=compound_atoms, highlightColor=(1, 0, 0))
    probe_img = Draw.MolToImage(probe_2d, size=(300, 300), highlightAtoms=probe_atoms, highlightColor=(0, 0, 1))
    
    # 縦方向に並べた画像を作成
    width = max(compound_img.width, probe_img.width)
    height = compound_img.height + probe_img.height + 40  # 40ピクセルの余白
    
    comparison_img = Image.new('RGB', (width, height), (255, 255, 255))
    
    # 化合物の画像を上部に配置
    comparison_img.paste(compound_img, ((width - compound_img.width) // 2, 0))
    
    # プローブの画像を下部に配置
    comparison_img.paste(probe_img, ((width - probe_img.width) // 2, compound_img.height + 40))
    
    # ラベルを追加
    draw = ImageDraw.Draw(comparison_img)
    draw.text((10, 10), "Compound", fill=(0, 0, 0))
    draw.text((10, compound_img.height + 20), "Probe", fill=(0, 0, 0))
    
    # 画像を保存
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    comparison_img.save(output_path)
    
    return comparison_img


def superimpose_molecules(compound_file, probe_file, atom_mapping):
    """
    化合物とプローブの分子を重ね合わせる

    Parameters
    ----------
    compound_file : str
        化合物のSDFファイルパス
    probe_file : str
        プローブのSDFファイルパス
    atom_mapping : list of tuple
        (化合物原子インデックス, プローブ原子インデックス)のリスト

    Returns
    -------
    rdkit.Chem.rdchem.Mol
        重ね合わせ後のプローブ分子
    """
    # 分子の読み込み
    compound = Chem.SDMolSupplier(compound_file)[0]
    probe = Chem.SDMolSupplier(probe_file)[0]
    
    # 3D座標の取得
    compound_conf = compound.GetConformer()
    probe_conf = probe.GetConformer()
    
    # 対応する原子の座標を抽出
    compound_coords = []
    probe_coords = []
    
    for compound_atom_idx, probe_atom_idx in atom_mapping:
        compound_coords.append(compound_conf.GetAtomPosition(compound_atom_idx))
        probe_coords.append(probe_conf.GetAtomPosition(probe_atom_idx))
    
    # NumPy配列に変換
    compound_coords = np.array([(p.x, p.y, p.z) for p in compound_coords])
    probe_coords = np.array([(p.x, p.y, p.z) for p in probe_coords])
    
    # Kabsch algorithmを使用して回転行列と平行移動ベクトルを計算
    rotation_matrix, translation_vector = calculate_optimal_transform(
        probe_coords, compound_coords
    )
    
    # プローブ分子の全原子に変換を適用
    transformed_probe = Chem.Mol(probe)
    transformed_conf = transformed_probe.GetConformer()
    
    for atom_idx in range(transformed_probe.GetNumAtoms()):
        pos = transformed_conf.GetAtomPosition(atom_idx)
        pos_array = np.array([pos.x, pos.y, pos.z])
        
        # 回転と平行移動を適用
        new_pos = np.dot(rotation_matrix, pos_array) + translation_vector
        
        # 新しい座標を設定
        transformed_conf.SetAtomPosition(atom_idx, 
                                        Chem.rdGeometry.Point3D(new_pos[0], new_pos[1], new_pos[2]))
    
    # 重ね合わせ結果の情報を表示
    print("\n重ね合わせ結果:")
    print(f"  化合物の原子数: {compound.GetNumAtoms()}")
    print(f"  プローブの原子数: {probe.GetNumAtoms()}")
    print(f"  対応関係の数: {len(atom_mapping)}")
    print(f"  回転行列:\n{rotation_matrix}")
    print(f"  平行移動ベクトル: {translation_vector}")
    
    # 重ね合わせ結果を可視化
    try:
        # 縦方向に並べた比較画像を作成
        output_path = "./output/superimposed.png"
        create_vertical_comparison_image(compound, transformed_probe, atom_mapping, output_path)
        
        print("\n重ね合わせ結果の可視化画像を保存しました:")
        print(f"  {output_path}")
    except Exception as e:
        print(f"可視化中にエラーが発生しました: {e}")
    
    return transformed_probe


def main():
    """
    メイン関数
    """
    import argparse
    
    parser = argparse.ArgumentParser(description='分子の重ね合わせ自動化ツール')
    parser.add_argument('--compound', required=True, help='化合物のSDFファイルパス')
    parser.add_argument('--probe', required=True, help='プローブのSDFファイルパス')
    parser.add_argument('--output', required=True, help='出力SDFファイルパス')
    args = parser.parse_args()
    
    # 原子対応関係の自動生成
    print(f"化合物 {args.compound} とプローブ {args.probe} の原子対応関係を自動生成中...")
    atom_mapping = generate_atom_mapping(args.compound, args.probe)
    
    # 原子対応関係の検証と修正
    atom_mapping = validate_and_refine_atom_mapping(args.compound, args.probe, atom_mapping)
    
    # 分子の重ね合わせ
    print(f"分子を重ね合わせ中...")
    transformed_probe = superimpose_molecules(args.compound, args.probe, atom_mapping)
    
    # 結果を保存
    writer = Chem.SDWriter(args.output)
    writer.write(transformed_probe)
    writer.close()
    
    print(f"重ね合わせ後のプローブ分子を {args.output} に保存しました")


if __name__ == "__main__":
    main()