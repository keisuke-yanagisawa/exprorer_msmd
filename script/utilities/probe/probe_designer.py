#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
プローブ設計機能を提供するモジュール
"""

import glob
import os
import sys

# スクリプトのディレクトリをPythonパスに追加
script_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if script_dir not in sys.path:
    sys.path.append(script_dir)

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs


def design_probe(non_common_part, probe_library_dir):
    """
    非共通部分に合うプローブ分子を設計する

    Parameters
    ----------
    non_common_part : dict
        非共通部分の情報 {'mol': rdkit.Chem.rdchem.Mol, 'atom_indices': list}
    probe_library_dir : str
        プローブライブラリのディレクトリパス

    Returns
    -------
    rdkit.Chem.rdchem.Mol
        設計されたプローブ分子
    """
    # 非共通部分の特徴を抽出
    non_common_mol = non_common_part['mol']
    
    # 非共通部分が空の場合、デフォルトのプローブを返す
    if non_common_mol.GetNumAtoms() == 0:
        print("警告: 非共通部分が空です。デフォルトのプローブを使用します。")
        # プローブライブラリから最初のプローブを取得
        probe_files = glob.glob(os.path.join(probe_library_dir, "*.sdf"))
        if not probe_files:
            raise ValueError(f"プローブライブラリディレクトリ {probe_library_dir} にSDFファイルが見つかりません")
        
        default_probe = Chem.SDMolSupplier(probe_files[0])[0]
        if default_probe is None:
            raise ValueError(f"デフォルトプローブ {probe_files[0]} を読み込めませんでした")
        
        return default_probe
    
    # 非共通部分のフィンガープリントを計算
    fp = AllChem.GetMorganFingerprintAsBitVect(non_common_mol, 2, nBits=1024)
    
    # プローブライブラリから類似したプローブを検索
    best_similarity = 0
    best_probe = None
    
    for probe_file in glob.glob(os.path.join(probe_library_dir, "*.sdf")):
        probe = Chem.SDMolSupplier(probe_file)[0]
        if probe is None:
            continue
        
        # プローブのフィンガープリントを計算
        probe_fp = AllChem.GetMorganFingerprintAsBitVect(probe, 2, nBits=1024)
        
        # 類似度を計算
        similarity = DataStructs.TanimotoSimilarity(fp, probe_fp)
        
        if similarity > best_similarity:
            best_similarity = similarity
            best_probe = probe
    
    # 最も類似度の高いプローブを選択
    if best_probe is None or best_similarity < 0.5:
        # 類似度が低い場合は、非共通部分自体をプローブとして使用
        designed_probe = Chem.Mol(non_common_mol)
        # 必要に応じて水素を追加
        designed_probe = Chem.AddHs(designed_probe)
        
        # 3D座標を生成
        AllChem.EmbedMolecule(designed_probe)
        AllChem.UFFOptimizeMolecule(designed_probe)
    else:
        designed_probe = Chem.Mol(best_probe)
    
    # プローブの最適化
    designed_probe = optimize_probe(designed_probe, non_common_mol)
    
    return designed_probe


def optimize_probe(probe, target_structure):
    """
    プローブをターゲット構造により適合させるための最適化

    Parameters
    ----------
    probe : rdkit.Chem.rdchem.Mol
        最適化するプローブ分子
    target_structure : rdkit.Chem.rdchem.Mol
        ターゲット構造（非共通部分）

    Returns
    -------
    rdkit.Chem.rdchem.Mol
        最適化されたプローブ分子
    """
    # 現在の実装では、プローブをそのまま返す
    # 実際の実装では、以下のような最適化が必要:
    # - 結合長や角度の調整
    # - 官能基の追加・削除・置換
    # - 立体配置の最適化
    # - エネルギー最小化
    
    # 3D座標が無い場合は生成
    if not probe.GetNumConformers():
        probe = Chem.AddHs(probe)
        AllChem.EmbedMolecule(probe)
        AllChem.UFFOptimizeMolecule(probe)
    
    return probe


def visualize_probe(probe, output_file):
    """
    プローブ分子を可視化する

    Parameters
    ----------
    probe : rdkit.Chem.rdchem.Mol
        プローブ分子
    output_file : str
        出力ファイルパス
    """
    # プローブをSDFファイルに保存
    writer = Chem.SDWriter(output_file)
    writer.write(probe)
    writer.close()
    
    print(f"プローブ分子を {output_file} に保存しました")


def main():
    """
    メイン関数
    """
    import argparse

    # 絶対インポートを使用
    from utilities.molecular.mcs_extractor import extract_non_common_parts
    
    parser = argparse.ArgumentParser(description='プローブ設計ツール')
    parser.add_argument('--mol1', required=True, help='1つ目の化合物のSDFファイルパス')
    parser.add_argument('--mol2', required=True, help='2つ目の化合物のSDFファイルパス')
    parser.add_argument('--probe-library', required=True, help='プローブライブラリのディレクトリパス')
    parser.add_argument('--output-dir', default=".", help='出力ディレクトリ')
    args = parser.parse_args()
    
    # 非共通部分を抽出
    non_common_parts = extract_non_common_parts(args.mol1, args.mol2)
    
    # 各非共通部分に対してプローブを設計
    probe1 = design_probe(non_common_parts['mol1'], args.probe_library)
    probe2 = design_probe(non_common_parts['mol2'], args.probe_library)
    
    # 出力ディレクトリを作成
    os.makedirs(args.output_dir, exist_ok=True)
    
    # プローブを可視化
    visualize_probe(probe1, os.path.join(args.output_dir, "probe1.sdf"))
    visualize_probe(probe2, os.path.join(args.output_dir, "probe2.sdf"))


if __name__ == "__main__":
    main()