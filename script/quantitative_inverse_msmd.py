#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
半自動Quantitative Inverse MSMDを実現するメインスクリプト
"""

import argparse
import os
import sys

import numpy as np
from rdkit import Chem

from script.utilities.molecular.mcs_extractor import extract_non_common_parts
from script.utilities.molecular.molecule_superimposer import (
    generate_atom_mapping,
    superimpose_molecules,
    validate_and_refine_atom_mapping,
)
from script.utilities.probe.probe_designer import design_probe


def setup_msmd_simulation(probe_file, protein_file, output_dir):
    """
    MSMDシミュレーションの設定を作成する

    Parameters
    ----------
    probe_file : str
        プローブのSDFファイルパス
    protein_file : str
        タンパク質構造のPDBファイルパス
    output_dir : str
        出力ディレクトリ

    Returns
    -------
    dict
        MSMDシミュレーション設定
    """
    # ここでは簡略化のため、設定を辞書として返す
    # 実際の実装では、YAMLファイルなどを生成する
    
    # プローブのIDを取得（ファイル名から拡張子を除いたもの）
    probe_id = os.path.splitext(os.path.basename(probe_file))[0]
    
    # 設定を作成
    settings = {
        "general": {
            "name": f"msmd_{probe_id}",
            "workdir": os.path.join(output_dir, probe_id),
            "iter_index": "1-3",  # 3回の独立試行
            "multiprocessing": True
        },
        "input": {
            "protein": {
                "pdb": protein_file
            },
            "probe": {
                "sdf": probe_file,
                "cid": probe_id
            }
        },
        "probe_profile": {
            "map": "nV",  # 体積マップ
            "profiles": [
                {"name": "anion", "atoms": [" OD1", " OD2", " OE1", " OE2"]},
                {"name": "cation", "atoms": [" NZ ", " NH1", " NH2"]},
                {"name": "donor", "atoms": [" N  "]},
                {"name": "acceptor", "atoms": [" O  "]},
                {"name": "hydrophobic", "atoms": [" CB "]}
            ]
        }
    }
    
    return settings


def run_msmd_simulation(settings):
    """
    MSMDシミュレーションを実行する

    Parameters
    ----------
    settings : dict
        MSMDシミュレーション設定

    Returns
    -------
    str
        シミュレーション結果のディレクトリパス
    """
    # ここでは実際のMSMDシミュレーションは実行せず、
    # 設定を表示するだけにする
    
    print("MSMDシミュレーション設定:")
    print(f"  プロジェクト名: {settings['general']['name']}")
    print(f"  作業ディレクトリ: {settings['general']['workdir']}")
    print(f"  タンパク質: {settings['input']['protein']['pdb']}")
    print(f"  プローブ: {settings['input']['probe']['sdf']}")
    
    # 実際の実装では、以下のようなコマンドを実行する
    # command = f"./exprorer_msmd {settings_file}"
    # subprocess.run(command, shell=True, check=True)
    
    # シミュレーション結果のディレクトリパスを返す
    return settings['general']['workdir']


def extract_residue_environment(simulation_dir):
    """
    プローブ周辺の残基環境を抽出する

    Parameters
    ----------
    simulation_dir : str
        シミュレーション結果のディレクトリパス

    Returns
    -------
    dict
        残基環境の情報
    """
    # ここでは実際の残基環境抽出は行わず、
    # ダミーデータを返す
    
    print(f"残基環境を抽出中: {simulation_dir}")
    
    # 実際の実装では、以下のようなコマンドを実行する
    # command = f"./probe_profile {settings_file}"
    # subprocess.run(command, shell=True, check=True)
    
    # ダミーデータを返す
    return {
        "anion": f"{simulation_dir}/anion_profile.dx",
        "cation": f"{simulation_dir}/cation_profile.dx",
        "donor": f"{simulation_dir}/donor_profile.dx",
        "acceptor": f"{simulation_dir}/acceptor_profile.dx",
        "hydrophobic": f"{simulation_dir}/hydrophobic_profile.dx"
    }


def create_odds_ratio_profile(residue_environment):
    """
    オッズ比に基づくプロファイルを作成する

    Parameters
    ----------
    residue_environment : dict
        残基環境の情報

    Returns
    -------
    dict
        オッズ比プロファイル
    """
    # ここでは実際のオッズ比計算は行わず、
    # 入力をそのまま返す
    
    print("オッズ比プロファイルを作成中")
    
    # 実際の実装では、各残基タイプについてオッズ比を計算する
    
    return residue_environment


def calculate_matching_score(protein_file, profile):
    """
    プロファイルとタンパク質構造の合致度を計算する

    Parameters
    ----------
    protein_file : str
        タンパク質構造のPDBファイルパス
    profile : dict
        プロファイル情報

    Returns
    -------
    float
        合致度スコア
    """
    # ここでは実際の合致度計算は行わず、
    # ダミーデータを返す
    
    print(f"合致度を計算中: {protein_file}")
    
    # 実際の実装では、タンパク質構造とプロファイルの合致度を計算する
    
    # ダミースコアを返す
    import random
    return random.uniform(0, 10)


def estimate_binding_strength(score1, score2):
    """
    合致度スコアの差から結合強度の差を推定する

    Parameters
    ----------
    score1 : float
        化合物1の合致度スコア
    score2 : float
        化合物2の合致度スコア

    Returns
    -------
    float
        結合強度差
    """
    # 合致度スコアの差を計算
    strength_diff = score1 - score2
    
    return strength_diff


def semi_automatic_quantitative_inverse_msmd(compound1_file, compound2_file, protein_file, probe_library_dir, output_dir, only_non_common_rings=False, separate_rings=False):
    """
    半自動Quantitative Inverse MSMDを実行する

    Parameters
    ----------
    compound1_file : str
        1つ目の化合物のSDFファイルパス
    compound2_file : str
        2つ目の化合物のSDFファイルパス
    protein_file : str
        タンパク質構造のPDBファイルパス
    probe_library_dir : str
        プローブライブラリのディレクトリパス
    output_dir : str
        出力ディレクトリ
    only_non_common_rings : bool, optional
        Trueの場合、複合環のうち非共通な環のみをプローブ化する, by default False
    separate_rings : bool, optional
        Trueの場合、非共通な環を別々のプローブとして抽出する, by default False

    Returns
    -------
    dict
        結果情報
    """
    # 出力ディレクトリを作成
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. 2つの化合物間の非共通部分を抽出
    print("ステップ1: 非共通部分の抽出")
    non_common_parts = extract_non_common_parts(compound1_file, compound2_file, only_non_common_rings, separate_rings)
    
    # 2. 各非共通部分に対してプローブを設計
    print("ステップ2: プローブの設計")
    probe1 = design_probe(non_common_parts['mol1'], probe_library_dir)
    probe2 = design_probe(non_common_parts['mol2'], probe_library_dir)
    
    # 3. プローブをファイルに保存
    probe1_file = os.path.join(output_dir, "probe1.sdf")
    probe2_file = os.path.join(output_dir, "probe2.sdf")
    
    writer1 = Chem.SDWriter(probe1_file)
    writer1.write(probe1)
    writer1.close()
    
    writer2 = Chem.SDWriter(probe2_file)
    writer2.write(probe2)
    writer2.close()
    
    # 4. ユーザーによる量子化学計算（手動ステップ）
    print("ステップ3: プローブの量子化学計算")
    print("  プローブ1とプローブ2の量子化学計算を実行してください")
    print("  計算が完了したら、任意のキーを押して続行してください")
    input()
    
    # 5. 化合物とプローブの重ね合わせ
    print("ステップ4: 分子の重ね合わせ")
    
    # 化合物1とプローブ1の原子対応関係を自動生成
    print("  化合物1とプローブ1の原子対応関係を自動生成中...")
    atom_mapping1 = generate_atom_mapping(compound1_file, probe1_file)
    atom_mapping1 = validate_and_refine_atom_mapping(compound1_file, probe1_file, atom_mapping1)
    
    # 化合物1とプローブ1の重ね合わせ
    transformed_probe1 = superimpose_molecules(compound1_file, probe1_file, atom_mapping1)
    transformed_probe1_file = os.path.join(output_dir, "transformed_probe1.sdf")
    
    writer1 = Chem.SDWriter(transformed_probe1_file)
    writer1.write(transformed_probe1)
    writer1.close()
    
    # 化合物2とプローブ2の原子対応関係を自動生成
    print("  化合物2とプローブ2の原子対応関係を自動生成中...")
    atom_mapping2 = generate_atom_mapping(compound2_file, probe2_file)
    atom_mapping2 = validate_and_refine_atom_mapping(compound2_file, probe2_file, atom_mapping2)
    
    # 化合物2とプローブ2の重ね合わせ
    transformed_probe2 = superimpose_molecules(compound2_file, probe2_file, atom_mapping2)
    transformed_probe2_file = os.path.join(output_dir, "transformed_probe2.sdf")
    
    writer2 = Chem.SDWriter(transformed_probe2_file)
    writer2.write(transformed_probe2)
    writer2.close()
    
    # 6. inverse MSMDシミュレーションの実行
    print("ステップ5: inverse MSMDシミュレーションの実行")
    
    # プローブ1のinverse MSMDシミュレーション
    print("  プローブ1のinverse MSMDシミュレーションを実行中...")
    probe1_settings = setup_msmd_simulation(transformed_probe1_file, protein_file, output_dir)
    probe1_simulation_dir = run_msmd_simulation(probe1_settings)
    
    # プローブ2のinverse MSMDシミュレーション
    print("  プローブ2のinverse MSMDシミュレーションを実行中...")
    probe2_settings = setup_msmd_simulation(transformed_probe2_file, protein_file, output_dir)
    probe2_simulation_dir = run_msmd_simulation(probe2_settings)
    
    # 7. プロファイル作成
    print("ステップ6: プロファイルの作成")
    
    # プローブ1のプロファイル作成
    probe1_residue_env = extract_residue_environment(probe1_simulation_dir)
    profile1 = create_odds_ratio_profile(probe1_residue_env)
    
    # プローブ2のプロファイル作成
    probe2_residue_env = extract_residue_environment(probe2_simulation_dir)
    profile2 = create_odds_ratio_profile(probe2_residue_env)
    
    # 8. 合致度計算
    print("ステップ7: 合致度計算と結合強度推定")
    
    # 各プロファイルの合致度を計算
    score1 = calculate_matching_score(protein_file, profile1)
    score2 = calculate_matching_score(protein_file, profile2)
    
    # 9. 結合強度の推定
    strength_diff = estimate_binding_strength(score1, score2)
    
    # 10. 結果の出力
    print("ステップ8: 結果出力")
    print(f"  化合物1の合致度スコア: {score1}")
    print(f"  化合物2の合致度スコア: {score2}")
    print(f"  結合強度差: {strength_diff}")
    
    if strength_diff > 0:
        print("  結論: 化合物1の方が強く結合すると予測されます")
    elif strength_diff < 0:
        print("  結論: 化合物2の方が強く結合すると予測されます")
    else:
        print("  結論: 両化合物の結合強度は同程度と予測されます")
    
    # 結果をファイルに保存
    result_file = os.path.join(output_dir, "quantitative_inverse_msmd_result.txt")
    with open(result_file, "w") as f:
        f.write(f"化合物1: {compound1_file}\n")
        f.write(f"化合物2: {compound2_file}\n")
        f.write(f"タンパク質: {protein_file}\n")
        f.write(f"プローブ1: {transformed_probe1_file}\n")
        f.write(f"プローブ2: {transformed_probe2_file}\n")
        f.write(f"化合物1の合致度スコア: {score1}\n")
        f.write(f"化合物2の合致度スコア: {score2}\n")
        f.write(f"結合強度差: {strength_diff}\n")
        
        if strength_diff > 0:
            f.write("結論: 化合物1の方が強く結合すると予測されます\n")
        elif strength_diff < 0:
            f.write("結論: 化合物2の方が強く結合すると予測されます\n")
        else:
            f.write("結論: 両化合物の結合強度は同程度と予測されます\n")
    
    print(f"結果を {result_file} に保存しました")
    
    return {
        'probe1': probe1,
        'probe2': probe2,
        'profile1': profile1,
        'profile2': profile2,
        'score1': score1,
        'score2': score2,
        'strength_diff': strength_diff
    }


def main():
    """
    メイン関数
    """
    parser = argparse.ArgumentParser(description='半自動Quantitative Inverse MSMD')
    parser.add_argument('--compound1', required=True, help='1つ目の化合物のSDFファイルパス')
    parser.add_argument('--compound2', required=True, help='2つ目の化合物のSDFファイルパス')
    parser.add_argument('--protein', required=True, help='タンパク質構造のPDBファイルパス')
    parser.add_argument('--probe-library', required=True, help='プローブライブラリのディレクトリパス')
    parser.add_argument('--output-dir', default="./output", help='出力ディレクトリ')
    parser.add_argument('--only-non-common-rings', action='store_true',
                        help='複合環のうち非共通な環のみをプローブ化する')
    parser.add_argument('--separate-rings', action='store_true',
                        help='非共通な環を別々のプローブとして抽出する')
    args = parser.parse_args()
    
    # 半自動Quantitative Inverse MSMDを実行
    semi_automatic_quantitative_inverse_msmd(
        args.compound1,
        args.compound2,
        args.protein,
        args.probe_library,
        args.output_dir,
        args.only_non_common_rings,
        args.separate_rings
    )


if __name__ == "__main__":
    main()