"""
quantitative inverse MSMD のテスト

1. data/4HW3A_lig.sdf と data/4HW2A_lig.sdf を参照して、どの部分構造をプローブとすべきかを判定する
  1-1. 4WH2_ligに対しては O=C([O-])c1cc2ccc(Cl)cc2[nH]1 がプローブとして選択される
  1-2. 4WH3_ligに対しては O=C([O-])c1cc2ccccc2s1 がプローブとして選択される
"""

import os
import sys

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

# プロジェクトのルートディレクトリをパスに追加
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../')))

from mcs_extractor import extract_non_common_parts
from probe_designer import design_probe


# テストデータのパスを設定
@pytest.fixture
def data_paths():
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    return {
        'hw2_lig_file': os.path.join(data_dir, '4HW2A_lig.sdf'),
        'hw3_lig_file': os.path.join(data_dir, '4HW3A_lig.sdf'),
        'probe_library_dir': data_dir
    }

def test_probe_selection(data_paths):
    """プローブ選択のテスト"""
    # 1. 非共通部分の抽出
    non_common_parts = extract_non_common_parts(
        data_paths['hw2_lig_file'], 
        data_paths['hw3_lig_file']
    )
    
    # 2. 各非共通部分に対してプローブを設計
    probe_hw2 = design_probe(non_common_parts['mol1'], data_paths['probe_library_dir'])
    probe_hw3 = design_probe(non_common_parts['mol2'], data_paths['probe_library_dir'])
    
    # 3. プローブのSMILES表現を取得
    probe_hw2_smiles = Chem.MolToSmiles(probe_hw2)
    probe_hw3_smiles = Chem.MolToSmiles(probe_hw3)
    
    # 4. 期待値との比較
    expected_hw2_smiles = "O=C(O)c1cc2ccc(Cl)cc2[nH]1"
    expected_hw3_smiles = "O=C(O)c1cc2ccccc2s1"
    
    # 正規化して比較（SMILESの表現が異なる場合があるため）
    probe_hw2_mol = Chem.MolFromSmiles(probe_hw2_smiles)
    probe_hw3_mol = Chem.MolFromSmiles(probe_hw3_smiles)
    expected_hw2_mol = Chem.MolFromSmiles(expected_hw2_smiles)
    expected_hw3_mol = Chem.MolFromSmiles(expected_hw3_smiles)
    
    probe_hw2_smiles_canonical = Chem.MolToSmiles(probe_hw2_mol)
    probe_hw3_smiles_canonical = Chem.MolToSmiles(probe_hw3_mol)
    expected_hw2_smiles_canonical = Chem.MolToSmiles(expected_hw2_mol)
    expected_hw3_smiles_canonical = Chem.MolToSmiles(expected_hw3_mol)
    
    # 結果を出力
    print(f"4WH2_ligのプローブ: {probe_hw2_smiles_canonical}")
    print(f"4WH3_ligのプローブ: {probe_hw3_smiles_canonical}")
    
    # アサーション
    assert probe_hw2_smiles_canonical == expected_hw2_smiles_canonical, \
        f"4WH2_ligのプローブが期待値と異なります。\n期待値: {expected_hw2_smiles_canonical}\n実際: {probe_hw2_smiles_canonical}"
    assert probe_hw3_smiles_canonical == expected_hw3_smiles_canonical, \
        f"4WH3_ligのプローブが期待値と異なります。\n期待値: {expected_hw3_smiles_canonical}\n実際: {probe_hw3_smiles_canonical}"
