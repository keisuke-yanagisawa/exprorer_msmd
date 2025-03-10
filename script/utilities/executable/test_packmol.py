import warnings
from pathlib import Path

import pytest
from script.utilities.executable.packmol import Packmol
from subprocess import CalledProcessError

@pytest.fixture
def test_data_dir():
    """Fixture providing the path to test data directory"""
    return Path("script/utilities/executable/test_data")

@pytest.fixture
def test_files(test_data_dir):
    """Fixture providing file paths required for testing"""
    return {
        'protein_pdb': test_data_dir / "tripeptide.pdb",
        'cosolv_pdb': test_data_dir / "A11.pdb",
        'cosolv_mol2': test_data_dir / "A11.mol2"
    }

def test_run_packmol(test_files):
    """Test for normal packing process"""
    packmol = Packmol()
    packmol.set(
        protein_pdb=test_files['protein_pdb'],
        cosolv_pdb=test_files['cosolv_pdb'],
        box_size=20,
        molar=0.5
    )
    packmol.run(seed=1)

@pytest.mark.parametrize("molar,expected_warning", [
    (0, RuntimeWarning),
    (1e-10, RuntimeWarning),
    (10, RuntimeWarning)
])
def test_warning_cases(test_files, molar, expected_warning):
    """Test for warnings with abnormal molar concentrations"""
    packmol = Packmol()
    with pytest.warns(expected_warning):
        packmol.set(
            protein_pdb=test_files['protein_pdb'],
            cosolv_pdb=test_files['cosolv_pdb'],
            box_size=20,
            molar=molar
        )
        packmol.run(seed=1)

def test_extremely_high_molar(test_files):
    """Test for error with extremely high molar concentration"""
    packmol = Packmol()
    with pytest.raises((RuntimeError, CalledProcessError)):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            packmol.set(
                protein_pdb=test_files['protein_pdb'],
                cosolv_pdb=test_files['cosolv_pdb'],
                box_size=20,
                molar=1e10
            )
        packmol.run(seed=1)

@pytest.mark.parametrize("test_case", [
    {
        'desc': "Non-existent protein file",
        'protein_pdb': Path("NOTHING.pdb"),
        'cosolv_pdb': 'cosolv_pdb',
        'expected_error': FileNotFoundError
    },
    {
        'desc': "Non-existent cosolvent file",
        'protein_pdb': 'protein_pdb',
        'cosolv_pdb': Path("NOTHING.pdb"),
        'expected_error': FileNotFoundError
    }
])
def test_missing_files(test_files, test_case):
    """Test for missing file cases"""
    packmol = Packmol()
    cosolv_pdb = test_files[test_case['cosolv_pdb']] if isinstance(test_case['cosolv_pdb'], str) else test_case['cosolv_pdb']
    protein_pdb = test_files[test_case['protein_pdb']] if isinstance(test_case['protein_pdb'], str) else test_case['protein_pdb']
    
    with pytest.raises(test_case['expected_error']):
        packmol.set(
            protein_pdb=protein_pdb,
            cosolv_pdb=cosolv_pdb,
            box_size=20,
            molar=0.1
        )
        packmol.run(seed=1)

@pytest.mark.parametrize("test_case", [
    {
        'desc': "Invalid cosolvent file format",
        'protein_pdb': 'protein_pdb',
        'cosolv_pdb': 'cosolv_mol2',
        'expected_error': ValueError
    },
    {
        'desc': "Cosolvent file contains protein",
        'protein_pdb': 'protein_pdb',
        'cosolv_pdb': 'protein_pdb',
        'expected_error': RuntimeError
    },
    {
        'desc': "Protein file contains cosolvent",
        'protein_pdb': 'cosolv_pdb',
        'cosolv_pdb': 'cosolv_pdb',
        'expected_error': RuntimeError
    }
])
def test_invalid_file_contents(test_files, test_case):
    """Test for invalid file content cases"""
    packmol = Packmol()
    with pytest.raises(test_case['expected_error']):
        packmol.set(
            protein_pdb=test_files[test_case['protein_pdb']],
            cosolv_pdb=test_files[test_case['cosolv_pdb']],
            box_size=20,
            molar=0.1
        )
        packmol.run(seed=1)
