import pytest
import numpy as np
from pathlib import Path

from script.utilities import gromacs
from script.addvirtatom2gro import center_of_mass, addvirtatom2gro


@pytest.fixture
def single_carbon_atom():
    """Test data for a single carbon atom"""
    atom = gromacs.GroAtom()
    atom.point = np.array([1.0, 1.0, 1.0])
    atom.atomic_mass = 12.0  # Carbon atom
    return atom


@pytest.fixture
def water_molecule():
    """Test data for a water molecule (O + 2H)"""
    oxygen = gromacs.GroAtom()
    oxygen.point = np.array([0.0, 0.0, 0.0])
    oxygen.atomic_mass = 16.0  # Oxygen atom

    hydrogen1 = gromacs.GroAtom()
    hydrogen1.point = np.array([0.1, 0.0, 0.0])
    hydrogen1.atomic_mass = 1.0  # Hydrogen atom

    hydrogen2 = gromacs.GroAtom()
    hydrogen2.point = np.array([0.0, 0.1, 0.0])
    hydrogen2.atomic_mass = 1.0  # Hydrogen atom

    return [oxygen, hydrogen1, hydrogen2]


@pytest.fixture
def co_molecule():
    """Test data for a carbon monoxide molecule (C + O)"""
    carbon = gromacs.GroAtom()
    carbon.point = np.array([0.0, 0.0, 0.0])
    carbon.atomic_mass = 12.0  # Carbon atom

    oxygen = gromacs.GroAtom()
    oxygen.point = np.array([1.0, 0.0, 0.0])
    oxygen.atomic_mass = 16.0  # Oxygen atom

    return [carbon, oxygen]


def test_single_atom(single_carbon_atom):
    """Test for center of mass calculation of a single atom"""
    result = center_of_mass([single_carbon_atom])
    np.testing.assert_array_almost_equal(result, single_carbon_atom.point)


def test_multiple_atoms(water_molecule):
    """Test for center of mass calculation of multiple atoms (water molecule)"""
    result = center_of_mass(water_molecule)
    # Manual calculation of water molecule's center of mass
    expected = (16.0 * water_molecule[0].point +
               1.0 * water_molecule[1].point +
               1.0 * water_molecule[2].point) / 18.0
    np.testing.assert_array_almost_equal(result, expected)


def test_different_masses(co_molecule):
    """Test for center of mass calculation of atoms with different masses"""
    result = center_of_mass(co_molecule)
    # Center of mass is closer to the heavier atom (oxygen)
    expected = np.array([16.0 / (12.0 + 16.0), 0.0, 0.0])
    np.testing.assert_array_almost_equal(result, expected)


@pytest.fixture
def basic_gro():
    """GRO file containing a single water molecule"""
    return """Simple water system
3
    1WAT     OW    1   0.000   0.000   0.000
    1WAT    HW1    2   0.100   0.000   0.000
    1WAT    HW2    3   0.000   0.100   0.000
   5.0   5.0   5.0
"""


@pytest.fixture
def multi_probe_gro():
    """GRO file containing multiple water molecules"""
    return """System with multiple probes
6
    1WAT     OW    1   0.000   0.000   0.000
    1WAT    HW1    2   0.100   0.000   0.000
    1WAT    HW2    3   0.000   0.100   0.000
    2WAT     OW    4   1.000   1.000   1.000
    2WAT    HW1    5   1.100   1.000   1.000
    2WAT    HW2    6   1.000   1.100   1.000
   5.0   5.0   5.0
"""


@pytest.fixture
def empty_gro():
    """Empty GRO file"""
    return """Empty system
0
5.00000     5.00000     5.00000
"""


def test_basic_conversion(basic_gro):
    """Test for basic GRO file conversion"""
    result = addvirtatom2gro(basic_gro, "WAT")
    
    # Verify virtual atom is added
    assert "VIS" in result
    # Verify virtual atom coordinates are near center of mass
    assert "0.006   0.006   0.000" in result  # Water molecule's center of mass


def test_multiple_probes(multi_probe_gro):
    """Test for GRO file conversion with multiple probe molecules"""
    result = addvirtatom2gro(multi_probe_gro, "WAT")
    
    # Verify two virtual atoms are added
    vis_count = result.count("VIS")
    assert vis_count == 2
    
    # Verify virtual atoms are placed at center of mass of both water molecules
    assert "0.006   0.006   0.000" in result  # First water molecule
    assert "1.006   1.006   1.000" in result  # Second water molecule


def test_no_probe(basic_gro):
    """Test for when probe does not exist"""
    result = addvirtatom2gro(basic_gro, "ETH")  # Non-existent probe ID
    
    # Verify no virtual atoms are added
    assert "VIS" not in result
    # Verify original GRO file content is preserved
    assert len(result.split("\n")) == len(basic_gro.split("\n"))


def test_empty_gro(empty_gro):
    """Test for processing empty GRO file"""
    result = addvirtatom2gro(empty_gro, "WAT")
    
    # Compare normalized expected and actual output
    def normalize_gro(gro_str):
        lines = gro_str.strip().split("\n")
        return "\n".join([
            lines[0],  # System name
            lines[1].strip(),  # Number of atoms
            lines[2].strip()  # Box size
        ])
    assert normalize_gro(result) == normalize_gro(empty_gro)