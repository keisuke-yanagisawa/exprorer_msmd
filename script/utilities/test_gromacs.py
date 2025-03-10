import pytest
import numpy as np
import warnings
from pathlib import Path

from script.utilities.gromacs import GroAtom, Gro, ATOMIC_NUMBER, ATOMIC_WEIGHT
from scipy import constants

class TestGroAtom:
    @pytest.fixture
    def empty_atom(self):
        """Fixture providing an empty GroAtom instance"""
        return GroAtom()

    @pytest.fixture
    def normal_atom(self):
        """Fixture providing a GroAtom instance with normal atom information"""
        return GroAtom("    1WAT     OW    1   0.000   0.000   0.000")

    def test_empty_initialization(self, empty_atom):
        """Test for creating an empty GroAtom instance

        Expected behavior:
        - All attributes should be initialized with appropriate default values
        """
        assert empty_atom.resi == -1
        assert empty_atom.resn == ""
        assert empty_atom.atomtype == ""
        assert empty_atom.atom_id == -1
        assert empty_atom.point == pytest.approx(np.zeros(3))
        assert empty_atom.velocity == pytest.approx(np.zeros(3))
        assert empty_atom.comment == ""
        assert empty_atom.atomic_mass == pytest.approx(0.0)

    @pytest.mark.parametrize("atom_str, expected", [
        (
            "    1WAT     OW    1   0.000   0.000   0.000",
            {
                "resi": 1,
                "resn": "WAT",
                "atomtype": "OW",
                "atom_id": 1,
                "point": np.array([0.0, 0.0, 0.0]),
                "velocity": np.array([0.0, 0.0, 0.0]),
                "comment": "",
                "atomic_mass": ATOMIC_WEIGHT[ATOMIC_NUMBER["O"]]
            }
        ),
        (
            "    1WAT     OW    1   0.000   0.000   0.000 ; water oxygen",
            {
                "resi": 1,
                "resn": "WAT",
                "atomtype": "OW",
                "atom_id": 1,
                "point": np.array([0.0, 0.0, 0.0]),
                "velocity": np.array([0.0, 0.0, 0.0]),
                "comment": "water oxygen",
                "atomic_mass": ATOMIC_WEIGHT[ATOMIC_NUMBER["O"]]
            }
        ),
        (
            "    1WAT     OW    1   0.000   0.000   0.000   1.000   2.000   3.000",
            {
                "resi": 1,
                "resn": "WAT",
                "atomtype": "OW",
                "atom_id": 1,
                "point": np.array([0.0, 0.0, 0.0]),
                "velocity": np.array([1.0, 2.0, 3.0]),
                "comment": "",
                "atomic_mass": ATOMIC_WEIGHT[ATOMIC_NUMBER["O"]]
            }
        )
    ])
    def test_parse_atom_variations(self, atom_str, expected):
        """Test for parsing various formats of atom information

        Parameters:
        - atom_str: Atom information string to parse
        - expected: Dictionary of expected attribute values

        Expected behavior:
        - All attributes should be correctly set to their expected values
        """
        atom = GroAtom(atom_str)
        assert atom.resi == expected["resi"]
        assert atom.resn == expected["resn"]
        assert atom.atomtype == expected["atomtype"]
        assert atom.atom_id == expected["atom_id"]
        np.testing.assert_array_almost_equal(atom.point, expected["point"])
        np.testing.assert_array_almost_equal(atom.velocity, expected["velocity"])
        assert atom.comment == expected["comment"]
        assert atom.atomic_mass == pytest.approx(expected["atomic_mass"])

    @pytest.mark.parametrize("invalid_str", [
        "    1WAT     OW    1   0.000   0.000",  # Missing coordinates
        "    1WAT     OW    1",  # No coordinates
    ])
    def test_parse_invalid_coordinates(self, invalid_str):
        """Test for parsing invalid coordinate information

        Parameters:
        - invalid_str: Invalid atom information string

        Expected behavior:
        - Should raise RuntimeError
        """
        with pytest.raises(RuntimeError, match="the dimension of atom coordinates/velocities are wrong"):
            GroAtom(invalid_str)

    def test_unknown_atomtype(self):
        """Test for handling unknown atom types

        Expected behavior:
        - Should raise RuntimeWarning
        - atomic_mass should be set to default value
        """
        warning_message = "atomtype XX is not matched to any atom names"
        with pytest.warns(RuntimeWarning, match=warning_message):
            atom = GroAtom("    1WAT     XX    1   0.000   0.000   0.000")
            assert atom.atomic_mass == pytest.approx(ATOMIC_WEIGHT[-1])

    def test_string_representation(self, normal_atom):
        """Test for string representation generation

        Expected behavior:
        - Should generate correctly formatted string
        """
        expected = "    1  WAT   OW    1   0.000   0.000   0.000"
        assert str(normal_atom) == expected

class TestGro:
    @pytest.fixture
    def gro_content(self):
        """Fixture providing GRO file content for testing"""
        return """Simple water system
3
    1WAT     OW    1   0.000   0.000   0.000
    1WAT    HW1    2   0.100   0.000   0.000
    1WAT    HW2    3   0.000   0.100   0.000
   5.0   5.0   5.0
"""

    @pytest.fixture
    def mock_gro_file(self, gro_content, tmp_path):
        """Fixture creating a temporary GRO file"""
        gro_file = tmp_path / "test.gro"
        gro_file.write_text(gro_content)
        return gro_file

    def test_empty_initialization(self):
        """Test for empty Gro instance creation"""
        gro = Gro()
        assert gro.description == ""
        assert gro.natoms == 0
        assert gro.box_size == [0, 0, 0]
        assert len(gro.atoms) == 0

    def test_file_parsing(self, mock_gro_file):
        """Test for GRO file parsing"""
        gro = Gro(str(mock_gro_file))
        assert gro.description == "Simple water system"
        assert gro.natoms == 3
        assert gro.box_size == [5.0, 5.0, 5.0]
        assert len(gro.atoms) == 3
        assert [atom.atomtype for atom in gro.atoms] == ["OW", "HW1", "HW2"]

    @pytest.mark.parametrize("search_params, expected_count", [
        ({"resi": 1}, 3),
        ({"resn": "WAT"}, 3),
        ({"atomtype": "OW"}, 1),
        ({"atom_id": 2}, 1),
        ({"resi": 1, "atomtype": "OW"}, 1),
    ])
    def test_get_atoms(self, mock_gro_file, search_params, expected_count):
        """Test for atom search functionality"""
        gro = Gro(str(mock_gro_file))
        atoms = gro.get_atoms(**search_params)
        assert len(atoms) == expected_count

    def test_add_atom(self, mock_gro_file):
        """Test for atom addition functionality"""
        gro = Gro(str(mock_gro_file))
        new_atom = GroAtom("    2  WAT   OW    4   1.000   1.000   1.000")
        gro.add_atom(new_atom)
        
        assert gro.natoms == 4
        assert gro.atoms[-1].resi == 2
        assert gro.atoms[-1].atom_id == 4

    def test_add_invalid_atom(self):
        """Test for invalid atom addition"""
        gro = Gro()
        with pytest.raises(TypeError, match="the input is NON-GRO_ATOM"):
            gro.add_atom("not a GroAtom")

    def test_molar_concentration(self, mock_gro_file):
        """Test for molar concentration calculation"""
        gro = Gro(str(mock_gro_file))
        molar = gro.molar("WAT")
        # Box size 5.0 nm × 5.0 nm × 5.0 nm = 125 nm3
        expected_molar = (1 / constants.N_A) / ((5*constants.nano)**3 * 10**3)  # nm3 -> cm3 * 10^3 = L
        assert molar == pytest.approx(expected_molar, rel=1e-5)

    def test_string_representation(self, mock_gro_file):
        """Test for string representation generation"""
        gro = Gro(str(mock_gro_file))
        expected = """Simple water system
     3
    1  WAT   OW    1   0.000   0.000   0.000
    1  WAT  HW1    2   0.100   0.000   0.000
    1  WAT  HW2    3   0.000   0.100   0.000
   5.00000     5.00000     5.00000
"""
        assert str(gro) == expected