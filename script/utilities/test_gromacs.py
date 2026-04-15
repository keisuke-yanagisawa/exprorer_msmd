import numpy as np
import pytest
from scipy import constants

from script.utilities.gromacs import Gro, GroAtom


class TestGroAtom:
    def test_default_initialization(self):
        """Test for creating an empty GroAtom instance"""
        atom = GroAtom()
        assert atom.resi == -1
        assert atom.resn == ""
        assert atom.atomtype == ""
        assert atom.atom_id == -1
        assert atom.point == pytest.approx(np.zeros(3))
        assert atom.velocity == pytest.approx(np.zeros(3))
        assert atom.comment == ""
        assert atom.atomic_mass == pytest.approx(0.0)

    def test_set_attributes(self):
        """Test for setting GroAtom attributes"""
        atom = GroAtom()
        atom.resi = 1
        atom.resn = "WAT"
        atom.atomtype = "OW"
        atom.atom_id = 1
        atom.point = np.array([0.5, 0.5, 0.5])
        atom.atomic_mass = 16.0

        assert atom.resi == 1
        assert atom.resn == "WAT"
        assert atom.atomtype == "OW"
        assert atom.atom_id == 1
        np.testing.assert_array_almost_equal(atom.point, [0.5, 0.5, 0.5])
        assert atom.atomic_mass == pytest.approx(16.0)

    def test_from_parmed(self):
        """Test creating GroAtom from ParmEd atom (via Gro file parsing)"""
        gro_content = "Test\n1\n    1WAT     OW    1   0.500   0.600   0.700\n   5.0   5.0   5.0\n"
        gro_file = _write_tmp_gro(gro_content)
        gro = Gro(str(gro_file))
        atoms = gro.atoms

        assert len(atoms) == 1
        atom = atoms[0]
        assert atom.resi == 1
        assert atom.resn == "WAT"
        assert atom.atomtype == "OW"
        assert atom.atom_id == 1
        np.testing.assert_array_almost_equal(atom.point, [0.5, 0.6, 0.7], decimal=3)
        assert atom.atomic_mass == pytest.approx(15.999, rel=1e-2)


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
        assert gro.box_size == [0.0, 0.0, 0.0]
        assert len(gro.atoms) == 0

    def test_file_parsing(self, mock_gro_file):
        """Test for GRO file parsing"""
        gro = Gro(str(mock_gro_file))
        assert gro.description == "Simple water system"
        assert gro.natoms == 3
        assert gro.box_size == pytest.approx([5.0, 5.0, 5.0])
        assert len(gro.atoms) == 3
        assert [atom.atomtype for atom in gro.atoms] == ["OW", "HW1", "HW2"]

    @pytest.mark.parametrize(
        "search_params, expected_count",
        [
            ({"resi": 1}, 3),
            ({"resn": "WAT"}, 3),
            ({"atomtype": "OW"}, 1),
            ({"atom_id": 2}, 1),
        ],
    )
    def test_get_atoms(self, mock_gro_file, search_params, expected_count):
        """Test for atom search functionality"""
        gro = Gro(str(mock_gro_file))
        atoms = gro.get_atoms(**search_params)
        assert len(atoms) == expected_count

    def test_add_atom(self, mock_gro_file):
        """Test for atom addition functionality"""
        gro = Gro(str(mock_gro_file))
        new_atom = GroAtom()
        new_atom.resi = 2
        new_atom.resn = "WAT"
        new_atom.atomtype = "OW"
        new_atom.point = np.array([1.0, 1.0, 1.0])
        new_atom.atomic_mass = 16.0
        gro.add_atom(new_atom)

        assert gro.natoms == 4
        # Verify the atom was added to the correct residue
        res2_atoms = gro.get_atoms(resi=2)
        assert len(res2_atoms) == 1
        assert res2_atoms[0].atomtype == "OW"

    def test_add_invalid_atom(self):
        """Test for invalid atom addition"""
        gro = Gro()
        with pytest.raises(TypeError, match="the input is NON-GRO_ATOM"):
            gro.add_atom("not a GroAtom")

    def test_molar_concentration(self, mock_gro_file):
        """Test for molar concentration calculation"""
        gro = Gro(str(mock_gro_file))
        molar = gro.molar("WAT")
        # Box size 5.0 nm x 5.0 nm x 5.0 nm = 125 nm3
        expected_molar = (1 / constants.N_A) / ((5 * constants.nano) ** 3 * 10**3)  # nm3 -> cm3 * 10^3 = L
        assert molar == pytest.approx(expected_molar, rel=1e-5)

    def test_string_representation(self, mock_gro_file):
        """Test for GRO format output via ParmEd"""
        gro = Gro(str(mock_gro_file))
        result = str(gro)
        lines = result.strip().split("\n")

        # Description preserved
        assert lines[0] == "Simple water system"
        # Atom count
        assert int(lines[1].strip()) == 3
        # Verify atoms are present (ParmEd standard GRO format)
        assert "OW" in lines[2]
        assert "HW1" in lines[3]
        assert "HW2" in lines[4]
        # Box values
        box_vals = [float(s) for s in lines[5].split()]
        assert box_vals == pytest.approx([5.0, 5.0, 5.0])


def _write_tmp_gro(content: str):
    """Helper to write GRO content to a temp file and return path."""
    import tempfile
    from pathlib import Path

    f = tempfile.NamedTemporaryFile(suffix=".gro", mode="w", delete=False)
    f.write(content)
    f.close()
    return Path(f.name)
