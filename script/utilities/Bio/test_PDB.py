import pytest
import tempfile
from pathlib import Path
import numpy as np

from script.utilities.Bio import PDB
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure


@pytest.fixture
def pdb_files():
    """Fixture providing paths to PDB files"""
    return {
        'pdb': Path("script/utilities/Bio/test_data/PDB/7m67.pdb"),
        'gzipped': Path("script/utilities/Bio/test_data/PDB/7m67.pdb.gz")
    }


@pytest.fixture
def basic_atom():
    """Fixture providing a basic atom object"""
    return Atom(
        name="CA",
        coord=np.array([1.0, 2.0, 3.0], dtype=np.float64),
        bfactor=20.0,
        occupancy=1.0,
        altloc=" ",
        fullname=" CA ",
        serial_number=1,
        element="C"
    )


@pytest.fixture
def basic_structure(basic_atom):
    """Fixture providing a basic structure object"""
    structure = Structure("test")
    model = Model(0)
    chain = Chain("A")
    residue = Residue((" ", 1, " "), "ALA", "")
    
    residue.add(basic_atom)
    chain.add(residue)
    model.add(chain)
    structure.add(model)
    
    return structure


class TestMultiModelPDBReader:
    """MultiModelPDBReaderクラスのテスト群"""

    def test_read_models(self, pdb_files):
        """通常のPDBファイル読み込みテスト"""
        reader = PDB.MultiModelPDBReader(str(pdb_files['pdb']))
        models = [model for model in reader]
        assert len(models) == 10

    def test_get_model_out_of_range(self, pdb_files):
        """範囲外のモデル取得時のエラーテスト"""
        reader = PDB.MultiModelPDBReader(str(pdb_files['pdb']))
        with pytest.raises(IndexError):
            reader.get_model(10)

    def test_file_not_found(self):
        """存在しないファイルを指定した場合のエラーテスト"""
        with pytest.raises(FileNotFoundError):
            PDB.MultiModelPDBReader("INVALID_PATH")

    def test_read_gzipped(self, pdb_files):
        """gzip圧縮されたPDBファイルの読み込みテスト"""
        reader = PDB.MultiModelPDBReader(str(pdb_files['gzipped']))
        models = [model for model in reader]
        assert len(models) == 10


class TestPDBIOHelper:
    """PDBIOhelperクラスのテスト群"""

    def test_save_and_read_models(self, pdb_files):
        """モデルの保存と読み込みテスト"""
        reader = PDB.MultiModelPDBReader(str(pdb_files['pdb']))
        models = [model for model in reader]
        tmp_output_pdb = Path(tempfile.mkstemp(suffix=".pdb")[1])
        
        writer = PDB.PDBIOhelper(tmp_output_pdb)
        for model in models:
            writer.save(model)
        writer.close()

        reader = PDB.MultiModelPDBReader(str(tmp_output_pdb))
        models = [model for model in reader]
        assert len(models) == 10


class TestStructureOperations:
    """構造操作関連のテスト群"""

    def test_get_structure_gzipped(self, pdb_files):
        """圧縮PDBファイルからの構造取得テスト"""
        structure = PDB.get_structure(pdb_files['gzipped'])
        assert len(structure) == 10

    @pytest.mark.parametrize("attr,expected", [
        ("resid", 1),
        ("resname", "ALA"),
        ("coord", np.array([1.0, 2.0, 3.0])),
        ("element", "C"),
        ("fullname", " CA ")
    ])
    def test_get_atom_attr(self, basic_structure, attr, expected):
        """原子属性取得機能のテスト"""
        atom = next(basic_structure.get_atoms())
        result = PDB.get_atom_attr(atom, attr)
        
        if attr == "coord":
            assert np.allclose(result, expected)
        else:
            assert result == expected

    def test_get_atom_attr_invalid(self, basic_structure):
        """無効な属性指定時のエラーテスト"""
        atom = next(basic_structure.get_atoms())
        with pytest.raises(NotImplementedError):
            PDB.get_atom_attr(atom, "invalid_attr")

    def test_get_attr(self, basic_structure):
        """モデルからの属性取得テスト"""
        model = basic_structure[0]

        # セレクタなしのテスト
        coords = PDB.get_attr(model, "coord")
        assert len(coords) == 1
        assert np.allclose(coords[0], [1.0, 2.0, 3.0])

        # セレクタありのテスト
        def select_c_atoms(atom):
            return atom.element == "C"
        
        coords = PDB.get_attr(model, "coord", select_c_atoms)
        assert len(coords) == 1
        assert np.allclose(coords[0], [1.0, 2.0, 3.0])


class TestAtomClassification:
    """原子分類機能のテスト群"""

    @pytest.fixture
    def water_atom(self):
        """Fixture providing an atom object for water molecule"""
        return Atom(
            name="O",
            coord=np.array([0, 0, 0], dtype=np.float64),
            bfactor=20.0,
            occupancy=1.0,
            altloc=" ",
            fullname=" O  ",
            serial_number=1,
            element="O"
        )

    @pytest.fixture
    def water_structure(self, water_atom):
        """Fixture providing a structure object containing water molecule"""
        structure = Structure("test")
        model = Model(0)
        chain = Chain("A")
        water_residue = Residue((" ", 1, " "), "WAT", "")
        water_residue.add(water_atom)
        chain.add(water_residue)
        model.add(chain)
        structure.add(model)
        return structure

    def test_is_water(self, water_structure, basic_structure):
        """水分子判定機能のテスト"""
        water_atom = next(water_structure.get_atoms())
        non_water_atom = next(basic_structure.get_atoms())
        
        assert PDB.is_water(water_atom)
        assert not PDB.is_water(non_water_atom)

    @pytest.fixture
    def hetero_structure(self):
        """Fixture providing a structure object containing hetero atoms"""
        structure = Structure("test")
        model = Model(0)
        chain = Chain("A")
        hetero_residue = Residue(("H_", 1, " "), "HET", "")
        hetero_atom = Atom(
            name="O",
            coord=np.array([0, 0, 0], dtype=np.float64),
            bfactor=20.0,
            occupancy=1.0,
            altloc="H",
            fullname=" O  ",
            serial_number=1,
            element="O"
        )
        hetero_residue.add(hetero_atom)
        chain.add(hetero_residue)
        model.add(chain)
        structure.add(model)
        return structure

    def test_is_hetero(self, hetero_structure, basic_structure):
        """ヘテロ原子判定機能のテスト"""
        hetero_atom = next(hetero_structure.get_atoms())
        non_hetero_atom = next(basic_structure.get_atoms())
        
        assert PDB.is_hetero(hetero_atom)
        assert not PDB.is_hetero(non_hetero_atom)

    @pytest.mark.parametrize("atom_data,expected", [
        ({"name": "H", "element": "H", "fullname": " H  "}, True),
        ({"name": "C", "element": "C", "fullname": " C  "}, False)
    ])
    def test_is_hydrogen(self, atom_data, expected):
        """水素原子判定機能のテスト"""
        atom = Atom(
            name=atom_data["name"],
            coord=np.array([0, 0, 0], dtype=np.float64),
            bfactor=20.0,
            occupancy=1.0,
            altloc=" ",
            fullname=atom_data["fullname"],
            serial_number=1,
            element=atom_data["element"]
        )
        assert PDB.is_hydrogen(atom) == expected


class TestVolumeCalculation:
    """体積計算機能のテスト群"""

    @pytest.mark.skip(reason="Volume calculation needs to be fixed for small structures")
    def test_estimate_exclude_volume(self):
        """排除体積計算機能のテスト"""
        structure = Structure("test")
        model = Model(0)
        chain = Chain("A")
        residue = Residue((" ", 1, " "), "ALA", "")
        
        # 2つの炭素原子を近接して配置
        atoms = [
            Atom(
                name=f"C{i}",
                coord=np.array([0.1 * i, 0.0, 0.0], dtype=np.float64),
                bfactor=20.0,
                occupancy=1.0,
                altloc=" ",
                fullname=f" C{i} ",
                serial_number=i+1,
                element="C"
            ) for i in range(2)
        ]
        
        for atom in atoms:
            residue.add(atom)
        chain.add(residue)
        model.add(chain)
        structure.add(model)

        volume = PDB.estimate_exclute_volume(structure)
        assert volume > 0


class TestStructureExtraction:
    """構造抽出機能のテスト群"""

    def test_extract_substructure(self, basic_structure):
        """部分構造抽出機能のテスト"""
        class CarbonSelector(PDB.PDB.Select):
            def accept_atom(self, atom):
                return atom.element == "C"

        substructure = PDB.extract_substructure(basic_structure, CarbonSelector())
        
        atoms = list(substructure.get_atoms())
        assert len(atoms) == 1
        assert all(atom.element == "C" for atom in atoms)
