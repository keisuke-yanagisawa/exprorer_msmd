import glob
import json
import os
import tempfile
from pathlib import Path
from unittest import TestCase

import pytest

from script.utilities.executable.parmchk import Parmchk

_TESTDATA_DIR = Path("script/utilities/executable/test_data")


def _get_ambertools_version() -> str:
    """Get the installed AmberTools version from conda metadata."""
    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    meta_files = glob.glob(f"{conda_prefix}/conda-meta/ambertools-*.json")
    if meta_files:
        with open(meta_files[0]) as f:
            return json.load(f).get("version", "unknown")
    return "unknown"


AMBERTOOLS_VERSION = _get_ambertools_version()
AMBERTOOLS_EXPECTED_DIR = _TESTDATA_DIR / f"ambertools_{AMBERTOOLS_VERSION}"


class TestParmchk(TestCase):
    def __init__(self, *args, **kwargs):
        super(TestParmchk, self).__init__(*args, **kwargs)
        self.probe_mol2 = _TESTDATA_DIR / "A11.mol2"
        self.expected_probe_frcmod = AMBERTOOLS_EXPECTED_DIR / "A11.frcmod"

    @pytest.mark.skipif(
        not AMBERTOOLS_EXPECTED_DIR.is_dir(),
        reason=f"No expected data for AmberTools {AMBERTOOLS_VERSION} (missing {AMBERTOOLS_EXPECTED_DIR})"
    )
    def test_parmchk(self):
        tmp_frcmod = Path(tempfile.mkstemp(suffix=".frcmod")[1])
        parmchk = Parmchk()
        parmchk.set(mol2=self.probe_mol2, at="gaff2")
        parmchk.run(frcmod=tmp_frcmod)

        expected = open(self.expected_probe_frcmod, "r").read()
        actual = open(tmp_frcmod, "r").read()
        self.assertEqual(expected, actual)

    def test_file_does_not_exist(self):
        parmchk = Parmchk()
        with self.assertRaises(FileNotFoundError):
            parmchk.set(mol2=Path("NOTHING.mol2"), at="gaff2")
            parmchk.run()

    def test_at_is_not_supported(self):
        parmchk = Parmchk()
        with self.assertRaises(ValueError):
            parmchk.set(mol2=self.probe_mol2, at="NOT_SUPPORTED_ID")  # type: ignore
            parmchk.run()
