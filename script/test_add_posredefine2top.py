import numpy as np
import pytest

from .add_posredefine2top import embed_posre

__PREFIX = "POSRE_"
__SINGLE_ATOM_IDS = np.array([1])


@pytest.fixture
def basic_topology():
    """Fixture providing basic topology string"""
    return """[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
1               2               yes             0.5     0.8333

[ moleculetype ]
; Name            nrexcl
Protein             3

[ atoms ]
;   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB
     1         HC      1    ACE    HH3      1     0.1123      1.008   ; qtot 0.1123
     2         CT      1    ACE    CH3      2    -0.3662      12.01   ; qtot -0.2539"""


@pytest.fixture
def multi_molecule_topology():
    """Fixture providing topology string containing multiple molecule types"""
    return """[ moleculetype ]
; Name            nrexcl
Water               3

[ atoms ]
     1         OW      1    SOL     OW      1    -0.834      16.00

[ moleculetype ]
; Name            nrexcl
Protein             3

[ atoms ]
     1         HC      1    ACE    HH3      1     0.1123      1.008"""


@pytest.fixture
def simple_topology():
    """Fixture providing simple topology string"""
    return """[ moleculetype ]
; Name            nrexcl
Protein             3

[ atoms ]
     1         HC      1    ACE    HH3      1     0.1123      1.008"""



def test_basic_embed_posre(basic_topology):
    """Test for basic position restraint embedding"""
    strength = [1000]

    expected_output = """[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
1               2               yes             0.5     0.8333

[ moleculetype ]
; Name            nrexcl
Protein             3

[ atoms ]
;   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB
     1         HC      1    ACE    HH3      1     0.1123      1.008   ; qtot 0.1123
     2         CT      1    ACE    CH3      2    -0.3662      12.01   ; qtot -0.2539

; Position restraints
; Position restraints
#ifdef POSRE_1000
[ position_restraints ]
; atom  type      fx      fy      fz

     1     1  4184  4184  4184
     2     1  4184  4184  4184
#endif"""

    result = embed_posre(basic_topology, np.array([1, 2]), __PREFIX, strength)
    assert normalize_string(result) == normalize_string(expected_output)

def test_basic_embed_posre_with_single_atom(basic_topology):
    """Test for when multiple atoms are available but only one atom is specified"""
    strength = [1000]

    expected_output = """[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
1               2               yes             0.5     0.8333

[ moleculetype ]
; Name            nrexcl
Protein             3

[ atoms ]
;   nr       type  resnr residue  atom   cgnr     charge       mass  typeB    chargeB      massB
     1         HC      1    ACE    HH3      1     0.1123      1.008   ; qtot 0.1123
     2         CT      1    ACE    CH3      2    -0.3662      12.01   ; qtot -0.2539

; Position restraints
; Position restraints
#ifdef POSRE_1000
[ position_restraints ]
; atom  type      fx      fy      fz

     1     1  4184  4184  4184
#endif"""
    
    result = embed_posre(basic_topology, __SINGLE_ATOM_IDS, __PREFIX, strength)
    assert normalize_string(result) == normalize_string(expected_output)


def test_multiple_moleculetypes(multi_molecule_topology):
    """Test for multiple molecule types"""
    strength = [500]

    expected_output = """[ moleculetype ]
; Name            nrexcl
Water               3

[ atoms ]
     1         OW      1    SOL     OW      1    -0.834      16.00


; Position restraints
; Position restraints
#ifdef POSRE_500
[ position_restraints ]
; atom  type      fx      fy      fz

     1     1  2092  2092  2092
#endif

[ moleculetype ]
; Name            nrexcl
Protein             3

[ atoms ]
     1         HC      1    ACE    HH3      1     0.1123      1.008"""

    result = embed_posre(multi_molecule_topology, __SINGLE_ATOM_IDS, __PREFIX, strength)
    assert normalize_string(result) == normalize_string(expected_output)


@pytest.mark.parametrize("strength_list,expected_values", [
    ([1000, 500, 100], [4184, 2092, 418]),
    ([1000], [4184]),
    ([100, 50], [418, 209])
])
def test_multiple_strengths(simple_topology, strength_list, expected_values):
    """Test for multiple strength values"""
    posre_blocks = []
    for strength, value in zip(strength_list, expected_values):
        posre_blocks.append(f"""; Position restraints
#ifdef POSRE_{strength}
[ position_restraints ]
; atom  type      fx      fy      fz

     1     1{value:6d}{value:6d}{value:6d}
#endif""")

    expected_output = f"""[ moleculetype ]
; Name            nrexcl
Protein             3

[ atoms ]
     1         HC      1    ACE    HH3      1     0.1123      1.008

; Position restraints
{chr(10).join(posre_blocks)}""" # chr(10) == '\n'

    result = embed_posre(simple_topology, __SINGLE_ATOM_IDS, __PREFIX, strength_list)
    assert normalize_string(result) == normalize_string(expected_output)


def test_empty_strength_list(simple_topology):
    """Test for empty strength list"""
    strength = []
    result = embed_posre(simple_topology, __SINGLE_ATOM_IDS, __PREFIX, strength)
    assert normalize_string(result) == normalize_string(simple_topology)


def normalize_string(s: str) -> str:
    """Normalize string (unify newlines, remove extra whitespace)"""
    return "\n".join(line.rstrip() for line in s.strip().replace('\r\n', '\n').split('\n'))
