import pytest
from script.addvirtatom2top import addvirtatom2top, VIS_INFO


@pytest.fixture
def basic_top():
    """基本的なTOPファイルの内容を準備"""
    return """
[ defaults ]
1   1   no

[ atomtypes ]
; name mass charge ptype sigma epsilon
C1    12.01   0.0000  A   3.50000e-01  2.76144e-01

[ moleculetype ]
; name    nrexcl
MOL      3

[ atoms ]
; nr type resnr residu atom cgnr charge mass
1   C1   1     MOL    C1   1    0.0000  12.0100
2   C1   1     MOL    C2   2    0.0000  12.0100

[ bonds ]
1  2  1  0.1400  392459.2
"""


@pytest.fixture
def expected_vis_info():
    """期待されるVIS_INFO文字列を準備"""
    return VIS_INFO.format(sigma=2, epsilon=4.184e-6)


def test_basic_conversion(basic_top, expected_vis_info):
    """基本的な変換のテスト"""
    result = addvirtatom2top(basic_top, ["MOL"])
    
    # VIS_INFOが追加されていることを確認
    assert expected_vis_info in result
    
    # 仮想原子が追加されていることを確認
    assert "3        VIS      1    MOL    VIS" in result
    
    # 仮想サイトの定義が追加されていることを確認
    assert "[ virtual_sitesn ]" in result
    assert "3   2  1 2" in result


def test_multiple_molecules():
    """複数の分子タイプを含むTOPファイルのテスト"""
    top_with_multiple = """
[ defaults ]
1   1   no

[ atomtypes ]
C1    12.01   0.0000  A   3.50000e-01  2.76144e-01

[ moleculetype ]
MOL1     3

[ atoms ]
1   C1   1     MOL1   C1   1    0.0000  12.0100

[ moleculetype ]
MOL2     3

[ atoms ]
1   C1   1     MOL2   C1   1    0.0000  12.0100
"""
    result = addvirtatom2top(top_with_multiple, ["MOL1"])
    
    # MOL1にのみ仮想原子が追加されていることを確認
    assert "2        VIS      1    MOL1    VIS" in result
    assert "VIS      1    MOL2    VIS" not in result


def test_custom_parameters(basic_top):
    """カスタムパラメータでの変換テスト"""
    custom_sigma = 3.0
    custom_epsilon = 1.0e-5
    
    result = addvirtatom2top(basic_top, ["MOL"], sigma=custom_sigma, epsilon=custom_epsilon)
    
    expected_vis_info = VIS_INFO.format(sigma=custom_sigma, epsilon=custom_epsilon)
    assert expected_vis_info in result


def test_no_probe_names(basic_top, expected_vis_info):
    """プローブ名が空リストの場合のテスト"""
    result = addvirtatom2top(basic_top, [])
    
    # VIS_INFOは追加されるが、仮想原子は追加されないことを確認
    assert expected_vis_info in result
    assert "VIS      1    MOL    VIS" not in result
    assert "[ virtual_sitesn ]" not in result


def test_probe_name_not_in_top(basic_top, expected_vis_info):
    """存在しない分子名を指定した場合のテスト"""
    result = addvirtatom2top(basic_top, ["NONEXISTENT"])
    
    # VIS_INFOは追加されるが、仮想原子は追加されないことを確認
    assert expected_vis_info in result
    assert "VIS      1    NONEXISTENT    VIS" not in result


def test_empty_input():
    """空の入力文字列のテスト"""
    result = addvirtatom2top("", ["MOL"])
    assert result == ""


def test_malformed_input():
    """不正な形式の入力のテスト"""
    malformed_top = """
[ invalid_section
missing_bracket
"""
    result = addvirtatom2top(malformed_top, ["MOL"])
    
    # エラーを発生させずに処理できることを確認
    assert "[ invalid_section" in result
    assert "missing_bracket" in result


def test_preserve_comments():
    """コメントが保持されることを確認するテスト"""
    top_with_comments = """
; This is a comment
[ atomtypes ] ; Section comment
; name mass charge ptype sigma epsilon
C1    12.01   0.0000  A   3.50000e-01  2.76144e-01 ; Atom comment
"""
    result = addvirtatom2top(top_with_comments, ["MOL"])
    
    # コメントが保持されていることを確認
    assert "; This is a comment" in result
    assert "[ atomtypes ] ; Section comment" in result
    assert "; name mass charge ptype sigma epsilon" in result
    assert "; Atom comment" in result