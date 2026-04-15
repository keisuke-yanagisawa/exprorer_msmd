import os
import tempfile
from typing import Dict, List

import parmed as pmd

VERSION = "1.0.0"

VIS_INFO = """
[ atomtypes ]
VIS      VIS          0.00000  0.00000   V     0.00000e+00   0.00000e+00 ; virtual interaction site

[ nonbond_params ]
; i j func sigma epsilon
VIS   VIS    1  {sigma:1.6e}   {epsilon:1.6e}
"""


def _get_molecule_atom_counts(top_string: str) -> Dict[str, int]:
    """Parse topology with ParmEd to extract molecule name → atom count mapping.

    Returns an empty dict if ParmEd cannot parse the topology (e.g. malformed input).
    """
    with tempfile.NamedTemporaryFile(suffix=".top", mode="w", delete=False) as f:
        f.write(top_string)
        tmppath = f.name
    try:
        struct = pmd.load_file(tmppath)
        result = {}
        for parm, _ in struct.split():
            if parm.residues:
                mol_name = parm.residues[0].name
                result[mol_name] = len(parm.atoms)
        return result
    except Exception:
        return {}
    finally:
        os.unlink(tmppath)


def _generate_virtual_atom_entry(atom_count: int, molecule_name: str) -> str:
    """Generate virtual atom and virtual_sitesn section text."""
    vis_id = atom_count + 1
    atom_ids = " ".join(str(x) for x in range(1, atom_count + 1))
    return (
        f"\n"
        f"                    {vis_id: 5d}        VIS      1    {molecule_name}    VIS  {vis_id: 5d} 0.00000000   0.000000\n"
        f"                    [ virtual_sitesn ]\n"
        f"                    {vis_id: 5d}   2  {atom_ids}\n"
        f"                    "
    )


def addvirtatom2top(top_string: str, probe_names: List[str], sigma: float = 2, epsilon: float = 4.184e-6) -> str:
    """TOPファイルに仮想原子の定義を追加する

    ParmEd でトポロジー構造をパースし、各分子の原子数を取得する。
    GROMACS 固有の virtual_sitesn / nonbond_params はテキスト挿入で対応。

    Args:
        top_string: 入力TOPファイルの内容
        probe_names: 仮想原子を追加する分子名のリスト
        sigma: VIS-VIS相互作用のシグマパラメータ
        epsilon: VIS-VIS相互作用のイプシロンパラメータ

    Returns:
        str: 仮想原子が追加されたTOPファイルの内容
    """
    if not top_string:
        return ""

    # ParmEd で分子名→原子数のマッピングを取得
    mol_atom_counts = _get_molecule_atom_counts(top_string)

    output_lines = []
    current_section = None
    current_molecule = None
    # ParmEd パース失敗時のフォールバック用
    manual_atom_count = 0

    for line in top_string.split("\n"):
        # コメントを分離してセクション検出
        content = line.split(";", 1)[0].strip() if ";" in line else line.strip()

        if content.startswith("[") and "]" in content:
            prev_section = current_section

            # atomtypes セクション終了後に VIS 定義を挿入
            if prev_section == "atomtypes":
                output_lines.append(VIS_INFO.format(sigma=sigma, epsilon=epsilon))

            # atoms セクション終了後にプローブ分子に仮想原子を挿入
            elif prev_section == "atoms" and current_molecule in probe_names:
                atom_count = mol_atom_counts.get(current_molecule, manual_atom_count)
                if atom_count > 0:
                    output_lines.append(_generate_virtual_atom_entry(atom_count, current_molecule))
                manual_atom_count = 0

            current_section = content[content.find("[") + 1 : content.find("]")].strip()
            if current_section == "moleculetype":
                current_molecule = None

        elif current_section == "atoms" and content:
            manual_atom_count += 1
        elif current_section == "moleculetype" and current_molecule is None and content:
            current_molecule = content.split()[0].strip()

        # 元の行を保持（空行も含む）
        if content or ";" in line:
            output_lines.append(line)
        else:
            output_lines.append("")

    return "\n".join(output_lines)
