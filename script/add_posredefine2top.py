import os
from typing import List

import jinja2
import numpy.typing as npt
from scipy import constants

VERSION = "2.0.0"


def _position_restraint(atom_id_list: npt.ArrayLike, prefix: str, weight) -> str:
    """Generate a position restraint #ifdef block using Jinja2 template."""
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(os.path.dirname(__file__)))
    template = env.get_template("./template/position_restraints")
    return template.render(
        {
            "define_name": f"{prefix}{weight}",
            "weight": weight,
            "weight_in_calorie": weight * constants.calorie,
            "atom_id_list": atom_id_list,
        }
    )


def _insert_posre_block(lines: List[str], atom_id_list: npt.ArrayLike, prefix: str, strength: list[int]):
    """Append position restraint blocks to the output lines."""
    lines.append("")
    lines.append("; Position restraints")
    for s in strength:
        lines.append(_position_restraint(atom_id_list, prefix, s))


def embed_posre(top_string: str, atom_id_list: npt.ArrayLike, prefix: str, strength: list[int]) -> str:
    """Embed position restraint records into a given topology string.

    #ifdef ブロックは GROMACS 固有のプリプロセッサ指令のためテキスト挿入で対応。
    Position restraints are only inserted into the first molecule type.
    """
    ret = []
    curr_section = None
    mol_count = 0
    in_first_molecule = False

    for line in top_string.split("\n"):
        content = line.split(";", 1)[0].strip() if ";" in line else line.strip()

        if content.startswith("[") and "]" in content:
            # セクション遷移時に位置拘束を挿入
            if curr_section == "atoms" and in_first_molecule and strength:
                _insert_posre_block(ret, atom_id_list, prefix, strength)
                ret.append("")

            curr_section = content[content.find("[") + 1 : content.find("]")].strip()
            if curr_section == "moleculetype":
                mol_count += 1
                in_first_molecule = mol_count == 1

        ret.append(line)

    # ファイルの最後が atoms セクションの場合の処理
    if curr_section == "atoms" and in_first_molecule and strength:
        _insert_posre_block(ret, atom_id_list, prefix, strength)

    return "\n".join(ret)
