from typing import List

VERSION = "1.0.0"

VIS_INFO = """
[ atomtypes ]
VIS      VIS          0.00000  0.00000   V     0.00000e+00   0.00000e+00 ; virtual interaction site

[ nonbond_params ]
; i j func sigma epsilon
VIS   VIS    1  {sigma:1.6e}   {epsilon:1.6e}
"""


def addvirtatom2top(top_string: str, probe_names: List[str], sigma: float = 2, epsilon: float = 4.184e-6) -> str:
    """TOPファイルに仮想原子の定義を追加する

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

    output_lines = []
    current_section = None
    current_molecule = None
    atom_count = 0

    for line in top_string.split("\n"):
        # コメントを保持しながら処理
        if ";" in line:
            content, comment = line.split(";", 1)
            content = content.strip()
            comment = ";" + comment
        else:
            content = line.strip()
            comment = ""

        # セクション開始の処理
        if content.startswith("["):
            prev_section = current_section
            if prev_section == "atomtypes":
                output_lines.append(VIS_INFO.format(sigma=sigma, epsilon=epsilon))
            elif prev_section == "atoms":
                if current_molecule in probe_names:
                    output_lines.append(
                        f"""
                    {atom_count+1: 5d}        VIS      1    {current_molecule}    VIS  {atom_count+1: 5d} 0.00000000   0.000000
                    [ virtual_sitesn ]
                    {atom_count+1: 5d}   2  {' '.join([str(x) for x in range(1, atom_count+1)])}
                    """
                    )
                atom_count = 0

            current_section = content[content.find("[") + 1 : content.find("]")].strip()
            if current_section == "moleculetype":
                current_molecule = None

        # セクション内容の処理
        elif current_section == "atoms" and content:
            atom_count += 1
        elif current_section == "moleculetype" and current_molecule is None and content:
            current_molecule = content.split()[0].strip()

        # 元の行を保持（空行も含む）
        if content or comment:
            output_lines.append(line)
        else:
            output_lines.append("")

    return "\n".join(output_lines)
