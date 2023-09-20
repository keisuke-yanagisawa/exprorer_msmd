from typing import List

VERSION = "1.0.0"

VIS_INFO = """
[ atomtypes ]
VIS      VIS          0.00000  0.00000   V     0.00000e+00   0.00000e+00 ; virtual interaction site

[ nonbond_params ]
; i j func sigma epsilon
VIS   VIS    1  {sigma:1.6e}   {epsilon:1.6e}
"""


def addvirtatom2top(top_string: str,
                    probe_names: List[str],
                    sigma: float = 2,
                    epsilon: float = 4.184e-6) -> str:
    """
    Add definition of a virtual atom to a top file
    Pseudo repulsion term (VIS-VIS nonbond LJ parameter) is added to the top file.
    """

    ret = []
    curr_section = None
    now_mol = None
    natoms = 0
    for line in top_string.split("\n"):
        line = line.split(";")[0].strip()
        if line.startswith("["):
            prev_section = curr_section
            if prev_section == "atomtypes":
                ret.append(
                    VIS_INFO.format(sigma=sigma, epsilon=epsilon)
                )
            elif prev_section == "atoms":
                if now_mol in probe_names:  # TODO
                    ret.append(f"""
                    {natoms+1: 5d}        VIS      1    {now_mol}    VIS  {natoms+1: 5d} 0.00000000   0.000000
                    [ virtual_sitesn ]
                    {natoms+1: 5d}   2  {' '.join([str(x) for x in range(1, natoms+1)])}
                    """)
                natoms = 0
            curr_section = line[line.find("[") + 1:line.find("]")].strip()
            if curr_section == "moleculetype":
                now_mol = None
        elif curr_section == "atoms" and line != "":
            natoms += 1
        elif curr_section == "moleculetype" and now_mol is None and line != "":
            now_mol = line.split()[0].strip()

        ret.append(line)
    ret = "\n".join(ret)
    return ret
