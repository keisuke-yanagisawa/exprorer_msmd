import argparse

VERSION = "1.0.0"

VIS_INFO = """
[ atomtypes ]
VIS      VIS          0.00000  0.00000   V     0.00000e+00   0.00000e+00 ; virtual interaction site

[ nonbond_params ]
; i j func sigma epsilon
VIS   VIS    1  {sigma:5e}   {epsilon:6e}
"""


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="modify topology records")
    parser.add_argument("-v,--version", action="version", version=VERSION)
    parser.add_argument("-i", required=True, help="input topology file")
    parser.add_argument("-o", required=True, help="output topology file")
    parser.add_argument("-ovis", required=True, help="output topology file for virtual repulsion")
    parser.add_argument("-cname", required=True, nargs="+", help="cosolvent name")
    parser.add_argument("-sigma", type=float, default=2,
                        help="sigma for virtual repulsion")
    parser.add_argument("-epsilon", type=float, default=4.184e-6,
                        help="epsilon for virtual repulsion")
    args = parser.parse_args()

    cosolvent_names = args.cname
    POSITION_RESTRAINT = True

    # modify topology file
    with open(args.ovis, "w") as fout:
        fout.write(VIS_INFO.format(sigma=args.sigma, epsilon=args.epsilon))

    with open(args.i) as fin:
        with open(args.o, "w") as fout:
            curr_section = None
            now_mol = None
            natoms = 0
            for line in fin:
                l = line.split(";")[0].strip()
                if l.startswith("["):
                    prev_section = curr_section
                    if prev_section == "atomtypes":
                        fout.write('#include "{}"\n\n'.format(args.ovis))
                    elif prev_section == "atoms":
                        if now_mol in cosolvent_names:  # TODO
                            fout.write("""{aid: 5d}        VIS      1    {cosolvname}    VIS  {aid: 5d} 0.00000000   0.000000

[ virtual_sitesn ]
{aid: 5d}   2  {lst}
""".format(aid=natoms+1, cosolvname=now_mol, lst=" ".join([str(x) for x in range(1, natoms+1)])))
                        natoms = 0
                    curr_section = l[l.find("[")+1:l.find("]")].strip()
                    if curr_section == "moleculetype":
                        now_mol = None
                elif curr_section == "atoms" and l != "":
                    natoms += 1
                elif curr_section == "moleculetype" and now_mol is None and l != "":
                    now_mol = l.split()[0].strip()

                fout.write(line)
