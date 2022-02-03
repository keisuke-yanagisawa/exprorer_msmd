import argparse
import parmed as pmd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert Amber files to Gromacs files")
    parser.add_argument("-iprefix", required=True, help="prefix of input files")
    parser.add_argument("-oprefix", required=True, help="prefix of output files")
    args = parser.parse_args()

    amber = pmd.load_file(f"{args.iprefix}.parm7", xyz=f"{args.iprefix}.rst7")
    amber.save(f"{args.oprefix}.top", overwrite=True)
    amber.save(f"{args.oprefix}.gro", overwrite=True)
