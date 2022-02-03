import numpy as np
from scipy import constants
from . import logger
aLogger = logger.Logger()
ATOMIC_NUMBER = {"H": 1,  "C": 6,  "N": 7,  "O": 8,  "F": 9, "Na": 11,
                 "Si": 14,  "P": 15,  "S": 16, "Cl": 17,  "K": 19, "Ca": 20,
                 "Br": 35,  "I": 53}
ATOMIC_STR_LEN = {1: 1,  6: 1,  7: 1,  8: 1,  9: 1, 11: 2, 14: 2, 15: 1, 16: 1,
                  17: 2, 19: 1, 20: 2, 35: 2, 53: 2}
ATOMIC_WEIGHT = {1: 1.008,  6: 12.01,  7: 14.01,  8: 16.00,  9: 19.00,
                 11: 22.99,  14: 28.09, 15: 30.97, 16: 32.07, 17: 35.45,
                 19: 39.10,  20: 40.08, 35: 79.90, 53: 126.9,
                 -1: 0.000}  # for dummy atoms


class Gro_atom():
    def __init__(self, string=""):
        self.resi = -1
        self.resn = ""
        self.atomtype = ""
        self.atom_id = -1
        self.point = np.array([])
        self.velocity = np.array([])
        self.comment = ""
        self.atomic_mass = 0

        if string != "":
            self.parse(string)

    def parse(self, string):
        if len(string.split(";")) > 1:  # there is a comment
            info = string.split(";")
            string = info[0].rstrip()
            self.comment = ";".join(info[1:]).strip()

        self.resi = int(string[0:5])
        self.resn = string[5:10].strip()
        self.atomtype = string[10:15].strip()
        self.atom_id = int(string[15:20])
        tmp = np.array([float(s) for s in string[20:].split()])
        if(len(tmp) == 3):
            self.point = tmp
            self.velocity = np.array([0, 0, 0])
        elif(len(tmp) == 6):
            self.point = tmp[:3]
            self.velocity = tmp[3:]
        else:
            aLogger.warn("the dimension of atom coordinates/velocities are wrong: {}".format(tmp))

        temp = [i for atype, i in ATOMIC_NUMBER.items() if self.atomtype.startswith(atype)]
        if(len(temp) == 1):
            self.atomic_num = temp[0]
        elif(len(temp) >= 2):
            self.atomic_num = temp[np.argmax([ATOMIC_STR_LEN[i] for i in temp])]
        else:  # len == 0
            aLogger.warn("atomtype {} is not matched to any atom names".format(self.atomtype))
            aLogger.warn("assume it is a kind of pseudo atom")
            self.atomic_num = -1

        self.atomic_mass = ATOMIC_WEIGHT[self.atomic_num]

    def __repr__(self):
        ret_str = "{:>5}{:>5}{:>5}{:>5}{:8.3f}{:8.3f}{:8.3f}" \
            .format(self.resi, self.resn, self.atomtype, self.atom_id,
                    *self.point)
        if(self.comment != ""):
            ret_str += " ; {}".format(self.comment)
        return ret_str


class Gro():
    def __init__(self, path=""):
        self.description = ""
        self.natoms = 0
        self.box_size = [0, 0, 0]
        self.atoms = []
        if path != "":
            self.parse(path)

    def parse(self, path):
        with open(path) as fin:
            lines = [l.rstrip() for l in fin.readlines()]
        self.description = lines[0]
        self.natoms = int(lines[1])
        self.box_size = [float(s) for s in lines[-1].split()]
        for line in lines[2:-1]:
            self.atoms.append(Gro_atom(line))

    def __update_atomid(self):
        for i in range(len(self.atoms)):
            self.atoms[i].atom_id = i+1

    def get_atoms(self, resi=-1, resn="", atomtype="", atom_id=-1, atomic_num=-1):
        ret_atoms = []
        for atom in self.atoms:
            if(resi != -1 and atom.resi != resi):
                continue
            elif(resn != "" and atom.resn != resn):
                continue
            elif(atomtype != "" and atom.atomtype != atomtype):
                continue
            elif(atom_id != -1 and atom.atom_id != atom_id):
                continue
            elif(atomic_num != -1 and atom.atomic_num != atomic_num):
                continue
            else:
                ret_atoms.append(atom)
        return ret_atoms

    def add_atom(self, atom):
        if not isinstance(atom, Gro_atom):
            aLogger.error("the input is NON-GRO_ATOM")
            return
        atom.atom_id = max([a.atom_id for a in self.atoms])+1
        self.atoms.append(atom)
        self.natoms = len(self.atoms)
        self.__sort_atoms()
        self.__update_atomid()

    def __sort_atoms(self):
        self.atoms = sorted(self.atoms, key=lambda a: a.resi)

    def __repr__(self):
        ret_str = "{}\n".format(self.description)
        ret_str += "{:>6}\n".format(self.natoms)
        for a in self.atoms:
            ret_str += "{}\n".format(a)
        ret_str += "{: 10.5f}  {: 10.5f}  {: 10.5f}\n".format(*self.box_size)
        return ret_str

    def molar(self, resn):
        focused_atoms = self.get_atoms(resn=resn)
        resis = {a.resi for a in focused_atoms}
        volume = self.box_size[0] * self.box_size[1] * self.box_size[2]  # in nanometer
        return (len(resis) / constants.N_A) / (volume * 1E-24)  # nm^3 -> cm^3
