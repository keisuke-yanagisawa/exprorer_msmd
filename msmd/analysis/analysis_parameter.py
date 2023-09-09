from typing import Final, Literal
from ..variable import VariableInterface


class CpptrajSelector(VariableInterface):
    @staticmethod
    def _validation(selector: str) -> None:
        """selectorの妥当性をチェックする"""
        # TODO どうやってチェックするのか？
        pass

    def __init__(self, selector: str):
        self.__selector: Final[str] = selector
        self._validation(self.__selector)

    def get(self) -> str:
        return self.__selector

    def __str__(self) -> str:
        return f"CpptrajSelector('{self.__selector}')"


class PMAPNormalization(VariableInterface):
    @staticmethod
    def _validation(normalization: Literal["total", "snapshot", "gfe"]) -> None:
        if normalization not in ["total", "snapshot", "gfe"]:
            raise ValueError(f"The normalization: {normalization} is not supported")

    def __init__(self, normalization: Literal["total", "snapshot", "gfe"]):
        self.__normalization: Final[Literal["total", "snapshot", "gfe"]] = normalization
        self._validation(self.__normalization)

    def get(self) -> Literal["total", "snapshot", "gfe"]:
        return self.__normalization

    def __str__(self):
        return f"PMAPNormalization('{self.__normalization}')"


# map: # settings to create a spatial probability distribution map (PMAP)
#   type: pmap
#   # valid_dist: 5            # PMAP valid distance from protein atoms
#   snapshot: 2001-4001:1 # steps / "NAME" trajectory will be used with "NAME|"
#   maps: # list of maps to be created
#     - suffix: nVH # suffix of the map file name
#       selector: (!@VIS)&(!@H*) # atoms to be used to create the PMAP
#     - suffix: nV
#       selector: (!@VIS)
#   map_size: 80 # 80 A * 80 A * 80 A
#   normalization: total # total or snapshot

# probe_profile: # settings to create a residue interaction profile (inverse MSMD)
#   resenv: # Extract residue environments around probe molecules
#     map: nVH # map suffix
#     threshold: 0.001 # probability threshold of "preferred" residue environment
#     # env_dist: 4    # residue environment distance from a probe molecule
#     # align: [" C1 ", " C2 ", " C3 ", " O1 "] # all heavy atoms are used if "align" is not defined

#   profile: # list of residue interaction profiles to be created
#     types:
#       - name: anion # name of the profile
#         atoms: # atoms to be used to create the profile
#           - ["ASP", " CB "]
#           - ["GLU", " CB "]
#       - name: cation
#         atoms:
#           - ["LYS", " CB "]
#           - ["HIE", " CB "]
#           - ["ARG", " CB "]
#       - name: aromatic
#         atoms:
#           - ["TYR", " CB "]
#           - ["TRP", " CB "]
#           - ["PHE", " CB "]
#       - name: hydrophilic
#         atoms:
#           - ["ASN", " CB "]
#           - ["GLN", " CB "]
#           - ["SER", " CB "]
#           - ["THR", " CB "]
#           - ["CYS", " CB "]
#       - name: hydrophobic
#         atoms:
#           - ["ILE", " CB "]
#           - ["LEU", " CB "]
#           - ["VAL", " CB "]
#           - ["ALA", " CB "]
#           - ["PRO", " CB "]
#           - ["MET", " CB "]
#       - name: neutral
#         atoms:
#           - ["GLY", " CB "]
#       - name: gly
#         atoms:
#           - ["GLY", " CB "]
#       - name: met
#         atoms:
#           - ["MET", " CB "]
