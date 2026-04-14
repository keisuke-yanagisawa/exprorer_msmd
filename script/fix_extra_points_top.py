"""Fix Extra Point (EP) atoms in GROMACS topology for GROMACS 2026+ compatibility.

ParmEd converts AMBER EP atoms as regular atoms (particle type "A", mass=0).
GROMACS 2026+ requires mass-0 atoms to be virtual sites (type "V") with proper
[ virtual_sites2 ] definitions.

This module post-processes the topology string to:
1. Change EP atomtype particle type from "A" to "V"
2. Remove EP bonds, angles, and dihedrals (replaced by virtual site definition)
3. Add [ virtual_sites2 ] definitions for each EP atom
4. Exclude EP atoms from VIS virtual_sitesn constructions
"""

import re


def fix_extra_points_top(top_string: str) -> str:
    """Convert EP atoms to proper GROMACS virtual sites in a topology string.

    EP atoms are identified by their atom type "EP" in the [ atoms ] section.
    For each EP, the parent atom (bonded to EP) and reference atom (bonded to
    parent, forming the linear arrangement) are determined from the [ bonds ]
    section, and a [ virtual_sites2 ] definition is generated.

    Parameters
    ----------
    top_string : str
        GROMACS topology file content.

    Returns
    -------
    str
        Modified topology with EP atoms as virtual sites.
    """
    if not re.search(r"^\s*EP\s+", top_string, re.MULTILINE):
        return top_string

    lines = top_string.split("\n")
    result = []

    section = None
    ep_indices = set()
    bonds_info = {}  # atom_idx -> [(other_idx, bond_length_nm)]
    ep_bonds = {}    # ep_idx -> (parent_idx, bond_length_nm)
    pending_vs2 = []

    def _compute_vs2():
        """Compute virtual_sites2 definitions from collected bond data."""
        vs2 = []
        for ep_idx in sorted(ep_indices):
            if ep_idx not in ep_bonds:
                continue
            parent_idx, ep_dist = ep_bonds[ep_idx]
            if parent_idx not in bonds_info:
                continue
            for other_idx, other_dist in bonds_info[parent_idx]:
                if other_idx not in ep_indices:
                    a = -(ep_dist / other_dist)
                    vs2.append(
                        f"{ep_idx:5d} {parent_idx:4d} {other_idx:4d}    1   {a:.6f}"
                    )
                    break
        return vs2

    def _emit_pending():
        """Insert pending virtual_sites2 lines into result."""
        nonlocal pending_vs2
        if pending_vs2:
            result.append("")
            result.append("[ virtual_sites2 ]")
            result.append("; Site  from    funct  a")
            for vs_line in pending_vs2:
                result.append(vs_line)
            result.append("")
            pending_vs2 = []

    def _reset():
        nonlocal ep_indices, bonds_info, ep_bonds, pending_vs2
        ep_indices = set()
        bonds_info = {}
        ep_bonds = {}
        pending_vs2 = []

    for line in lines:
        stripped = line.strip()

        # --- Section headers ---
        if stripped.startswith("["):
            new_section = stripped[stripped.find("[") + 1 : stripped.find("]")].strip()

            if new_section in ("moleculetype", "system") and ep_indices:
                pending_vs2 = _compute_vs2()
                _emit_pending()
                _reset()

            section = new_section
            result.append(line)
            continue

        # Empty / comment lines
        if not stripped or stripped.startswith(";"):
            result.append(line)
            continue

        # --- [ atomtypes ]: fix EP particle type A → V ---
        if section == "atomtypes":
            if re.match(r"EP\s+", stripped):
                line = re.sub(r"(\s+)A(\s+)", r"\1V\2", line, count=1)
            result.append(line)
            continue

        # --- [ atoms ]: identify EP atoms ---
        if section == "atoms":
            parts = stripped.split()
            if len(parts) >= 2:
                try:
                    idx = int(parts[0])
                    atype = parts[1]
                    if atype == "EP":
                        ep_indices.add(idx)
                except (ValueError, IndexError):
                    pass
            result.append(line)
            continue

        # --- [ virtual_sitesn ]: exclude EP atoms from VIS construction ---
        if section == "virtual_sitesn" and ep_indices:
            parts = stripped.split()
            if len(parts) >= 3:
                try:
                    # Format: site_idx funct atom1 atom2 ...
                    site_idx = int(parts[0])
                    funct = int(parts[1])
                    atom_refs = [int(p) for p in parts[2:]]
                    filtered = [a for a in atom_refs if a not in ep_indices]
                    if len(filtered) != len(atom_refs):
                        line = f"                       {site_idx}   {funct}  {' '.join(str(a) for a in filtered)}"
                except (ValueError, IndexError):
                    pass
            result.append(line)
            continue

        # --- [ bonds ]: filter EP bonds, collect bond data ---
        if section == "bonds" and ep_indices:
            parts = stripped.split()
            if len(parts) >= 2:
                try:
                    ai, aj = int(parts[0]), int(parts[1])
                    dist = float(parts[3]) if len(parts) >= 4 else None

                    if ai in ep_indices or aj in ep_indices:
                        ep_idx = ai if ai in ep_indices else aj
                        par_idx = aj if ai in ep_indices else ai
                        if dist is not None:
                            ep_bonds[ep_idx] = (par_idx, dist)
                        continue  # remove EP bond

                    if dist is not None:
                        bonds_info.setdefault(ai, []).append((aj, dist))
                        bonds_info.setdefault(aj, []).append((ai, dist))
                except (ValueError, IndexError):
                    pass
            result.append(line)
            continue

        # --- [ angles ]: filter EP angles ---
        if section == "angles" and ep_indices:
            parts = stripped.split()
            if len(parts) >= 3:
                try:
                    indices = [int(parts[j]) for j in range(3)]
                    if any(i in ep_indices for i in indices):
                        continue
                except (ValueError, IndexError):
                    pass
            result.append(line)
            continue

        # --- [ dihedrals ]: filter EP dihedrals ---
        if section == "dihedrals" and ep_indices:
            parts = stripped.split()
            if len(parts) >= 4:
                try:
                    indices = [int(parts[j]) for j in range(4)]
                    if any(i in ep_indices for i in indices):
                        continue
                except (ValueError, IndexError):
                    pass
            result.append(line)
            continue

        # --- [ pairs ]: filter EP pairs ---
        if section == "pairs" and ep_indices:
            parts = stripped.split()
            if len(parts) >= 2:
                try:
                    ai, aj = int(parts[0]), int(parts[1])
                    if ai in ep_indices or aj in ep_indices:
                        continue
                except (ValueError, IndexError):
                    pass
            result.append(line)
            continue

        result.append(line)

    # Final flush
    if ep_indices:
        pending_vs2 = _compute_vs2()
        _emit_pending()

    return "\n".join(result)
