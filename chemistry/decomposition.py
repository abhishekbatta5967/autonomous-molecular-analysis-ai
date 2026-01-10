from rdkit import Chem

def ring_based_decomposition(mol):
    """
    Ring-based structural decomposition with aromatic safety.
    """

    # Full sanitize original molecule
    Chem.SanitizeMol(mol)

    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    if not atom_rings:
        return [], mol

    # Collect ring atom indices
    ring_atoms = set()
    for ring in atom_rings:
        ring_atoms.update(ring)

    # ---- Ring submolecule ----
    ring_mol = Chem.PathToSubmol(mol, list(ring_atoms))
    Chem.SanitizeMol(
        ring_mol,
        sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^
                    Chem.SanitizeFlags.SANITIZE_KEKULIZE
    )

    # ---- Chain submolecule ----
    chain_atoms = [
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetIdx() not in ring_atoms
    ]

    chain_mol = None
    if chain_atoms:
        chain_mol = Chem.PathToSubmol(mol, chain_atoms)
        Chem.SanitizeMol(
            chain_mol,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^
                        Chem.SanitizeFlags.SANITIZE_KEKULIZE
        )

    return [ring_mol], chain_mol
