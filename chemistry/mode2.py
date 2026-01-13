from rdkit import Chem
import random

def combine_molecules(mol1, mol2):
    """
    Prototype-level molecule combination:
    - No chemical reaction
    - Graph-level merge
    """
    try:
        combined = Chem.CombineMols(mol1, mol2)
        Chem.SanitizeMol(combined)
        return combined, None
    except Exception as e:
        return None, str(e)
    

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


def single_displacement(mol_a, mol_bc):
    """
    Relaxed structural single displacement:
    Replaces a weakly connected atom in BC
    """

    mol_bc = Chem.AddHs(Chem.Mol(mol_bc))
    rw = Chem.RWMol(mol_bc)

    # Collect heavy atoms only
    heavy_atoms = [a for a in rw.GetAtoms() if a.GetSymbol() != "H"]

    if len(heavy_atoms) < 2:
        return None, None, "Molecule too simple for displacement"

    # 1️⃣ Try terminal heavy atoms
    candidates = [a for a in heavy_atoms if a.GetDegree() == 1]

    # 2️⃣ Fallback: lowest-degree heavy atoms
    if not candidates:
        min_deg = min(a.GetDegree() for a in heavy_atoms)
        candidates = [a for a in heavy_atoms if a.GetDegree() == min_deg]

    # 3️⃣ Final fallback: any heavy atom except first
    if not candidates:
        candidates = heavy_atoms[1:]

    displaced_atom = candidates[-1]
    displaced_symbol = displaced_atom.GetSymbol()
    displaced_idx = displaced_atom.GetIdx()

    # Remove atom
    rw.RemoveAtom(displaced_idx)
    product_ac = rw.GetMol()

    product_ac = Chem.RemoveHs(product_ac)

    try:
        Chem.SanitizeMol(
            product_ac,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^
                        Chem.SanitizeFlags.SANITIZE_KEKULIZE
        )
    except:
        return None, None, "Displacement produced unstable structure"

    # Displaced atom as separate product
    product_b = Chem.MolFromSmiles(displaced_symbol)

    return product_ac, product_b, None



def double_displacement(mol_ab, mol_cd):
    """
    Structural double displacement:
    AB + CD → AD + CB
    """

    frags_ab = Chem.GetMolFrags(mol_ab, asMols=True)
    frags_cd = Chem.GetMolFrags(mol_cd, asMols=True)

    if len(frags_ab) < 1 or len(frags_cd) < 1:
        return None, None, "Invalid molecules"

    # Take first fragment of each
    a = frags_ab[0]
    c = frags_cd[0]

    product_ad = Chem.CombineMols(a, c)
    product_cb = Chem.CombineMols(frags_cd[-1], frags_ab[-1])

    try:
        Chem.SanitizeMol(product_ad, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^
                                        Chem.SanitizeFlags.SANITIZE_KEKULIZE)
        Chem.SanitizeMol(product_cb, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^
                                        Chem.SanitizeFlags.SANITIZE_KEKULIZE)
    except:
        return None, None, "Invalid displacement products"

    return product_ad, product_cb, None


def simple_rearrangement(mol):
    """
    Structural rearrangement by bond shuffling.
    Keeps atom count and formula same.
    """

    mol = Chem.Mol(mol)
    rw = Chem.RWMol(mol)

    bonds = list(rw.GetBonds())
    atoms = list(rw.GetAtoms())

    # Need at least 3 atoms to rearrange
    if len(atoms) < 3 or len(bonds) < 2:
        return None, "Molecule too small to rearrange"

    # Pick a random bond to break
    bond = random.choice(bonds)
    a1 = bond.GetBeginAtomIdx()
    a2 = bond.GetEndAtomIdx()
    rw.RemoveBond(a1, a2)

    # Pick two different atoms to connect
    possible_atoms = [a.GetIdx() for a in atoms if a.GetIdx() not in (a1, a2)]
    if len(possible_atoms) < 2:
        return None, "No valid rearrangement found"

    new_a1, new_a2 = random.sample(possible_atoms, 2)
    rw.AddBond(new_a1, new_a2, Chem.BondType.SINGLE)

    new_mol = rw.GetMol()

    try:
        Chem.SanitizeMol(
            new_mol,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^
                        Chem.SanitizeFlags.SANITIZE_KEKULIZE
        )
    except:
        return None, "Rearrangement produced invalid structure"

    return new_mol, None