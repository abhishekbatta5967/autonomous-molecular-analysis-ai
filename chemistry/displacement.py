from rdkit import Chem

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
