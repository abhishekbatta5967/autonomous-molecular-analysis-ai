from rdkit import Chem
import random

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
