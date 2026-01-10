from rdkit import Chem

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