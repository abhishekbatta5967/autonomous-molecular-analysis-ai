from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors

def identity_properties(name, formula, smiles):
    return {
        "Compound Name": name,
        "Molecular Formula": formula,
        "SMILES Representation": smiles
    }



def atomic_properties(mol):
    return {
        "Molecular Weight": f"{round(Descriptors.MolWt(mol), 2)} g/mol",
        "Num Atoms": mol.GetNumAtoms(),
        "Num Bonds": mol.GetNumBonds(),
        "Num Electrons": sum(atom.GetAtomicNum() for atom in mol.GetAtoms()),
        "Formal Charge": f"{Descriptors.NumValenceElectrons(mol)} e"
    }

def structural_properties(mol):
    try:
        return {
            "Rings": Lipinski.RingCount(mol),
            "Aromatic Rings": Lipinski.NumAromaticRings(mol),
            "H-Bond Donors": Lipinski.NumHDonors(mol),
            "H-Bond Acceptors": Lipinski.NumHAcceptors(mol),
            "Rotatable Bonds": Lipinski.NumRotatableBonds(mol),
            "TPSA": round(Descriptors.TPSA(mol), 2),
            "logP": round(Descriptors.MolLogP(mol), 2)
        }
    except Exception:
        return {
            "Rings": "N/A",
            "Aromatic Rings": "N/A",
            "H-Bond Donors": Lipinski.NumHDonors(mol),
            "H-Bond Acceptors": Lipinski.NumHAcceptors(mol),
            "Rotatable Bonds": "N/A",
            "TPSA": "N/A",
            "logP": "N/A"
        }

def predicted_properties(mol):
    return {
        "Predicted Solubility": "Coming Soon",
        "Predicted Toxicity": "Coming Soon",
        "Drug-likeness": "Coming Soon"
    }

