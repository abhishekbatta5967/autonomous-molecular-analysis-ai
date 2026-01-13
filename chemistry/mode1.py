from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
from rdkit import Chem

PT = Chem.GetPeriodicTable()

def identity_properties(name, formula, smiles):
    return {
        "Compound Name": name,
        "Molecular Formula": formula
    }

ELEMENT_DATA = {
    1: {"config": "1s¹", "valence": 1},
    6: {"config": "[He] 2s² 2p²", "valence": 4},
    7: {"config": "[He] 2s² 2p³", "valence": 5},
    8: {"config": "[He] 2s² 2p⁴", "valence": 6},
    # extend as needed
}


def atomic_properties(mol):
    atom = mol.GetAtomWithIdx(0)  # Example: properties of the first atom
    Z = atom.GetAtomicNum()
    fc = atom.GetFormalCharge()

    return {
        "Atomic Numbers": Z,
        "Atomic Mass": PT.GetAtomicWeight(Z),
        "Electron Configuration": ELEMENT_DATA.get(Z, {}).get("config", "Unknown"),
        "Valence Electrons": ELEMENT_DATA.get(Z, {}).get("valence", "Unknown"),
        "Valency": atom.GetTotalValence(),
        "Formal Charge": fc,
        "Oxidation State (est.)": fc,
        "Ionic Radius": "N/A" if fc == 0 else "Estimated",
        "Ionization Energy": "Element-based",
        "Electron Affinity": "Element-based",
        "Nuclear Charge": f"+{Z}",
        "Hybridization": atom.GetHybridization().name,
        "Coordination Number": len(atom.GetNeighbors()),
        "Bonding Capacity": atom.GetTotalValence(),
        "Spin State": "Estimated",
        "Polarizability": "Estimated",
        "Atomic Symmetry Contribution": atom.GetDegree(),
        "Quantum Numbers": "Valence shell (theoretical)",
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

