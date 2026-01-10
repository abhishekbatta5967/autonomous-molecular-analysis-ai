import py3Dmol
from rdkit import Chem
from rdkit.Chem import AllChem

def render_molecule_3d(mol):
    mol = Chem.Mol(mol)
    rw = Chem.RWMol(mol)

    # 🔥 FORCE-REMOVE AROMATICITY (atoms)
    for atom in rw.GetAtoms():
        atom.SetIsAromatic(False)

    # 🔥 FORCE-REMOVE AROMATICITY (bonds)
    for bond in rw.GetBonds():
        bond.SetIsAromatic(False)
        bond.SetBondType(Chem.BondType.SINGLE)

    mol = rw.GetMol()

    # Re-sanitize safely (no kekulization issues now)
    Chem.SanitizeMol(
        mol,
        sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^
                    Chem.SanitizeFlags.SANITIZE_KEKULIZE
    )

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol)

    block = Chem.MolToMolBlock(mol)

    view = py3Dmol.view(width=1500, height=900)
    view.addModel(block, "mol")

    view.setStyle({
        "stick": {"radius": 0.15},
        "sphere": {"scale": 0.25, "colorscheme": "Jmol"},
        "bond": {"order": True}
    })


    # TRUE transparent background
    view.setBackgroundColor("transparent")

    # Center & fit molecule
    view.zoomTo()

    view.setClickable(True)
    view.setHoverable(True)
    view.setViewStyle({"style": "outline"})

    # Smooth rotation animation
    #view.spin(True)

    return view
