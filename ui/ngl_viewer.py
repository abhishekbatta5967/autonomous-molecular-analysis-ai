from rdkit import Chem
from rdkit.Chem import AllChem

def mol_to_pdb_block(mol):
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol)
    return Chem.MolToPDBBlock(mol)

def render_ngl_viewer(pdb_block, show_labels=True, height=700):

    label_block = """
        comp.addRepresentation("label", {
            labelType: "atomname",
            labelSize: 0.55,
            color: "white"
        });
    """ if show_labels else ""

    html = f"""
<div id="ngl_container" style="width:100%; height:{height}px;"></div>

<script src="https://unpkg.com/ngl@latest/dist/ngl.js"></script>

<script>
(function() {{
    const container = document.getElementById("ngl_container");
    if (!container) return;

    container.innerHTML = "";

    const stage = new NGL.Stage("ngl_container", {{
        backgroundColor: "transparent"
    }});

    stage.mouseControls.remove("scroll");
    stage.mouseControls.remove("pan");

    const pdbData = `{pdb_block}`;

    stage.loadFile(
        new Blob([pdbData], {{ type: "text/plain" }}),
        {{ ext: "pdb" }}
    ).then(function(comp) {{

        comp.addRepresentation("ball+stick", {{
            scale: 3.0,
            multipleBond: true,
            bondScale: 0.3,
            atomRadius: 0.9
        }});

        {label_block}

        stage.autoView();
        

        stage.signals.clicked.add(function(picking) {{
            if (!picking) return;

            if (picking.atom) {{
                window.parent.postMessage({{
                    type: "atom",
                    element: picking.atom.element,
                    index: picking.atom.index,
                    name: picking.atom.atomname
                }}, "*");
            }}

            if (picking.bond) {{
                window.parent.postMessage({{
                    type: "bond",
                    order: picking.bond.bondOrder
                }}, "*");
            }}
        }});
    }});
}})();
</script>
"""
    return html
