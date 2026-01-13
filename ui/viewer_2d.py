from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdDepictor
import streamlit as st
from io import BytesIO

rdDepictor.SetPreferCoordGen(True)

def render_2d_structure(mol):
    img = Draw.MolToImage(mol, size=(800, 550))
    buf = BytesIO()
    img.save(buf, format="PNG")
    st.image(buf.getvalue())