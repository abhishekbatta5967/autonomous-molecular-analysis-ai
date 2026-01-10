import pubchempy as pcp
from rdkit import Chem
import streamlit as st


def load_molecule(user_input: str):
    query = user_input.strip()

    compounds = []

    # Try name search
    try:
        compounds = pcp.get_compounds(query, "name")
    except:
        pass

    # If name fails, try formula search
    if not compounds:
        try:
            compounds = pcp.get_compounds(query, "formula")
        except:
            pass

    if not compounds:
        return None, None, None, None, "Compound not found in PubChem"

    comp = compounds[0]

    smiles = comp.connectivity_smiles or comp.isomeric_smiles
    if not smiles:
        return None, None, None, None, "SMILES not available"

    formula = comp.molecular_formula or "Unknown"
    name = comp.iupac_name or comp.synonyms[0] if comp.synonyms else "Unknown"

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None, None, None, "Invalid molecule structure"

    Chem.SanitizeMol(mol)
    return mol, name, formula, smiles, None

@st.cache_data(show_spinner=False)
def load_molecule_cached(user_input):
    return load_molecule(user_input)
