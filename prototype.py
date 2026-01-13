import streamlit as st
from rdkit import Chem
import json
from chemistry.molecule import load_molecule_cached as load_molecule
from chemistry.mode1 import atomic_properties, structural_properties, identity_properties, predicted_properties
from chemistry.mode2 import combine_molecules, ring_based_decomposition, single_displacement, double_displacement, simple_rearrangement
from chemistry.mode3 import suggest_modifications
from ui.ngl_viewer import mol_to_pdb_block, render_ngl_viewer
from ui.viewer_2d import render_2d_structure
from llm.description import ai_description, mol_to_text, mol_to_smiles



st.set_page_config(layout="wide")
st.title("🧬 Autonomous Molecular Analysis System")

mode_tabs = st.tabs([
    "Mode 1: Analyze Compound",
    "Mode 2: Compounds Reactions",
    "Mode 3: Optimize Properties"
])

with mode_tabs[0]:
    st.subheader("Compound Analysis")
    st.info("🔍 Mode 1 analyzes molecular structure and properties from a single compound.")


    user_input = st.text_input("Enter compound name or molecular formula:")

    if user_input:
        with st.spinner("⚗️ Resolving molecular identity and computing properties..."):
            mol, name, formula, smiles, error = load_molecule(user_input)

            if error:
                st.error(error)
            else:
                st.subheader("Molecular Properties")
                tabs = st.tabs(["Atomic", "Structural", "Predicted"])

                with tabs[0]:
                    for k, v in identity_properties(name, formula, smiles).items():
                        st.write(f"**{k}:** {v}")
                    for k, v in atomic_properties(mol).items():
                        st.write(f"**{k}:** {v}")

                with tabs[1]:
                    for k, v in structural_properties(mol).items():
                        st.write(f"**{k}:** {v}")

                with tabs[2]:
                    st.info("Predicted properties will be available in the full model")
                
                col1, col2 = st.columns(2, gap="large")

                with col1:
                    st.subheader("2D Molecular Structure")
                    render_2d_structure(mol)

                    st.download_button("Download 2D Img", data=Chem.MolToSmiles(mol), file_name= f"{user_input}.smiles")

                with col2:
                    col1, col2 = st.columns(2, gap="large")
                    with col1:
                        st.subheader("3D Molecular Structure")
                    with col2:
                        show_labels = st.toggle("Show Atom Labels", value=True)
                    pdb_block = mol_to_pdb_block(mol)
                    html = render_ngl_viewer(pdb_block, show_labels=show_labels)
                    st.components.v1.html(html, height=550)
                    st.download_button("Download 3D PDB", data=pdb_block, file_name= f"{user_input}.pdb")

            st.subheader("🧠 Compound Description (AI Generated)")

            with st.spinner("🧪 Interpreting chemical structure..."):
                input_list = [name, formula, smiles]
                description = ai_description("Compound Analysis", input_list, None)

            st.write(description)

with mode_tabs[1]:
    st.subheader("🧪 Compound Reactions")
    st.info("🔁 Mode 2 simulates structural reaction categories for educational purposes.")


    reaction_tabs = st.tabs([
        "🔗 Combination",
        "🧯 Decomposition",
        "🔁 Single Displacement",
        "🔄 Double Displacement",
        "🔃 Rearrangement"
    ])

    with reaction_tabs[0]:
        st.markdown("### 🔗 Combination Reaction (A + B → AB)")

        col1, col2 = st.columns(2, gap="large")

        with col1:
            compound_a = st.text_input(
                "Compound A (name or formula)",
                key="compound_a"
            )

        with col2:
            compound_b = st.text_input(
                "Compound B (name or formula)",
                key="compound_b"
            )

        combine_btn = st.button("🔗 Combine Compounds")

        if combine_btn:
            if not compound_a or not compound_b:
                st.error("please provide both compounds.")
            else:
                with st.spinner("Resolving compunds and combining structures..."):
                    mol_a, name_a, formula_a, smiles_a, err_a = load_molecule(compound_a)
                    mol_b, name_b, formula_b, smiles_b, err_b = load_molecule(compound_b)

                    if err_a:
                        st.error(f"Compound A error: {err_a}")
                    elif err_b:
                        st.error(f"Compound B error: {err_b}")
                    else:
                        combined_mol, err_c = combine_molecules(mol_a, mol_b)

                        if err_c:
                            st.error(f"Combination falied: {err_c}")
                        else:
                            st.success("Compounds combined successfully!")

                            # ---------- DISPLAY BASIC INFO ----------
                            st.markdown("### 🧬 Compound Overview")
                            st.write(f"**A:** {name_a} ({formula_a})")
                            st.write(f"**B:** {name_b} ({formula_b})")

                            # ---------- PROPERTIES ----------
                            st.markdown("### 📊 Properties Comparison")

                            colA, colB, colC = st.columns(3, gap="large")

                            with colA:
                                st.markdown("#### 🅰 Compound A")
                                st.write(f"**Name:** {name_a}")
                                st.write(f"**Formula:** {formula_a}")
                                for k, v in atomic_properties(mol_a).items():
                                    st.write(f"**{k}:** {v}")
                                for k, v in structural_properties(mol_a).items():
                                    st.write(f"**{k}:** {v}")

                            with colB:
                                st.markdown("#### 🅱 Compound B")
                                st.write(f"**Name:** {name_b}")
                                st.write(f"**Formula:** {formula_b}")
                                for k, v in atomic_properties(mol_b).items():
                                    st.write(f"**{k}:** {v}")
                                for k, v in structural_properties(mol_b).items():
                                    st.write(f"**{k}:** {v}")
                            
                            with colC:
                                st.markdown("#### 🔗 Combined")
                                st.write(f"**Source:** {name_a} + {name_b}")
                                atomic_properties(combined_mol)
                                structural_properties(combined_mol)


                            # ---------- 3D VIEW ----------
                            col1, col2 = st.columns(2, gap="large")

                            with col1:
                                st.subheader("2D Molecular Structure")
                                render_2d_structure(mol)

                                st.download_button("Download 2D Img", data=Chem.MolToSmiles(mol), file_name= user_input + ".smiles")

                            with col2:
                                col1, col2 = st.columns(2, gap="large")
                                with col1:
                                    st.subheader("3D Molecular Structure")
                                with col2:
                                    show_labels = st.toggle("Show Atom Labels", value=True)
                            pdb_block = mol_to_pdb_block(combined_mol)
                            html = render_ngl_viewer(pdb_block)
                            st.components.v1.html(html, height=600)

                            # ---------- AI EXPLANATION ----------

                            # Prepare dynamic values
                            input_list = [
                                mol_to_text(mol_a, name_a),
                                mol_to_text(mol_b, name_b)
                            ]

                            output_list = None
                            output_list = [
                                mol_to_smiles(combined_mol)
                            ]

                            # LLM explanation
                            with st.markdown("🧠 AI Explanation"):
                                explanation = ai_description(
                                    "Combination",
                                    input_list,
                                    output_list,
                                )
                                st.write(explanation)


    with reaction_tabs[1]:
        st.markdown("### 🧯  Ring-Based Decomposition Reaction (AB → A + B)")

        decomp_input = st.text_input(
            "Compound to Decompose (name or formula)",
            key="decomp_input"
        )

        if st.button("Run Decomposition Reaction"):
            if not decomp_input:
                st.error("Please enter a compound.")
            else:
                with st.spinner("Performing ring-based structural decomposition..."):
                    mol, name, formula, smiles, err = load_molecule(decomp_input)

                    if err:
                        st.error(err)
                    else:
                        ring_mols, chain_mol = ring_based_decomposition(mol)

                        st.success("Decomposition completed (structural).")

                        # ---- ORIGINAL MOLECULE ----
                        st.markdown("## 🧬 Original Molecule")
                        pdb_block = mol_to_pdb_block(mol)
                        render_ngl_viewer(pdb_block, height=900)

                        # ---- RING FRAGMENTS ----
                        if ring_mols:
                            st.markdown("## 🔵 Ring Fragment(s)")
                            for i, r in enumerate(ring_mols, 1):
                                st.markdown(f"### Ring {i}")
                                atomic_properties(r)
                                structural_properties(r)
                                pdb_block = mol_to_pdb_block(r)
                                html = render_ngl_viewer(pdb_block)
                                st.components.v1.html(html, height=600)
                        else:
                            st.info("No ring structures detected.")

                        # ---- CHAIN FRAGMENT ----
                        if chain_mol:
                            st.markdown("## 🟢 Chain / Side-Chain Fragment")
                            atomic_properties(chain_mol)
                            structural_properties(chain_mol)
                            pdb_block = mol_to_pdb_block(chain_mol)
                            html = render_ngl_viewer(pdb_block)
                            st.components.v1.html(html, height=600)

                        # ---------- AI EXPLANATION ----------

                            # Prepare dynamic values
                            input_list = [
                                mol_to_text(mol)
                            ]

                            output_list = None
                            output_list = [
                                (mol_to_smiles(r) for r in ring_mols),
                                mol_to_smiles(chain_mol)
                            ]

                            # LLM explanation
                            with st.subheader("🧠 AI Explanation"):
                                explanation = ai_description(
                                    "Decomposition",
                                    input_list,
                                    output_list
                                )
                                st.write(explanation)


    with reaction_tabs[2]:
        st.markdown("### 🔁 Single Displacement (A + BC → AC + B)")

        col1, col2 = st.columns(2)
        with col1:
            single_a = st.text_input("Compound A", key="single_a")
        with col2:
            single_bc = st.text_input("Compound BC", key="single_bc")

        if st.button("Run Single Displacement Reaction"):
            with st.spinner("Simulating single displacement..."):
                mol_a, _, _, _, err1 = load_molecule(single_a)
                mol_bc, _, _, _, err2 = load_molecule(single_bc)

                if err1 or err2:
                    st.error("Failed to load compounds.")
                else:
                    prod_ac, prod_b, err = single_displacement(mol_a, mol_bc)

                    if err:
                        st.error(err)
                    else:
                        st.success("Single displacement completed.")

                        col1, col2 = st.columns(2, gap="large")
                        with col1:
                            st.markdown("### Molecule A")
                            pdb_block = mol_to_pdb_block(mol_a)
                            html = render_ngl_viewer(pdb_block)
                            st.components.v1.html(html, height=600)

                        with col2:
                            st.markdown("### Molecule BC")
                            pdb_block = mol_to_pdb_block(mol_bc)
                            html = render_ngl_viewer(pdb_block)
                            st.components.v1.html(html, height=600)

                        col3, col4 = st.columns(2)
                        with col3:
                            st.markdown("### Product AC")
                            pdb_block = mol_to_pdb_block(prod_ac)
                            html = render_ngl_viewer(pdb_block)
                            st.components.v1.html(html, height=600)

                        with col4:
                            st.markdown("### Product B")
                            pdb_block = mol_to_pdb_block(prod_b)
                            html = render_ngl_viewer(pdb_block)
                            st.components.v1.html(html, height=600)

                            # ---------- AI EXPLANATION ----------
                            name_a = ""
                            name_bc = ""

                            # Prepare dynamic values
                            input_list = [
                                mol_to_text(mol_a, name_a),
                                mol_to_text(mol_bc, name_bc)
                            ]


                            output_list = None
                            output_list = [
                                mol_to_smiles(prod_ac),
                                mol_to_smiles(prod_b)
                            ]

                            # LLM explanation
                            with st.subheader("🧠 AI Explanation"):
                                explanation = ai_description(
                                    "Single Displacement",
                                    input_list,
                                    output_list
                                )
                                st.write(explanation)


    with reaction_tabs[3]:
        st.markdown("### 🔄 Double Displacement (AB + CD → AD + CB)")

        col1, col2 = st.columns(2)
        with col1:
            double_ab = st.text_input("Compound AB", key="double_ab")
        with col2:
            double_cd = st.text_input("Compound CD", key="double_cd")

        if st.button("Run Double Displacement Reaction"):
            with st.spinner("Simulating double displacement..."):
                mol_ab, _, _, _, err1 = load_molecule(double_ab)
                mol_cd, _, _, _, err2 = load_molecule(double_cd)

                if err1 or err2:
                    st.error("Failed to load compounds.")
                else:
                    prod_ad, prod_cb, err = double_displacement(mol_ab, mol_cd)

                    if err:
                        st.error(err)
                    else:
                        st.success("Double displacement completed.")

                        col1, col2 = st.columns(2, gap="large")
                        with col1:
                            st.markdown("### Molecule AB")
                            pdb_block = mol_to_pdb_block(mol_ab)
                            html = render_ngl_viewer(pdb_block)
                            st.components.v1.html(html, height=600)

                        with col2:
                            st.markdown("### Molecule CD")
                            pdb_block = mol_to_pdb_block(mol_cd)
                            html = render_ngl_viewer(pdb_block)
                            st.components.v1.html(html, height=600)

                        col3, col4 = st.columns(2, gap="large")
                        with col3:
                            st.markdown("### Product AD")
                            pdb_block = mol_to_pdb_block(prod_ad)
                            html = render_ngl_viewer(pdb_block)
                            st.components.v1.html(html, height=600)

                        with col4:
                            st.markdown("### Product CB")
                            pdb_block = mol_to_pdb_block(prod_cb)
                            html = render_ngl_viewer(pdb_block)
                            st.components.v1.html(html, height=600)

                        # ---------- AI EXPLANATION ----------
                            name_ab = ""
                            name_cd = ""

                            # Prepare dynamic values
                            input_list = [
                                mol_to_text(mol_ab, name_ab),
                                mol_to_text(mol_cd, name_cd)
                            ]

                            output_list = None
                            output_list = [
                                mol_to_smiles(prod_ad),
                                mol_to_smiles(prod_cb)
                            ]

                            # LLM explanation
                            with st.subheader("🧠 AI Explanation"):
                                explanation = ai_description(
                                    "Double Displacement",
                                    input_list,
                                    output_list,
                                )
                                st.write(explanation)

    with reaction_tabs[4]:
        st.markdown("### 🔃 Structural Rearrangement")

        rearr_input = st.text_input(
            "Compound for Rearrangement (name or formula)",
            key="rearr_input"
        )

        if st.button("Run Rearrangement"):
            if not rearr_input:
                st.error("Please enter a compound.")
            else:
                with st.spinner("Performing structural rearrangement..."):
                    mol, name, formula, smiles, err = load_molecule(rearr_input)

                    if err:
                        st.error(err)
                    else:
                        new_mol, err2 = simple_rearrangement(mol)

                        if err2:
                            st.error(err2)
                        else:
                            st.success("Rearrangement completed.")

                            col1, col2 = st.columns(2, gap="large")

                            with col1:
                                st.markdown("#### Original Structure")
                                atomic_properties(mol)
                                structural_properties(mol)
                                pdb_block = mol_to_pdb_block(mol)
                                html = render_ngl_viewer(pdb_block)
                                st.components.v1.html(html, height=600)

                            with col2:
                                st.markdown("#### Rearranged Structure")
                                atomic_properties(new_mol)
                                structural_properties(new_mol)
                                pdb_block = mol_to_pdb_block(new_mol)
                                html = render_ngl_viewer(pdb_block)
                                st.components.v1.html(html, height=600)

                            # ---------- AI EXPLANATION ----------

                            # Prepare dynamic values
                            input_list = [
                                mol_to_text(mol)
                            ]

                            output_list = None
                            output_list = [
                                mol_to_smiles(new_mol)
                            ]

                            # LLM explanation
                            with st.subheader("🧠 AI Explanation"):
                                explanation = ai_description(
                                   "Rearrangement",
                                    input_list,
                                    output_list
                                )
                                st.write(explanation)



with mode_tabs[2]:
    st.markdown("## 🎯 Mode 3: Property Optimization")
    st.info("🎯 Mode 3 suggests property-driven molecular modifications using heuristic rules.")


    compound_input = st.text_input(
        "Enter compound (name or formula)",
        key="mode3_compound"
    )

    goal = st.selectbox(
        "Desired property improvement",
        [
            "Increase solubility",
            "Increase lipophilicity",
            "Reduce molecular weight",
            "Increase stability",
            "Reduce polarity"
        ]
    )

    if st.button("Suggest Modifications"):
        with st.spinner("Analyzing molecule..."):
            mol, name, formula, smiles, err = load_molecule(compound_input)

            if err:
                st.error(err)
            else:
                atomic = atomic_properties(mol)
                structural = structural_properties(mol)

                col1, col2 = st.columns(2, gap="large")

                with col1:
                    st.subheader("📊 Current Properties")
                    st.json({**atomic, **structural})

                    suggestions = suggest_modifications(goal)

                    st.subheader("🧪 Suggested Structural Changes")
                    for i, s in enumerate(suggestions, 1):
                        st.markdown(f"**Step {i}:** {s}")

                with col2:
                    st.subheader("3D Molecular Structure")
                    pdb_block = mol_to_pdb_block(mol)
                    html = render_ngl_viewer(pdb_block)
                    st.components.v1.html(html, height=600)

            with st.expander("🧠 AI Explanation"):
                explanation = ai_description(
                    "Property Optimization",
                    [name, smiles],
                    [goal]
                )
                st.markdown(explanation)


st.markdown("""
---
⚠️ **Disclaimer**  
This system performs *structural and conceptual simulations* of molecules and reactions.  
It does **not** predict real chemical feasibility, reaction kinetics, or laboratory outcomes.
""")



