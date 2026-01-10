import streamlit as st
from chemistry.molecule import load_molecule_cached as load_molecule
from chemistry.properties import (
    atomic_properties,
    structural_properties,
    identity_properties,
    predicted_properties
)
from ui.viewer import render_molecule_3d
from llm.description import cached_reaction_description, mol_to_text, mol_to_smiles
from chemistry.combine import combine_molecules
from chemistry.decomposition import ring_based_decomposition
from chemistry.rearrangement import simple_rearrangement
from chemistry.displacement import single_displacement, double_displacement
from chemistry.optimization import suggest_modifications



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
                col1, col2 = st.columns(2, gap="large")

                with col1:
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

                with col2:
                    st.subheader("3D Molecular Structure")
                    view = render_molecule_3d(mol)
                    st.components.v1.html(view._make_html(), height=540)

            st.subheader("🧠 Compound Description (AI Generated)")

            with st.spinner("🧪 Interpreting chemical structure..."):
                input_list = [name, formula, smiles]
                description = cached_reaction_description("Compound Analysis", input_list, None)

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
                            st.markdown("### 🧊 3D Combined Structure")
                            view = render_molecule_3d(combined_mol)
                            st.components.v1.html(
                                view._make_html(),
                                height=520
                            )

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
                                explanation = cached_reaction_description(
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
                        view = render_molecule_3d(mol)
                        st.components.v1.html(view._make_html(), height=900)

                        # ---- RING FRAGMENTS ----
                        if ring_mols:
                            st.markdown("## 🔵 Ring Fragment(s)")
                            for i, r in enumerate(ring_mols, 1):
                                st.markdown(f"### Ring {i}")
                                atomic_properties(r)
                                structural_properties(r)
                                view = render_molecule_3d(r)
                                st.components.v1.html(view._make_html(), height=900)
                        else:
                            st.info("No ring structures detected.")

                        # ---- CHAIN FRAGMENT ----
                        if chain_mol:
                            st.markdown("## 🟢 Chain / Side-Chain Fragment")
                            atomic_properties(chain_mol)
                            structural_properties(chain_mol)
                            view = render_molecule_3d(chain_mol)
                            st.components.v1.html(view._make_html(), height=900)

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
                                explanation = cached_reaction_description(
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
                            view = render_molecule_3d(mol_a)
                            st.components.v1.html(view._make_html(), height=800)

                        with col2:
                            st.markdown("### Molecule BC")
                            view = render_molecule_3d(mol_bc)
                            st.components.v1.html(view._make_html(), height=800)

                        col3, col4 = st.columns(2)
                        with col3:
                            st.markdown("### Product AC")
                            view = render_molecule_3d(prod_ac)
                            st.components.v1.html(view._make_html(), height=800)

                        with col4:
                            st.markdown("### Product B")
                            view = render_molecule_3d(prod_b)
                            st.components.v1.html(view._make_html(), height=800)

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
                                explanation = cached_reaction_description(
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
                            view = render_molecule_3d(mol_ab)
                            st.components.v1.html(view._make_html(), height=800)

                        with col2:
                            st.markdown("### Molecule CD")
                            view = render_molecule_3d(mol_cd)
                            st.components.v1.html(view._make_html(), height=800)

                        col3, col4 = st.columns(2, gap="large")
                        with col3:
                            st.markdown("### Product AD")
                            view = render_molecule_3d(prod_ad)
                            st.components.v1.html(view._make_html(), height=800)

                        with col4:
                            st.markdown("### Product CB")
                            view = render_molecule_3d(prod_cb)
                            st.components.v1.html(view._make_html(), height=800)

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
                                explanation = cached_reaction_description(
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
                                view = render_molecule_3d(mol)
                                st.components.v1.html(view._make_html(), height=900)

                            with col2:
                                st.markdown("#### Rearranged Structure")
                                atomic_properties(new_mol)
                                structural_properties(new_mol)
                                view = render_molecule_3d(new_mol)
                                st.components.v1.html(view._make_html(), height=900)

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
                                explanation = cached_reaction_description(
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
                    view = render_molecule_3d(mol)
                    st.components.v1.html(view._make_html(), height=540)

            with st.expander("🧠 AI Explanation"):
                explanation = cached_reaction_description(
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



