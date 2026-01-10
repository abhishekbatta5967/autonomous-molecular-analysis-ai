from groq import Groq
from dotenv import load_dotenv
import os
import streamlit as st
from rdkit import Chem

# Load environment variables
load_dotenv()

# Initialize Groq client using environment variable
client = Groq(api_key=os.getenv("GROQ_API_KEY"))


# ---------------------------------------------------------
# PROMPT BUILDER (Unified for all modes)
# ---------------------------------------------------------
def build_chemistry_prompt(
    mode: str,
    reaction_type: str | None,
    inputs: list[str],
    outputs: list[str] | None
):
    """
    Builds a chemistry-focused prompt for compound analysis,
    reaction interpretation, or property optimization.
    """

    prompt = f"""
You are a computational chemistry assistant analyzing molecular structures
and chemical transformations at a conceptual level.

MODE:
{mode}

REACTION TYPE:
{reaction_type if reaction_type else "Not applicable"}

INPUT COMPOUNDS / CONTEXT:
{chr(10).join(f"- {i}" for i in inputs)}

OUTPUT COMPOUNDS / RESULTS:
{chr(10).join(f"- {o}" for o in outputs) if outputs else "None"}

TASKS:

1. Structural Analysis
   - Identify functional groups, ring systems, and bonding patterns.
   - Comment on aromaticity, heteroatoms, and molecular complexity.

2. Reaction Interpretation (if applicable)
   - Classify the reaction correctly.
   - Describe the transformation in terms of bond rearrangement or exchange.
   - Avoid mechanisms, kinetics, or synthesis details.

3. Property Reasoning
   - Discuss qualitative trends in:
     • molecular weight
     • polarity
     • hydrogen bonding capability
     • rigidity vs flexibility
   - Focus on reasoning, not numeric prediction.

4. Optimization Insight
   - Suggest conceptual structural modifications that could improve:
     • solubility
     • stability
     • reactivity control
   - No synthesis routes or experimental conditions.

5. Scientific Clarity
   - Use precise chemistry terminology.
   - Keep explanation concise, educational, and well-structured.
   - Limit response to 180 words.

SAFETY CONSTRAINTS:
- No synthesis instructions
- No lab procedures
- No safety or handling advice
"""

    return prompt.strip()


# ---------------------------------------------------------
# CACHED LLM CALL (Used by ALL modes)
# ---------------------------------------------------------
@st.cache_data(show_spinner=False)
def cached_reaction_description(
    task_type: str,
    _inputs,
    _outputs
):
    """
    Cached LLM explanation generator.

    task_type : str
        Describes the current mode (Compound Analysis,
        Reaction Analysis, Property Optimization, etc.)

    _inputs, _outputs :
        Ignored by Streamlit hashing (can be lists or complex objects)
    """

    prompt = build_chemistry_prompt(
        mode=task_type,
        reaction_type=task_type if "Reaction" in task_type else None,
        inputs=_inputs,
        outputs=_outputs
    )

    response = client.chat.completions.create(
        model="llama-3.1-8b-instant",
        messages=[
            {"role": "system", "content": "You are a chemistry teaching assistant."},
            {"role": "user", "content": prompt}
        ],
        temperature=0.35,
        max_tokens=300
    )

    return response.choices[0].message.content.strip()




def mol_to_text(mol, fallback_name=None):
    if fallback_name:
        return fallback_name
    return Chem.MolToSmiles(mol)

def mol_to_smiles(mol):
    return Chem.MolToSmiles(mol)

    