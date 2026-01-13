def suggest_modifications(goal):
    rules = {
        "Increase solubility": [
            "Add a hydroxyl (–OH) group",
            "Introduce an amine functional group",
            "Reduce hydrophobic aromatic rings"
        ],
        "Increase lipophilicity": [
            "Add alkyl side chains",
            "Replace polar groups with non-polar groups",
            "Increase aromatic content"
        ],
        "Reduce molecular weight": [
            "Remove long alkyl chains",
            "Simplify ring systems",
            "Remove bulky substituents"
        ],
        "Increase stability": [
            "Reduce reactive functional groups",
            "Add aromatic stabilization",
            "Remove strained ring systems"
        ],
        "Reduce polarity": [
            "Replace –OH with –CH₃",
            "Reduce heteroatoms",
            "Increase hydrocarbon content"
        ]
    }

    return rules.get(goal, ["No suggestions available"])
