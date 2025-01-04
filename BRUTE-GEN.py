from rdkit import Chem

# Original SMILES
original_smiles = "C1=CC=C2C(=C1)C(=C(N2)O)C3=NC4=CC=CC=C4C3=O"

# Load the molecule
mol = Chem.MolFromSmiles(original_smiles)

# SMARTS pattern for carbonyl oxygen
carbonyl_smarts = "[OH]"  # Correct SMARTS pattern to match carbonyl groups
carbonyl_pattern = Chem.MolFromSmarts(carbonyl_smarts)

# Fragments to replace the carbonyl oxygen
fragments_smiles = [
    "N",        # Amine
    "CC",       # Ethyl group
    "O",        # Hydroxyl group
    "C(=O)",    # Carbonyl group (for demonstration)
    "C#N"       # Nitrile group
]

new_smiles_list = []

if carbonyl_pattern:
    for fragment_smiles in fragments_smiles:
        # Create a molecule from the fragment SMILES
        fragment_mol = Chem.MolFromSmiles(fragment_smiles)
        if fragment_mol:
            try:
                # Replace the carbonyl oxygen with the fragment
                new_mol = Chem.ReplaceSubstructs(mol, carbonyl_pattern, fragment_mol, replaceAll=True)[0]
                # Generate SMILES for the new molecule
                new_smiles = Chem.MolToSmiles(new_mol)
                new_smiles_list.append(new_smiles)
            except Exception as e:
                print(f"Error processing fragment {fragment_smiles}: {e}")
        else:
            print(f"Failed to parse fragment SMILES: {fragment_smiles}")
else:
    print("Failed to parse SMARTS pattern for carbonyl oxygen.")

# Print all generated SMILES
for smi in new_smiles_list:
    print(smi)
