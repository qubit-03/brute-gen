from rdkit import Chem

original_smiles = "C1=CC=C2C(=C1)C(=C(N2)O)C3=NC4=CC=CC=C4C3=O"
mol = Chem.MolFromSmiles(original_smiles)

# SMARTS pattern for carbonyl oxygen
carbonyl_smarts = "[OH]"  #Change it with any other group X to replace it
carbonyl_pattern = Chem.MolFromSmarts(carbonyl_smarts)

# Example fragments for demonstration
fragments_smiles = [
    "N",        
    "CC",       
    "O",       
    "C(=O)",    
    "C#N"       
]

new_smiles_list = []

if carbonyl_pattern:
    for fragment_smiles in fragments_smiles:
        
        fragment_mol = Chem.MolFromSmiles(fragment_smiles)
        if fragment_mol:
            try:
            
                new_mol = Chem.ReplaceSubstructs(mol, carbonyl_pattern, fragment_mol, replaceAll=True)[0]
            
                new_smiles = Chem.MolToSmiles(new_mol)
                new_smiles_list.append(new_smiles)
            except Exception as e:
                print(f"Error processing fragment {fragment_smiles}: {e}")
        else:
            print(f"Failed to parse fragment SMILES: {fragment_smiles}")
else:
    print("Failed to parse SMARTS pattern for carbonyl oxygen.")

for smi in new_smiles_list:
    print(smi)

