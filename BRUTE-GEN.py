from rdkit import Chem
from itertools import product
from rdkit import RDLogger
import mols2grid
RDLogger.DisableLog('rdApp.*')

fragments_smiles = [
    "N",        
    "CC",       
    "O",       
    "C(=O)",    
    "C#N"       
]

# This is the core with '*' in replacement point
core_smiles = "CC1(C)CCC(CC1)NC[C@@H](O)-*"
core = Chem.MolFromSmiles(core_smiles)

# Defining '*' as smarts
rgroup_smarts = Chem.MolFromSmarts("[#0]")


# rgroup_smiles indicated python list of smiles.

rgroup_smiles = fragments_smiles
rgroup_mols = [Chem.MolFromSmiles(smi) for smi in rgroup_smiles]

# Create a empty list products where smiles will be appended after replacement. 
# You can make it as generator that yields one smiles at a time when called.
products = []

for r_mol in rgroup_mols:
    if r_mol is None:
        continue

    replaced = Chem.ReplaceSubstructs(
        core,
        rgroup_smarts,
        r_mol,
        replaceAll=True
    )

    if replaced and '.' not in Chem.MolToSmiles(replaced[0]):
        try:
            prod = replaced[0]
            Chem.SanitizeMol(prod)
            products.append(Chem.MolToSmiles(prod))
        except:
            pass


# Use mols2grid to display the replaced structure
products = sorted(set(products))
smi = [Chem.MolFromSmiles(i) for i in products]
mols2grid.display(smi)
