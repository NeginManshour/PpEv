from pathlib import Path
from Bio.PDB import PDBParser, PDBIO
from argparse import ArgumentParser

# Args: path of pdb file, Chain 1 ID, Chain 2 ID
arg_parser=ArgumentParser(description="Clean native PDB")
arg_parser.add_argument('pdb_path',metavar='<pdb_path>',type=str,nargs=1,help='path to pdb file')
arg_parser.add_argument('chain1',metavar='<chain1>',type=str,nargs=1,help='ID of chain 1')
arg_parser.add_argument('chain2',metavar='<chain2>',type=str,nargs=1,help='ID of chain 2')
args = arg_parser.parse_args()

# Load the reference PDB file
parser = PDBParser()

pdb_path = Path(args.pdb_path[0])
pdb_file = open(pdb_path, 'r')
ref_structure = parser.get_structure("ref", pdb_file)
ref_model = ref_structure[0]
pdb_file.close()

# IDs for the 2 chains to keep 
ref_chain1_id = args.chain1[0]
ref_chain2_id = args.chain2[0]

# Find all extra chains
extra_chain_ids = []
for chain in ref_model.get_chains():
    chain_id = chain.get_id()
    if (chain_id != ref_chain1_id) and (chain_id != ref_chain2_id):
        extra_chain_ids.append(chain_id)

# Remove extra chains
for id in extra_chain_ids:
    ref_model.detach_child(id)

# Get the chains of ref_structure
ref_chains = ref_structure.get_chains()
ref_chain1 = next(ref_chains)
ref_chain2 = next(ref_chains)

# Change chain IDs to "A" and "B"
if (ref_chain2.id != "A"):    
    ref_chain1.id = "A"
    ref_chain2.id = "B"
else:
    ref_chain2.id = ref_chain1.id + "1"
    ref_chain1.id = "A"
    ref_chain2.id = "B"

# Renumber atoms in a chain starting with 1
# Also remove hetatm
def renumber_atoms(chain):
    i = 1
    het_ids = [] # List of IDs for hetatm

    for ref_residue in chain:
        # Remove hetatm
        if ref_residue.get_id()[0] != ' ': 
            het_ids.append(ref_residue.id)
        else:
            for atom in ref_residue:
                atom.serial_number = i
                i += 1
    
    # Remove hetatm
    for id in het_ids:
        chain.detach_child(id)

renumber_atoms(ref_chain1)
renumber_atoms(ref_chain2)

# Write out the new PDB file
out_path = pdb_path.parent / (pdb_path.stem + "_" + ref_chain1_id + ref_chain2_id + ".pdb")
#out_path = os.path.join(pdb_path[:-4] + "_clean.pdb")
io = PDBIO()
io.set_structure(ref_structure)
io.save(str(out_path))
