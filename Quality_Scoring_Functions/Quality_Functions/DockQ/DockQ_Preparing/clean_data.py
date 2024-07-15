import os
from Bio.PDB import PDBParser, PDBIO
from argparse import ArgumentParser
from pathlib import Path
from tqdm import tqdm

re_numbering = 0


### HELPER METHODS ###

# Get the structure of a PDB given its path (Path object)
def load_pdb(pdb_path, id):
    parser = PDBParser()
    pdb = open(pdb_path, 'r')
    structure = parser.get_structure(id, pdb)
    pdb.close()
    return structure

# Change chain names of predicted structures (OPTIONAL)
def change_chain_names(out_chain1, out_chain2, ref_chain1_id, ref_chain2_id):
    # Change chain IDs
    if (out_chain2.id != ref_chain1_id):    
        out_chain1.id = ref_chain1_id
        out_chain2.id = ref_chain2_id
    else:
        out_chain2.id = ref_chain1_id + "1"
        out_chain1.id = ref_chain1_id
        out_chain2.id = ref_chain2_id

# Clean extra residues from the compare structrue, return the final cleaned structure
def clean_structure(compare_structure, ref_structure, ref_chain1_id, ref_chain2_id, re_numbering):
        
    # Select the chains to compare
    ref_chain1 = ref_structure[0][ref_chain1_id]
    ref_chain2 = ref_structure[0][ref_chain2_id]
    compare_chain1 = compare_structure[0]["B"]
    compare_chain2 = compare_structure[0]["C"]

    # Create a dictionary of residue names for each reference chain
    ref_resnames1 = {}
    for ref_residue in ref_chain1:
        ref_resnames1[ref_residue.id[1]] = ref_residue.resname

    ref_resnames2 = {}
    for ref_residue in ref_chain2:
        ref_resnames2[ref_residue.id[1]] = ref_residue.resname
        #print("ref_resnames2 =", ref_resnames2)

    # Create lists of matching residues for each compare chain
    matching_residues1 = []
    for compare_residue in compare_chain1:
        if compare_residue.id[1] in ref_resnames1 and compare_residue.resname == ref_resnames1[compare_residue.id[1]]:
            matching_residues1.append(compare_residue)

    matching_residues2 = []
    for compare_residue in compare_chain2:
        if (compare_residue.id[1] + re_numbering) in ref_resnames2 and compare_residue.resname == ref_resnames2[(compare_residue.id[1] + re_numbering)]:
            matching_residues2.append(compare_residue)

    # Create new structures with only the matching residues
    out_structure = compare_structure.copy()
    out_model = out_structure[0]
    out_chain1 = out_model["B"]
    out_chain1.child_list = matching_residues1
    out_chain2 = out_model["C"]
    out_chain2.child_list = matching_residues2

    # Change predicted structure chain IDs if desired
    change_chain_names(out_chain1, out_chain2, ref_chain1_id, ref_chain2_id)

    return out_structure


######################


# Args: path of comparison pdb file, path of reference pdb file, Chain 1 ID, Chain 2 ID
arg_parser=ArgumentParser(description="Clean predicted PDB")
arg_parser.add_argument('compdir_path',metavar='<compdir_path>',type=str,nargs=1,help='path to directory containing comparison pdb files')
arg_parser.add_argument('ref_path',metavar='<ref_path>',type=str,nargs=1,help='path to reference pdb file')
arg_parser.add_argument('chain1',metavar='<chain1>',type=str,nargs=1,help='ID of chain 1 in reference structure')
arg_parser.add_argument('chain2',metavar='<chain2>',type=str,nargs=1,help='ID of chain 2 in reference structure')
args = arg_parser.parse_args()

# Load the reference PDB
ref_path = Path(args.ref_path[0])
ref_structure = load_pdb(ref_path, "ref")
ref_model = ref_structure[0]

# IDs for the 2 chains to keep 
ref_chain1_id = args.chain1[0]
ref_chain2_id = args.chain2[0]

# Set the directory containing the PDB files to compare
pdb_dir = Path(args.compdir_path[0])

# Output directory located in same directory as pdb_dir
output_dir = pdb_dir.parent / ("cleaned_" + pdb_dir.stem) 
# Make the output directory
try:
    output_dir.mkdir()
# Do nothing if directory already exists
except FileExistsError:
    pass

# Iterate over the pdb files in the directory
pdb_list = list(pdb_dir.iterdir())
for pdb_file in tqdm(pdb_list, desc="Cleaning PDB files"):
    if pdb_file.suffix == ".pdb":
        compare_structure = load_pdb(pdb_file, "comp")
        out_structure = clean_structure(compare_structure, ref_structure, ref_chain1_id, ref_chain2_id, re_numbering)

        # Write out the new PDB file
        out_path = output_dir / (pdb_file.stem + "_clean.pdb")
        # out_path = os.path.join(pdb_path[:-4] + "_clean.pdb")
        io = PDBIO()
        io.set_structure(out_structure)
        io.save(str(out_path))


