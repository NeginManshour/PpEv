"""
This code is for cleaning a directory of predicted PDB structures. 
"""

from pathlib import Path
import os
from argparse import ArgumentParser


def is_fixed_pdb(pdb_file):
    """
    Given the path to a PDB file, check if it is a "fixed" PDB,
    i.e. it has been run with fix_numbering.pl
    """
    # IMPORTANT: I changed the fix_numbering.pl code to output files in the following format:
    # example.pdb.pdb
    if pdb_file.suffix == ".pdb" and pdb_file.stem[-4:] == ".pdb":
        return True
    else:
        return False


def clean(pdb_dir, ref_path):
    """
    Clean a directory of predicted PDB structures according to a provided
    reference PDB. Runs `fix_numbering.pl` and `clean_data.py`. 
    :param pdb_dir: (Path or str) path to directory of PDBs to clean
    :param ref_path: (Path or str) path to the reference PDB file
    :return cleaned_dir: (Path) path to the newly created directory of cleaned PDBs
    """

    # Parent directory of this script
    parent_dir = Path(__file__).parent 

    clean_path = parent_dir / "clean_data.py"

    # Iterate over the PDB files in the directory
    # Run fix_numbering.pl on each one
    pdb_list = list(pdb_dir.iterdir())
    fix_path = parent_dir / "fix_numbering.pl"
    fix_count = 0
    for pdb_file in pdb_list:
        if pdb_file.suffix == ".pdb":
            os.system(str(fix_path) + " " + str(pdb_file) + " " + str(ref_path))
            fix_count += 1
    print()
    print("Ran fix_numbering.pl on", fix_count, "files")

    # Create a new directory to hold the fixed PDBs
    # This will be the compare directory for clean_data
    compare_dir = pdb_dir.parent / (pdb_dir.stem + "_compare") 
    # Make the output directory
    try:
        compare_dir.mkdir()
    # Do nothing if directory already exists
    except FileExistsError:
        pass

    # Move all the fixed PDB files into compare_dir
    pdb_list = list(pdb_dir.iterdir())
    move_count = 0
    for pdb_file in pdb_list:
        if is_fixed_pdb(pdb_file):
            os.system("mv -i " + str(pdb_file) + " " + str(compare_dir))
            move_count += 1
    print()
    print("Moved", move_count, "files into compare directory")
    print()

    # Run clean_data.py on each PDB file in compare_dir
    os.system("python " + str(clean_path) + " " + str(compare_dir) + " " + str(ref_path) + " A B")

    # Directory containing cleaned PDBs
    cleaned_dir = compare_dir.parent / ("cleaned_" + compare_dir.stem) 

    print("Cleaned files in", str(cleaned_dir))

    return cleaned_dir


if __name__ == '__main__':

    # Args: path to PDB directory, path to reference PDB
    arg_parser=ArgumentParser(description="Clean a directory of PDBs")
    arg_parser.add_argument('pdb_dir',metavar='<pdb_dir>',type=str,nargs=1,help='path to directory containing predicted PDB files')
    arg_parser.add_argument('ref_path',metavar='<ref_path>',type=str,nargs=1,help='path to reference pdb file')

    args = arg_parser.parse_args()
    pdb_directory = Path(args.pdb_dir[0])
    reference_pdb = Path(args.ref_path[0])

    clean(pdb_directory, reference_pdb)