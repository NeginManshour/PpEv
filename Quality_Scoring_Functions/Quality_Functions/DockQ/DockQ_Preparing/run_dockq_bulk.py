"""
Program for processing a large number of predicted protein-peptide complex structures
Given a directory full of predicted PDBs, clean all the PDB files and run DockQ on the cleaned directory.
Write the DockQ results to an Excel file. 
"""

from pathlib import Path
import os
from argparse import ArgumentParser
import pandas as pd
from DockQ.DockQ import load_PDB, run_on_all_native_interfaces
from clean_predicted import clean
from tqdm import tqdm


def get_pdb_rank(file):
    """
    Get the rank of a PDB from its filename. The PDB name is assumed to 
    be in the following format: {string}_{#}.pdb
    :param file: (str) name of PDB file
    :return rank: (int) rank of PDB
    """
    file_split = file.split("_")
    file_split1 = file_split[1].split(".")
    rank = int(file_split1[0]) 
    return rank


def dockq_info_to_list(file, alphafold_rank, results):
    """
    Collect the important DockQ informatation for one PDB in a list.
    :param file: (str) name of PDB file
    :param alphafold_rank: (int) alphafold rank of PDB file
    :param results: (dict) dictionary of DockQ results
    :return info_list: (list) list of DockQ information
    """
    info_list = [file, alphafold_rank, 
                 results["len1"], results["len2"], results["fnat"], results["fnonnat"],
                 results["irms"], results["Lrms"], results["DockQ_F1"], results["DockQ"]]
    return info_list


def sort_key(info_list):
    """
    Specify line sorting criterion for DockQ scores.
    :param info_list: (list) list of DockQ information for one PDB
    :return: (float) the item by which to sort 
    """
    # The DockQ score should be the last element of the list
    return(float(info_list[-1]))


if __name__ == '__main__':

    # Parse arguments
    arg_parser=ArgumentParser(description="Clean and run DockQ on many PDBs")
    arg_parser.add_argument('pdb_dir',metavar='<pdb_dir>',type=str,nargs=1,help='path to directory containing predicted PDB files')
    arg_parser.add_argument('ref_path',metavar='<ref_path>',type=str,nargs=1,help='path to reference pdb file')
    arg_parser.add_argument('-xl',metavar='xl', default='', type=str,nargs=1, help='path to Excel file to save data')
    arg_parser.add_argument('-clean', action='store_true', help='clean data before running DockQ')

    args = arg_parser.parse_args()
    pdb_dir = Path(args.pdb_dir[0])
    ref_path = Path(args.ref_path[0])
    xl_path = None
    if args.xl != "":
        xl_path = Path(args.xl[0])

    if args.clean:
        pdb_dir = clean(pdb_dir, ref_path)

    # Load the native PDB as a Biopython model
    native = load_PDB(str(ref_path))

    # Run DockQ for each PDB file in model_dir
    pdb_list = list(pdb_dir.iterdir())
    dockq_data = [] # List of DockQ data
    for pdb_file in tqdm(pdb_list, desc="Running DockQ on PDB files"):

        if pdb_file.suffix == ".pdb":
            model = load_PDB(str(pdb_file)) # Load predicted model

            # Parse AlphaFold rank from filename
            filename = str(pdb_file.stem)
            af_rank = get_pdb_rank(filename) # Get AlphaFold rank of PDB

            # Run DockQ
            dockq_results, dockq_score = run_on_all_native_interfaces(model, native)
            result_dict = dockq_results[("A", "B")]

            # Extract necessary results 
            dockq_data.append(dockq_info_to_list(filename, af_rank, result_dict))
    
    # Sort data by DockQ score
    dockq_data.sort(reverse=True, key=sort_key)
    for i in range(len(dockq_data)):
        # Add DockQ rank to each sample
        dockq_data[i].append(i)

    # Store the data in a Pandas DataFrame
    columns = ["File Name", "Alphafold Rank", "len1", "len2",
            "Fnat", "Fnonnat", "iRMS", "LRMS", "DockQ_F1", "DockQ", "DockQ Rank"]
    dockq_df = pd.DataFrame(data=dockq_data, columns=columns)

    # Save the data to an Excel file
    if xl_path is None:
        protein_name = (ref_path.stem).split("_")[0]
        xl_path = pdb_dir.parent / (protein_name + "_DockQ_data.xlsx")
    
    dockq_df.to_excel(str(xl_path), sheet_name="DockQ", index=False)
    print("Saved DockQ results to ", str(xl_path))