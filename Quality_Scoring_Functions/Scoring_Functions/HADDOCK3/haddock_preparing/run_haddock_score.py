from argparse import ArgumentParser
import print_pdb_list
from pathlib import Path
import os


"""
Run HADDOCK emscoring-mdscoring on a directory of PDBs (should be 1000 PDB files).
Requires the PDBs to follow a certain name format. 
Args: PDB directory path, total number of PDB files, number of PDB files to include in a batch
1) Get a list of n PDB files from the directory
2) Create and run a cfg file for HADDOCK scoring
3) Repeat until all PDBs in the directory have been scored
4) Read the tsv files containing the output scores, compile all the data into a single file
"""

def create_cfg(output_dir, pdb_dir, start, stop, type):
    """
    Create a CFG file to run HADDOCK emscoring-mdscoring on a list of PDBs.
    Returns a pair: (path to output cfg file, path to HADDOCK output directory)
    output_dir = (str) path to directory where the CFG file will be saved
    pdb_dir = (str) path to PDB directory
    start = (int) index at which to start listing PDBs
    stop = (int) index at which to stop listing PDBs (exclusive)
    type = (str) type of scoring: either "em" or "md"
    """

    # Get a list of PDBs as a string
    pdb_list = print_pdb_list.get_list(pdb_dir, start, stop)

    filename = type + "_" + str(start) + "-" + str(stop - 1)
    out_path = Path(output_dir) / (filename + ".cfg")
    out_file = open(out_path, 'w') # Path of output cfg file

    # Write the file
    run_dir_name = "run-" + filename # Name of the HADDOCK run directory
    out_file.write('run_dir = "' + run_dir_name + '"' + '\n')
    out_file.write("molecules = " + pdb_list + "\n" + "\n")
    out_file.write("[topoaa]" + "\n" + "autohis = true" + "\n" + "\n")
    out_file.write("[" + type + "scoring]" + "\n" + "tolerance = 20" + "\n" + "\n")
    out_file.close()

    run_dir_path = Path(output_dir) / run_dir_name
    return (out_path, run_dir_path)


def parse_data(tsv_file, data_list):
    """
    Parse scoring data from a tsv output file and add the data to a list.
    """
    # Read each line of tsv file, add to data list
    line_str = tsv_file.readline() # Skip the first line (column labels)
    line_str = tsv_file.readline() 
    while(line_str != ''):
        line_list = line_str.split("\t") # Convert line to a list using a tab as separator
        data_list.append(line_list)
        line_str = tsv_file.readline()
    return data_list


def sort_key(line):
    """
    Specify line sorting criterion
    """
    # 4th column = HADDOCK score
    return float(line[3])


def write_data(tsv_file, data_list):
    """
    Write scoring data from a list to a tsv file.
    """

    # First line = column labels
    tsv_file.write("#" + "\t" "old_rank" + "\t" + "structure_name" + "\t" + "score" + "\n")
    
    # Write the data
    for i in range(len(data_list)):
        line = data_list[i]
        # Start each line with the index number
        tsv_file.write(str(i) + "\t")

        # Get original (AlphaFold) rank
        file_split = line[1].split("_")
        file_split1 = file_split[1].split(".")
        af_rank = int(file_split1[0]) # AlphaFold rank
        tsv_file.write(str(af_rank) + "\t")

        # Write content of each line, skip the "structure" and "md5" columns
        line = data_list[i]
        tsv_file.write(line[1] + "\t") # Original structure name
        tsv_file.write(line[3]) # Score


def main():

    arg_parser=ArgumentParser(description="Run HADDOCK scoring on a PDB directory")
    arg_parser.add_argument('dir_path',metavar='<dir_path>',type=str,nargs=1,help='path to pdb directory')
    arg_parser.add_argument('num_files',metavar='<num_files>',type=str,nargs=1,help='total number of PDB files in the directory')
    arg_parser.add_argument('num_batch',metavar='<num_batch>',type=str,nargs=1,help='number of PDBs to score at a time')
    args = arg_parser.parse_args()

    dir_path = args.dir_path[0]
    num_files = int(args.num_files[0])
    num_batch = int(args.num_batch[0])

    cwd = Path().absolute() # Current working directory
    output_dir = cwd / "haddock_scoring" # Output directory
    os.system("mkdir " + str(output_dir))

    emscoring_data = [] # List of emscoring data from output
    mdscoring_data = [] # List of mdscoring data from output

    # Loop through PDB file numbers, create cfg files and run HADDOCK
    for i in range(0, num_files, num_batch):
        # Get start and stop indices
        start = i
        stop = i + num_batch
        # Create emscoring cfg file and get the path to the file and output dir
        cfg = create_cfg(output_dir, dir_path, start, stop, "em")

        os.chdir(str(output_dir)) # Enter the output directory
        os.system("haddock3 " + str(cfg[0])) # Run HADDOCK on cfg file
        
        em_path = cfg[1] / "1_emscoring/emscoring.tsv"
        em_file = open(em_path, 'r') # Open emscoring output file
        emscoring_data = parse_data(em_file, emscoring_data) # Parse the data
        em_file.close()

        # ---------------------------------------------------------------------
        # Create mdscoring cfg file and get the path to the file and output dir
        cfg = create_cfg(output_dir, dir_path, start, stop, "md")

        os.system("haddock3 " + str(cfg[0])) # Run HADDOCK on cfg file

        md_path = cfg[1] / "1_mdscoring/mdscoring.tsv"
        md_file = open(md_path, 'r')
        mdscoring_data = parse_data(md_file, mdscoring_data)
        md_file.close()
    
    # Sort the data lists
    emscoring_data.sort(reverse=False, key=sort_key) 
    mdscoring_data.sort(reverse=False, key=sort_key)

    # Write emscoring data
    em_out_path = output_dir / "emscoring_total.tsv"
    emscoring_file = open(em_out_path, 'w')
    write_data(emscoring_file, emscoring_data)
    emscoring_file.close()
    print("emscoring data written")

    # Write mdscoring data
    md_out_path = output_dir / "mdscoring_total.tsv"
    mdscoring_file = open(md_out_path, 'w')
    write_data(mdscoring_file, mdscoring_data)
    mdscoring_file.close()
    print("mdscoring data written")



if __name__ == '__main__':
  main() 