# import pyperclip
from argparse import ArgumentParser
from pathlib import Path


"""
Print out a list of all the PDB file names in a directory, and copy the list to clipboard (optional).
For use with HADDOCK scoring on Lewis. Requires the PDBs to follow a certain name format. 
Args: PDB directory path. start index, stop index
"""

# dir_path = "/home/renjz/data/alphafold/unique_out_TB/ranked-34_8ahs/"

def main():

    arg_parser=ArgumentParser(description="Print and copy PDB list")
    arg_parser.add_argument('dir_path',metavar='<dir_path>',type=str,nargs=1,help='path to pdb directory')
    arg_parser.add_argument('start',metavar='<start>',type=str,nargs=1,help='number at which to start')
    arg_parser.add_argument('stop',metavar='<stop>',type=str,nargs=1,help='number at which to stop (exclusive)')
    args = arg_parser.parse_args()

    dir_path = args.dir_path[0]
    start = int(args.start[0])
    stop = int(args.stop[0])

    list = get_list(dir_path, start, stop)
    # pyperclip.copy(list)
    # print("paste=", pyperclip.paste())
    print(list)


# Convert a list to a string
def list_to_str(list):
    str = "["
    for i in range(len(list)):
        str += "'" + list[i] + "'"
        if i < len(list) - 1:
            str += ", "
    str += "]"
    return str


def get_list(dir_path, start, stop):

    # Name format: ranked_[#].pdb_clean.pdb

    pdb_list = [] # List of PDB files
    for i in range(start, stop):
        name = "ranked_" + str(i) + ".pdb_clean.pdb"
        path = dir_path + name
        pdb_list.append(path)

    list_str = list_to_str(pdb_list)
    return list_str


if __name__ == '__main__':
  main() 