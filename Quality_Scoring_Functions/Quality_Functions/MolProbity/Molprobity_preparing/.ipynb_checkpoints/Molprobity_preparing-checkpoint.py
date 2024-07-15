from os import path
from optparse import OptionParser
import uuid
from typing import Optional, Dict
from Bio.PDB import PDBParser, SASA
from Bio import PDB
import shutil
import csv
import os
import pandas as pd
import json
import numpy as np
import subprocess


def binding_sites(cutoff, pdbfile, protein_chain ='A', peptide_chain = 'B', CA_only =False):

    parser = PDBParser()

    structure = parser.get_structure('pre_ranked', pdbfile)
    model = structure[0]
    chain_protein = model[protein_chain]
    chain_peptide = model[peptide_chain]

    peptide_len = len(chain_peptide)
    cutoff = float(cutoff)
    peptide_bind_sites = {}
    protein_bind_sites = {}

    for residue1 in chain_protein:
        protein_position = residue1.id[1]
        #print("protein_position = ", protein_position)


        for residue2 in chain_peptide:
            peptide_position = residue2.id[1]
            #print("peptide_position =", peptide_position)

            if CA_only:

                try:
                    distance = residue2['CA'] - residue1['CA']

                except KeyError:
                    print("No CA")
                    continue

                if distance <= cutoff:
                    peptide_bind_sites[peptide_position] = distance
                    protein_bind_sites[protein_position] = distance
            else:
                #all atoms
                for atom1 in residue1.get_atoms():
                    for atom2 in residue2.get_atoms():
                            distance = atom1 - atom2
                            if distance <= cutoff:
                                if((residue1.resname != 'HOH') and (residue2.resname != 'HOH')):
                                    protein_bind_sites[protein_position] = residue1.resname
                                    peptide_bind_sites[peptide_position] = residue2.resname
                                    break

    num_of_protein_binding_sites = len(protein_bind_sites)
    num_of_peptide_binding_sites = len(peptide_bind_sites)
    #print("protein_bind_sites =", protein_bind_sites)
    #print("peptide_bind_sites =", peptide_bind_sites)


    return protein_bind_sites, peptide_bind_sites, num_of_protein_binding_sites, num_of_peptide_binding_sites

############################################################################################################

def Calculate_protein_IOU(protein_pre, protein_nat):

    intersection = 0
    difference = 0

    for key in protein_pre.keys():
        if key not in protein_nat.keys():
            difference += 1
        elif protein_pre[key] != protein_nat[key]:
            difference += 1
        else:
            intersection += 1

    #print("intersection, difference =", intersection, difference)

    protein_IOU = (intersection / (intersection + difference)) * 100

    return protein_IOU

def Calculate_peptide_IOU(peptide_pre, peptide_nat, key_offset):

    intersection = 0
    difference = 0

    for key in peptide_pre.keys():
        shifted_key = key + key_offset
        if shifted_key not in peptide_nat.keys():
            difference += 1
        elif peptide_pre[key] != peptide_nat[shifted_key]:
            difference += 1
        else:
            intersection += 1

    #print("intersection, difference =", intersection, difference)

    peptide_IOU = (intersection / (intersection + difference)) * 100

    return peptide_IOU
############################################################################################################
############################################################################################################

import os
import re
import subprocess

def molprobity(pdbfile):
    # Path to the phenix.molprobity executables
    molprobity_reduce_exe = "phenix.reduce"
    molprobity_exe = "phenix.molprobity"

    # Print the operation being performed
    print('molprobity_score ' + pdbfile)

    # Generate the output file name from the pdbfile
    pdb_reduce_output = pdbfile[:-4] + "FH.pdb"
    pdb_options = f" -FLIP -BUILD {pdbfile} &> {pdb_reduce_output}"

    # Run phenix.reduce
    os.system(f"{molprobity_reduce_exe} {pdb_options}")

    # Run phenix.molprobity and capture the output
    proc = subprocess.Popen([molprobity_exe, pdb_reduce_output], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = proc.communicate()

    # Regular expressions for data extraction
    patterns = {
        "molprobity_score": re.compile(r"MolProbity score\s+=\s+(\d+\.\d+)"),
        "Clashscore": re.compile(r"Clashscore\s*=\s*(\d+\.\d+)"),
        "Cis_nonPro": re.compile(r"Cis-general\s+:\s+(\d+\.\d+)"),
        "Twisted_peptide": re.compile(r"Twisted General\s+:\s+(\d+\.\d+)"),
        "Rama_Z_whole": re.compile(r"whole:\s+([-\d\.]+)"),
        "Rama_Z_helix": re.compile(r"helix:\s+([-\d\.]+)"),
        "Rama_Z_sheet": re.compile(r"sheet:\s+([-\d\.]+)"),
        "Rama_Z_loop": re.compile(r"loop\s*:\s+([-\d\.]+)")
    }

    # Dictionary to hold the results
    results = {}

    try:
        # Read the MolProbity score from the output file
        with open('molprobity.out', 'r') as f:
            output = f.read()

        # Search the file for each pattern and add the results to the dictionary
        for key, pattern in patterns.items():
            match = pattern.search(output)
            if match:
                # Adding "%" to the values for Cis_nonPro and Twisted_peptide
                #results[key] = f"{match.group(1)}%" if key in ["Cis_nonPro", "Twisted_peptide"] else match.group(1)
                results[key] = match.group(1)
            else:
                print(f"Could not find {key} in the output file.")
                results[key] = "Not Found"

    except FileNotFoundError:
        print(f"The molprobity.out file does not exist.")
        results = {key: "Not Found" for key in patterns}

    return results

############################################################################################################
def main():

    os.chdir('/home/nmn5x/data/IOU/Molprobity/TF/7zx4_Mol/7zx4')
    nat_pdb = "7ZX4_BD.pdb"
    pdb_dir = "./"


    #pdbfile_energy = "peptide.pdb"
    #pdbfile_molprobity = "./ranked_0_out.pdb"
    protein_dict_nat, peptide_dict_nat, bind_protein_nat , bind_peptide_nat = binding_sites(5, nat_pdb, 'A', 'B', False)
    #print("bind_protein_native , bind_peptide_native =", bind_protein_nat, bind_peptide_nat)

    #molprobity_score, Clashscore, Cis_nonPro, Twisted_peptide= molprobity(nat_pdb)
    results = molprobity(nat_pdb)
    # Extract each value from the dictionary
    molprobity_score = results.get("molprobity_score", "Not Found")
    Clashscore = results.get("Clashscore", "Not Found")
    Cis_nonPro = results.get("Cis_nonPro", "Not Found")
    Twisted_peptide = results.get("Twisted_peptide", "Not Found")
    Rama_Z_whole = results.get("Rama_Z_whole", "Not Found")
    Rama_Z_helix = results.get("Rama_Z_helix", "Not Found")
    Rama_Z_sheet = results.get("Rama_Z_sheet", "Not Found")
    Rama_Z_loop = results.get("Rama_Z_loop", "Not Found")
    current_directory = os.getcwd()
    print("Current directory_native:", current_directory)
    #Stability = caculate_energy(nat_pdb)

    # Create an empty DataFrame to store the binding site data
    data = pd.DataFrame(columns=["bind_protein_pre", "bind_peptide_pre", "Molprobity_score", "Clashscore", "Cis_nonPro", "Twisted_peptide", "Rama_Z_whole","Rama_Z_helix", "Rama_Z_sheet", "Rama_Z_loop"])

    # Add a row for the native data
    data.loc["native"] = {"bind_protein_pre": bind_protein_nat,"bind_peptide_pre": bind_peptide_nat, "Molprobity_score": molprobity_score, "Clashscore": Clashscore, "Cis_nonPro": Cis_nonPro, "Twisted_peptide": Twisted_peptide, "Rama_Z_whole": Rama_Z_whole, "Rama_Z_helix": Rama_Z_helix, "Rama_Z_sheet": Rama_Z_sheet, "Rama_Z_loop": Rama_Z_loop}
    
    for pdb_file in os.listdir(pdb_dir):
        #if pdb_file.endswith(".pdb") and pdb_file.startswith("ranked"):
        if pdb_file.endswith("_clean.pdb") and pdb_file.startswith("ranked"):
            print("pdb_file =", pdb_file)
            pdbfile = os.path.join(pdb_dir, pdb_file)
            print("pdbfile directory=", pdbfile)
            #os.chdir('/home/nmn5x/data/IOU')
            current_directory = os.getcwd()

            # Print the current directory
            print("Current directory:", current_directory)
            protein_dict_pre, peptide_dict_pre, bind_protein_pre , bind_peptide_pre = binding_sites(5, pdbfile, 'A', 'B', False)
           #print("bind_protein_pre , bind_peptide_pre =", bind_protein_pre, bind_peptide_pre)

            #protein_IOU = Calculate_protein_IOU(protein_dict_pre, protein_dict_nat)
            #protein_IOU = 1
            #peptide_IOU = Calculate_peptide_IOU(peptide_dict_pre, peptide_dict_nat, 368)
            #peptide_IOU = 1
            #molprobity_score, Clashscore, Cis_nonPro, Twisted_peptide= molprobity(pdbfile)
            #molprobity_score = 1
            results = molprobity(pdbfile)
            # Extract each value from the dictionary
            molprobity_score = results.get("molprobity_score", "Not Found")
            Clashscore = results.get("Clashscore", "Not Found")
            Cis_nonPro = results.get("Cis_nonPro", "Not Found")
            Twisted_peptide = results.get("Twisted_peptide", "Not Found")
            Rama_Z_whole = results.get("Rama_Z_whole", "Not Found")
            Rama_Z_helix = results.get("Rama_Z_helix", "Not Found")
            Rama_Z_sheet = results.get("Rama_Z_sheet", "Not Found")
            Rama_Z_loop = results.get("Rama_Z_loop", "Not Found")
            #Stability = caculate_energy(pdb_file)
            # Print results or handle them as needed
            print(f"Results for {pdbfile}:")
            print(f"MolProbity Score: {molprobity_score}")
            print(f"Clashscore: {Clashscore}")
            print(f"Cis-nonPro: {Cis_nonPro}")
            print(f"Twisted Peptide: {Twisted_peptide}")
            print(f"Rama-Z Whole: {Rama_Z_whole}")
            print(f"Rama-Z Helix: {Rama_Z_helix}")
            print(f"Rama-Z Sheet: {Rama_Z_sheet}")
            print(f"Rama-Z Loop: {Rama_Z_loop}")
            # Add a row to the DataFrame with the binding site data for this PDB file
            data.loc[pdb_file] = {"bind_protein_pre": bind_protein_pre, "bind_peptide_pre": bind_peptide_pre, "Molprobity_score": molprobity_score, "Clashscore": Clashscore, "Cis_nonPro": Cis_nonPro, "Twisted_peptide": Twisted_peptide, "Rama_Z_whole": Rama_Z_whole, "Rama_Z_helix": Rama_Z_helix, "Rama_Z_sheet": Rama_Z_sheet, "Rama_Z_loop": Rama_Z_loop}

    # Write the DataFrame to an Excel file
    os.chdir('/home/nmn5x/data/IOU/Molprobity/TF/7zx4_Mol')
    data.to_excel("Binidng_Molprobity_TF_7zx4.xlsx")
    #energy = caculate_energy(pdbfile_energy)
    #Molprobity_score = molprobity(pdbfile_molprobity)


if __name__ == "__main__":
    main()

