#!/usr/bin/python3

import numpy as np
from vina import Vina
import sys

data_directory = str(sys.argv[1])
output_file = str(sys.argv[2]) 
complex_id = str(sys.argv[3])
#complex_id = "7prx"

#Writing header to output csv data file
output = open(output_file, mode = "w")
output.write("complex_id,file,scoring_function,total,lig_inter,flex_inter,other_inter,flex_intra,lig_intra,torsions,best_pose_lig_intra\n")
output.close()

#Opening file list to read in file names

for i in range(9):
    receptor_file = data_directory + "/ranked_" + str(i) + ".pdb_clean_receptor.pdbqt"
    ligand_file = data_directory + "/ranked_" + str(i) + ".pdb_clean_ligand.pdbqt"	

    #Calculating Box Center and Size
    #Read in ligand information
    x = np.array([])
    y = np.array([])
    z = np.array([])

    data_file = open(ligand_file, mode = 'r')
    for line in data_file:
        data = line.split()

        if data[0] == "ATOM":
            x_index = 6 # Index of x-coord in pdbqt file

            # Adjust x_index if there is no space between chain ID and residue number
            if not data[4].isalpha():
                x_index = 5

            x = np.append(x, float(data[x_index]))
            y = np.append(y, float(data[x_index + 1]))
            z = np.append(z, float(data[x_index + 2])) 
            
    data_file.close()
    #Calculate box size (5 A larger than the ligand in each direction)
    centers = []
    max_width = 0
    for array in [x,y,z]:
            max = np.max(array)
            center = (max + np.min(array))/2
            width = max - center
            centers.append(center)
            if width > max_width:
                max_width = width
    box_width = (max_width + 5) * 2
    box_size = [box_width, box_width, box_width]

    for sf in ['vina', 'vinardo']:
        #Scoring using vina
        v = Vina(sf_name=sf)
        v.set_receptor(receptor_file)
        v.set_ligand_from_file(ligand_file)

        v.compute_vina_maps(center=centers, box_size=box_size)
        energy = v.score()

        #Writing output to file as csv data
        output = open(output_file, mode = "a")
        output.write(complex_id + "," + "ranked_" + str(i) + ".pdb" + "," + sf + "," + np.array2string(np.array(energy), separator=",", precision=3, floatmode="fixed")[2:-1] + "\n")
        output.close()
