#!/bin/bash
#SBATCH --partition=hpc4
#SBATCH --time=2-00:00:00
#SBATCH --mem=100G
#SBATCH --job-name=Vina_dock

#reading in input variables

#data_directory="/home/nmn5x/data/IOU/Vina/Data_Vina/7prx_test"
#output_file="/home/nmn5x/data/IOU/Vina/Data_Vina"

#reading in input variables
if ! [ -d "$1" ]
then
    echo "Error: location of data (parent directory of template_free and template_based) needed as first argument"
    echo "Exiting"
    exit
else
    data_directory=$1
fi

#Getting pdb ID from data directory based on the expected file
pdb_id=$(basename "$data_directory")

#optional file name for data output
if [ "$2" != "" ]
then
    output_file=$2
else
    output_file=${SLURM_SUBMIT_DIR}/vina_output_${pdb_id}
fi


#Activating enviroment to produce pdbqt files and adding ADFRSuite to path  
module load miniconda3
source activate auto_dock_vina #path to vina enviroment goes here

export PATH=/home/nmn5x/data/IOU/Vina/ADFRsuite_x86_64Linux_1.0/bin:$PATH #Path to ADFRsuite executable goes here

directories=(${data_directory}/template_free ${data_directory}/template_based)
for d in ${directories[@]}
do 
    cd $d 
    for i in $(seq 0 999) #preparing files for vina
    do
        file_name=${d}/ranked_${i}.pdb_clean
        

        # Check if receptor pdbqt file already exists
        if [[ ! -f "${file_name}_receptor.pdbqt" ]]; then
            if [[ -f "${file_name}_receptor.pdb" ]]; then
                echo "${file_name}_receptor.pdb already exists, skipping generation."
            else
                grep -E "ATOM +[0-9]+ +[A-Z]+ +[A-Z]{3} A[0-9]+" "${file_name}.pdb" > "${file_name}_receptor.pdb"
            fi
            reduce "${file_name}_receptor.pdb"
            prepare_receptor -r "${file_name}_receptor.pdb" -o "${file_name}_receptor.pdbqt"
        else
            echo "${file_name}_receptor.pdbqt already exists, skipping generation."
        fi

        # Check if ligand pdbqt file already exists
        if [[ ! -f "${file_name}_ligand.pdbqt" ]]; then
            if [[ -f "${file_name}_ligand.pdb" ]]; then
                echo "${file_name}_ligand.pdb already exists, skipping generation."
            else
                grep -E "ATOM +[0-9]+ +[A-Z]+ +[A-Z]{3} B[0-9]+" "${file_name}.pdb" > "${file_name}_ligand.pdb"
            fi
            reduce "${file_name}_ligand.pdb"
            prepare_ligand -l "${file_name}_ligand.pdb" -o "${file_name}_ligand.pdbqt"
        else
            echo "${file_name}_ligand.pdbqt already exists, skipping generation."
        fi
    done
        
    

    rm ${d}/ranked_*_ligand.pdb
    rm ${d}/ranked_*_receptor.pdb
    template_state=$(basename "$d")
    python3 ${SLURM_SUBMIT_DIR}/vina_score_and_output_5A.py $d ${output_file}_${template_state}.csv $pdb_id #Scoring with vina
done
