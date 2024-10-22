# Comprehensive Evaluation of AlphaFold-Multimer, AlphaFold3 and ColabFold, and Scoring Functions in Predicting Protein-Peptide Complex Structures

## Abstract
This project focuses on evaluating three state-of-the-art tools—AlphaFold-Multimer, ColabFold, and AlphaFold3—for predicting protein-peptide complex structures using Template-Based (TB) and Template-Free (TF) methods. AlphaFold-Multimer excels in TB predictions, while AlphaFold3 produces superior protein-like structures. ColabFold shows versatility in both settings. The study also assesses various scoring functions, revealing that while AlphaFold's built-in scoring is the best, combining multiple scoring functions can enhance the accuracy of predictions. The findings highlight the potential for improving protein-peptide docking predictions by leveraging these tools and scoring strategies.

![Fig_1](https://github.com/user-attachments/assets/02b69ce0-9fb7-472b-b089-82456aaaa6f1)

## Preparing the Predicted Structures for the Evaluation Process

### Installation
To preprocess the predicted structures, the environment below should be installed:
```commandline
conda env create -f evaluation.yaml
```
Each qualifying and scoring function has its corresponding environment, which is provided in its respective directory.

### Cleaning Data
Before evaluating the predicted structures, it is important to clean the PDB files (i.e. remove unnecessary lines, fix residue numbering, etc.). 

First, clean the native PDB structure by running `Quality_Scoring_Functions/Quality_Functions/DockQ/DockQ_Preparing/clean_native.py` using the following command:
```commandline
python clean_native.py path/to/pdb chain_id1 chain_id2
```
You must provide as arguments the path to the native structure PDB file as well as the IDs of the two chains you wish to keep. All other chains will be removed.

Next, clean the directory of predicted PDB structures by running the following command: 
```commandline
python clean_predicted.py path/to/pdb_directory path/to/native_pdb
```

`clean_predicted.py` runs `fix_numbering.pl` and `clean_data.py`. `fix_numbering.pl` requires you to have 'needle' installed as part of the EMBOSS package: http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/needle.html.

The cleaned PDB files will be located in a directory called `cleaned_{original directory name}_compare`. 

### DockQ

Install DockQ v2.1.1 from this GitHub repository: https://github.com/bjornwallner/DockQ.

Run `Quality_Scoring_Functions/Quality_Functions/DockQ/DockQ_Preparing/run_dockq_bulk.py` using the following command:

```commandline
python run_dockq_bulk.py {pdb_dir} {native_pdb}
```

`pdb_dir` is the directory of PDB files, and `native_pdb` is the PDB of the native structure. The DockQ results will be saved in an Excel (.xlsx) file located in the parent directory `pdb_dir`.

Optional arguments:
- `-xl`: Specify a path to an existing Excel file to save the data
- `-clean`: Clean the PDB files by running `fix_numbering.pl` and `clean_data.py` before running DockQ

## Datasets
The cleaned datasets used for this study can be accessed from the following links:

- [AlphaFold-Multimer, Template-based Dataset(AFM-TB)](https://drive.google.com/file/d/1p2cHTfgjrTj1wHCPjcFqs-7zPzNysxZq/view?usp=drive_link)
- [AlphaFold-Multimer, Template-free Dataset(AFM-TF)](https://drive.google.com/file/d/1ATmbF25mEcMMFyEGb02PRFv-aPLl8H63/view?usp=drive_link)
- [ColabFold, Template-Based Dataset(CF-TB)](https://drive.google.com/file/d/1tCRil26Ggz8YIaRq3LCk-OFHdJhw-PVM/view?usp=drive_link)
- [ColabFold, Template-Free Dataset(CF-TF)](https://drive.google.com/file/d/1GUYmu-BZ9bWJzK_yGUgNn_miAd128kDA/view?usp=drive_link)
- [AlphaFold3 Dataset(AF3)](https://drive.google.com/file/d/1eESMSK3aLyFf7y5JbI_eYO29j7xAkVdP/view?usp=drive_link)
