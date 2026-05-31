#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import json

# Load configuration and define a global variable
def load_config():
    global config
    with open('config_Haddock.json') as config_file:
        config = json.load(config_file)

# Call the function to load the configuration
load_config()


# In[ ]:


import os
import pandas as pd
import json

# Load configuration and define a global variable
def load_config():
    global config
    with open('config_Haddock.json') as config_file:
        config = json.load(config_file)

# Call the function to load the configuration
load_config()

def convert_tsv_to_excel(source_directory, target_directory, subdirectory_name):
    """
    Converts all TSV files found in the source_directory into Excel files
    and saves them into the target_directory along with a specified subdirectory.

    Args:
    - source_directory (str): The path to the directory containing TSV files.
    - target_directory (str): The path to the directory where Excel files will be saved.
    - subdirectory_name (str): The name of the subdirectory where Excel files will be placed.
    """
    # Append the subdirectory to the target directory
    full_target_directory = os.path.join(target_directory, subdirectory_name)

    # Ensure the target directory exists
    if not os.path.exists(full_target_directory):
        os.makedirs(full_target_directory)

    # Iterate over all files in the source directory
    for filename in os.listdir(source_directory):
        if filename.endswith('.tsv'):
            # Construct the full file paths
            source_file_path = os.path.join(source_directory, filename)
            # Change the file extension from .tsv to .xlsx for the output file
            target_file_path = os.path.join(full_target_directory, filename.replace('.tsv', '.xlsx'))

            # Load the TSV file
            df = pd.read_csv(source_file_path, sep='\t')
            # Save the dataframe to an Excel file
            df.to_excel(target_file_path, index=False)

            # Optionally print a confirmation message
            #print(f'Converted {filename} to Excel and saved as {os.path.basename(target_file_path)}')

# Example usage
source_directory = config['Haddock_em_directory']
target_directory = config['output_directories']
subdirectory_name = 'haddock_em_TF_ex'

convert_tsv_to_excel(source_directory, target_directory, subdirectory_name)


# In[ ]:


import os
import pandas as pd
import json


def convert_tsv_to_excel(source_directory, target_directory, subdirectory_name):
    """
    Converts all TSV files found in the source_directory into Excel files
    and saves them into the target_directory along with a specified subdirectory.

    Args:
    - source_directory (str): The path to the directory containing TSV files.
    - target_directory (str): The path to the directory where Excel files will be saved.
    - subdirectory_name (str): The name of the subdirectory where Excel files will be placed.
    """
    # Append the subdirectory to the target directory
    full_target_directory = os.path.join(target_directory, subdirectory_name)

    # Ensure the target directory exists
    if not os.path.exists(full_target_directory):
        os.makedirs(full_target_directory)

    # Iterate over all files in the source directory
    for filename in os.listdir(source_directory):
        if filename.endswith('.tsv'):
            # Construct the full file paths
            source_file_path = os.path.join(source_directory, filename)
            # Change the file extension from .tsv to .xlsx for the output file
            target_file_path = os.path.join(full_target_directory, filename.replace('.tsv', '.xlsx'))

            # Load the TSV file
            df = pd.read_csv(source_file_path, sep='\t')
            # Save the dataframe to an Excel file
            df.to_excel(target_file_path, index=False)

            # Optionally print a confirmation message
            #print(f'Converted {filename} to Excel and saved as {os.path.basename(target_file_path)}')

# Example usage
source_directory = config['Haddock_md_directory']
target_directory = config['output_directories']
subdirectory_name = 'haddock_md_TF_ex'

convert_tsv_to_excel(source_directory, target_directory, subdirectory_name)



# In[ ]:


import os
import shutil
import pandas as pd
from openpyxl import load_workbook
import json


def extract_pdb_id_from_dockq(filename):
    """Extracts the PDB ID from a DockQ filename."""
    parts = filename.split('_')
    return parts[0]  # PDB ID is the first part

def extract_pdb_id_from_haddock(filename):
    """Extracts the PDB ID from a Haddock filename."""
    parts = filename.split('_')
    return parts[2]  # PDB ID is after 'emscoring_total'

def is_valid_excel_file(filepath):
    """Checks if the file is a valid Excel file."""
    try:
        pd.read_excel(filepath, engine='openpyxl')
        return True
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return False

def copy_haddock_sheets_to_dockq(dockq_dir, haddock_dir, output_dir, subdirectory_name):
    """
    Copies sheets from Haddock Excel files to corresponding DockQ Excel files
    into an output directory that includes a specified subdirectory.

    Args:
    - dockq_dir (str): Directory containing DockQ files.
    - haddock_dir (str): Directory containing Haddock files.
    - output_dir (str): Base directory where output files will be saved.
    - subdirectory_name (str): The name of the subdirectory where files will be placed.
    """
    full_output_dir = os.path.join(output_dir, subdirectory_name)
    if not os.path.exists(full_output_dir):
        os.makedirs(full_output_dir)

    dockq_files = os.listdir(dockq_dir)
    haddock_files = os.listdir(haddock_dir)

    # Map Haddock files to their PDB IDs
    haddock_map = {extract_pdb_id_from_haddock(f).upper(): f for f in haddock_files if 'emscoring_total' in f and is_valid_excel_file(os.path.join(haddock_dir, f))}

    for dockq_file in dockq_files:
        if dockq_file.endswith('_TF.xlsx') and 'DockQ_data' in dockq_file:
            pdb_id = extract_pdb_id_from_dockq(dockq_file).upper()
            if pdb_id in haddock_map:
                dockq_path = os.path.join(dockq_dir, dockq_file)
                haddock_path = os.path.join(haddock_dir, haddock_map[pdb_id])
                output_path = os.path.join(full_output_dir, dockq_file)

                # Copy the DockQ file to the output directory if it's not already there
                if not os.path.exists(output_path):
                    shutil.copyfile(dockq_path, output_path)

                try:
                    haddock_df = pd.read_excel(haddock_path, sheet_name='Sheet1', engine='openpyxl')
                except Exception as e:
                    print(f"Failed to read {haddock_path}: {e}")
                    continue

                # Load the existing workbook
                book = load_workbook(output_path)

                # Remove the existing sheet if it exists
                if 'Hadd_em' in book.sheetnames:
                    std = book['Hadd_em']
                    book.remove(std)

                # Save the updated workbook
                book.save(output_path)

                # Add the new sheet to the workbook
                with pd.ExcelWriter(output_path, engine='openpyxl', mode='a') as writer:
                    haddock_df.to_excel(writer, sheet_name='Hadd_em', index=False)

                print(f'Updated {dockq_file} with Hadd_em data for PDB ID {pdb_id}')
            else:
                print(f"No matching Hadd_em file found for {dockq_file} with PDB ID {pdb_id}")

# Use the paths from the already loaded config dictionary
dockq_dir = config['TF_uni']
haddock_dir = os.path.join(config['output_directories'], 'haddock_em_TF_ex')
output_dir = config['output_directories']
subdirectory_name = 'haddock_em_TF_DockQ'

copy_haddock_sheets_to_dockq(dockq_dir, haddock_dir, output_dir, subdirectory_name)


# In[ ]:


import os
import shutil
import pandas as pd
from openpyxl import load_workbook
import json


def extract_pdb_id_from_dockq(filename):
    """Extracts the PDB ID from a DockQ filename."""
    parts = filename.split('_')
    return parts[0]  # PDB ID is the first part

def extract_pdb_id_from_haddock(filename):
    """Extracts the PDB ID from a Haddock filename."""
    parts = filename.split('_')
    return parts[2]  # PDB ID is after 'emscoring_total'

def is_valid_excel_file(filepath):
    """Checks if the file is a valid Excel file."""
    try:
        pd.read_excel(filepath, engine='openpyxl')
        return True
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return False

def copy_haddock_sheets_to_dockq(dockq_dir, haddock_dir, output_dir, subdirectory_name):
    """
    Copies sheets from Haddock Excel files to corresponding DockQ Excel files
    into an output directory that includes a specified subdirectory.

    Args:
    - dockq_dir (str): Directory containing DockQ files.
    - haddock_dir (str): Directory containing Haddock files.
    - output_dir (str): Base directory where output files will be saved.
    - subdirectory_name (str): The name of the subdirectory where files will be placed.
    """
    full_output_dir = os.path.join(output_dir, subdirectory_name)
    if not os.path.exists(full_output_dir):
        os.makedirs(full_output_dir)

    dockq_files = os.listdir(dockq_dir)
    haddock_files = os.listdir(haddock_dir)

    # Map Haddock files to their PDB IDs
    haddock_map = {extract_pdb_id_from_haddock(f).upper(): f for f in haddock_files if 'mdscoring_total' in f and is_valid_excel_file(os.path.join(haddock_dir, f))}

    for dockq_file in dockq_files:
        if dockq_file.endswith('_TF.xlsx') and 'DockQ_data' in dockq_file:
            pdb_id = extract_pdb_id_from_dockq(dockq_file).upper()
            if pdb_id in haddock_map:
                dockq_path = os.path.join(dockq_dir, dockq_file)
                haddock_path = os.path.join(haddock_dir, haddock_map[pdb_id])
                output_path = os.path.join(full_output_dir, dockq_file)

                # Copy the DockQ file to the output directory if it's not already there
                if not os.path.exists(output_path):
                    shutil.copyfile(dockq_path, output_path)

                try:
                    haddock_df = pd.read_excel(haddock_path, sheet_name='Sheet1', engine='openpyxl')
                except Exception as e:
                    print(f"Failed to read {haddock_path}: {e}")
                    continue

                # Load the existing workbook
                book = load_workbook(output_path)

                # Remove the existing sheet if it exists
                if 'Hadd_md' in book.sheetnames:
                    std = book['Hadd_md']
                    book.remove(std)

                # Save the updated workbook
                book.save(output_path)

                # Add the new sheet to the workbook
                with pd.ExcelWriter(output_path, engine='openpyxl', mode='a') as writer:
                    haddock_df.to_excel(writer, sheet_name='Hadd_md', index=False)

                print(f'Updated {dockq_file} with Hadd_md data for PDB ID {pdb_id}')
            else:
                print(f"No matching Hadd_md file found for {dockq_file} with PDB ID {pdb_id}")

# Use the paths from the already loaded config dictionary
dockq_dir = config['TF_uni']
haddock_dir = os.path.join(config['output_directories'], 'haddock_md_TF_ex')
output_dir = config['output_directories']
subdirectory_name = 'haddock_md_TF_DockQ'

copy_haddock_sheets_to_dockq(dockq_dir, haddock_dir, output_dir, subdirectory_name)


# In[ ]:


import os
import pandas as pd
from scipy.stats import spearmanr

def calculate_spearman_for_directory(directory_path, output_file):
    # List to store results
    results = []
    missing_sheet_files = []

    # Iterate through all files in the directory
    for filename in os.listdir(directory_path):
        if filename.endswith(".xlsx"):
            file_path = os.path.join(directory_path, filename)

            try:
                # Load the Excel file
                excel_data = pd.ExcelFile(file_path)
            except Exception as e:
                print(f'An error occurred while loading the file "{filename}": {e}')
                continue

            # Check if the expected sheets exist
            if 'Sheet' in excel_data.sheet_names and 'Hadd_em' in excel_data.sheet_names:
                # Load the data from each sheet
                sheet_data = pd.read_excel(file_path, sheet_name='Sheet')
                hadd_data = pd.read_excel(file_path, sheet_name='Hadd_em')

                # Extract relevant columns and rename them for clarity
                sheet_relevant = sheet_data[['AlphaFold Rank', 'DockQ Rank']]
                hadd_relevant = hadd_data[['old_rank', '#']].copy()

                # Create a mapping from 'AlphaFold Rank' to 'DockQ Rank'
                alpha_to_dockq_map = sheet_relevant.set_index('AlphaFold Rank')['DockQ Rank'].to_dict()

                # Map 'old_rank' in hadd_relevant to 'DockQ Rank' using the created mapping
                hadd_relevant.loc[:, 'DockQ Rank'] = hadd_relevant['old_rank'].map(alpha_to_dockq_map)

                # Drop rows with NaN values in 'DockQ Rank'
                filtered_hadd_df = hadd_relevant.dropna(subset=['DockQ Rank'])

                # Calculate the Spearman correlation
                if not filtered_hadd_df.empty:
                    spearman_corr, _ = spearmanr(filtered_hadd_df['#'], filtered_hadd_df['DockQ Rank'])

                    # Append the result to the list
                    results.append({'File Name': filename, 'Spearman Correlation': spearman_corr})
                else:
                    results.append({'File Name': filename, 'Spearman Correlation': None})
            else:
                missing_sheet_files.append(filename)

    # Create a DataFrame from the results
    results_df = pd.DataFrame(results)

    # Save the results to an Excel file
    results_df.to_excel(output_file, index=False)

    # Print files missing the required sheets
    if missing_sheet_files:
        print("Files missing the 'Sheet' or 'Hadd_em' sheet:")
        for file in missing_sheet_files:
            print(file)

# Usage
directory_path = os.path.join(config['output_directories'], 'haddock_em_TF_DockQ')
output_directory = config['Spearman_Correlation_directory_TF']
output_file = os.path.join(output_directory, "correlations_hadd_em_TF.xlsx")

calculate_spearman_for_directory(directory_path, output_file)


# In[ ]:


import pandas as pd

def calculate_positive_negative_percentages_and_average(file_path):
    # Load the Excel file
    data = pd.read_excel(file_path)

    # Count the total number of entries
    total_entries = len(data['Spearman Correlation'])

    # Count the number of positive and negative Hadd_em_Spearman Correlation values
    positive_count = data[data['Spearman Correlation'] > 0].shape[0]
    negative_count = data[data['Spearman Correlation'] < 0].shape[0]

    # Calculate percentages
    positive_percentage = (positive_count / total_entries) * 100
    negative_percentage = (negative_count / total_entries) * 100

    # Calculate the average of Hadd_em_Spearman Correlation values
    average_correlation = data['Spearman Correlation'].mean()

    return positive_percentage, negative_percentage, average_correlation

# Example usage
file_path = f"{config['Spearman_Correlation_directory_TF']}/correlations_hadd_em_TF.xlsx"
positive_percentage, negative_percentage, average_correlation = calculate_positive_negative_percentages_and_average(file_path)
print(f"Positive Hadd_em_Spearman Correlation: {positive_percentage:.2f}%")
print(f"Negative Hadd_em_Spearman Correlation: {negative_percentage:.2f}%")
print(f"Average Hadd_em_Spearman Correlation: {average_correlation:.2f}")


# In[ ]:


import os
import pandas as pd

# Directory containing the Excel files
directory_path = f"{config['output_directories']}/haddock_em_TF_DockQ"

# Initialize a list to store the results
results = []

# List all Excel files in the directory
excel_files = [file for file in os.listdir(directory_path) if file.endswith('.xlsx')]

for file_name in excel_files:
    file_path = os.path.join(directory_path, file_name)
    xls = pd.ExcelFile(file_path)

    # Extract pdb_id from file name
    pdb_id = file_name.split('_')[0]

    # Load data from "Foldx" sheet
    foldx_df = pd.read_excel(file_path, sheet_name='Hadd_em')
    structure_name_for_rank_zero = foldx_df[foldx_df['#'] == 0]['structure_name'].iloc[0].replace('_clean_haddock.pdb', '')

    # Load data from "Sheet" sheet
    sheet_df = pd.read_excel(file_path, sheet_name='Sheet')

    # Highest DockQ score
    highest_dockq_score = sheet_df['DockQ'].max()

    # DockQ score for the structure from Foldx
    highest_foldx_dockq_score = sheet_df[sheet_df['File Name'].str.contains(structure_name_for_rank_zero)]['DockQ'].max()

    # Calculate Loss
    score_Loss = highest_dockq_score - highest_foldx_dockq_score

    # Append results
    results.append([pdb_id, highest_dockq_score, highest_foldx_dockq_score, score_Loss])

# Create a DataFrame from the results
results_df = pd.DataFrame(results, columns=['File Name', 'DockQ', 'Haddock_emscore ranked', 'Loss'])

# Save the results into a new Excel file
output_directory = config["DockQ_Loss_directory_TF"]
output_path = os.path.join(output_directory, 'Haddock_emscore_TF.xlsx')
results_df.to_excel(output_path, index=False)

print(f"Results have been saved to {output_path}")



# In[ ]:


import os
import pandas as pd
from scipy.stats import spearmanr

def calculate_spearman_for_directory(directory_path, output_file):
    # List to store results
    results = []
    missing_sheet_files = []

    # Iterate through all files in the directory
    for filename in os.listdir(directory_path):
        if filename.endswith(".xlsx"):
            file_path = os.path.join(directory_path, filename)

            try:
                # Load the Excel file
                excel_data = pd.ExcelFile(file_path)
            except Exception as e:
                print(f'An error occurred while loading the file "{filename}": {e}')
                continue

            # Check if the expected sheets exist
            if 'Sheet' in excel_data.sheet_names and 'Hadd_md' in excel_data.sheet_names:
                # Load the data from each sheet
                sheet_data = pd.read_excel(file_path, sheet_name='Sheet')
                hadd_data = pd.read_excel(file_path, sheet_name='Hadd_md')

                # Extract relevant columns and rename them for clarity
                sheet_relevant = sheet_data[['AlphaFold Rank', 'DockQ Rank']]
                hadd_relevant = hadd_data[['old_rank', '#']].copy()

                # Create a mapping from 'AlphaFold Rank' to 'DockQ Rank'
                alpha_to_dockq_map = sheet_relevant.set_index('AlphaFold Rank')['DockQ Rank'].to_dict()

                # Map 'old_rank' in hadd_relevant to 'DockQ Rank' using the created mapping
                hadd_relevant.loc[:, 'DockQ Rank'] = hadd_relevant['old_rank'].map(alpha_to_dockq_map)

                # Drop rows with NaN values in 'DockQ Rank'
                filtered_hadd_df = hadd_relevant.dropna(subset=['DockQ Rank'])

                # Calculate the Spearman correlation
                if not filtered_hadd_df.empty:
                    spearman_corr, _ = spearmanr(filtered_hadd_df['#'], filtered_hadd_df['DockQ Rank'])

                    # Append the result to the list
                    results.append({'File Name': filename, 'Spearman Correlation': spearman_corr})
                else:
                    results.append({'File Name': filename, 'Spearman Correlation': None})
            else:
                missing_sheet_files.append(filename)

    # Create a DataFrame from the results
    results_df = pd.DataFrame(results)

    # Save the results to an Excel file
    results_df.to_excel(output_file, index=False)

    # Print files missing the required sheets
    if missing_sheet_files:
        print("Files missing the 'Sheet' or 'Hadd_md' sheet:")
        for file in missing_sheet_files:
            print(file)

directory_path = os.path.join(config['output_directories'], 'haddock_md_TF_DockQ')
output_directory = config['Spearman_Correlation_directory_TF']
output_file = os.path.join(output_directory, "correlations_hadd_md_TF.xlsx")

calculate_spearman_for_directory(directory_path, output_file)



# In[ ]:


import pandas as pd

def calculate_positive_negative_percentages_and_average(file_path):
    # Load the Excel file
    data = pd.read_excel(file_path)

    # Count the total number of entries
    total_entries = len(data['Spearman Correlation'])

    # Count the number of positive and negative Hadd_em_Spearman Correlation values
    positive_count = data[data['Spearman Correlation'] > 0].shape[0]
    negative_count = data[data['Spearman Correlation'] < 0].shape[0]

    # Calculate percentages
    positive_percentage = (positive_count / total_entries) * 100
    negative_percentage = (negative_count / total_entries) * 100

    # Calculate the average of Hadd_em_Spearman Correlation values
    average_correlation = data['Spearman Correlation'].mean()

    return positive_percentage, negative_percentage, average_correlation

# Example usage
file_path = f"{config['Spearman_Correlation_directory_TF']}/correlations_hadd_md_TF.xlsx"
positive_percentage, negative_percentage, average_correlation = calculate_positive_negative_percentages_and_average(file_path)
print(f"Positive Hadd_em_Spearman Correlation: {positive_percentage:.2f}%")
print(f"Negative Hadd_em_Spearman Correlation: {negative_percentage:.2f}%")
print(f"Average Hadd_em_Spearman Correlation: {average_correlation:.2f}")



# In[ ]:


import os
import pandas as pd

# Directory containing the Excel files
#directory_path = '/Users/neginmanshour/Desktop/Haddock_Evaluation/Data/TB/haddock_md_TB_DockQ'
directory_path = f"{config['output_directories']}/haddock_md_TF_DockQ"

# Initialize a list to store the results
results = []

# List all Excel files in the directory
excel_files = [file for file in os.listdir(directory_path) if file.endswith('.xlsx')]

for file_name in excel_files:
    file_path = os.path.join(directory_path, file_name)
    xls = pd.ExcelFile(file_path)

    # Extract pdb_id from file name
    pdb_id = file_name.split('_')[0]

    # Load data from "Foldx" sheet
    foldx_df = pd.read_excel(file_path, sheet_name='Hadd_md')
    structure_name_for_rank_zero = foldx_df[foldx_df['#'] == 0]['structure_name'].iloc[0].replace('_clean_haddock.pdb', '')

    # Load data from "Sheet" sheet
    sheet_df = pd.read_excel(file_path, sheet_name='Sheet')

    # Highest DockQ score
    highest_dockq_score = sheet_df['DockQ'].max()

    # DockQ score for the structure from Foldx
    highest_foldx_dockq_score = sheet_df[sheet_df['File Name'].str.contains(structure_name_for_rank_zero)]['DockQ'].max()

    # Calculate Loss
    score_Loss = highest_dockq_score - highest_foldx_dockq_score

    # Append results
    results.append([pdb_id, highest_dockq_score, highest_foldx_dockq_score, score_Loss])

# Create a DataFrame from the results
results_df = pd.DataFrame(results, columns=['File Name', 'DockQ', 'Haddock_emscore ranked', 'Loss'])

# Save the results into a new Excel file
output_directory = config["DockQ_Loss_directory_TF"]
output_path = os.path.join(output_directory, 'Haddock_mdscore_TF.xlsx')
results_df.to_excel(output_path, index=False)

print(f"Results have been saved to {output_path}")



# In[ ]:




