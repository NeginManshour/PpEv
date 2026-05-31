#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import pandas as pd

# Directories where the original Excel files are stored
directories = {
    'AlphaFold_Multimer_All_TB': 'path to /AlphaFold_Multimer_All_TB',
    'Pyrosetta_DockQ_TB': 'path to /Pyrosetta_DockQ_TB',
    'haddock_em_TB_DockQ': 'path to /haddock_em_TB_DockQ',
    'haddock_md_TB_DockQ': 'path to /haddock_md_TB_DockQ',
    'Deep_TB_DockQ': 'path to /Deep_DockQ_TB',
    'Dove_TB_DockQ': 'path to /Dove_DockQ_TB',
    'Vina_TB_DockQ': 'path to /vina_DockQ_TB',
    'Vinardo_TB_DockQ': 'path to /vinardo_DockQ_TB',
    'Foldx_S_TB_DockQ': 'path to /stability_TB',
    'Foldx_int_TB_DockQ': 'path to /interaction_TB'
}

# Output directory to save the merged Excel files
output_directory = 'path to /Heatmaps, Spearman_Correlation/Data/Merged_TB'

# Ensure output directory exists
os.makedirs(output_directory, exist_ok=True)

# Specify which sheet to extract from each Excel file (if they are all named consistently as the first sheet, this could be simplified)
sheet_names = {
    'Sheet': 'AlphaFold_Multimer_All_TB',  # Assuming the default sheet name or modify as needed
    'Pyrosetta': 'Pyrosetta_DockQ_TB',
    'Hadd_em': 'haddock_em_TB_DockQ',
    'Hadd_md': 'haddock_md_TB_DockQ',
    'Deep_GNN': 'Deep_TB_DockQ',
    'gnn_dove': 'Dove_TB_DockQ',
    'Vina': 'Vina_TB_DockQ',
    'Vinardo': 'Vinardo_TB_DockQ',
    'Foldx': 'Foldx_S_TB_DockQ',
    'Foldx_Int': 'Foldx_int_TB_DockQ'

}

# Process each protein sample
for file_name in os.listdir(directories['AlphaFold_Multimer_All_TB']):
    if file_name.endswith('.xlsx'):
        pdb_id = file_name.split('_')[0]
        writer = pd.ExcelWriter(os.path.join(output_directory, f'{pdb_id}_merged.xlsx'), engine='xlsxwriter')

        for sheet, dir_key in sheet_names.items():
            file_path = os.path.join(directories[dir_key], file_name)
            if os.path.exists(file_path):
                data = pd.read_excel(file_path, sheet_name=sheet)  # Load data from the specified sheet
                data.to_excel(writer, sheet_name=sheet, index=False)

        writer.save()

print("All files have been processed and saved.")



# In[ ]:


import os
import pandas as pd

# Directories where the original Excel files are stored
directories = {
    'AlphaFold_Multimer_All_TF': 'path to /AlphaFold_Multimer_All_TF',
    'Pyrosetta_DockQ_TF': 'path to /Pyrosetta_DockQ_TF',
    'haddock_em_TF_DockQ': 'path to /haddock_em_TF_DockQ',
    'haddock_md_TF_DockQ': 'path to /haddock_md_TF_DockQ',
    'Deep_TF_DockQ': 'path to /Deep_DockQ_TF',
    'Dove_TF_DockQ': 'path to /Dove_DockQ_TF',
    'Vina_TF_DockQ': 'path to /vina_DockQ_TF',
    'Vinardo_TF_DockQ': 'path to /vinardo_DockQ_TF',
    'Foldx_S_TF_DockQ': 'path to /stability_TF',
    'Foldx_int_TF_DockQ': 'path to /interaction_TF'
}

# Output directory to save the merged Excel files
output_directory = 'path to /Heatmaps, Spearman_Correlation/Data/Merged_TF'

# Ensure output directory exists
os.makedirs(output_directory, exist_ok=True)

# Specify which sheet to extract from each Excel file (if they are all named consistently as the first sheet, this could be simplified)
sheet_names = {
    'Sheet': 'AlphaFold_Multimer_All_TF',  # Assuming the default sheet name or modify as needed
    'Pyrosetta': 'Pyrosetta_DockQ_TF',
    'Hadd_em': 'haddock_em_TF_DockQ',
    'Hadd_md': 'haddock_md_TF_DockQ',
    'Deep_GNN': 'Deep_TF_DockQ',
    'gnn_dove': 'Dove_TF_DockQ',
    'Vina': 'Vina_TF_DockQ',
    'Vinardo': 'Vinardo_TF_DockQ',
    'Foldx': 'Foldx_S_TF_DockQ',
    'Foldx_Int': 'Foldx_int_TF_DockQ'
}

# Process each protein sample
for file_name in os.listdir(directories['AlphaFold_Multimer_All_TF']):
    if file_name.endswith('.xlsx'):
        pdb_id = file_name.split('_')[0]
        writer = pd.ExcelWriter(os.path.join(output_directory, f'{pdb_id}_merged.xlsx'), engine='xlsxwriter')

        for sheet, dir_key in sheet_names.items():
            file_path = os.path.join(directories[dir_key], file_name)
            if os.path.exists(file_path):
                data = pd.read_excel(file_path, sheet_name=sheet)  # Load data from the specified sheet
                data.to_excel(writer, sheet_name=sheet, index=False)

        writer.save()

print("All files have been processed and saved.")




# In[ ]:


import pandas as pd
import os

def process_excel_files(input_directory, output_directory):
    error_files = []  # List to store names of problematic files

    # Loop through all Excel files in the input directory
    for filename in os.listdir(input_directory):
        if filename.endswith(".xlsx"):
            # Extract PDB ID from the filename (assuming it's before the first underscore)
            pdb_id = filename.split('_')[0]
            input_path = os.path.join(input_directory, filename)
            output_path = os.path.join(output_directory, f"{pdb_id}_Ranked_TB.xlsx")

            try:
                # Process each file
                create_combined_excel(input_path, output_path)
                #print(f"Processed {filename} into {output_path}")
            except Exception as e:
                print(f"Error processing {filename}: {e}")
                error_files.append(filename)

    # Print out the problematic files
    if error_files:
        print("Files that had issues:")
        for file in error_files:
            print(file)

def create_combined_excel(input_path, output_path):
    try:
        # Load data from different sheets
        data = pd.read_excel(input_path, sheet_name='Pyrosetta', engine='openpyxl')
        foldx_data = pd.read_excel(input_path, sheet_name='Foldx', engine='openpyxl')
        foldx_int_data = pd.read_excel(input_path, sheet_name='Foldx_Int', engine='openpyxl')
        hadd_em_data = pd.read_excel(input_path, sheet_name='Hadd_em', engine='openpyxl')
        hadd_md_data = pd.read_excel(input_path, sheet_name='Hadd_md', engine='openpyxl')

        # Optional sheets
        try:
            vina_data = pd.read_excel(input_path, sheet_name='Vina', engine='openpyxl')
            vina_rank = vina_data["alphafold rank"]
        except Exception as e:
            print(f"Vina sheet not found in {input_path}: {e}")
            vina_rank = None
        try:
            vinardo_data = pd.read_excel(input_path, sheet_name='Vinardo', engine='openpyxl')
            vinardo_rank = vinardo_data["alphafold rank"]
        except Exception as e:
            print(f"Vinardo sheet not found in {input_path}: {e}")
            vinardo_rank = None
        try:
            gnn_dove_data = pd.read_excel(input_path, sheet_name='gnn_dove', engine='openpyxl')
            gnn_dove_rank = gnn_dove_data["AlphaFold Rank"]
        except Exception as e:
            print(f"GNN Dove sheet not found in {input_path}: {e}")
            gnn_dove_rank = None
        try:
            deep_gnn_data = pd.read_excel(input_path, sheet_name='Deep_GNN', engine='openpyxl')
            deep_gnn_rank = deep_gnn_data["AlphaRank"]
        except Exception as e:
            print(f"Deep GNN sheet not found in {input_path}: {e}")
            deep_gnn_rank = None

        # Rename and select columns, make a copy to avoid SettingWithCopyWarning
        modified_data = data[["Pyrosetta_rank", "alphafold_rank"]].copy()
        modified_data.rename(columns={"Pyrosetta_rank": "AlphaFold Rank", "alphafold_rank": "Pyrosetta Rank"}, inplace=True)

        # Add and rename columns from other sheets
        modified_data["Fold_S Rank"] = foldx_data["alphafold rank"]
        modified_data["Foldx_Int Rank"] = foldx_int_data["alphafold rank"]
        modified_data["Hadd_em Rank"] = hadd_em_data["old_rank"]
        modified_data["Hadd_md Rank"] = hadd_md_data["old_rank"]
        if vina_rank is not None:
            modified_data["Vina Rank"] = vina_rank
        if vinardo_rank is not None:
            modified_data["Vinardo Rank"] = vinardo_rank
        if gnn_dove_rank is not None:
            modified_data["GNN_Dove Rank"] = gnn_dove_rank
        if deep_gnn_rank is not None:
            modified_data["Deep_GNN Rank"] = deep_gnn_rank

        # Save the final DataFrame to a new Excel file
        modified_data.to_excel(output_path, index=False)
    except Exception as e:
        raise Exception(f"Failed to process {input_path}: {str(e)}")

# Usage example:
# Replace the paths with the actual paths you intend to use
process_excel_files('path to/Merged_TB', 'output to/Merged/Rank_scoring_TB_vina')



# In[ ]:


import pandas as pd
import os

def process_excel_files(input_directory, output_directory):
    error_files = []  # List to store names of problematic files

    # Loop through all Excel files in the input directory
    for filename in os.listdir(input_directory):
        if filename.endswith(".xlsx"):
            # Extract PDB ID from the filename (assuming it's before the first underscore)
            pdb_id = filename.split('_')[0]
            input_path = os.path.join(input_directory, filename)
            output_path = os.path.join(output_directory, f"{pdb_id}_Ranked_TF.xlsx")

            try:
                # Process each file
                create_combined_excel(input_path, output_path)
                #print(f"Processed {filename} into {output_path}")
            except Exception as e:
                print(f"Error processing {filename}: {e}")
                error_files.append(filename)

    # Print out the problematic files
    if error_files:
        print("Files that had issues:")
        for file in error_files:
            print(file)

def create_combined_excel(input_path, output_path):
    try:
        # Load data from different sheets
        data = pd.read_excel(input_path, sheet_name='Pyrosetta', engine='openpyxl')
        foldx_data = pd.read_excel(input_path, sheet_name='Foldx', engine='openpyxl')
        foldx_int_data = pd.read_excel(input_path, sheet_name='Foldx_Int', engine='openpyxl')
        hadd_em_data = pd.read_excel(input_path, sheet_name='Hadd_em', engine='openpyxl')
        hadd_md_data = pd.read_excel(input_path, sheet_name='Hadd_md', engine='openpyxl')

        # Optional sheets
        try:
            vina_data = pd.read_excel(input_path, sheet_name='Vina', engine='openpyxl')
            vina_rank = vina_data["alphafold rank"]
        except Exception as e:
            print(f"Vina sheet not found in {input_path}: {e}")
            vina_rank = None
        try:
            vinardo_data = pd.read_excel(input_path, sheet_name='Vinardo', engine='openpyxl')
            vinardo_rank = vinardo_data["alphafold rank"]
        except Exception as e:
            print(f"Vinardo sheet not found in {input_path}: {e}")
            vinardo_rank = None
        try:
            gnn_dove_data = pd.read_excel(input_path, sheet_name='gnn_dove', engine='openpyxl')
            gnn_dove_rank = gnn_dove_data["AlphaFold Rank"]
        except Exception as e:
            print(f"GNN Dove sheet not found in {input_path}: {e}")
            gnn_dove_rank = None
        try:
            deep_gnn_data = pd.read_excel(input_path, sheet_name='Deep_GNN', engine='openpyxl')
            deep_gnn_rank = deep_gnn_data["AlphaRank"]
        except Exception as e:
            print(f"Deep GNN sheet not found in {input_path}: {e}")
            deep_gnn_rank = None

        # Rename and select columns, make a copy to avoid SettingWithCopyWarning
        modified_data = data[["Pyrosetta_rank", "alphafold_rank"]].copy()
        modified_data.rename(columns={"Pyrosetta_rank": "AlphaFold Rank", "alphafold_rank": "Pyrosetta Rank"}, inplace=True)

        # Add and rename columns from other sheets
        modified_data["Fold_S Rank"] = foldx_data["alphafold rank"]
        modified_data["Foldx_Int Rank"] = foldx_int_data["alphafold rank"]
        modified_data["Hadd_em Rank"] = hadd_em_data["old_rank"]
        modified_data["Hadd_md Rank"] = hadd_md_data["old_rank"]
        if vina_rank is not None:
            modified_data["Vina Rank"] = vina_rank
        if vinardo_rank is not None:
            modified_data["Vinardo Rank"] = vinardo_rank
        if gnn_dove_rank is not None:
            modified_data["GNN_Dove Rank"] = gnn_dove_rank
        if deep_gnn_rank is not None:
            modified_data["Deep_GNN Rank"] = deep_gnn_rank

        # Save the final DataFrame to a new Excel file
        modified_data.to_excel(output_path, index=False)
    except Exception as e:
        raise Exception(f"Failed to process {input_path}: {str(e)}")

# Usage example:
# Replace the paths with the actual paths you intend to use
process_excel_files('path to/Merged/Merged_TF', 'output to/Merged/Rank_scoring_TF_vina')



# In[ ]:


import os
import pandas as pd

# Function to process and modify Excel files
def add_tb_to_numbers(input_dir, output_dir):
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Iterate over all files in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith(".xlsx") or filename.endswith(".xls"):
            file_path = os.path.join(input_dir, filename)
            df = pd.read_excel(file_path, sheet_name=None)

            # Add '_TB' to the end of all numbers in the DataFrame
            for sheet_name, sheet_data in df.items():
                df[sheet_name] = sheet_data.applymap(lambda x: f"{x}_TB" if isinstance(x, (int, float)) else x)

            # Save the modified DataFrame to a new Excel file in the output directory
            modified_file_path = os.path.join(output_dir, filename)
            with pd.ExcelWriter(modified_file_path) as writer:
                for sheet_name, sheet_data in df.items():
                    sheet_data.to_excel(writer, sheet_name=sheet_name, index=False)

# Define your input and output directories
input_directory = 'path to/Merged/Rank_scoring_TB_vina'
output_directory = 'output to/Merged/Rank_scoring_TB_TB_md'

# Process the files
add_tb_to_numbers(input_directory, output_directory)



# In[ ]:


import os
import pandas as pd

# Function to process and modify Excel files
def add_tf_to_numbers(input_dir, output_dir):
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Iterate over all files in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith(".xlsx") or filename.endswith(".xls"):
            file_path = os.path.join(input_dir, filename)
            df = pd.read_excel(file_path, sheet_name=None)

            # Add '_TF' to the end of all numbers in the DataFrame
            for sheet_name, sheet_data in df.items():
                df[sheet_name] = sheet_data.applymap(lambda x: f"{x}_TF" if isinstance(x, (int, float)) else x)

            # Save the modified DataFrame to a new Excel file in the output directory
            modified_file_path = os.path.join(output_dir, filename)
            with pd.ExcelWriter(modified_file_path) as writer:
                for sheet_name, sheet_data in df.items():
                    sheet_data.to_excel(writer, sheet_name=sheet_name, index=False)

# Define your input and output directories
input_directory = 'path to/Merged/Rank_scoring_TF_vina'
output_directory = 'output to/Merged/Rank_scoring_TF_TF_md'

# Process the files
add_tf_to_numbers(input_directory, output_directory)


# In[ ]:


import os
import pandas as pd

def update_cell_name(cell_value, pdb_id):
    if isinstance(cell_value, str) and cell_value.endswith('_TB'):
        base_value = cell_value[:-3]
        return f"{base_value}_{pdb_id}_TB"
    return cell_value

def update_excel_files(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for filename in os.listdir(input_dir):
        if filename.endswith('.xlsx'):
            # Extract pdb_id from the filename (assuming it's the part before the first underscore or period)
            pdb_id = filename.split('_')[0]

            file_path = os.path.join(input_dir, filename)
            df = pd.read_excel(file_path)
            df_updated = df.applymap(lambda x: update_cell_name(x, pdb_id))

            output_file_path = os.path.join(output_dir, filename)
            df_updated.to_excel(output_file_path, index=False)
            #print(f"Updated file saved: {output_file_path}")

# Example usage
input_directory = 'path to/Merged/Rank_scoring_TB_TB_md'  # Replace with the path to your input directory
output_directory = 'output to/Merged/Rank_scoring_TB_TB_pdb_md'  # Replace with the path to your output directory

update_excel_files(input_directory, output_directory)



# In[ ]:


import os
import pandas as pd

def update_cell_name(cell_value, pdb_id):
    if isinstance(cell_value, str) and cell_value.endswith('_TF'):
        base_value = cell_value[:-3]
        return f"{base_value}_{pdb_id}_TF"
    return cell_value

def update_excel_files(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for filename in os.listdir(input_dir):
        if filename.endswith('.xlsx'):
            # Extract pdb_id from the filename (assuming it's the part before the first underscore or period)
            pdb_id = filename.split('_')[0]

            file_path = os.path.join(input_dir, filename)
            df = pd.read_excel(file_path)
            df_updated = df.applymap(lambda x: update_cell_name(x, pdb_id))

            output_file_path = os.path.join(output_dir, filename)
            df_updated.to_excel(output_file_path, index=False)
            #print(f"Updated file saved: {output_file_path}")

# Example usage
input_directory = 'path to/Merged/Rank_scoring_TF_TF_md'  # Replace with the path to your input directory
output_directory = 'output to/Merged/Rank_scoring_TF_TF_pdb_md'  # Replace with the path to your output directory

update_excel_files(input_directory, output_directory)



# In[ ]:


import os
import pandas as pd

# Define the directory containing the Excel files and the output directory
input_directory = 'path to/Merged/Rank_scoring_TB_TB_pdb_md'
output_directory = 'output to/Merged/Rank_scoring_TB_TB_pdb_three_md'

# Define the columns to keep
columns_to_keep = ['AlphaFold Rank', 'Fold_S Rank', 'Hadd_md Rank']

# Ensure the output directory exists
os.makedirs(output_directory, exist_ok=True)

# Process each Excel file in the input directory
for filename in os.listdir(input_directory):
    if filename.endswith('.xlsx'):
        # Construct the full file path for input and output
        input_file_path = os.path.join(input_directory, filename)
        output_file_path = os.path.join(output_directory, filename)

        # Load the Excel file
        df = pd.read_excel(input_file_path)

        # Keep only the specified columns
        filtered_df = df[columns_to_keep]

        # Save the filtered DataFrame to a new Excel file in the output directory
        filtered_df.to_excel(output_file_path, index=False)

        #print(f'Processed {filename} and saved to {output_file_path}')


# In[ ]:


import os
import pandas as pd

# Define the directory containing the Excel files and the output directory
input_directory = 'path to/Merged/Rank_scoring_TF_TF_pdb_md'
output_directory = 'output to/Merged/Rank_scoring_TF_TF_pdb_three_md'

# Define the columns to keep
columns_to_keep = ['AlphaFold Rank', 'Fold_S Rank', 'Hadd_md Rank']

# Ensure the output directory exists
os.makedirs(output_directory, exist_ok=True)

# Process each Excel file in the input directory
for filename in os.listdir(input_directory):
    if filename.endswith('.xlsx'):
        # Construct the full file path for input and output
        input_file_path = os.path.join(input_directory, filename)
        output_file_path = os.path.join(output_directory, filename)

        # Load the Excel file
        df = pd.read_excel(input_file_path)

        # Keep only the specified columns
        filtered_df = df[columns_to_keep]

        # Save the filtered DataFrame to a new Excel file in the output directory
        filtered_df.to_excel(output_file_path, index=False)

        #print(f'Processed {filename} and saved to {output_file_path}')




# In[ ]:


import os
import pandas as pd

# Directories
tb_dir = 'path to/Merged/Rank_scoring_TB_TB_pdb_three_md'
tf_dir = 'path to/Merged/Rank_scoring_TF_TF_pdb_three_md'
output_dir = 'output to/Merged/combined_TB_TF_md'

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Extract the numerical parts and pdb_id from the rank column
def extract_rank_and_pdb_id(rank):
    """Extracts the numerical part and pdb_id of the rank ignoring the '_TB' or '_TF' suffix."""
    if pd.isna(rank):
        return None, None
    parts = rank.split('_')
    try:
        # Convert the rank part to an integer (handles float conversion)
        rank_number = int(float(parts[0]))
    except ValueError:
        return None, None  # Return a flag value if conversion fails
    return rank_number, parts[1]

# Apply the extraction to specified columns and find common numbers
def process_column(tb_df, tf_df, column_name, target_count=10):
    tb_df[[f'{column_name} Number', 'pdb_id']] = tb_df[column_name].apply(lambda x: pd.Series(extract_rank_and_pdb_id(x)))
    tf_df[[f'{column_name} Number', 'pdb_id']] = tf_df[column_name].apply(lambda x: pd.Series(extract_rank_and_pdb_id(x)))

    # Check for NaN values or flag values indicating conversion issues
    tb_nan = tb_df[tb_df[[f'{column_name} Number', 'pdb_id']].isnull().any(axis=1)]
    tf_nan = tf_df[tf_df[[f'{column_name} Number', 'pdb_id']].isnull().any(axis=1)]

    tb_df.dropna(subset=[f'{column_name} Number', 'pdb_id'], inplace=True)
    tf_df.dropna(subset=[f'{column_name} Number', 'pdb_id'], inplace=True)

    def find_multiple_similar_numbers(tb_df, tf_df, target_count=10):
        window_size = 1
        common_numbers_positions = []

        while len(common_numbers_positions) < target_count:
            tb_top_ranks = tb_df.head(window_size)
            tf_top_ranks = tf_df.head(window_size)

            # Find the intersection of the two series
            common_numbers = set(tb_top_ranks[f'{column_name} Number']).intersection(set(tf_top_ranks[f'{column_name} Number']))

            for number in common_numbers:
                if all(number != x[0] for x in common_numbers_positions):
                    tb_row = tb_df[tb_df[f'{column_name} Number'] == number].iloc[0]
                    tf_row = tf_df[tf_df[f'{column_name} Number'] == number].iloc[0]
                    tb_position = tb_row.name + 1
                    tf_position = tf_row.name + 1
                    pdb_id = tb_row['pdb_id']
                    common_numbers_positions.append((number, pdb_id, tb_position, tf_position))

            window_size += 1

            # To avoid infinite loop in case there are not enough common numbers
            if window_size > max(len(tb_df), len(tf_df)):
                break

        # Sort the common numbers based on the positions (TB position + TF position)
        sorted_common_numbers_positions = sorted(common_numbers_positions, key=lambda x: x[2] + x[3])

        return sorted_common_numbers_positions, window_size

    sorted_common_numbers_positions, window_size = find_multiple_similar_numbers(tb_df, tf_df, target_count)

    # Combine the lists and sort by TB and TF positions
    combined_list = []
    for num, pdb_id, tb_pos, tf_pos in sorted_common_numbers_positions:
        combined_list.append((f"{num}_{pdb_id}_TB", tb_pos))
        combined_list.append((f"{num}_{pdb_id}_TF", tf_pos))

    # Sort by position
    sorted_combined_list = sorted(combined_list, key=lambda x: x[1])

    return pd.DataFrame(sorted_combined_list, columns=[f'{column_name}', 'Position']), tb_nan, tf_nan

# Track files with NaN values
files_with_nan = []

# Process each file in the tb_dir
for tb_file_name in os.listdir(tb_dir):
    if tb_file_name.endswith('.xlsx'):
        tb_file_path = os.path.join(tb_dir, tb_file_name)

        # Extract the PDB ID from the tb_file_name
        pdb_id = tb_file_name.split('_')[0]

        # Construct the corresponding tf_file_name
        tf_file_name = f'{pdb_id}_Ranked_TF.xlsx'
        tf_file_path = os.path.join(tf_dir, tf_file_name)

        if os.path.exists(tf_file_path):
            tb_df = pd.read_excel(tb_file_path)
            tf_df = pd.read_excel(tf_file_path)

            columns_to_process = ["AlphaFold Rank", "Fold_S Rank", "Hadd_md Rank"]
            results = []
            nan_info = {'tb_file': tb_file_name, 'tf_file': tf_file_name, 'columns': {}}

            for column in columns_to_process:
                result_df, tb_nan, tf_nan = process_column(tb_df, tf_df, column, target_count=10)
                results.append(result_df)

                if not tb_nan.empty or not tf_nan.empty:
                    nan_info['columns'][column] = {'tb_nan': tb_nan, 'tf_nan': tf_nan}

            # Save NaN information if there are any NaN values found
            if nan_info['columns']:
                files_with_nan.append(nan_info)

            # Merge all results into one DataFrame
            final_df = pd.concat(results, axis=1)

            # Remove the 'Position' columns
            position_columns = [col for col in final_df.columns if 'Position' in col]
            final_df = final_df.drop(columns=position_columns)

            # Save the final DataFrame to an Excel file
            output_file_path = os.path.join(output_dir, f'{pdb_id}_combined.xlsx')
            final_df.to_excel(output_file_path, index=False)

            #print(f"Results saved to {output_file_path}")

# Print files with NaN values
for nan_file in files_with_nan:
    print(f"File {nan_file['tb_file']} and {nan_file['tf_file']} have NaN values in the following columns:")
    for column, data in nan_file['columns'].items():
        print(f"  Column: {column}")
        print(f"    TB NaN values:\n{data['tb_nan']}")
        print(f"    TF NaN values:\n{data['tf_nan']}")


# In[ ]:


import os
import pandas as pd

# Define the directory containing the Excel files
input_directory = 'path to/Merged/combined_TB_TF_md'
output_directory = 'output to/Merged/combined_TB_TF_top20_md'

# Ensure the output directory exists
os.makedirs(output_directory, exist_ok=True)

# Process each Excel file in the input directory
for filename in os.listdir(input_directory):
    if filename.endswith('.xlsx'):
        file_path = os.path.join(input_directory, filename)

        # Load the Excel file
        excel_data = pd.ExcelFile(file_path)

        # Assume the sheet we need is the first one
        sheet_name = excel_data.sheet_names[0]

        # Load the data from the sheet
        df = pd.read_excel(file_path, sheet_name=sheet_name)

        # Keep only the top 10 rows
        df_top_10 = df.head(20)

        # Save the modified dataframe to a new Excel file in the output directory
        output_path = os.path.join(output_directory, filename)
        df_top_10.to_excel(output_path, index=False)

print("Processing complete. Files saved to", output_directory)



# In[ ]:


import os
import pandas as pd
from matplotlib_venn import venn3
import matplotlib.pyplot as plt

# Define the path to the directory containing the Excel files
directory_path = 'path to/Merged/combined_TB_TF_top20_md'

# Define the scoring function column names
scoring_functions = ['AlphaFold Rank', 'Fold_S Rank', 'Hadd_md Rank']

# Initialize a dictionary to store models for each scoring function
all_models = {func: set() for func in scoring_functions}

# Function to add models to respective sets
def merge_models(df):
    for func in scoring_functions:
        all_models[func].update(df[func].dropna().unique())

# Iterate over all Excel files in the directory
for filename in os.listdir(directory_path):
    if filename.endswith('.xlsx'):
        file_path = os.path.join(directory_path, filename)
        df = pd.read_excel(file_path)
        merge_models(df)

# Calculate common models among all three scoring functions
common_all_three = set.intersection(*[all_models[func] for func in scoring_functions])

# Plot the Venn diagram
plt.figure(figsize=(8, 8))
venn3([all_models[func] for func in scoring_functions], set_labels=scoring_functions)
plt.title('Common Models Among All Scoring Functions')
plt.show()

# Print the common models
print("Common models among all three scoring functions:")
print(common_all_three)




# In[ ]:


import pandas as pd
import os

# Load the common models file
common_models_path = 'path to/Merged/Output_TB_TF/common_models_md_20.xlsx'
common_models_df = pd.read_excel(common_models_path)

# Print the first few rows of the DataFrame to check its structure
print("Common models DataFrame head:")
print(common_models_df.head())

# Check the number of columns in the DataFrame
print("Number of columns in the common models DataFrame:", len(common_models_df.columns))

# Add a new column for Combined DockQ at index 1
common_models_df.insert(1, 'Combined DockQ', '')

# Define the directories for TB and TF files
tb_directory = 'path to /AlphaFold_Multimer_All_TB'
tf_directory = 'path to/AlphaFold_Multimer_All_TF'

# Function to get DockQ value from the corresponding file
def get_dockq_value(model_name):
    parts = model_name.split('_')
    first_number = parts[0]
    pdb_id = parts[1]
    file_type = parts[-1]

    directory = tb_directory if file_type == 'TB' else tf_directory
    file_path = os.path.join(directory, f'{pdb_id}_DockQ_data_{file_type}.xlsx')

    if os.path.exists(file_path):
        dockq_df = pd.read_excel(file_path)
        file_name = f'ranked_{first_number}.pdb_clean'
        dockq_value = dockq_df.loc[dockq_df['File Name'] == file_name, 'DockQ'].values
        if len(dockq_value) > 0:
            return dockq_value[0]
    return None

# Extract DockQ values for each model in the common models file
common_models_df['Combined DockQ'] = common_models_df['File name'].apply(get_dockq_value)

# Save the updated DataFrame to a new Excel file
output_file_path = 'output to/Merged/Output_TB_TF/updated_common_models_md_20.xlsx'
common_models_df.to_excel(output_file_path, index=False)

print(f"Updated common models file has been saved to {output_file_path}")



