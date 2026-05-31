#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import os

# Define the paths to the directories containing the Excel files
directory_tb = 'path to/sum_files_TB'  # Path to the directory with 7PRX_TB.xlsx files
directory_dockq = 'path to/AlphaFold_Multimer_All_TB'  # Path to the directory with 7PRX_DockQ_data_TB.xlsx files

# List all Excel files in both directories
files_tb = [f for f in os.listdir(directory_tb) if f.endswith('.xlsx')]
files_dockq = [f for f in os.listdir(directory_dockq) if f.endswith('.xlsx')]

# Create an empty DataFrame to store the final results
final_data = pd.DataFrame(columns=['File Name', 'DockQ'])

# Process each file in the tb directory
for file_tb in files_tb:
    pdb_id = file_tb.split('_')[0]  # Extract pdb_id from file name like '7PRX_TB.xlsx'
    corresponding_file = f"{pdb_id}_DockQ_data_TB.xlsx"  # Construct the corresponding DockQ file name

    # Check if the corresponding DockQ file exists in the dockq directory
    if corresponding_file in files_dockq:
        # Load data from both files
        data_tb = pd.read_excel(os.path.join(directory_tb, file_tb))
        data_dockq = pd.read_excel(os.path.join(directory_dockq, corresponding_file))

        # Get the file number from the first file and construct the new file name
        file_number = data_tb.iloc[0]['File Name'].split('.')[0]
        new_file_name = f"{file_number}_{pdb_id}_TB"

        # Filter the second DataFrame for the new file name and extract the DockQ value
        filtered_row = data_dockq[data_dockq['File Name'].str.contains(f"ranked_{file_number}.pdb_clean", na=False)]
        dockq_value = filtered_row['DockQ'].values[0] if not filtered_row.empty else None

        # Append the results to the final DataFrame
        final_data = pd.concat([final_data, pd.DataFrame({'File Name': [new_file_name], 'DockQ': [dockq_value]})], ignore_index=True)

# Save the final DataFrame to an Excel file
output_file_path = 'output to /dockq_values_TB.xlsx'
final_data.to_excel(output_file_path, index=False)

print(f"Final data saved to {output_file_path}")



# In[ ]:


import pandas as pd

# Load the Excel files
updated_common_models = pd.read_excel('path to/dockq_values_TB.xlsx')
dockq_tb = pd.read_excel('path to/DockQ_TB.xlsx')
dockq_tf = pd.read_excel('path to/DockQ_TF.xlsx')

# Function to fetch the DockQ score based on the model type and pdb_id
def fetch_dockq(row):
    pdb_id, model_type = row['File Name'].rsplit('_', 1)
    # Extract only the pdb_id without the prefix
    pdb_id = pdb_id.split('_')[1]
    # Determine the appropriate DockQ file based on the type
    if model_type == 'TB':
        dockq_value = dockq_tb[dockq_tb['File Name'].str.upper() == pdb_id]['DockQ']
    else:
        dockq_value = dockq_tf[dockq_tf['File Name'].str.upper() == pdb_id]['DockQ']
    # Return the first value found
    return dockq_value.iloc[0] if not dockq_value.empty else None

# Apply the function to each row in the dataframe
updated_common_models['DockQ'] = updated_common_models.apply(fetch_dockq, axis=1)

# Calculate the Loss value as the difference between DockQ and Combined DockQ
updated_common_models['Loss'] = updated_common_models['DockQ'] - updated_common_models['Combined DockQ']

# Reorder the columns including the new Loss column
updated_common_models = updated_common_models[['File Name', 'DockQ', 'Combined DockQ', 'Loss']]

# Save the updated DataFrame to an Excel file
output_path = 'output to/dockq_values_TB_dockq.xlsx'
updated_common_models.to_excel(output_path, index=False)



# In[ ]:


import pandas as pd
import os

# Define the paths to the directories containing the Excel files
directory_tb = 'path to/sum_files_TF'  # Path to the directory with 7PRX_TB.xlsx files
directory_dockq = 'path to/Data/AlphaFold_Multimer_All_TF'  # Path to the directory with 7PRX_DockQ_data_TB.xlsx files

# List all Excel files in both directories
files_tb = [f for f in os.listdir(directory_tb) if f.endswith('.xlsx')]
files_dockq = [f for f in os.listdir(directory_dockq) if f.endswith('.xlsx')]

# Create an empty DataFrame to store the final results
final_data = pd.DataFrame(columns=['File Name', 'DockQ'])

# Process each file in the tb directory
for file_tb in files_tb:
    pdb_id = file_tb.split('_')[0]  # Extract pdb_id from file name like '7PRX_TB.xlsx'
    corresponding_file = f"{pdb_id}_DockQ_data_TF.xlsx"  # Construct the corresponding DockQ file name

    # Check if the corresponding DockQ file exists in the dockq directory
    if corresponding_file in files_dockq:
        # Load data from both files
        data_tb = pd.read_excel(os.path.join(directory_tb, file_tb))
        data_dockq = pd.read_excel(os.path.join(directory_dockq, corresponding_file))

        # Get the file number from the first file and construct the new file name
        file_number = data_tb.iloc[0]['File Name'].split('.')[0]
        new_file_name = f"{file_number}_{pdb_id}_TF"

        # Filter the second DataFrame for the new file name and extract the DockQ value
        filtered_row = data_dockq[data_dockq['File Name'].str.contains(f"ranked_{file_number}.pdb_clean", na=False)]
        dockq_value = filtered_row['DockQ'].values[0] if not filtered_row.empty else None

        # Append the results to the final DataFrame
        final_data = pd.concat([final_data, pd.DataFrame({'File Name': [new_file_name], 'DockQ': [dockq_value]})], ignore_index=True)

# Save the final DataFrame to an Excel file
output_file_path = 'output to/dockq_values_TF.xlsx'
final_data.to_excel(output_file_path, index=False)

print(f"Final data saved to {output_file_path}")


# In[ ]:


import pandas as pd

# Load the Excel files
updated_common_models = pd.read_excel('path to/dockq_values_TF.xlsx')
dockq_tb = pd.read_excel('path to/Output_TB_TF/DockQ_TB.xlsx')
dockq_tf = pd.read_excel('output to/Output_TB_TF/DockQ_TF.xlsx')

# Function to fetch the DockQ score based on the model type and pdb_id
def fetch_dockq(row):
    pdb_id, model_type = row['File Name'].rsplit('_', 1)
    # Extract only the pdb_id without the prefix
    pdb_id = pdb_id.split('_')[1]
    # Determine the appropriate DockQ file based on the type
    if model_type == 'TB':
        dockq_value = dockq_tb[dockq_tb['File Name'].str.upper() == pdb_id]['DockQ']
    else:
        dockq_value = dockq_tf[dockq_tf['File Name'].str.upper() == pdb_id]['DockQ']
    # Return the first value found
    return dockq_value.iloc[0] if not dockq_value.empty else None

# Apply the function to each row in the dataframe
updated_common_models['DockQ'] = updated_common_models.apply(fetch_dockq, axis=1)

# Calculate the Loss value as the difference between DockQ and Combined DockQ
updated_common_models['Loss'] = updated_common_models['DockQ'] - updated_common_models['Combined DockQ']

# Reorder the columns including the new Loss column
updated_common_models = updated_common_models[['File Name', 'DockQ', 'Combined DockQ', 'Loss']]

# Save the updated DataFrame to an Excel file
output_path = 'output to/dockq_values_TF_new.xlsx'
updated_common_models.to_excel(output_path, index=False)

