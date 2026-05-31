#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import json

# Load configuration and define a global variable
def load_config():
    global config
    with open('config_DockQ.json') as config_file:
        config = json.load(config_file)

# Call the function to load the configuration
load_config()


# In[ ]:


import pandas as pd
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Load configuration and define a global variable
def load_config():
    global config
    with open('config_DockQ.json') as config_file:
        config = json.load(config_file)

# Call the function to load the configuration
load_config()

# Directory where your Excel files are stored
input_directory = config['AlphaFold_directory']
output_directory = config['Spearman_Correlation_directory_TF']  # Specify the output directory here

# Initialize a list to store the results
results = []

# Loop through each Excel file in the directory
for filename in os.listdir(input_directory):
    if filename.endswith('.xlsx') and not filename.startswith('~$'):  # Skip temporary files
        file_path = os.path.join(input_directory, filename)
        try:
            # Load the data
            data = pd.read_excel(file_path)
            # Calculate Spearman's rank correlation
            spearman_corr, _ = spearmanr(data['AlphaFold Rank'], data['DockQ Rank'])
            # Append the results
            results.append((filename, spearman_corr))
            # Plotting
            sns.scatterplot(x='DockQ Rank', y='AlphaFold Rank', data=data)
            plt.title(f"{filename} - Spearman's Correlation: {spearman_corr:.3f}")
            plt.xlabel('DockQ Rank')
            plt.ylabel('AlphaFold Rank')
            plt.show()
        except Exception as e:
            print(f"Error processing file {filename}: {e}")

# Save the Spearman's correlation coefficients to a new Excel file
output_file_path = os.path.join(output_directory, "spearman_correlation_AF_TF.xlsx")
results_df = pd.DataFrame(results, columns=['File Name', 'Spearman Correlation'])
results_df.to_excel(output_file_path, index=False)

print(f"Process completed. Spearman correlation results saved to '{output_file_path}'.")


# In[ ]:


import pandas as pd

# Load the Excel file
# Make sure to include the file extension, such as '.xlsx', in the file path
file_path = f"{config['Spearman_Correlation_directory_TF']}/spearman_correlation_AF_TF.xlsx" # Adjusted file path with extension
try:
    df = pd.read_excel(file_path)

    # Calculate the percentages of positive and negative Spearman correlation values
    positive_percentage = (df['Spearman Correlation'] > 0).mean() * 100
    negative_percentage = (df['Spearman Correlation'] < 0).mean() * 100

    # Calculate the average of Spearman Correlation values
    average_correlation = df['Spearman Correlation'].mean()

    # Print the results
    print(f"Percentage of positive correlations: {positive_percentage:.2f}%")
    print(f"Percentage of negative correlations: {negative_percentage:.2f}%")
    print(f"Average Spearman Correlation: {average_correlation:.2f}")
except FileNotFoundError:
    print(f"File not found. Please check the file path: {file_path}")


# In[ ]:


import os
import pandas as pd


def process_excel_files(input_directory, output_directory):
    # Define the output DataFrame
    output_df = pd.DataFrame(columns=['File Name', 'DockQ', 'Alphafold ranked', 'Loss'])

    # Iterate over each file in the input directory
    for file in os.listdir(input_directory):
        # Skip temporary or hidden files
        if file.endswith(".xlsx") and not file.startswith("~$"):
            # Extract the pdb_id from the file name
            pdb_id = file.split('_')[0]

            # Load the Excel file with error handling
            file_path = os.path.join(input_directory, file)
            try:
                data = pd.read_excel(file_path, engine='openpyxl')
            except Exception as e:
                print(f"Error reading {file}: {e}")
                continue

            # Find the highest DockQ score
            highest_dockq_score = data['DockQ'].max()

            # Find the DockQ score for "ranked_0.pdb_clean"
            dockq_score_for_ranked_0 = data.loc[data['File Name'] == 'ranked_0.pdb_clean', 'DockQ'].values[0] if not data.loc[data['File Name'] == 'ranked_0.pdb_clean', 'DockQ'].empty else None

            # Calculate the Loss
            Loss = highest_dockq_score - dockq_score_for_ranked_0 if dockq_score_for_ranked_0 is not None else None

            # Append the results to the output DataFrame using pd.concat
            new_row = pd.DataFrame({'File Name': [pdb_id], 'DockQ': [highest_dockq_score], 'Alphafold ranked': [dockq_score_for_ranked_0], 'Loss': [Loss]})
            output_df = pd.concat([output_df, new_row], ignore_index=True)

    # Save the output DataFrame to an Excel file in the output directory
    output_file_path = os.path.join(output_directory, 'AlphaFold_TF.xlsx')
    output_df.to_excel(output_file_path, index=False)

# Example usage
input_directory = config["AlphaFold_directory"]
output_directory = config["DockQ_Loss_directory_TF"]
process_excel_files(input_directory, output_directory)


# In[ ]:


########################### All Models DockQ #############################

import os
import pandas as pd
import matplotlib.pyplot as plt

def categorize_quality(dockq_value):
    if dockq_value > 0.8:
        return 'High quality'
    elif 0.5 <= dockq_value <= 0.8:
        return 'Medium quality'
    elif 0.2 <= dockq_value < 0.5:
        return 'Acceptable quality'
    elif dockq_value < 0.2:
        return 'Incorrect'
    else:
        return 'Undefined'

directory_path = config['AlphaFold_directory']

quality_counts = {
    'High quality': 0,
    'Medium quality': 0,
    'Acceptable quality': 0,
    'Incorrect': 0
}

total_count = 0

for file_name in os.listdir(directory_path):
    if file_name.endswith('.xlsx') or file_name.endswith('.xls'):
        if file_name.startswith('~$'):
            continue

        file_path = os.path.join(directory_path, file_name)

        try:
            df = pd.read_excel(file_path)
            if df.empty:
                continue
        except Exception as e:
            continue

        if 'DockQ' not in df.columns:
            continue

        df['Quality Category'] = df['DockQ'].apply(categorize_quality)

        # Summarize counts for each category
        quality_counts_update = df['Quality Category'].value_counts()
        for quality, count in quality_counts_update.items():
            quality_counts[quality] += count
        total_count += len(df)

# Calculate overall percentages for each quality level
if total_count > 0:
    quality_percentages = {quality: (count / total_count) * 100 for quality, count in quality_counts.items()}

    # Data to plot
    labels = list(quality_percentages.keys())
    sizes = list(quality_percentages.values())
    # Update the colors as per your choice
    colors = ['blue', 'orange', 'green', 'red']
    explode = (0.1, 0, 0, 0)  # explode the first slice (High quality)

    # Plot
    plt.figure(figsize=(8, 8))
    plt.pie(sizes, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', 
            shadow=True, startangle=140, textprops={'fontsize': 16})
    plt.title('Overall DockQ Quality Distribution', fontsize=20)
    plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.show()
else:
    print("No data to plot.")


# In[ ]:


########################### FIRST RANKED DockQ #############################

import pandas as pd

def categorize_quality(dockq_value):
    if dockq_value > 0.8:
        return 'High quality'
    elif 0.5 <= dockq_value <= 0.8:
        return 'Medium quality'
    elif 0.2 <= dockq_value < 0.5:
        return 'Acceptable quality'
    elif dockq_value < 0.2:
        return 'Incorrect'
    else:
        return 'Undefined'

# Load the Excel file
file_path = f"{config['DockQ_Loss_directory_TF']}/AlphaFold_TF.xlsx"
data = pd.read_excel(file_path)

# Apply the categorization function to the DockQ column
data['AlphaFold-Quality Category'] = data['Alphafold ranked'].apply(categorize_quality)

# Calculate the percentage of each category
category_counts = data['AlphaFold-Quality Category'].value_counts(normalize=True) * 100

# Print the results
print(category_counts)



# ColabFold Analysis 
# 
# For Template_Based and Template_Free models, with analysis their First-Ranked predicted structures.

# In[ ]:


import os
import pandas as pd
import json

# Load configuration and define a global variable
def load_config():
    global config
    with open('config_DockQ.json') as config_file:
        config = json.load(config_file)

# Call the function to load the configuration
load_config()

def process_excel_files(input_directory, output_directory):
    # Define the output DataFrame
    output_df = pd.DataFrame(columns=['File Name', 'DockQ', 'ColabFold ranked', 'Loss'])

    # Iterate over each file in the input directory
    for file in os.listdir(input_directory):
        if file.endswith(".xlsx") and not file.startswith('~$'):  # Skip temporary files
            # Extract the pdb_id from the file name
            pdb_id = file.split('_')[0]

            # Load the Excel file
            file_path = os.path.join(input_directory, file)
            try:
                data = pd.read_excel(file_path, engine='openpyxl')

                # Find the highest DockQ score
                highest_dockq_score = data['DockQ'].max()

                # Find the DockQ score for files ending with "_rank_1.pdb_clean"
                rank_1_row = data[data['File Name'].str.endswith("_rank_1.pdb_clean")]
                dockq_score_for_rank_1 = rank_1_row['DockQ'].values[0] if not rank_1_row.empty else None

                # Calculate the Loss
                Loss = highest_dockq_score - dockq_score_for_rank_1 if dockq_score_for_rank_1 is not None else None

                # Append the results to the output DataFrame
                output_df = pd.concat([output_df, pd.DataFrame([{
                    'File Name': pdb_id,
                    'DockQ': highest_dockq_score,
                    'ColabFold ranked': dockq_score_for_rank_1,
                    'Loss': Loss
                }])], ignore_index=True)

            except Exception as e:
                print(f"Error processing file {file}: {e}")

    # Save the output DataFrame to an Excel file in the output directory
    output_file_path = os.path.join(output_directory, 'ColabFold_TF.xlsx')
    output_df.to_excel(output_file_path, index=False, engine='openpyxl')
    print(f"All data extracted and written to {output_file_path}")

# Example usage
input_directory = config["ColabFold_directory"]  # Update the path as per your requirement
output_directory = config["DockQ_Loss_directory_TF"]  # Update the path as per your requirement

# Ensure the output directory exists
os.makedirs(output_directory, exist_ok=True)

process_excel_files(input_directory, output_directory)


# In[ ]:


import os
import pandas as pd
import matplotlib.pyplot as plt

def categorize_quality(dockq_value):
    if dockq_value > 0.8:
        return 'High quality'
    elif 0.5 <= dockq_value <= 0.8:
        return 'Medium quality'
    elif 0.2 <= dockq_value < 0.5:
        return 'Acceptable quality'
    elif dockq_value < 0.2:
        return 'Incorrect'
    else:
        return 'Undefined'

directory_path = config['ColabFold_directory']

quality_counts = {
    'High quality': 0,
    'Medium quality': 0,
    'Acceptable quality': 0,
    'Incorrect': 0
}

total_count = 0

for file_name in os.listdir(directory_path):
    if file_name.endswith('.xlsx') or file_name.endswith('.xls'):
        if file_name.startswith('~$'):
            continue

        file_path = os.path.join(directory_path, file_name)

        try:
            df = pd.read_excel(file_path)
            if df.empty:
                continue
        except Exception as e:
            continue

        if 'DockQ' not in df.columns:
            continue

        df['Quality Category'] = df['DockQ'].apply(categorize_quality)

        # Summarize counts for each category
        quality_counts_update = df['Quality Category'].value_counts()
        for quality, count in quality_counts_update.items():
            quality_counts[quality] += count
        total_count += len(df)

# Calculate overall percentages for each quality level
if total_count > 0:
    quality_percentages = {quality: (count / total_count) * 100 for quality, count in quality_counts.items()}

    # Data to plot
    labels = list(quality_percentages.keys())
    sizes = list(quality_percentages.values())
    # Update the colors as per your choice
    colors = ['blue', 'orange', 'green', 'red']
    explode = (0.1, 0, 0, 0)  # explode the first slice (High quality)

    # Plot
    plt.figure(figsize=(8, 8))
    plt.pie(sizes, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', 
            shadow=True, startangle=140, textprops={'fontsize': 16})
    plt.title('Overall DockQ Quality Distribution', fontsize=20)
    plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.show()
else:
    print("No data to plot.")


# In[ ]:


import pandas as pd

def categorize_quality(dockq_value):
    if dockq_value > 0.8:
        return 'High quality'
    elif 0.5 <= dockq_value <= 0.8:
        return 'Medium quality'
    elif 0.2 <= dockq_value < 0.5:
        return 'Acceptable quality'
    elif dockq_value < 0.2:
        return 'Incorrect'
    else:
        return 'Undefined'

# Load the Excel file
file_path = f"{config['DockQ_Loss_directory_TF']}/ColabFold_TF.xlsx"  # Update this to your file path
data = pd.read_excel(file_path)

# Apply the categorization function to the DockQ column
data['ColabFold-Quality Category'] = data['ColabFold ranked'].apply(categorize_quality)

# Calculate the percentage of each category
category_counts = data['ColabFold-Quality Category'].value_counts(normalize=True) * 100

# Print the results
print(category_counts)


# In[ ]:


import os
import pandas as pd
import matplotlib.pyplot as plt

def categorize_quality(dockq_value):
    if dockq_value > 0.8:
        return 'High quality'
    elif 0.5 <= dockq_value <= 0.8:
        return 'Medium quality'
    elif 0.2 <= dockq_value < 0.5:
        return 'Acceptable quality'
    elif dockq_value < 0.2:
        return 'Incorrect'
    else:
        return 'Undefined'

directory_path = config['AlphaFold3_directory']

quality_counts = {
    'High quality': 0,
    'Medium quality': 0,
    'Acceptable quality': 0,
    'Incorrect': 0
}

total_count = 0

for file_name in os.listdir(directory_path):
    if file_name.endswith('.xlsx') or file_name.endswith('.xls'):
        if file_name.startswith('~$'):
            continue

        file_path = os.path.join(directory_path, file_name)

        try:
            df = pd.read_excel(file_path)
            if df.empty:
                continue
        except Exception as e:
            continue

        if 'DockQ' not in df.columns:
            continue

        df['Quality Category'] = df['DockQ'].apply(categorize_quality)

        # Summarize counts for each category
        quality_counts_update = df['Quality Category'].value_counts()
        for quality, count in quality_counts_update.items():
            quality_counts[quality] += count
        total_count += len(df)

# Calculate overall percentages for each quality level
if total_count > 0:
    quality_percentages = {quality: (count / total_count) * 100 for quality, count in quality_counts.items()}

    # Data to plot
    labels = list(quality_percentages.keys())
    sizes = list(quality_percentages.values())
    # Update the colors as per your choice
    colors = ['blue', 'orange', 'green', 'red']
    explode = (0.1, 0, 0, 0)  # explode the first slice (High quality)

    # Plot
    plt.figure(figsize=(8, 8))
    plt.pie(sizes, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', 
            shadow=True, startangle=140, textprops={'fontsize': 16})
    plt.title('Overall DockQ Quality Distribution', fontsize=20)
    plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    plt.show()
else:
    print("No data to plot.")


# In[ ]:


import os
import pandas as pd

def process_excel_files(input_directory, output_directory):
    # Define the output DataFrame
    output_df = pd.DataFrame(columns=['File Name', 'DockQ', 'AF3 ranked', 'Loss'])

    # Iterate over each file in the input directory
    for file in os.listdir(input_directory):
        if file.endswith(".xlsx"):
            file_path = os.path.join(input_directory, file)
            try:
                # Load the Excel file
                data = pd.read_excel(file_path, engine='openpyxl')

                # Find the highest DockQ score
                highest_dockq_score = data['DockQ'].max()

                # Find the DockQ score for files ending with "_model_0.pdb_clean"
                rank_1_row = data[data['File Name'].str.endswith("_model_0.pdb_clean")]
                dockq_score_for_rank_1 = rank_1_row['DockQ'].values[0] if not rank_1_row.empty else None

                # Calculate the Loss
                loss = highest_dockq_score - dockq_score_for_rank_1 if dockq_score_for_rank_1 is not None else None

                # Append the results to the output DataFrame using pd.concat
                new_row = pd.DataFrame([{
                    'File Name': file.split('_')[0],
                    'DockQ': highest_dockq_score,
                    'AF3 ranked': dockq_score_for_rank_1,
                    'Loss': loss
                }])
                output_df = pd.concat([output_df, new_row], ignore_index=True)

            except Exception as e:
                print(f"Failed to process {file}: {e}")

    # Save the output DataFrame to an Excel file in the output directory
    output_file_path = os.path.join(output_directory, 'AF3_Loss.xlsx')
    output_df.to_excel(output_file_path, index=False, engine='openpyxl')
    print(f"All data extracted and written to {output_file_path}")

# Update the paths as per your requirement
input_directory = config['AlphaFold3_directory']
output_directory = config['DockQ_Loss_directory_TB']
process_excel_files(input_directory, output_directory)


# In[ ]:


import pandas as pd

def categorize_quality(dockq_value):
    if dockq_value > 0.8:
        return 'High quality'
    elif 0.5 <= dockq_value <= 0.8:
        return 'Medium quality'
    elif 0.2 <= dockq_value < 0.5:
        return 'Acceptable quality'
    elif dockq_value < 0.2:
        return 'Incorrect'
    else:
        return 'Undefined'

# Load the Excel file
file_path = f"{config['DockQ_Loss_directory_TB']}/AF3_Loss.xlsx"  # Update this to your file path
data = pd.read_excel(file_path)

# Apply the categorization function to the DockQ column
data['AlphaFold3-Quality Category'] = data['AF3 ranked'].apply(categorize_quality)

# Calculate the percentage of each category
category_counts = data['AlphaFold3-Quality Category'].value_counts(normalize=True) * 100

# Print the results
print(category_counts)



# In[ ]:


import matplotlib.pyplot as plt
import numpy as np
import json

# Extract the directory path from the config
figures_directory = config["Main_Figures_directory"]

# Data for TB and TF of AlphaFold and ColabFold
categories = ['High', 'Medium', 'Acceptable', 'Incorrect']
# AlphaFold values
af3_values = [34.9, 25.4, 18.3, 21.4]
af_tb_values = [26.5, 30.1, 25.3, 18.1]
af_tf_values = [22.5, 35.2, 24.9, 17.4]
# ColabFold values
cf_tb_values = [21.3, 30.6, 27.7, 20.3]
cf_tf_values = [20.3, 30.3, 27.3, 22]

# Positions of the bar-groups
barWidth = 0.15  # Narrower bar width to fit all bars
r1 = np.arange(len(categories))  # Positions for the first set of bars
r2 = [x + barWidth for x in r1]  # Positions for the second set
r3 = [x + barWidth for x in r2]  # Positions for the third set
r4 = [x + barWidth for x in r3]  # Positions for the fourth set
r5 = [x + barWidth for x in r4]  # Positions for the fifth set

plt.figure(figsize=(10, 8), facecolor='white')  # Set figure background to white
plt.bar(r1, af3_values, color='#FF8E44', width=barWidth, edgecolor='grey', label='AF3')
plt.bar(r2, af_tb_values, color='#006400', width=barWidth, edgecolor='grey', label='AFM-TB')
plt.bar(r3, af_tf_values, color='#b2d8b2', width=barWidth, edgecolor='grey', label='AFM-TF')  # Very light green
plt.bar(r4, cf_tb_values, color='darkblue', width=barWidth, edgecolor='grey', label='CF-TB')
plt.bar(r5, cf_tf_values, color='#bfbfff', width=barWidth, edgecolor='grey', label='CF-TF')

# Add xticks on the middle of the group bars, adjust font sizes
plt.xticks([r + 1.5*barWidth for r in range(len(categories))], categories, rotation=0, fontsize=22)
plt.ylabel('Percentage of DockQ scores (%)', fontsize=23)  # Larger font size for y-axis label
plt.xlabel('DockQ score categories', fontsize=23)

# Adjust y-axis font size
plt.tick_params(axis='y', labelsize=18)  # Adjusted y-axis font size

# Create legend & Show graphic with larger font
plt.title('AFM, CF, and AF3 DockQ scores (all models)', fontsize=23, fontweight='bold', pad=20)  # Added pad parameter


# Create legend to display the color code mapping
plt.legend(fontsize=18, loc='upper right')  # Adjust location and font size as needed


# Construct the output path
output_filename = 'Fig_2a.jpeg'
output_path = f"{figures_directory}/{output_filename}"

# Save the figure with a white background
plt.savefig(output_path, format='png', dpi=1000, facecolor='white')

plt.show()


# In[ ]:


import matplotlib.pyplot as plt
import numpy as np
import json

# Extract the directory path from the config
figures_directory = config["Main_Figures_directory"]

# Data for TB and TF of AlphaFold and ColabFold
categories = ['High', 'Medium', 'Acceptable', 'Incorrect']
# AlphaFold values
af3_TF_values = [34.6, 24.8, 20.4, 20.3]
af3_TB_values = [34.9, 28.9, 16.5, 19.7]
af_tb_values = [26.5, 30.1, 25.3, 18.1]
af_tf_values = [22.5, 35.2, 24.9, 17.4]

# Positions of the bar-groups
barWidth = 0.15  # Narrower bar width to fit all bars
r1 = np.arange(len(categories))  # Positions for AF3_TB
r2 = [x + barWidth for x in r1]  # Positions for AF3_TF
r3 = [x + barWidth for x in r2]  # Positions for AFM-TB
r4 = [x + barWidth for x in r3]  # Positions for AFM-TF

plt.figure(figsize=(10, 8), facecolor='white')  # Set figure background to white

# Plot AF3_TB (darker orange) and AF3_TF (lighter orange) side by side
# AF3 bars (TB brown, TF yellow)
plt.bar(r1, af3_TB_values, color='#A05200', width=barWidth, edgecolor='grey', label='AF3-L-TB')  # brown
plt.bar(r2, af3_TF_values, color='#FDD835', width=barWidth, edgecolor='grey', label='AF3-L-TF')  # yellow
# Plot the rest
plt.bar(r3, af_tb_values, color='#006400', width=barWidth, edgecolor='grey', label='AFM-TB')
plt.bar(r4, af_tf_values, color='#b2d8b2', width=barWidth, edgecolor='grey', label='AFM-TF')

# X-ticks in the center of the group
plt.xticks([r + 1.5*barWidth for r in range(len(categories))], categories, rotation=0, fontsize=22)
plt.ylabel('Percentage of DockQ scores (%)', fontsize=23)
plt.xlabel('DockQ score categories', fontsize=23)
plt.tick_params(axis='y', labelsize=18)
plt.title('AFM and AF3-L DockQ scores (1000 models)', fontsize=23, fontweight='bold', pad=20)
plt.legend(fontsize=18, loc='upper right')
output_filename = 'Fig_2a_AF3.jpeg'
output_path = f"{figures_directory}/{output_filename}"
plt.savefig(output_path, format='png', dpi=1000, facecolor='white')
plt.show()


# In[ ]:


import matplotlib.pyplot as plt
import numpy as np

figures_directory = config["Main_Figures_directory"]

# Data for the first-ranked TB and TF of AlphaFold_Multimer and ColabFold
categories = ['High', 'Medium', 'Acceptable', 'Incorrect']
# AlphaFold values
af3_values = [35.59, 25.42, 18.64, 20.33]
af_tb_values = [31.66, 33.33, 20.00, 15.00]
af_tf_values = [33.33, 30.00, 20.00, 16.66]
# ColabFold values
cf_tb_values = [23.33, 33.3, 30.0, 13.3]
cf_tf_values = [21.6, 33.3, 28.3, 16.6]

# Positions of the bar-groups
barWidth = 0.15  # Narrower bar width to fit all bars
r1 = np.arange(len(categories))  # Positions for the first set of bars
r2 = [x + barWidth for x in r1]  # Positions for the second set
r3 = [x + barWidth for x in r2]  # Positions for the third set
r4 = [x + barWidth for x in r3]  # Positions for the fourth set
r5 = [x + barWidth for x in r4]  # Positions for the fifth set

# Create the plot with adjusted figure size
plt.figure(figsize=(10, 8), facecolor='white')  # Set figure background to white
plt.bar(r1, af3_values, color='#FF8E44', width=barWidth, edgecolor='grey', label='AF3')
plt.bar(r2, af_tb_values, color='#006400', width=barWidth, edgecolor='grey', label='AFM-TB')
plt.bar(r3, af_tf_values, color='#b2d8b2', width=barWidth, edgecolor='grey', label='AFM-TF')
plt.bar(r4, cf_tb_values, color='darkblue', width=barWidth, edgecolor='grey', label='CF-TB')
plt.bar(r5, cf_tf_values, color='#bfbfff', width=barWidth, edgecolor='grey', label='CF-TF')

# Add xticks on the middle of the group bars, adjust font sizes
plt.xticks([r + 1.5*barWidth for r in range(len(categories))], categories, rotation=0, fontsize=22)
plt.ylabel('Percentage of DockQ scores (%)', fontsize=23)  # Larger font size for y-axis label
plt.xlabel('DockQ score categories', fontsize=23)

# Adjust y-axis font size
plt.tick_params(axis='y', labelsize=20)  # Adjusted y-axis font size

# Create legend to display the color code mapping
plt.legend(fontsize=18, loc='upper right')  # Adjust location and font size as needed


# Create legend & Show graphic with larger font
plt.title('AFM, CF, and AF3 DockQ scores(first-ranked)', fontsize=23, fontweight='bold', pad=20)  # Added pad parameter

# Construct the output path
output_filename = 'Fig_2b.jpeg'
output_path = f"{figures_directory}/{output_filename}"

# Save the figure with a white background
plt.savefig(output_path, format='png', dpi=1000, facecolor='white')

plt.show()


# In[ ]:


"""
Plot Top-1 (max) DockQ vs. sampling size for ONE protein, using config_DockQ.json.

- Looks for Excel files under: <config["High_sampling"]>/<PROTEIN_CODE>/*.xlsx
- Saves the figure to: <config[OUTPUT_DIR_KEY]>/<OUTPUT_BASENAME>
"""

import re
import json
from pathlib import Path
from typing import Union, Optional

import pandas as pd
import matplotlib.pyplot as plt

# ========== USER CONTROLS ==========
PROTEIN_CODE    = "8Q7K"                       # subfolder name under High_sampling
OUTPUT_DIR_KEY  = "Supp_Figures_directory"     # or "Main_Figures_directory"
CUSTOM_TITLE    = f"Effect of sampling size on top-1 DockQ ({PROTEIN_CODE})"
OUTPUT_BASENAME = f"{PROTEIN_CODE}_Sampling.png"
CONFIG_FILE     = "config_DockQ.json"
# ====================================

# ---- Plot appearance ----
FIGSIZE, DPI = (10, 6), 180
TITLE_FONTSIZE, LABEL_FONTSIZE, TICK_FONTSIZE, LEGEND_FONTSIZE = 22, 22, 17, 17
GRID_ALPHA = 0.35

COLOR_TOP1 = "#1f77b4"
LINEWIDTH_TOP1 = 4
MARKER_TOP1 = "o"
MARKERSIZE_TOP1 = 7

AUTO_ZOOM, Y_PAD = True, 0.02
FILE_GLOB = "*.xlsx"
SAMPLING_REGEX = re.compile(r"_([0-9]+)\.xlsx$", re.IGNORECASE)


def load_config(cfg_path: Union[str, Path] = CONFIG_FILE) -> dict:
    cfg_path = Path(cfg_path)
    with cfg_path.open("r", encoding="utf-8") as f:
        return json.load(f)


def sampling_from_name(fname: str) -> Optional[int]:
    m = SAMPLING_REGEX.search(fname)
    return int(m.group(1)) if m else None


def _find_dockq_series(df: pd.DataFrame) -> Optional[pd.Series]:
    # Try exact first, then case-insensitive fallback
    if "DockQ" in df.columns:
        return pd.to_numeric(df["DockQ"], errors="coerce")
    for c in df.columns:
        if str(c).strip().lower() == "dockq":
            return pd.to_numeric(df[c], errors="coerce")
    return None


def load_top1_by_sampling(folder: Path) -> pd.DataFrame:
    rows = []
    for fp in sorted(folder.glob(FILE_GLOB)):
        s = sampling_from_name(fp.name)
        if s is None:
            continue
        try:
            # Requires openpyxl installed in your env
            df = pd.read_excel(fp)
        except Exception as e:
            print(f"[WARN] Could not read {fp}: {e}")
            continue

        dockq_series = _find_dockq_series(df)
        if dockq_series is None:
            print(f"[WARN] No 'DockQ' column (case-insensitive) in {fp}. Skipping.")
            continue

        dockq = dockq_series.dropna()
        if dockq.empty:
            print(f"[WARN] No valid DockQ values in {fp}.")
            continue

        rows.append((s, float(dockq.max())))

    return pd.DataFrame(rows, columns=["Sampling", "Top1"]).sort_values("Sampling")


def plot_improvement(top1_df: pd.DataFrame, title: str, save_path: Optional[Path] = None):
    if top1_df.empty:
        print("[INFO] No data to plot.")
        return

    df = top1_df.copy().sort_values("Sampling")

    plt.figure(figsize=FIGSIZE, dpi=DPI)
    plt.plot(
        df["Sampling"], df["Top1"],
        marker=MARKER_TOP1, markersize=MARKERSIZE_TOP1,
        linewidth=LINEWIDTH_TOP1, color=COLOR_TOP1,
        label="Top-1 DockQ"
    )

    plt.xlabel("Number of Models Sampled", fontsize=LABEL_FONTSIZE)
    plt.ylabel("DockQ Score", fontsize=LABEL_FONTSIZE)
    plt.title(title, fontsize=TITLE_FONTSIZE, fontweight='bold')

    if AUTO_ZOOM:
        y_min = df["Top1"].min() - Y_PAD
        y_max = df["Top1"].max() + Y_PAD
        plt.ylim(max(0.0, y_min), min(1.0, y_max))
    else:
        plt.ylim(0.0, 1.0)

    plt.grid(True, linestyle="--", alpha=GRID_ALPHA)
    plt.xticks(fontsize=TICK_FONTSIZE)
    plt.yticks(fontsize=TICK_FONTSIZE)
    plt.legend(fontsize=LEGEND_FONTSIZE)
    plt.tight_layout()

    if save_path is not None:
        save_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, dpi=1000, format='png')
        print(f"[Saved] {save_path}")

    plt.show()
    plt.close()


def main():
    # Load config
    config = load_config(CONFIG_FILE)

    high_sampling_key = "High_sampling"
    if high_sampling_key not in config:
        raise SystemExit(
            f"Config missing '{high_sampling_key}'. Available keys: {list(config.keys())}"
        )

    # Input: <High_sampling>/<PROTEIN_CODE>
    high_sampling_root = Path(config[high_sampling_key])
    protein_dir = high_sampling_root / PROTEIN_CODE
    if not protein_dir.is_dir():
        raise SystemExit(f"Folder not found: {protein_dir}")

    # Output directory from config (Supp or Main)
    if OUTPUT_DIR_KEY not in config:
        raise SystemExit(
            f"Output dir key '{OUTPUT_DIR_KEY}' not found in config. "
            f"Available: {list(config.keys())}"
        )
    out_dir = Path(config[OUTPUT_DIR_KEY])
    out_dir.mkdir(parents=True, exist_ok=True)

    output_path = out_dir / OUTPUT_BASENAME

    print(f"[INFO] Reading Excel files from: {protein_dir}")
    print(f"[INFO] Saving figure to:        {output_path}")

    df = load_top1_by_sampling(protein_dir)
    plot_improvement(df, title=CUSTOM_TITLE, save_path=output_path)


if __name__ == "__main__":
    main()


# In[ ]:




