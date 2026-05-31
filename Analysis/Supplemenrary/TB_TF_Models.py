#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import json

# Load configuration and define a global variable
def load_config():
    global config
    with open('config_Distributions.json') as config_file:
        config = json.load(config_file)

# Call the function to load the configuration
load_config()


# In[ ]:


import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ttest_rel, combine_pvalues
import warnings

# Define the paths to the folders (replace these paths with actual folder paths)
folder_paths = {
    "AFM-TB":'path to /PpEv/Analysis/Supplemenrary/Data/AlphaFold_Multimer_All_TB_5',
    "AFM-TF": 'path to/PpEv/Analysis/Supplemenrary/Data/AlphaFold_Multimer_All_TF_5',
    "AF3": 'path to/PpEv/Analysis/Supplemenrary/Data/AlphaFold3_All',
    "CF-TB": 'path to/PpEv/Analysis/Supplemenrary/Data/ColabFold_All_TB',
    "CF-TF": 'path to/PpEv/Analysis/Supplemenrary/Data/ColabFold_All_TF',
}

# Initialize dictionary to store t-test results
summary_results = {}

# Create output directory for detailed results
output_dir = os.path.join(os.getcwd(), "t_test_detailed_results")
os.makedirs(output_dir, exist_ok=True)

# Compare each folder pair
folder_names = list(folder_paths.keys())
for i in range(len(folder_names)):
    for j in range(i + 1, len(folder_names)):
        folder1 = folder_names[i]
        folder2 = folder_names[j]
        comparison_key = f"{folder1}_vs_{folder2}"

        print(f"Comparing {folder1} vs {folder2}...")

        # Get files from both folders
        try:
            files1 = {f.split("_DockQ_data")[0]: os.path.join(folder_paths[folder1], f) 
                     for f in os.listdir(folder_paths[folder1]) if f.endswith(".xlsx")}
            files2 = {f.split("_DockQ_data")[0]: os.path.join(folder_paths[folder2], f) 
                     for f in os.listdir(folder_paths[folder2]) if f.endswith(".xlsx")}
        except FileNotFoundError as e:
            print(f"Error accessing folder: {e}")
            continue

        # Find common PDB IDs in both folders
        common_pdb_ids = set(files1.keys()).intersection(set(files2.keys()))
        if not common_pdb_ids:
            print(f"No common PDB IDs found between {folder1} and {folder2}")
            continue

        print(f"Found {len(common_pdb_ids)} common PDB IDs")

        valid_results_count = 0
        p_values = []
        pdb_detail_df = pd.DataFrame(columns=['PDB_ID', 't_statistic', 'p_value', 'sample_size'])

        for pdb_id in common_pdb_ids:
            try:
                df1 = pd.read_excel(files1[pdb_id])
                df2 = pd.read_excel(files2[pdb_id])

                required_cols = ['AlphaFold Rank', 'DockQ']
                if not all(col in df1.columns for col in required_cols) or \
                   not all(col in df2.columns for col in required_cols):
                    print(f"Warning: Missing required columns in {pdb_id}")
                    continue

                merged_df = pd.merge(df1[required_cols], 
                                     df2[required_cols], 
                                     on='AlphaFold Rank', 
                                     suffixes=('_1', '_2'))
                if len(merged_df) < 2:
                    print(f"Warning: Not enough data points for {pdb_id} (n={len(merged_df)})")
                    continue

                t_statistic, p_value = ttest_rel(merged_df['DockQ_1'], merged_df['DockQ_2'])
                result = {
                    'PDB_ID': pdb_id,
                    't_statistic': t_statistic,
                    'p_value': p_value,
                    'sample_size': len(merged_df)
                }
                pdb_detail_df = pd.concat([pdb_detail_df, pd.DataFrame([result])], ignore_index=True)

                if not np.isnan(t_statistic) and not np.isnan(p_value):
                    p_values.append(p_value)
                    valid_results_count += 1

            except Exception as e:
                print(f"Error processing {pdb_id}: {e}")
                continue

        pdb_detail_df.to_excel(os.path.join(output_dir, f"{comparison_key}_detailed.xlsx"), index=False)

        if valid_results_count > 0:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                try:
                    combined_stat, combined_p = combine_pvalues(p_values, method='fisher')
                except:
                    combined_stat, combined_p = np.nan, np.nan

            summary_results[comparison_key] = {
                'n_common_pdbs': len(common_pdb_ids),
                'n_valid_tests': valid_results_count,
                'n_significant': sum(1 for p in p_values if p < 0.05),
                'median_p': np.median(p_values),
                'combined_p_fisher': combined_p
            }
        else:
            print(f"No valid t-test results for {comparison_key}")
            summary_results[comparison_key] = {
                'n_common_pdbs': len(common_pdb_ids),
                'n_valid_tests': 0
            }

# Create p-value matrix for heatmap
p_value_matrix = pd.DataFrame(np.nan, index=folder_names, columns=folder_names)
for comparison, result in summary_results.items():
    if 'combined_p_fisher' in result:
        folder1, folder2 = comparison.split('_vs_')
        p_value_matrix.loc[folder1, folder2] = result.get('combined_p_fisher', np.nan)
        p_value_matrix.loc[folder2, folder1] = result.get('combined_p_fisher', np.nan)
for folder in folder_names:
    p_value_matrix.loc[folder, folder] = 1.0

# Set global font size and style for plots
plt.rcParams.update({
    'font.size': 16,
    'font.family': 'Arial',
    'axes.labelsize': 16,
    'axes.titlesize': 18,
    'xtick.labelsize': 17,
    'ytick.labelsize': 17
})

# Version 1: Standard heatmap with scientific notation for p-values
plt.figure(figsize=(12, 10))
log_p_matrix = -np.log10(p_value_matrix)

ax = sns.heatmap(
    log_p_matrix, 
    annot=p_value_matrix, 
    cmap="viridis", 
    fmt=".3e",
    cbar_kws={'label': '-log10(p-value)'}, 
    annot_kws={"size": 16, "color": "black", "weight": "bold"}
)

plt.title("Two-tailed Paired t-test (P-values)", fontsize=24, fontweight='bold', pad=20)
plt.xlabel("Method", fontsize=22, labelpad=15)
plt.ylabel("Method", fontsize=22, labelpad=15)
ax.set_xticklabels(ax.get_xticklabels(), fontweight='bold')
ax.set_yticklabels(ax.get_yticklabels(), fontweight='bold')
cbar = ax.collections[0].colorbar
cbar.ax.set_title('-log10(p-value)', fontsize=16, pad=10)
cbar.ax.axhline(y=-np.log10(0.05), color='red', linestyle='--')
cbar.ax.text(1.5, -np.log10(0.05), 'p=0.05', va='center', ha='left', color='red')

# Set diagonal cell annotations to white for better visibility
n_methods = len(folder_names)
for i in range(n_methods):
    diag_index = i * n_methods + i
    ax.texts[diag_index].set_color("white")

plt.tight_layout()
plt.savefig(os.path.join(os.getcwd(), "paired_ttest_pvalue_heatmap_scientific.png"), dpi=1000, bbox_inches='tight')
plt.close()

print("Analysis complete. Scientific p-value heatmap saved.")


# In[ ]:


import matplotlib.pyplot as plt
import numpy as np

# Data for each model/metric
categories = ['High', 'Medium', 'Acceptable', 'Incorrect']

# AlphaFold-related values
af3_values    = [38.0, 2.0, 12.0, 48.0]
af_tb_values  = [40.0, 10.0, 20.0, 30.0]
af_tf_values  = [40.0, 10.0, 20.0, 30.0]

# ColabFold-related values
cf_tb_values  = [20.0, 18.0, 24.0, 38.0]
cf_tf_values  = [18.0, 16.0, 30.0, 36.0]

# Additional metrics
GL_values     = [16.0, 14.0, 22.0, 48.0]
MD_values     = [0.0, 2.0, 34.0, 64.0]
PP_values     = [0.0, 0.0, 24.0, 76.0]
RB_values     = [0.0, 0.0, 0.0, 100.0]  # Fixed duplicate assignment

# Define labels and corresponding colors (order must match the data lists)
labels_list = ['AF3', 
               'AF_TB', 
               'AF_TF', 
               'CF_TB', 
               'CF_TF', 
               'GalaxyPepDock', 'MDockPeP2', 'pepATTRACT', 'Robetta']

colors_list = ['#FF8E44', '#006400', '#b2d8b2', 
               'darkblue', '#bfbfff', '#8A2BE2', '#00CED1', '#FF6347', 'brown']

# Group all data lists into one list (order must match labels_list)
data_values = [af3_values, af_tb_values, af_tf_values, 
               cf_tb_values, cf_tf_values, GL_values, MD_values, PP_values, RB_values]

# Directory to save figures (read from the configuration file)
figures_directory = config["Supp_Figures_directory"]

# For each category, create a separate bar chart
for i, category in enumerate(categories):
    # Create a new figure for this category
    plt.figure(figsize=(10, 6), facecolor='white')

    # X positions for the bars (one per model/metric)
    x_pos = np.arange(len(labels_list))

    # Extract the value for the current category from each data list
    values = [data[i] for data in data_values]

    # Plot the bars
    barWidth = 0.6  # Adjust this width if needed
    plt.bar(x_pos, values, color=colors_list, width=barWidth, edgecolor='grey')

    # Set the x-axis tick labels and rotate if needed
    plt.xticks(x_pos, labels_list, rotation=45, ha='right', fontsize=12)

    # Set the axis labels and title
    plt.ylabel('Percentage of DockQ scores (%)', fontsize=16)
    plt.xlabel('Prediction tools', fontsize=16)
    plt.title(f'DockQ scores for {category} predicted structures', fontsize=18, fontweight='bold', pad=15)

    # Optionally, adjust layout for a neat fit
    plt.tight_layout()

    # Save the figure with a unique filename for this category in the specified directory
    output_filename = f'Figure_S7_{category}.jpeg'
    output_path = f"{figures_directory}/{output_filename}"
    plt.savefig(output_path, format='jpeg', dpi=1000, facecolor='white')
    print(f"Figure saved as {output_path}")

    # Display and then close the figure before creating the next one
    plt.show()
    plt.close()


# In[ ]:


import pandas as pd
import matplotlib.pyplot as plt
import os

# Load configuration and define a global variable
def load_config():
    global config
    with open('config_Distributions.json') as config_file:
        config = json.load(config_file)

# Call the function to load the configuration
load_config()



def load_data(pdb_ids, tb_directory, tf_directory):
    data = {}
    for pdb_id in pdb_ids:
        # Construct file paths for TB and TF based on the given directories and pdb_id
        tb_file_path = f"{tb_directory}/{pdb_id}_DockQ_data_TB.xlsx"
        tf_file_path = f"{tf_directory}/{pdb_id}_DockQ_data_TF.xlsx"

        # Load data from Excel files
        tb_data = pd.read_excel(tb_file_path)['DockQ']
        tf_data = pd.read_excel(tf_file_path)['DockQ']

        data[pdb_id] = (tb_data, tf_data)
    return data

def print_boxplot_stats(bplot, labels):
    # Print statistics for each box in the boxplot
    for i, label in enumerate(labels):
        median = bplot['medians'][i].get_ydata()[0]
        q1 = bplot['boxes'][i].get_path().vertices[1, 1]
        q3 = bplot['boxes'][i].get_path().vertices[2, 1]
        whisker_low = bplot['whiskers'][2*i].get_ydata()[1]
        whisker_high = bplot['whiskers'][2*i+1].get_ydata()[1]
        cap_low = bplot['caps'][2*i].get_ydata()[1]
        cap_high = bplot['caps'][2*i+1].get_ydata()[1]
        mean = bplot['means'][i].get_ydata()[0]

        print(f"Statistics for {label}:")
        print(f"  Median: {median}")
        print(f"  Q1 (25th percentile): {q1}")
        print(f"  Q3 (75th percentile): {q3}")
        print(f"  Whisker low: {whisker_low}")
        print(f"  Whisker high: {whisker_high}")
        print(f"  Cap low: {cap_low}")
        print(f"  Cap high: {cap_high}")
        print(f"  Mean: {mean}\n")

def plot_dockq_boxplots(data, output_directory, dpi=1000):
    # Create a box plot
    plt.figure(figsize=(12, 8))
    all_data = []
    labels = []
    colors = []

    for pdb_id, (tb_data, tf_data) in data.items():
        all_data.extend([tb_data, tf_data])
        labels.extend([f'{pdb_id} TB', f'{pdb_id} TF'])
        colors.extend(['blue', 'orange'])

    bplot = plt.boxplot(all_data, patch_artist=True, notch=True, showmeans=True)

    # Color the boxes and set labels
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)

    plt.xticks(range(1, len(labels)+1), labels, rotation=45, ha="right", fontsize=16)
    plt.title('Evaluating TB and TF DockQ Box Plot for Multiple PDB IDs', fontsize=20)
    plt.ylabel('DockQ Score', fontsize=20)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()

    # Save the figure with the specified DPI and white background
    output_file_path = f"{output_directory}/Fig_S16.png"
    plt.savefig(output_file_path, dpi=dpi, facecolor='white')
    plt.show()

    # Print the boxplot statistics
    print_boxplot_stats(bplot, labels)

def main(pdb_ids, tb_directory, tf_directory, output_directory):
    data = load_data(pdb_ids, tb_directory, tf_directory)
    plot_dockq_boxplots(data, output_directory)

# Example usage:
tb_directory = "/Users/neginmanshour/Desktop/PpEv/Quality_Scoring_Functions/Quality_Functions/DockQ/Data/AlphaFold_Multimer_All_TB"  # Directory for template-based files
tf_directory = "/Users/neginmanshour/Desktop/PpEv/Quality_Scoring_Functions/Quality_Functions/DockQ/Data/AlphaFold_Multimer_All_TF"  # Directory for template-free files
output_directory = "/Users/neginmanshour/Desktop/PpEv/Figures/Supplementary"  # Directory to save the output image
pdb_ids = ["7SXJ", "7PRX", "7QOX", "8D7P", "8AHS"]  # Add your PDB IDs here
main(pdb_ids, tb_directory, tf_directory, output_directory)




# In[ ]:




