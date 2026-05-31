#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import json

# Load configuration and define a global variable
def load_config():
    global config
    with open('config_Heatmap.json') as config_file:
        config = json.load(config_file)

# Call the function to load the configuration
load_config()


# In[ ]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os

# Define the new correlation matrix data for template-based models
data = {
    "AFM-Score": [1, 0.28, 0.19, 0.23, 0.3, 0.15, 0.05, 0.12, 0.15, -0.02],
    "PyRosetta": [0.28, 1, 0.13, 0.13, 0.14, 0.06, 0.02, 0.12, 0.13, 0.00],
    "FoldX-Stability": [0.19, 0.13, 1, 0.16, 0.11, 0.04, 0.02, 0.07, 0.07, 0.00],
    "FoldX-Interaction": [0.23, 0.13, 0.16, 1, 0.12, 0.07, 0.02, 0.07, 0.08, 0.00],
    "HADDOCK_emscorem": [0.3, 0.14, 0.11, 0.12, 1, 0.1, 0.01, 0.06, 0.08, 0.00],
    "HADDOCK-mdscore": [0.15, 0.06, 0.11, 0.12, 0.1, 1, 0.02, 0.02, 0.02, 0.01],
    "GNN_DOVE": [0.05, 0.02, 0.02, 0.02, 0.01, 0.01, 1, 0.01, 0.00, 0.02],
    "Vina": [0.12, 0.12, 0.07, 0.07, 0.06, 0.02, 0.01, 1, 0.25, -0.02],
    "Vinardo": [0.15, 0.13, 0.07, 0.08, 0.08, 0.02, 0.00, 0.25, 1, -0.01],
    "DeepRank-GNN-esm": [-0.02, 0.00, 0.00, 0.00, 0.00, 0.00, 0.02, -0.02, -0.01, 1]
}

# Convert the data to a DataFrame
correlation_matrix = pd.DataFrame(
    data, 
    index=["AFM-Score", "PyRosetta", "FoldX-Stability", "FoldX-Interaction", "HADDOCK_emscore", "HADDOCK-mdscore", "GNN_DOVE", "Vina", "Vinardo", "DeepRank-GNN-esm"]
)

# Calculate the average excluding the diagonal values (1's)
average_values = []
for column in correlation_matrix:
    col_values = correlation_matrix[column].tolist()
    avg = np.mean([val for val in col_values if val != 1])
    average_values.append(avg)

# Add a new row for the averages
correlation_matrix.loc['Mean Value'] = average_values

# Create a heatmap
plt.figure(figsize=(12, 10), facecolor='white')  # Set facecolor to white for the figure
sns.set(font_scale=1.6)  # Adjust font scale for clearer annotations

# Pass a label for the colorbar via cbar_kws
heatmap = sns.heatmap(
    correlation_matrix, 
    annot=True, fmt=".2f", 
    cmap="RdBu_r", center=0, 
    annot_kws={"size": 18},
    cbar_kws={"label": "Correlation coefficient (r)"}
)
heatmap.set_title('Spearman Correlation matrix (TB)', fontsize=23, fontweight='bold', pad=20)

# Set the background color of the heatmap to white
heatmap.set_facecolor('white')

# Modify color bar tick label font size
colorbar = heatmap.collections[0].colorbar
colorbar.ax.tick_params(labelsize=20)

# Optionally, add a label to the colorbar (if not set by cbar_kws)
colorbar.set_label("Correlation coefficient (r)", fontsize=20, labelpad=15)

# Improve layout
plt.xticks(rotation=90, ha='right', fontsize=18)  # Rotate x labels for better readability
plt.yticks(rotation=0, fontsize=18)  # Ensure y labels are horizontal for clarity
plt.tight_layout()  # Adjust layout to make sure everything fits without overlap

# Define save path
save_path = f"{config['Main_Figures_directory']}/Fig_6a.jpeg"

# Check if directory exists and save the plot
if os.path.exists(os.path.dirname(save_path)):
    if os.access(os.path.dirname(save_path), os.W_OK):
        plt.savefig(save_path, dpi=1000, facecolor='white', bbox_inches='tight')  
        print(f"Plot successfully saved to {save_path}")
    else:
        print("Error: No write permission to the directory.")
else:
    print("Error: Directory does not exist.")

# Display the plot
plt.show()


# In[ ]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os

# Define the new correlation matrix data
data = {
    "AFM-Score": [1, 0.27, 0.22, 0.23, 0.29, 0.18, 0.08, 0.15, 0.15, -0.02],
    "PyRosetta": [0.27, 1, 0.13, 0.15, 0.19, 0.1, 0.08, 0.15, 0.14, -0.01],
    "FoldX-Stability": [0.22, 0.13, 1, 0.13, 0.13, 0.08, 0.06, 0.09, 0.08, 0.00],
    "FoldX-Interaction": [0.23, 0.15, 0.13, 1, 0.15, 0.09, 0.05, 0.12, 0.12, 0.00],
    "HADDOCK-emscore": [0.29, 0.19, 0.13, 0.15, 1, 0.12, 0.07, 0.12, 0.13, 0.00],
    "HADDOCK-mdscore": [0.18, 0.1, 0.08, 0.09, 0.12, 1, 0.05, 0.07, 0.06, 0.00],
    "GNN_DOVE": [0.08, 0.08, 0.06, 0.05, 0.07, 0.05, 1, 0.04, 0.03, 0.01],
    "Vina": [0.15, 0.15, 0.09, 0.12, 0.12, 0.07, 0.04, 1, 0.24, 0.00],
    "Vinardo": [0.15, 0.14, 0.08, 0.12, 0.13, 0.06, 0.03, 0.24, 1, 0.00],
    "DeepRank-GNN-esm": [-0.02, -0.01, 0.00, 0.00, 0.00, 0.00, 0.01, 0.00, 0.00, 1]
}

# Convert the data to a DataFrame
correlation_matrix = pd.DataFrame(
    data, 
    index=[
        "AFM-Score", "PyRosetta", "FoldX-Stability", "FoldX-Interaction", "HADDOCK-emscore",
        "HADDOCK-mdscore", "GNN_DOVE", "Vina", "Vinardo", "DeepRank-GNN-esm"
    ]
)

# Calculate the average excluding the diagonal values (1's)
average_values = []
for column in correlation_matrix:
    col_values = correlation_matrix[column].tolist()
    avg = np.mean([val for val in col_values if val != 1])
    average_values.append(avg)

# Add a new row for the averages
correlation_matrix.loc['Mean Value'] = average_values

# Create a heatmap with a colorbar label using cbar_kws
plt.figure(figsize=(12, 10), facecolor='white')
sns.set(font_scale=1.6)
heatmap = sns.heatmap(
    correlation_matrix, 
    annot=True, fmt=".2f", 
    cmap="RdBu_r", center=0, 
    annot_kws={"size": 18},
    cbar_kws={"label": "Correlation coefficient (r)"}
)
heatmap.set_title('Spearman Correlation matrix (TF)', fontsize=23, fontweight='bold', pad=20)

# Modify color bar font size and add label explicitly (if needed)
colorbar = heatmap.collections[0].colorbar
colorbar.ax.tick_params(labelsize=20)
colorbar.set_label("Correlation coefficient (r)", fontsize=20, labelpad=15)

plt.xticks(rotation=90, ha='right', fontsize=18)
plt.yticks(rotation=0, fontsize=18)
plt.tight_layout()

# Define save path
save_path = f"{config['Main_Figures_directory']}/Fig_6b.jpeg"

# Check if directory exists and save the plot
if os.path.exists(os.path.dirname(save_path)):
    if os.access(os.path.dirname(save_path), os.W_OK):
        plt.savefig(save_path, dpi=1000, facecolor='white', bbox_inches='tight')
        print(f"Plot successfully saved to {save_path}")
    else:
        print("Error: No write permission to the directory.")
else:
    print("Error: Directory does not exist.")

plt.show()


# In[ ]:


import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import json

# Load configuration and define a global variable
def load_config():
    global config
    with open('config_Heatmap.json') as config_file:
        config = json.load(config_file)

# Call the function to load the configuration
load_config()

# Set the directory path as specified in the configuration file
directory_path = config['DockQ_Loss_directory_TB']

# Construct file paths
file_paths = {
    'AFM-Score': os.path.join(directory_path, "AlphaFold_TB.xlsx"),
    'PyRosetta': os.path.join(directory_path, "Pyrosetta_TB.xlsx"),
    'FoldX-Stability': os.path.join(directory_path, "Foldx_stability_TB.xlsx"),
    'FoldX-Interaction': os.path.join(directory_path, "Foldx_Interaction_TB.xlsx"),
    'HADDOCK-emscore': os.path.join(directory_path, "Haddock_emscore_TB.xlsx"),
    'HADDOCK-mdscore': os.path.join(directory_path, "Haddock_mdscore_TB.xlsx"),
    'GNN_DOVE': os.path.join(directory_path, "gnn_dove_TB.xlsx"),
    'Vina': os.path.join(directory_path, "Vina_loss_TB.xlsx"),
    'Vinardo': os.path.join(directory_path, "Vinardo_loss_TB.xlsx"),
    'DeepRank-GNN-esm': os.path.join(directory_path, "Deep_GNN_TB.xlsx")
}

# Load the data from the Excel files into DataFrames
data_frames = {label: pd.read_excel(path) for label, path in file_paths.items()}

# Prepare the data for plotting
data_to_plot = [df['Loss'].dropna() for df in data_frames.values()]

# Calculate means and sort data by mean values
means = [np.mean(data) for data in data_to_plot]
sorted_indices = np.argsort(means)
data_to_plot_sorted = [data_to_plot[i] for i in sorted_indices]
labels_sorted = [list(data_frames.keys())[i] for i in sorted_indices]

# Create the boxplot with a white background
fig, ax = plt.subplots(figsize=(24, 14), facecolor='white')  # Set figure background to white
ax.set_facecolor('white')  # Set axes background to white

boxplot = ax.boxplot(data_to_plot_sorted, labels=labels_sorted, patch_artist=True,
                      boxprops=dict(linestyle='-', linewidth=4, color='darkblue'),
                      medianprops=dict(linestyle='-', linewidth=4, color='firebrick'),
                      whiskerprops=dict(linestyle='-', linewidth=4, color='darkblue'),
                      capprops=dict(linestyle='-', linewidth=4, color='darkblue'),
                      flierprops=dict(marker='o', markerfacecolor='darkblue', markersize=12, linestyle='none'))

# Set colors for each box
colors = plt.cm.Blues(np.linspace(0.3, 0.9, len(labels_sorted)))
for patch, color in zip(boxplot['boxes'], colors):
    patch.set_facecolor(color)

# Add a legend for the median line
median_legend = plt.Line2D([], [], color='firebrick', linewidth=4, label='Median value')
legend = ax.legend(handles=[median_legend], fontsize=32, loc='upper left')
legend.get_frame().set_facecolor('white')  # Set legend background to white
legend.get_frame().set_linewidth(0)  # Remove legend border

# Set titles and labels with increased font size
ax.set_title('Loss value comparison of scoring functions (TB)', fontsize=40, fontweight='bold', pad=20)
ax.set_ylabel('Loss', fontsize=50, fontweight='bold')
ax.set_xticklabels(labels_sorted, rotation=45, fontsize=32)
ax.tick_params(axis='y', labelsize=35)
ax.grid(True, linestyle='--', linewidth=0.5)

# Set axis lines and border
for spine in ax.spines.values():
    spine.set_linewidth(2)
    spine.set_color('black')

plt.tight_layout()

# Define save path and save the plot
save_path = f"{config['Main_Figures_directory']}/Fig_6c.jpeg"
plt.savefig(save_path, dpi=1000, facecolor='white')  # Ensure saved figure has a white background
plt.show()

# Print statistics for each box
print("Statistics for each box in the plot:")
for label, data in zip(labels_sorted, data_to_plot_sorted):
    mean = np.mean(data)
    median = np.median(data)
    q1 = np.percentile(data, 25)
    q3 = np.percentile(data, 75)
    iqr = q3 - q1
    minimum = data.min()
    maximum = data.max()
    spread = maximum - minimum
    print(f"\n{label}:")
    print(f"Mean: {mean}")
    print(f"Median: {median}")
    print(f"Q1 (25th percentile): {q1}")
    print(f"Q3 (75th percentile): {q3}")
    print(f"IQR: {iqr}")
    print(f"Minimum: {minimum}")
    print(f"Maximum: {maximum}")
    print(f"Spread: {spread}")

# Answering the questions
print("\nThe red line in each box plot represents the median value of the data.")
print("The plot is sorted based on the mean values of the 'Loss' data for each scoring function.")


# In[ ]:


import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import json

# Load configuration and define a global variable
def load_config():
    global config
    with open('config_Heatmap.json') as config_file:
        config = json.load(config_file)

# Call the function to load the configuration
load_config()

directory_path = config['DockQ_Loss_directory_TF']

# Define paths for each scoring function's Excel file
scoring_files = {
    'AFM-Score': os.path.join(directory_path, "AlphaFold_TF.xlsx"),
    'PyRosetta': os.path.join(directory_path, "Pyrosetta_TF.xlsx"),
    'FoldX-Stability': os.path.join(directory_path, "Foldx_stability_TF.xlsx"),
    'FoldX-Interaction': os.path.join(directory_path, "Foldx_Interaction_TF.xlsx"),
    'HADDOCK-emscore': os.path.join(directory_path, "Haddock_emscore_TF.xlsx"),
    'HADDOCK-mdscore': os.path.join(directory_path, "Haddock_mdscore_TF.xlsx"),
    'GNN_DOVE': os.path.join(directory_path, "gnn_dove_TF.xlsx"),
    'Vina': os.path.join(directory_path, "Vina_loss_TF.xlsx"),
    'Vinardo': os.path.join(directory_path, "Vinardo_loss_TF.xlsx"),
    'DeepRank-GNN-esm': os.path.join(directory_path, "Deep_GNN_TF.xlsx")
}

# Extract loss values from each Excel file
data_to_plot = []
labels = []
for label, path in scoring_files.items():
    df = pd.read_excel(path)
    loss_data = df['Loss'].dropna()
    data_to_plot.append(loss_data)
    labels.append(label)

# Calculate means and sort data by mean values
means = [np.mean(data) for data in data_to_plot]
sorted_indices = np.argsort(means)
data_to_plot_sorted = [data_to_plot[i] for i in sorted_indices]
labels_sorted = [labels[i] for i in sorted_indices]

# Plot configuration and generation with white background
fig, ax = plt.subplots(figsize=(24, 14), facecolor='white')  # Set figure background to white
ax.set_facecolor('white')  # Set axes background to white

colors = plt.cm.copper(np.linspace(0.9, 0.3, len(labels_sorted)))
boxprops = dict(linestyle='-', linewidth=4, color='#8B4513')
medianprops = dict(linestyle='-', linewidth=4, color='firebrick')
whiskerprops = dict(linestyle='-', linewidth=4, color='#8B4513')
capprops = dict(linestyle='-', linewidth=4, color='#8B4513')
flierprops = dict(marker='o', markerfacecolor='#8B4513', markersize=12, linestyle='none')

boxplot = ax.boxplot(data_to_plot_sorted, labels=labels_sorted, patch_artist=True,
                     boxprops=boxprops, medianprops=medianprops,
                     whiskerprops=whiskerprops, capprops=capprops, flierprops=flierprops)

for patch, color in zip(boxplot['boxes'], colors):
    patch.set_facecolor(color)

# Add a legend for the median line
median_legend = plt.Line2D([], [], color='firebrick', linewidth=4, label='Median value')
legend = ax.legend(handles=[median_legend], fontsize=32, loc='upper left')
legend.get_frame().set_facecolor('white')  # Set legend background to white
legend.get_frame().set_linewidth(0)  # Remove legend border

# Set titles and labels with increased font size
ax.set_title('Loss value comparison of scoring functions (TF)', fontsize=40, fontweight='bold', pad=20)
ax.set_ylabel('Loss', fontsize=50, fontweight='bold')
ax.set_xticklabels(labels_sorted, rotation=45, fontsize=30)
ax.tick_params(axis='y', labelsize=35)
ax.grid(True, linestyle='--', linewidth=0.5)

# Set axis lines and border
for spine in ax.spines.values():
    spine.set_linewidth(2)
    spine.set_color('black')

plt.tight_layout()

# Save the figure with a white background
save_path = f"{config['Main_Figures_directory']}/Fig_6d.jpeg"
plt.savefig(save_path, dpi=1000, facecolor='white')  # Ensure saved figure has a white background
plt.show()

# Print statistics for each box
print("Statistics for each box in the plot:")
for label, data in zip(labels_sorted, data_to_plot_sorted):
    mean = np.mean(data)
    median = np.median(data)
    q1 = np.percentile(data, 25)
    q3 = np.percentile(data, 75)
    iqr = q3 - q1
    minimum = data.min()
    maximum = data.max()
    spread = maximum - minimum
    print(f"\n{label}:")
    print(f"Mean: {mean}")
    print(f"Median: {median}")
    print(f"Q1 (25th percentile): {q1}")
    print(f"Q3 (75th percentile): {q3}")
    print(f"IQR: {iqr}")
    print(f"Minimum: {minimum}")
    print(f"Maximum: {maximum}")
    print(f"Spread: {spread}")


# In[ ]:


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import json  # Importing JSON for configuration

# Load configuration and define a global variable
def load_config():
    global config
    with open('config_Heatmap.json') as config_file:
        config = json.load(config_file)

# Call the function to load the configuration
load_config()

directory_path = config['Spearman_Correlation_directory_TB']

# Construct file paths
file_paths = {
    "spearman_correlation_AF_TB.xlsx": "AFM-Score",
    "spearman_correlation_Pyrosetta_TB.xlsx": "PyRosetta",
    "correlations_stability_TB.xlsx": "FoldX-Stability",
    "correlations_Interaction_TB.xlsx": "FoldX-Interaction",
    "correlations_hadd_em_TB.xlsx": "HADDOCK-emscore",
    "correlations_hadd_md_TB.xlsx": "HADDOCK-mdscore",
    "correlations_gnn_dove_TB.xlsx": "GNN_DOVE",
    "correlations_vina_TB.xlsx": "Vina",
    "correlations_vinardo_TB.xlsx": "Vinardo",
    "correlations_Deep_GNN_TB.xlsx": "DeepRank-GNN-esm"
}

# Load data and prepare for plotting
data_with_labels = []
for file_name, label in file_paths.items():
    file_path = os.path.join(directory_path, file_name)
    df = pd.read_excel(file_path)
    data_series = df['Spearman Correlation'].dropna()
    median = data_series.median()
    iqr = data_series.quantile(0.75) - data_series.quantile(0.25)
    data_with_labels.append((data_series, label, median, iqr))

# Sort data based on median and IQR
data_with_labels.sort(key=lambda x: (-x[2], x[3]))

# Unzip the sorted data and labels for plotting
data_to_plot, labels, _, _ = zip(*data_with_labels)

# Define a gradient of blue colors (darker for higher medians)
colors = plt.cm.Blues(np.linspace(0.3, 0.9, len(labels)))  # Reversed order for darker to lighter

# Create the figure with white background
plt.figure(figsize=(24, 14), facecolor='white')
ax = plt.gca()
ax.set_facecolor('white')  # Set axes background to white

# Add border around the entire plot
ax.spines['top'].set_visible(True)
ax.spines['right'].set_visible(True)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)
ax.spines['top'].set_linewidth(2)
ax.spines['right'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.spines['top'].set_color('black')
ax.spines['right'].set_color('black')
ax.spines['bottom'].set_color('black')
ax.spines['left'].set_color('black')

# Generate multiple boxplots with custom colors and properties
for i, data in enumerate(data_to_plot):
    box = plt.boxplot(
        data, 
        positions=[i + 1], 
        patch_artist=True, 
        widths=0.6,
        boxprops=dict(facecolor=colors[i], edgecolor='darkblue', linewidth=2),
        medianprops=dict(color='firebrick', linewidth=2),
        whiskerprops=dict(color='darkblue', linewidth=2),
        capprops=dict(color='darkblue', linewidth=2),
        flierprops=dict(marker='o', markerfacecolor='darkblue', markersize=8, markeredgecolor='darkblue'),
        labels=[labels[i]]
    )

# Customize x-axis to show labels for each position
plt.xticks(range(1, len(labels) + 1), labels, rotation=45, fontsize=32)
plt.yticks(fontsize=35)

# Add a legend with white background but no border
median_line = plt.Line2D([], [], color='firebrick', linewidth=4, label='Median value')
legend = plt.legend(handles=[median_line], fontsize=32, loc='upper right', 
                   frameon=True, facecolor='white', edgecolor='none')
legend.get_frame().set_linewidth(0)  # Remove border from legend

# Title and labels
plt.title('Spearman Correlation value of scoring functions (TB)', fontsize=37, fontweight='bold', pad=20)
plt.ylabel('Spearman Correlation', fontsize=35, fontweight='bold')
plt.grid(True, linestyle='--', linewidth=0.5)

plt.tight_layout()

# Define save path and save with white background
save_path = f"{config['Main_Figures_directory']}/Fig_6e.jpeg"
plt.savefig(save_path, dpi=1000, facecolor='white', edgecolor='black', 
           bbox_inches='tight', pad_inches=0.1)

plt.show()

# Print statistics for each box
print("Statistics for each box in the plot:")
for label, data in zip(labels, data_to_plot):
    median = np.median(data)
    q1 = np.percentile(data, 25)
    q3 = np.percentile(data, 75)
    iqr = q3 - q1
    minimum = data.min()
    maximum = data.max()
    spread = maximum - minimum
    print(f"\n{label}:")
    print(f"Median: {median}")
    print(f"Q1 (25th percentile): {q1}")
    print(f"Q3 (75th percentile): {q3}")
    print(f"IQR: {iqr}")
    print(f"Minimum: {minimum}")
    print(f"Maximum: {maximum}")
    print(f"Spread: {spread}")

# Answering the questions
print("\nThe red line in each box plot represents the median value of the data.")
print("The plot is sorted based on the median values of the 'Spearman Correlation' data for each scoring function.")


# In[ ]:


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np  # Ensure numpy is imported
import os  # For path operations
import json

# Load configuration and define a global variable
def load_config():
    global config
    with open('config_Heatmap.json') as config_file:
        config = json.load(config_file)

# Call the function to load the configuration
load_config()

directory_path = config['Spearman_Correlation_directory_TF']

# Define file paths and corresponding labels in a dictionary
file_paths = {
    "spearman_correlation_AF_TF.xlsx": "AFM-Score",
    "spearman_correlation_Pyrosetta_TF.xlsx": "PyRosetta",
    "correlations_stability_TF.xlsx": "FoldX-Stability",
    "correlations_Interaction_TF.xlsx": "FoldX-Interaction",
    "correlations_hadd_em_TF.xlsx": "HADDOCK-emscore",
    "correlations_hadd_md_TF.xlsx": "HADDOCK-mdscore",
    "correlations_gnn_dove_TF.xlsx": "GNN_DOVE",
    "correlations_vina_TF.xlsx": "Vina",
    "correlations_vinardo_TF.xlsx": "Vinardo",
    "correlations_Deep_GNN_TF.xlsx": "DeepRank-GNN-esm"
}

# Load data and prepare for plotting
data_with_labels = []
for file_name, label in file_paths.items():
    file_path = os.path.join(directory_path, file_name)
    df = pd.read_excel(file_path)
    data_series = df['Spearman Correlation'].dropna()
    median = data_series.median()
    iqr = data_series.quantile(0.75) - data_series.quantile(0.25)
    data_with_labels.append((data_series, label, median, iqr))

# Sort data based on median and then by IQR
data_with_labels.sort(key=lambda x: (-x[2], x[3]))  # Sort primarily by median (desc), secondarily by IQR (asc)

# Unzip the sorted data for plotting
data_to_plot, labels, _, _ = zip(*data_with_labels)

# Define a gradient of brown colors (lighter for higher medians)
colors = plt.cm.copper(np.linspace(0.9, 0.3, len(data_to_plot)))  # Light to dark brown

# Create the figure with white background
plt.figure(figsize=(24, 14), facecolor='white')
ax = plt.gca()
ax.set_facecolor('white')  # Set axes background to white

# Add border around the entire plot
ax.spines['top'].set_visible(True)
ax.spines['right'].set_visible(True)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)
ax.spines['top'].set_linewidth(2)
ax.spines['right'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.spines['top'].set_color('black')
ax.spines['right'].set_color('black')
ax.spines['bottom'].set_color('black')
ax.spines['left'].set_color('black')

# Generate boxplot with a gradient of colors
for i, data in enumerate(data_to_plot):
    plt.boxplot(data, positions=[i + 1], widths=0.6, patch_artist=True,
                boxprops=dict(facecolor=colors[i], edgecolor='darkblue', linewidth=2),
                medianprops=dict(color='firebrick', linewidth=2),
                whiskerprops=dict(color='darkblue', linewidth=2),
                capprops=dict(color='darkblue', linewidth=2),
                flierprops=dict(marker='o', markerfacecolor='darkblue', markersize=8, markeredgecolor='darkblue'),
                labels=[f"{labels[i]}"])

# Add a legend with white background but no border
median_line = plt.Line2D([], [], color='firebrick', linewidth=4, label='Median value')
legend = plt.legend(handles=[median_line], fontsize=32, loc='upper right', 
                   frameon=True, facecolor='white', edgecolor='none')
legend.get_frame().set_linewidth(0)  # Remove border from legend

# Title and labels
plt.title('Spearman Correlation value of scoring functions (TF)', fontsize=37, fontweight='bold', pad=20)
plt.ylabel('Spearman Correlation', fontsize=35, fontweight='bold')
plt.xticks(rotation=45, fontsize=32)
plt.yticks(fontsize=35)
plt.grid(True, linestyle='--', linewidth=0.5)

plt.tight_layout()

# Define save path and save with white background
save_path = f"{config['Main_Figures_directory']}/Fig_6f.jpeg"
plt.savefig(save_path, dpi=1000, facecolor='white', edgecolor='black', 
           bbox_inches='tight', pad_inches=0.1)

plt.show()

# Print median, IQR, min, max, and spread for each box
print("Statistics for each box in the plot:")
for label, data in zip(labels, data_to_plot):
    median = np.median(data)
    q1 = np.percentile(data, 25)
    q3 = np.percentile(data, 75)
    iqr = q3 - q1
    minimum = data.min()
    maximum = data.max()
    spread = maximum - minimum
    print(f"\n{label}:")
    print(f"Median: {median}")
    print(f"Q1 (25th percentile): {q1}")
    print(f"Q3 (75th percentile): {q3}")
    print(f"IQR: {iqr}")
    print(f"Minimum: {minimum}")
    print(f"Maximum: {maximum}")
    print(f"Spread: {spread}")

# Answering the questions
print("\nThe red line in each box plot represents the median value of the data.")
print("The plot is sorted based on the median values of the 'Spearman Correlation' data for each scoring function.")


# In[ ]:


import os
import pandas as pd
from matplotlib_venn import venn3
import matplotlib.pyplot as plt

# Define the path to the directory containing the Excel files
directory_path = '/Users/neginmanshour/Desktop/PpEv/Analysis/Heatmaps_Spearman_Correlation/Data/combined_TB_TF_top20_md'
output_directory = '/Users/neginmanshour/Desktop/PpEv/Analysis/Heatmaps_Spearman_Correlation/Data/Output_TB_TF_venn'
output_file = os.path.join(output_directory, 'common_models_md_20.xlsx')
image_output_path = os.path.join(output_directory, 'venn_diagram_fig_6_final.png')

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

# Save the common models to an Excel file
common_models_df = pd.DataFrame(list(common_all_three), columns=['File name'])
common_models_df.to_excel(output_file, index=False)

# Plot the Venn diagram with custom colors
plt.figure(figsize=(8, 8))
venn = venn3([all_models[func] for func in scoring_functions], set_labels=None)

# Set custom text labels and font size
labels = {
    '100': 'AFM-Score',
    '010': 'FoldX-Stability',
    '001': 'HADDOCK-mdscore',
    '110': 'AFM/Foldx',
    '101': 'AFM/MD',
    '011': 'Foldx/MD',
    '111': 'com-3'
}
for label_id, label_text in labels.items():
    venn.get_label_by_id(label_id).set_text(label_text)
    venn.get_label_by_id(label_id).set_fontsize(15)  # Change font size here

# Change the color of each circle
venn.get_patch_by_id('100').set_color('#ff9999')
venn.get_patch_by_id('010').set_color('#66b3ff')
venn.get_patch_by_id('001').set_color('#99ff99')

plt.title('Common Models Among three Scoring Functions', fontsize=18, fontweight='bold')

# Save the figure
plt.savefig(image_output_path, facecolor='white', bbox_inches='tight')

plt.show()

print(f"Common models among all three scoring functions have been saved to {output_file}")
print(f"Venn diagram image has been saved to {image_output_path}")



# In[ ]:


import os
import pandas as pd
from matplotlib_venn import venn3
import matplotlib.pyplot as plt

# Load configuration and define a global variable
def load_config():
    global config
    with open('config_Heatmap.json') as config_file:
        config = json.load(config_file)

# Call the function to load the configuration
load_config()

# Assume config is already loaded with necessary paths
# Define the path to the directory containing the Excel files
directory_path = config["combined_directory"]  # Use the loaded config for directory path
output_directory = config["common_outputs"]
output_file = os.path.join(output_directory, 'common_models_md_20.xlsx')
image_output_path = os.path.join(config["Main_Figures_directory"], 'Fig_7a.jpeg')

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

# Save the common models to an Excel file
common_models_df = pd.DataFrame(list(common_all_three), columns=['File name'])
common_models_df.to_excel(output_file, index=False)

# Plot the Venn diagram with custom colors
plt.figure(figsize=(8, 8))
venn = venn3([all_models[func] for func in scoring_functions], set_labels=None)

# Set custom text labels and font size
labels = {
    '100': 'AFM-Score',
    '010': 'FoldX-Stability',
    '001': 'HADDOCK-mdscore',
    '110': 'AFM/Foldx',
    '101': 'AFM/MD',
    '011': 'Foldx/MD',
    '111': 'com-3'
}
for label_id, label_text in labels.items():
    venn.get_label_by_id(label_id).set_text(label_text)
    venn.get_label_by_id(label_id).set_fontsize(15)  # Change font size here

# Change the color of each circle
venn.get_patch_by_id('100').set_color('#ff9999')
venn.get_patch_by_id('010').set_color('#66b3ff')
venn.get_patch_by_id('001').set_color('#99ff99')

plt.title('Common models among three scoring functions', fontsize=18, fontweight='bold')

# Save the figure
plt.savefig(image_output_path, facecolor='white', bbox_inches='tight')

plt.show()

print(f"Common models among all three scoring functions have been saved to {output_file}")
print(f"Venn diagram image has been saved to {image_output_path}")


# In[ ]:


import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import json

# Load configuration and define a global variable
def load_config():
    global config
    with open('config_Heatmap.json') as config_file:
        config = json.load(config_file)

# Call the function to load the configuration
load_config()

# Set the directory path as specified in the configuration file
directory_path = config['Loss_consensus_TB']

# Construct file paths
file_paths = {
    'AFM-Score': os.path.join(directory_path, "AlphaFold_TB.xlsx"),
    'FoldX-Stability': os.path.join(directory_path, "Foldx_stability_TB.xlsx"),
    'HADDOCK-mdscore': os.path.join(directory_path, "Haddock_mdscore_TB.xlsx"),
    'Com-score-3': os.path.join(directory_path, "updated_common_models_final_md_20.xlsx"),
    'AFM-Score-FoldX': os.path.join(directory_path, "updated_common_models_af_foldx_final.xlsx"),
    'AFM-Score-HADDOCK': os.path.join(directory_path, "updated_common_models_af_haddock_final.xlsx"),
    'FoldX-HADDOCK': os.path.join(directory_path, "updated_common_models_foldx_haddock_final.xlsx")
}

# Load the data from the Excel files into DataFrames
data_frames = {label: pd.read_excel(path) for label, path in file_paths.items()}

# Prepare the data for plotting
data_to_plot = [df['Loss'].dropna() for df in data_frames.values()]

# Calculate means and sort data by mean values
means = [np.mean(data) for data in data_to_plot]
sorted_indices = np.argsort(means)
data_to_plot_sorted = [data_to_plot[i] for i in sorted_indices]
labels_sorted = [list(data_frames.keys())[i] for i in sorted_indices]

# Create a mapping of labels to colors
color_mapping = {
    'AFM-Score': '#ff9999',  # Pink
    'FoldX-Stability': '#7FD4FF',  # Blue
    'HADDOCK-mdscore': '#99ff99',  # Green
    'Com-score-3': '#ffcc99',  # Orange
    'AFM-Score-FoldX': '#D2B48C',  # Light Purple
    'AFM-Score-HADDOCK': '#c2c2f0',  # Light Pink
    'FoldX-HADDOCK': '#66b3ff'  # Light Green
}

# Create the figure with white background
plt.figure(figsize=(20, 11), facecolor='white')
ax = plt.gca()
ax.set_facecolor('white')  # Set axes background to white

# Add border around the entire plot
ax.spines['top'].set_visible(True)
ax.spines['right'].set_visible(True)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)
ax.spines['top'].set_linewidth(2)
ax.spines['right'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.spines['top'].set_color('black')
ax.spines['right'].set_color('black')
ax.spines['bottom'].set_color('black')
ax.spines['left'].set_color('black')

# Create the boxplot
boxplot = plt.boxplot(data_to_plot_sorted, labels=labels_sorted, patch_artist=True,
                     showfliers=False,  # This will remove outliers from the plot
                     boxprops=dict(linestyle='-', linewidth=2, edgecolor='darkblue'),
                     medianprops=dict(linestyle='-', linewidth=2, color='firebrick'),
                     whiskerprops=dict(linestyle='-', linewidth=2, color='darkblue'),
                     capprops=dict(linestyle='-', linewidth=2, color='darkblue'),
                     flierprops=dict(marker='o', markerfacecolor='darkblue', markersize=8, linestyle='none'))

# Set colors for each box based on the sorted labels
for patch, label in zip(boxplot['boxes'], labels_sorted):
    patch.set_facecolor(color_mapping[label])
    patch.set_edgecolor('darkblue')
    patch.set_linewidth(2.0)

# Add legend with white background but no border
median_line = plt.Line2D([], [], color='firebrick', linewidth=4, label='Median value')
legend = plt.legend(handles=[median_line], loc='upper left', fontsize=25,
                   frameon=True, facecolor='white', edgecolor='none')
legend.get_frame().set_linewidth(0)  # Remove border from legend

# Set titles and labels with increased font size
plt.title('Loss value comparison of scoring functions (TB)', fontsize=38, fontweight='bold', pad=20)
plt.ylabel('Loss', fontsize=37, fontweight='bold')
plt.xticks(rotation=45, fontsize=30)
plt.yticks(fontsize=35)
plt.grid(True, linestyle='--', linewidth=0.5)
plt.tight_layout()

# Define save path and save the plot with white background and border
save_path = f"{config['Main_Figures_directory']}/Fig_7b.jpeg"
plt.savefig(save_path, dpi=1000, facecolor='white', edgecolor='black', 
           bbox_inches='tight', pad_inches=0.1)

plt.show()

# Print statistics for each box
print("Statistics for each box in the plot:")
for label, data in zip(labels_sorted, data_to_plot_sorted):
    mean = np.mean(data)
    median = np.median(data)
    q1 = np.percentile(data, 25)
    q3 = np.percentile(data, 75)
    iqr = q3 - q1
    minimum = data.min()
    maximum = data.max()
    spread = maximum - minimum
    print(f"\n{label}:")
    print(f"Mean: {mean}")
    print(f"Median: {median}")
    print(f"Q1 (25th percentile): {q1}")
    print(f"Q3 (75th percentile): {q3}")
    print(f"IQR: {iqr}")
    print(f"Minimum: {minimum}")
    print(f"Maximum: {maximum}")
    print(f"Spread: {spread}")

# Answering the questions
print("\nThe red line in each box plot represents the median value of the data.")
print("The plot is sorted based on the mean values of the 'Loss' data for each scoring function.")


# In[ ]:


import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import json

# Load configuration and define a global variable
def load_config():
    global config
    with open('config_Heatmap.json') as config_file:
        config = json.load(config_file)

# Call the function to load the configuration
load_config()

# Set the directory path as specified in the configuration file
directory_path = config['Loss_consensus_TF']

# Construct file paths
file_paths = {
    'AFM-Score': os.path.join(directory_path, "AlphaFold_TF.xlsx"),
    'FoldX-Stability': os.path.join(directory_path, "Foldx_stability_TF.xlsx"),
    'HADDOCK-mdscore': os.path.join(directory_path, "Haddock_mdscore_TF.xlsx"),
    'Com-score-3': os.path.join(directory_path, "updated_common_models_final_md_20.xlsx"),
    'AFM-Score-FoldX': os.path.join(directory_path, "updated_common_models_af_foldx_final.xlsx"),
    'AFM-Score-HADDOCK': os.path.join(directory_path, "updated_common_models_af_haddock_final.xlsx"),
    'FoldX-HADDOCK': os.path.join(directory_path, "updated_common_models_foldx_haddock_final.xlsx")
}

# Load the data from the Excel files into DataFrames
data_frames = {label: pd.read_excel(path) for label, path in file_paths.items()}

# Prepare the data for plotting
data_to_plot = [df['Loss'].dropna() for df in data_frames.values()]

# Calculate means and sort data by mean values
means = [np.mean(data) for data in data_to_plot]
sorted_indices = np.argsort(means)
data_to_plot_sorted = [data_to_plot[i] for i in sorted_indices]
labels_sorted = [list(data_frames.keys())[i] for i in sorted_indices]

# Create a mapping of labels to colors
color_mapping = {
    'AFM-Score': '#ff9999',  # Pink
    'FoldX-Stability': '#7FD4FF',  # Blue
    'HADDOCK-mdscore': '#99ff99',  # Green
    'Com-score-3': '#ffcc99',  # Orange
    'AFM-Score-FoldX': '#D2B48C',  # Light Purple
    'AFM-Score-HADDOCK': '#c2c2f0',  # Light Pink
    'FoldX-HADDOCK': '#66b3ff'  # Light Green
}

# Create the figure with white background
plt.figure(figsize=(20, 11), facecolor='white')
ax = plt.gca()
ax.set_facecolor('white')  # Set axes background to white

# Add border around the entire plot
ax.spines['top'].set_visible(True)
ax.spines['right'].set_visible(True)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)
ax.spines['top'].set_linewidth(2)
ax.spines['right'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.spines['top'].set_color('black')
ax.spines['right'].set_color('black')
ax.spines['bottom'].set_color('black')
ax.spines['left'].set_color('black')

# Create the boxplot
boxplot = plt.boxplot(data_to_plot_sorted, labels=labels_sorted, patch_artist=True,
                     showfliers=False,  # This will remove outliers from the plot
                     boxprops=dict(linestyle='-', linewidth=2, edgecolor='darkblue'),
                     medianprops=dict(linestyle='-', linewidth=2, color='firebrick'),
                     whiskerprops=dict(linestyle='-', linewidth=2, color='darkblue'),
                     capprops=dict(linestyle='-', linewidth=2, color='darkblue'),
                     flierprops=dict(marker='o', markerfacecolor='darkblue', markersize=8, linestyle='none'))

# Set colors for each box based on the sorted labels
for patch, label in zip(boxplot['boxes'], labels_sorted):
    patch.set_facecolor(color_mapping[label])
    patch.set_edgecolor('darkblue')
    patch.set_linewidth(2.0)

# Add legend with white background but no border
median_line = plt.Line2D([], [], color='firebrick', linewidth=4, label='Median value')
legend = plt.legend(handles=[median_line], loc='upper left', fontsize=25,
                   frameon=True, facecolor='white', edgecolor='none')
legend.get_frame().set_linewidth(0)  # Remove border from legend

# Set titles and labels with increased font size
plt.title('Loss value comparison of scoring functions (TF)', fontsize=38, fontweight='bold', pad=20)
plt.ylabel('Loss', fontsize=37, fontweight='bold')
plt.xticks(rotation=45, fontsize=30)
plt.yticks(fontsize=35)
plt.grid(True, linestyle='--', linewidth=0.5)
plt.tight_layout()

# Define save path and save the plot with white background and border
save_path = f"{config['Main_Figures_directory']}/Fig_7c.jpeg"
plt.savefig(save_path, dpi=1000, facecolor='white', edgecolor='black', 
           bbox_inches='tight', pad_inches=0.1)

plt.show()

# Print statistics for each box
print("Statistics for each box in the plot:")
for label, data in zip(labels_sorted, data_to_plot_sorted):
    mean = np.mean(data)
    median = np.median(data)
    q1 = np.percentile(data, 25)
    q3 = np.percentile(data, 75)
    iqr = q3 - q1
    minimum = data.min()
    maximum = data.max()
    spread = maximum - minimum
    print(f"\n{label}:")
    print(f"Mean: {mean}")
    print(f"Median: {median}")
    print(f"Q1 (25th percentile): {q1}")
    print(f"Q3 (75th percentile): {q3}")
    print(f"IQR: {iqr}")
    print(f"Minimum: {minimum}")
    print(f"Maximum: {maximum}")
    print(f"Spread: {spread}")

# Answering the questions
print("\nThe red line in each box plot represents the median value of the data.")
print("The plot is sorted based on the mean values of the 'Loss' data for each scoring function.")


# In[ ]:


import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import json
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

# Load configuration and define a global variable
def load_config():
    global config
    with open('config_Heatmap.json') as config_file:
        config = json.load(config_file)

# Call the function to load the configuration
load_config()

# Set the directory path as specified in the configuration file
directory_path = config['Loss_consensus_TB']

# Construct file paths
file_paths = {
    'AFM-Score': os.path.join(directory_path, "AlphaFold_TB.xlsx"),
    'PyRosetta': os.path.join(directory_path, "Pyrosetta_TB.xlsx"),
    'FoldX-Stability': os.path.join(directory_path, "Foldx_stability_TB.xlsx"),
    'FoldX-Interaction': os.path.join(directory_path, "Foldx_Interaction_TB.xlsx"),
    'HADDOCK-emscore': os.path.join(directory_path, "Haddock_emscore_TB.xlsx"),
    'HADDOCK-mdscore': os.path.join(directory_path, "Haddock_mdscore_TB.xlsx"),
    'GNN_DOVE': os.path.join(directory_path, "gnn_dove_TB.xlsx"),
    'Vina': os.path.join(directory_path, "Vina_loss_TB.xlsx"),
    'Vinardo': os.path.join(directory_path, "Vinardo_loss_TB.xlsx"),
    'DeepRank-GNN-esm': os.path.join(directory_path, "Deep_GNN_TB.xlsx"),
    'Com-score-3': os.path.join(directory_path, "updated_common_models_final_md_20.xlsx"),
    'Com-score-4': os.path.join(directory_path, "updated_common_models_final_pyrosetta_20.xlsx"),
    'Consensus-TB': os.path.join(directory_path, "dockq_values_TB_1000.xlsx"),
}

# Load the data from the Excel files into DataFrames
data_frames = {label: pd.read_excel(path) for label, path in file_paths.items()}

# Prepare the data for plotting, removing NaN values
data_to_plot = [df['Loss'].dropna() for df in data_frames.values()]

# Calculate means and sort data by mean values
means = [np.mean(data) for data in data_to_plot]
sorted_indices = np.argsort(means)
data_to_plot_sorted = [data_to_plot[i] for i in sorted_indices]
labels_sorted = [list(data_frames.keys())[i] for i in sorted_indices]
sorted_means = [means[i] for i in sorted_indices]  # Ensure mean values are sorted

# Create the figure with white background
plt.figure(figsize=(24, 11), facecolor='white')
ax = plt.gca()
ax.set_facecolor('white')  # Set axes background to white

# Add border around the entire plot
ax.spines['top'].set_visible(True)
ax.spines['right'].set_visible(True)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)
ax.spines['top'].set_linewidth(2)
ax.spines['right'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.spines['top'].set_color('black')
ax.spines['right'].set_color('black')
ax.spines['bottom'].set_color('black')
ax.spines['left'].set_color('black')

# Create the boxplot with specified styles and colors
boxplot = plt.boxplot(data_to_plot_sorted, labels=labels_sorted, patch_artist=True,
                     boxprops=dict(linestyle='-', linewidth=2, edgecolor='darkblue'),
                     medianprops=dict(linestyle='-', linewidth=2, color='firebrick'),
                     whiskerprops=dict(linestyle='-', linewidth=2, color='darkblue'),
                     capprops=dict(linestyle='-', linewidth=2, color='darkblue'),
                     flierprops=dict(marker='', alpha=0))  # Outliers not shown

# Set colors for each box
colors = plt.cm.Blues(np.linspace(0.3, 0.9, len(labels_sorted)))
for patch, color in zip(boxplot['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_edgecolor('darkblue')
    patch.set_linewidth(2.0)

# Get the positions of the upper whiskers
whiskers = [whisker.get_ydata()[1] for whisker in boxplot['whiskers'][1::2]]

# Annotate the mean values above the upper whisker of each box plot with a slight offset
offset = 0.02  # Adjust this value as needed
for i, mean in enumerate(sorted_means):
    plt.text(i + 1, whiskers[i] + offset, f'{mean:.3f}',
             horizontalalignment='center', color='black',
             fontsize=24, verticalalignment='bottom')

# --- New code to annotate boxes with star annotations ---

# Define which boxes have different sample sizes and the corresponding star annotations.
# For example:
#   Consensus-TB (8 samples)   -> 1 star
#   Com-score-3 (16 samples)    -> 2 stars
#   Com-score-4 (7 samples)     -> 3 stars
stars_dict = {
    'Consensus-TB': '*',
    'Com-score-3': '**',
    'Com-score-4': '***'
}

# Set a vertical offset for the star annotation. Adjust as needed.
star_offset = offset * 5

# Loop over sorted labels and, if a label is in stars_dict, annotate its box with the corresponding stars.
for i, label in enumerate(labels_sorted):
    if label in stars_dict:
        plt.text(i + 1, whiskers[i] + star_offset, stars_dict[label],
                 horizontalalignment='center', color='black',
                 fontsize=30, verticalalignment='bottom')

# Create custom legend handles
median_line = mlines.Line2D([], [], color='firebrick', linewidth=4, label='Median value')
mean_text = mpatches.Patch(color='none', label='Numbers indicate mean values')

# Option to change legend location
legend_location = 'upper left'  # Change this value to specify a different location

# Add the legend with white background but no border
legend = plt.legend(handles=[median_line, mean_text], fontsize=25, loc=legend_location,
                   handlelength=2.5, handletextpad=1.5, frameon=True,
                   facecolor='white', edgecolor='none')
legend.get_frame().set_linewidth(0)  # Remove border from legend

# Set titles and labels
plt.title('Loss value comparison of scoring functions (TB)', fontsize=40, fontweight='bold')
plt.ylabel('Loss', fontsize=45, fontweight='bold')
plt.xticks(rotation=45, fontsize=25)
plt.yticks(fontsize=35)
plt.grid(True, linestyle='--', linewidth=0.5)
plt.tight_layout()

# Define save path and save the plot with white background and border
save_path = f"{config['Main_Figures_directory']}/Fig_7d.jpeg"
plt.savefig(save_path, dpi=1000, facecolor='white', edgecolor='black', 
           bbox_inches='tight', pad_inches=0.1)

plt.show()

# Print statistics for each box
print("Statistics for each box in the plot:")
for label, data in zip(labels_sorted, data_to_plot_sorted):
    mean = np.mean(data)
    median = np.median(data)
    q1 = np.percentile(data, 25)
    q3 = np.percentile(data, 75)
    iqr = q3 - q1
    minimum = data.min()
    maximum = data.max()
    spread = maximum - minimum
    print(f"\n{label}:")
    print(f"Mean: {mean}")
    print(f"Median: {median}")
    print(f"Q1 (25th percentile): {q1}")
    print(f"Q3 (75th percentile): {q3}")
    print(f"IQR: {iqr}")
    print(f"Minimum: {minimum}")
    print(f"Maximum: {maximum}")
    print(f"Spread: {spread}")


# In[ ]:


import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import json
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

# Load configuration and define a global variable
def load_config():
    global config
    with open('config_Heatmap.json') as config_file:
        config = json.load(config_file)

# Call the function to load the configuration
load_config()

directory_path = config['Loss_consensus_TF']

# Define paths for each scoring function's Excel file
scoring_files = {
    'AFM-Score': os.path.join(directory_path, "AlphaFold_TF.xlsx"),
    'PyRosetta': os.path.join(directory_path, "Pyrosetta_TF.xlsx"),
    'FoldX-Stability': os.path.join(directory_path, "Foldx_stability_TF.xlsx"),
    'FoldX-Interaction': os.path.join(directory_path, "Foldx_Interaction_TF.xlsx"),
    'HADDOCK-emscore': os.path.join(directory_path, "Haddock_emscore_TF.xlsx"),
    'HADDOCK-mdscore': os.path.join(directory_path, "Haddock_mdscore_TF.xlsx"),
    'GNN_DOVE': os.path.join(directory_path, "gnn_dove_TF.xlsx"),
    'Vina': os.path.join(directory_path, "Vina_loss_TF.xlsx"),
    'Vinardo': os.path.join(directory_path, "Vinardo_loss_TF.xlsx"),
    'DeepRank-GNN-esm': os.path.join(directory_path, "Deep_GNN_TF.xlsx"),
    'Com-score-3': os.path.join(directory_path, "updated_common_models_final_md_20.xlsx"),
    'Com-score-4': os.path.join(directory_path, "updated_common_models_final_pyrosetta_20.xlsx"),
    'Consensus-TF': os.path.join(directory_path, "dockq_values_TF_1000.xlsx"),
}

# Extract loss values from each Excel file
data_to_plot = []
labels = []
for label, path in scoring_files.items():
    df = pd.read_excel(path)
    loss_data = df['Loss'].dropna()
    data_to_plot.append(loss_data)
    labels.append(label)

# Calculate means and sort data by mean values
means = [np.mean(data) for data in data_to_plot]
sorted_indices = np.argsort(means)
data_to_plot_sorted = [data_to_plot[i] for i in sorted_indices]
labels_sorted = [labels[i] for i in sorted_indices]
sorted_means = [means[i] for i in sorted_indices]  # Ensure mean values are sorted

# Plot configuration and generation
fig, ax = plt.subplots(figsize=(24, 11), facecolor='white')  # Set figure background to white
ax.set_facecolor('white')  # Set axes background to white

# Add border around the entire plot
ax.spines['top'].set_visible(True)
ax.spines['right'].set_visible(True)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)
ax.spines['top'].set_linewidth(2)
ax.spines['right'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.spines['top'].set_color('black')
ax.spines['right'].set_color('black')
ax.spines['bottom'].set_color('black')
ax.spines['left'].set_color('black')

colors = plt.cm.copper(np.linspace(0.9, 0.3, len(labels_sorted)))
boxprops = dict(linestyle='-', linewidth=2, color='#8B4513')
medianprops = dict(linestyle='-', linewidth=2, color='firebrick')
whiskerprops = dict(linestyle='-', linewidth=2, color='#8B4513')
capprops = dict(linestyle='-', linewidth=2, color='#8B4513')
flierprops = dict(marker='', alpha=0)  # Hide outliers

boxplot = ax.boxplot(data_to_plot_sorted, labels=labels_sorted, patch_artist=True,
                     boxprops=boxprops, medianprops=medianprops,
                     whiskerprops=whiskerprops, capprops=capprops, flierprops=flierprops)

for patch, color in zip(boxplot['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_edgecolor('#8B4513')
    patch.set_linewidth(2.0)

# Get the positions of the upper whiskers for each box
whiskers = [whisker.get_ydata()[1] for whisker in boxplot['whiskers'][1::2]]

# Annotate the mean values above the upper whisker of each box plot with a slight offset
offset = 0.02  # Adjust this value as needed
for i, mean in enumerate(sorted_means):
    ax.text(i + 1, whiskers[i] + offset, f'{mean:.3f}', 
            horizontalalignment='center', color='black', fontsize=24, verticalalignment='bottom')

# --- New code to annotate boxes with stars ---

# Define which boxes have different sample sizes and the corresponding star annotation.
stars_dict = {
    'Consensus-TF': '*',
    'Com-score-3': '**',
    'Com-score-4': '***'
}

# Set a vertical offset for the star annotation. Adjust as needed.
star_offset = offset * 5

# Loop over sorted labels; if a label is in our stars_dict, annotate its box with the corresponding stars.
for i, label in enumerate(labels_sorted):
    if label in stars_dict:
        ax.text(i + 1, whiskers[i] + star_offset, stars_dict[label],
                horizontalalignment='center', color='black', fontsize=30, verticalalignment='bottom')

# Create custom legend handles
median_line = mlines.Line2D([], [], color='firebrick', linewidth=4, label='Median value')
mean_text = mpatches.Patch(color='none', label='Numbers indicate mean values')

# Option to change legend location
legend_location = 'upper left'

# Add the legend with white background but no border
legend = ax.legend(handles=[median_line, mean_text], fontsize=25, loc=legend_location, 
                   handlelength=2.5, handletextpad=1.5, frameon=True,
                   facecolor='white', edgecolor='none')
legend.get_frame().set_linewidth(0)  # Remove border from legend

# Set titles and labels
ax.set_title('Loss value comparison of scoring functions (TF)', fontsize=40, fontweight='bold')
ax.set_ylabel('Loss', fontsize=45, fontweight='bold')
ax.set_xticklabels(labels_sorted, rotation=45, fontsize=28)
ax.tick_params(axis='y', labelsize=35)
ax.grid(True, linestyle='--', linewidth=0.5)

plt.tight_layout()

# Save the figure with white background and border
save_path = f"{config['Main_Figures_directory']}/Fig_7e.jpeg"
plt.savefig(save_path, dpi=1000, facecolor='white', edgecolor='black', 
           bbox_inches='tight', pad_inches=0.1)

plt.show()

# Print statistics for each box
print("Statistics for each box in the plot:")
for label, data in zip(labels_sorted, data_to_plot_sorted):
    mean = np.mean(data)
    median = np.median(data)
    q1 = np.percentile(data, 25)
    q3 = np.percentile(data, 75)
    iqr = q3 - q1
    minimum = data.min()
    maximum = data.max()
    spread = maximum - minimum
    print(f"\n{label}:")
    print(f"Mean: {mean}")
    print(f"Median: {median}")
    print(f"Q1 (25th percentile): {q1}")
    print(f"Q3 (75th percentile): {q3}")
    print(f"IQR: {iqr}")
    print(f"Minimum: {minimum}")
    print(f"Maximum: {maximum}")
    print(f"Spread: {spread}")


# In[ ]:




