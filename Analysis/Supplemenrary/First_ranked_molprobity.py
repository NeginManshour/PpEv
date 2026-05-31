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


import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import json

# Load configuration and define a global variable
def load_config():
    global config
    with open('config_Distributions.json') as config_file:
        config = json.load(config_file)

# Call the function to load the configuration
load_config()

# Define the directories for each category from the configuration file
directories = {
    'Molprobity_TB': os.path.join(config["First_ranked_MolProbity"], 'AFM-TB.xlsx'),
    'Molprobity_TF': os.path.join(config["First_ranked_MolProbity"], 'AFM-TF.xlsx'),
    'Colab_Mol_TB': os.path.join(config["First_ranked_MolProbity"], 'CF-TB.xlsx'),
    'Colab_Mol_TF': os.path.join(config["First_ranked_MolProbity"], 'CF-TF.xlsx'),
    'AF3': os.path.join(config["First_ranked_MolProbity"], 'AF3.xlsx'),
    'Native_Molprobity_TF': os.path.join(config["First_ranked_MolProbity"], 'GT.xlsx')
}

# List to store individual DataFrames for each category
category_dfs = {cat: [] for cat in directories}

# Loop through all files in each directory
for category, file_path in directories.items():
    df = pd.read_excel(file_path)
    df['Category'] = category
    category_dfs[category] = df

# Combine DataFrames for each category
combined_df = pd.concat(category_dfs.values(), ignore_index=True)

# Ensure the relevant columns are numeric
parameters = ['Rama_Z_whole', 'Rama_Z_helix', 'Rama_Z_sheet', 'Rama_Z_loop', 'Clashscore']
for param in parameters:
    combined_df[param] = pd.to_numeric(combined_df[param], errors='coerce')

# Define the new names for categories
category_order_renamed = ['AFM-TB', 'AFM-TF', 'CF-TB', 'CF-TF', 'AF3', 'GT']
category_mapping = {
    'Molprobity_TB': 'AFM-TB',
    'Molprobity_TF': 'AFM-TF',
    'Colab_Mol_TB': 'CF-TB',
    'Colab_Mol_TF': 'CF-TF',
    'AF3': 'AF3',
    'Native_Molprobity_TF': 'GT'
}

# Updating the category names in the DataFrame
combined_df['Category'] = combined_df['Category'].map(category_mapping)

# Set font sizes for each part of the plot
def plot_box_plots(output_directory, title_fontsize=22, xlabel_fontsize=20, ylabel_fontsize=20, xticks_fontsize=18, yticks_fontsize=20):
    plot_colors = ['lightgreen', 'lightblue', 'lightcoral', 'lightpink', 'lightyellow', 'lightgrey']

    # Create the output directory if it doesn't exist
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Plotting boxplots for each parameter
    for i, (param, color) in enumerate(zip(parameters, plot_colors), start=9):  # Start naming from Fig_S9
        plt.figure(figsize=(12, 8))

        # Drop rows with NaN values in the parameter column
        df_to_plot = combined_df.dropna(subset=[param])

        # Creating the boxplot with specified order and new names
        box = sns.boxplot(x='Category', y=param, data=df_to_plot, order=category_order_renamed, palette=[color])

        # Change the color of the median line to red
        for line in box.lines:
            # Lines 4n+2 (n=0,1,2,...) are the median lines
            if box.lines.index(line) % 6 == 4:
                line.set_color('red')
                line.set_linewidth(2)

        plt.title(f'Distribution of {param} (first-ranked models)', fontsize=title_fontsize)
        plt.xlabel('', fontsize=0)  # Hide the x-axis label
        plt.ylabel(param, fontsize=ylabel_fontsize)
        plt.xticks(fontsize=xticks_fontsize, rotation=45)
        plt.yticks(fontsize=yticks_fontsize)
        plt.grid(True)

        # Save the plot as a JPEG image with dpi 1000
        plot_path = os.path.join(output_directory, f'Figure_S{i}.jpeg')
        plt.savefig(plot_path, dpi=1000, format='jpeg', bbox_inches='tight')

        # Display the plot
        plt.show()

        # Print statistics for each category and parameter
        print(f"Statistics for {param}:")
        for category in category_order_renamed:
            data = df_to_plot[df_to_plot['Category'] == category][param].dropna()
            median = data.median()
            q1 = data.quantile(0.25)
            q3 = data.quantile(0.75)
            iqr = q3 - q1
            minimum = data.min()
            maximum = data.max()
            spread = maximum - minimum
            print(f"{category}:")
            print(f"  Median: {median}")
            print(f"  Q1 (25th percentile): {q1}")
            print(f"  Q3 (75th percentile): {q3}")
            print(f"  IQR: {iqr}")
            print(f"  Minimum: {minimum}")
            print(f"  Maximum: {maximum}")
            print(f"  Spread: {spread}")
        print("\n")

# Usage example
output_directory = config["Supp_Figures_directory"]
plot_box_plots(output_directory, title_fontsize=22, xlabel_fontsize=20, ylabel_fontsize=20, xticks_fontsize=18, yticks_fontsize=20)



# In[ ]:


import os
import pandas as pd
import plotly.graph_objects as go
import json
import numpy as np

# Load configuration and define a global variable
def load_config():
    global config
    with open('config_Distributions.json') as config_file:
        config = json.load(config_file)

# Call the function to load the configuration
load_config()

# Define the directory path where all the Excel files are located
directory_path = config["First_ranked_MolProbity"]

# Define the output directory to save the plot image
output_directory = config["Main_Figures_directory"]

# Define the file names
file_names = ['AFM-TB.xlsx', 'AFM-TF.xlsx', 'AF3.xlsx', 'CF-TB.xlsx', 'CF-TF.xlsx', 'GT.xlsx']

# Load the data from the first sheet of each file
dfs = {file_name: pd.read_excel(os.path.join(directory_path, file_name), sheet_name='Sheet1') for file_name in file_names}

# Select the columns for the radar plot
columns_for_radar = ['Clashscore', 'Rama_Z_whole', 'Rama_Z_helix', 'Rama_Z_sheet', 'Rama_Z_loop']

# Replace 'Not Found' with NaN and convert to numeric for all datasets
for df in dfs.values():
    for column in columns_for_radar:
        df[column] = pd.to_numeric(df[column], errors='coerce')

# Calculate the RMS values for the Rama_Z columns and median for Clashscore
rms_values = {file_name: {} for file_name in file_names}

for file_name, df in dfs.items():
    for column in columns_for_radar[1:]:  # Skipping the first column, which is 'Clashscore'
        rms_value = np.sqrt(np.mean(df[column] ** 2))
        rms_values[file_name][column] = rms_value
    # Add median Clashscore to the rms_values dictionary
    rms_values[file_name]['Clashscore'] = df['Clashscore'].median()

# Convert to DataFrame for easy handling
rms_df = pd.DataFrame(rms_values).T

# Add the new parameters: Twisted Peptides and Cis Non-Proline
twisted_peptides = {
    'GT.xlsx': 1.7,
    'AFM-TB.xlsx': 11.7,
    'AFM-TF.xlsx': 15.0,
    'CF-TB.xlsx': 16.7,
    'CF-TF.xlsx': 26.7,
    'AF3.xlsx': 0.0
}

cis_non_proline = {
    'GT.xlsx': 3.4,
    'AFM-TB.xlsx': 1.7,
    'AFM-TF.xlsx': 1.7,
    'CF-TB.xlsx': 3.3,
    'CF-TF.xlsx': 5.0,
    'AF3.xlsx': 6.7
}

rms_df['Twisted Peptides'] = pd.Series(twisted_peptides)
rms_df['cis non-proline'] = pd.Series(cis_non_proline)

# Print the Rama_Z scores and corresponding RMS values for each dataset
for column in columns_for_radar[1:]:
    print(f"\n{column}:")
    for file_name in file_names:
        print(f"{file_name.replace('.xlsx', '')}: {rms_df.at[file_name, column]:.4f}")

# Update columns for radar to include the new parameters with updated labels
columns_for_radar += ['Twisted Peptides', 'cis non-proline']

# Scale the RMS and new parameter values where lower is better
scaled_rms_df = pd.DataFrame()
scale_min = 0.25
scale_max = 1.0

for column in columns_for_radar:
    min_val = rms_df[column].min()
    max_val = rms_df[column].max()
    scaled_rms_df[column] = scale_max - (rms_df[column] - min_val) / (max_val - min_val) * (scale_max - scale_min)

# Reorder the datasets by their scaled values
reordered_data = {name: [0]*len(columns_for_radar) for name in rms_df.index}
for i, category in enumerate(columns_for_radar):
    sorted_datasets = scaled_rms_df[category].sort_values(ascending=False)
    for dataset_name in sorted_datasets.index:
        reordered_data[dataset_name][i] = sorted_datasets[dataset_name]

# Prepare the categories and values for the radar plot
categories = columns_for_radar

# Create the radar plot
fig = go.Figure()

# Define colors for each dataset using a high-contrast palette
colors = {
    'AFM-TB.xlsx': 'rgba(31, 119, 180, 1.0)',   # Blue for AFM-TB
    'AFM-TF.xlsx': 'rgba(255, 127, 14, 1.0)',    # Orange for AFM-TF
    'AF3.xlsx':    'rgba(44, 160, 44, 1.0)',      # Green for AF3
    'CF-TB.xlsx':  'rgba(214, 39, 40, 1.0)',       # Red for CF-TB
    'CF-TF.xlsx':  'rgba(148, 103, 189, 1.0)',     # Purple for CF-TF
    'GT.xlsx':     'rgba(0, 0, 0, 1.0)'            # Black for GT
}

# Add the datasets to the radar plot with solid lines
for dataset_name, values in reordered_data.items():
    fig.add_trace(go.Scatterpolar(
        r=values + [values[0]],  # Complete the loop
        theta=categories + [categories[0]],  # Complete the loop
        fill=None,  
        opacity=1.0,  
        name=dataset_name,
        line=dict(color=colors[dataset_name])
    ))

# Update the layout with larger and black angular axis label font
fig.update_layout(
    polar=dict(
        radialaxis=dict(
            visible=True,
            range=[0, 1],
            showticklabels=False
        ),
        angularaxis=dict(
            tickfont=dict(size=20, color='black')  # Set larger font size and color to black
        )
    ),
    showlegend=False,
    title=dict(
        text='Analysing the Parameters of MolProbity',
        font=dict(size=22, color='black'),
        xanchor='center',
        yanchor='top',
        x=0.5
    )
)

# Show the plot inline
fig.show()

# Save the plot as an image with high resolution
output_image_path = os.path.join(output_directory, 'Fig_3b.jpeg')
fig.write_image(output_image_path, format='png', scale=10, engine='kaleido')

print(f'Radar plot without legend saved at {output_image_path}')


# In[ ]:


import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

# Define directories manually for input files (colors remain unchanged)
directories = {
    'AF3': ('/Users/neginmanshour/Desktop/PpEv/Analysis/Supplemenrary/Data/AF3', '#FF8E44'),
    'AFM-TB': ('/Users/neginmanshour/Desktop/PpEv/Analysis/Supplemenrary/Data/AF_TB_5', 'darkgreen'),
    'AFM-TF': ('/Users/neginmanshour/Desktop/PpEv/Analysis/Supplemenrary/Data/AF_TF_5', '#b2d8b2'),
    'Relaxed-TB': ('/Users/neginmanshour/Desktop/PpEv/Analysis/Supplemenrary/Data/relaxed_TB', '#800080'),
    'Relaxed-TF': ('/Users/neginmanshour/Desktop/PpEv/Analysis/Supplemenrary/Data/relaxed_TF', '#D8BFD8')
}

# Set global font sizes and styles for the plot
plt.rcParams.update({
    'font.size': 15, 
    'axes.labelsize': 20, 
    'axes.titlesize': 32, 
    'xtick.labelsize': 26, 
    'ytick.labelsize': 26, 
    'legend.fontsize': 25
})

# Function to process each directory and return the MolProbity score distributions
def process_directory(directory_path):
    # Create detailed ranges: 0.0-0.2, 0.2-0.4, etc.
    range_counts = {f"{i/10:.1f}-{(i/10 + 0.2):.1f}": 0 for i in range(0, 30, 2)}

    total_files = 0
    for file in os.listdir(directory_path):
        if file.endswith(".xlsx") and not file.startswith("~$"):
            file_path = os.path.join(directory_path, file)
            try:
                df = pd.read_excel(file_path)
                average_score = df['Molprobity_score'].mean()
                index = int(average_score * 10 // 2) * 2

                # Ensure index is within our range bounds
                if 0 <= index < 30:
                    range_key = f"{index/10:.1f}-{(index/10 + 0.2):.1f}"
                    range_counts[range_key] += 1

                total_files += 1
            except Exception as e:
                print(f"Error processing file {file}: {e}")

    # Calculate percentages
    range_percentages = {k: (v / total_files * 100) if total_files else 0 for k, v in range_counts.items()}

    # Sort range keys to ensure correct display order
    sorted_keys = sorted(range_percentages.keys(), key=lambda x: float(x.split('-')[0]))
    sorted_percentages = [range_percentages[k] for k in sorted_keys]

    return sorted_keys, sorted_percentages

# Set up the figure
fig, ax = plt.subplots(figsize=(22, 12))
fig.patch.set_facecolor('white')
ax.set_facecolor('white')
bar_width = 0.17

# Plotting
for i, (category, (directory_path, color)) in enumerate(directories.items()):
    labels, percentages = process_directory(directory_path)
    bar_positions = [x + (bar_width * i) for x in range(len(labels))]
    ax.bar(bar_positions, percentages, width=bar_width, color=color, label=category, edgecolor='grey')

ax.set_ylabel('Average percentage of MolProbity scores(%)', fontsize=30, fontweight='bold')
ax.set_xlabel('MolProbity score ranges', fontsize=30, fontweight='bold')
ax.set_xticks([x + bar_width / 2 * (len(directories) - 1) for x in range(len(labels))])
ax.set_xticklabels(labels, rotation=45)
plt.title('Distribution of average MolProbity scores (all models)', fontweight='bold', pad=25)
plt.legend()
plt.tight_layout()

# Save the figure using the directory specified in the config file
figures_directory = config["Supp_Figures_directory"]
output_filename = "Figure_S14_b.jpeg"
output_path = os.path.join(figures_directory, output_filename)
plt.savefig(output_path, dpi=1000, facecolor='white')
plt.show()
print(f"Figure saved at: {output_path}")


# In[ ]:


import matplotlib.pyplot as plt
import numpy as np
import os

# Data for the plot
categories = ['GT', 'AFM-TB', 'AFM-TF', 'Relaxed-TB', 'Relaxed-TF', 'AF3']
non_zero_percentages = [20.0, 10.0, 10.0, 10.0, 16.0, 30.0]
colors_zero = ['#00441b', '#006d2c', '#238b45', '#41ab5d', '#74c476', '#9be699']

# Setup for plotting
fig, ax = plt.subplots(figsize=(14, 16))  # Adjusted figure size for tighter bars
bar_width = 0.5  # Increase bar width to fill more space

# Calculate positions to bring bars closer together
spacing = 0.1  # Decrease spacing between bars
non_zero_bar_positions = np.arange(len(categories)) * (bar_width + spacing)

# Plot bars
bars_non_zero = ax.bar(non_zero_bar_positions, non_zero_percentages, color=colors_zero, edgecolor='black', width=bar_width, label='Non-Zero')

# Labels above bars with larger font size
for bar, percentage in zip(bars_non_zero, non_zero_percentages):
    ax.text(bar.get_x() + bar.get_width() / 2., bar.get_height(), f'{percentage:.1f}%', ha='center', va='bottom', fontsize=35)

# Set labels and legend
ax.set_ylabel('Percentage of cis non-proline samples(%)', fontsize=38, fontweight='bold')
ax.set_xlabel('Prediction methods', fontsize=38, fontweight='bold')
ax.set_xticks(non_zero_bar_positions)
ax.set_xticklabels(categories, rotation=45, fontsize=30, fontweight='bold')
ax.set_title('cis non-proline (top 5 models)', fontsize=40, fontweight='bold', pad=20)
plt.yticks(fontsize=33)
plt.tight_layout()

# Save the figure using the directory specified in the config file
output_directory = config["Supp_Figures_directory"]
output_filename = 'Figure_S14_C_1.jpeg'
output_path = os.path.join(output_directory, output_filename)
plt.savefig(output_path, format='jpeg', dpi=1000)
print(f"Figure saved to: {output_path}")

# Display the plot
plt.show()


# In[ ]:


import matplotlib.pyplot as plt
import numpy as np
import os

# Twisted-peptide on All samples
# Data
categories = ['GT', 'AFM-TB', 'AFM-TF', 'Relaxed-TB', 'Relaxed-TF', 'AF3']
non_zero_percentages = [0.0, 10.0, 0.0, 16.0, 14.0, 0.0]
colors_zero = ['#08306b', '#08519c', '#2171b5', '#4292c6', '#6baed6', '#9ecae1']

# Plotting setup
fig, ax = plt.subplots(figsize=(14, 16))  # Adjusted figure size for tighter bars
bar_width = 0.5  # Increase bar width to fill more space

# Calculate positions to bring bars closer together
spacing = 0.1  # Decrease spacing between bars
non_zero_bar_positions = np.arange(len(categories)) * (bar_width + spacing)

# Plot "Non-Zero" bars using the specified colors
bars_non_zero = ax.bar(non_zero_bar_positions, non_zero_percentages, color=colors_zero, 
                       edgecolor='black', width=bar_width, label='Non-Zero')

# Add percentages above bars with larger font size
for bar, percentage in zip(bars_non_zero, non_zero_percentages):
    ax.text(bar.get_x() + bar.get_width() / 2., bar.get_height(), f'{percentage:.1f}%', 
            ha='center', va='bottom', fontsize=35)

# Set labels, ticks, and title
ax.set_ylabel('Percentage of Twisted Peptides samples(%)', fontsize=38, fontweight='bold')
ax.set_xlabel('Prediction methods', fontsize=38, fontweight='bold')
ax.set_xticks(non_zero_bar_positions)
ax.set_xticklabels(categories, rotation=45, fontsize=30, fontweight='bold')
ax.set_title('Twisted Peptide (top 5 models)', fontsize=40, fontweight='bold', pad=20)
plt.yticks(fontsize=33)
plt.tight_layout()

# Save the figure using the directory specified in the config file
output_directory = config["Supp_Figures_directory"]
output_filename = 'Figure_S14_C_2.jpeg'
output_path = os.path.join(output_directory, output_filename)
plt.savefig(output_path, format='jpeg', dpi=1000)
print(f"Figure saved to: {output_path}")

# Display the plot
plt.show()


# In[ ]:


import matplotlib.pyplot as plt
import numpy as np
import os

# Twisted-peptide on All samples
# Data
categories = ['GT', 'AFM-TB', 'AFM-TF', 'Relaxed-TB', 'Relaxed-TF', 'AF3']
non_zero_percentages = [3.66, 13.78, 16.32, 1.39, 2.10, 4.32]
colors_zero = ['#4a1486', '#6a1b9a', '#8e24aa', '#ab47bc', '#ce93d8', '#f3e5f5']  # Using "Zero" colors for simplicity

# Plotting setup
fig, ax = plt.subplots(figsize=(14, 16))  # Adjusted figure size for tighter bars
bar_width = 0.5  # Increase bar width to fill more space

# Calculate positions to bring bars closer together
spacing = 0.1  # Decrease spacing between bars
non_zero_bar_positions = np.arange(len(categories)) * (bar_width + spacing)

# Plot "Non-Zero" bars using the specified colors
bars_non_zero = ax.bar(non_zero_bar_positions, non_zero_percentages, color=colors_zero, 
                       edgecolor='black', width=bar_width, label='Non-Zero')

# Add percentages above bars with larger font size
for bar, percentage in zip(bars_non_zero, non_zero_percentages):
    ax.text(bar.get_x() + bar.get_width() / 2., bar.get_height(), f'{percentage:.1f}%', 
            ha='center', va='bottom', fontsize=35)

# Set labels, ticks, and title
ax.set_ylabel('Percentage of Clashscores samples(%)', fontsize=38, fontweight='bold')
ax.set_xlabel('Prediction methods', fontsize=38, fontweight='bold')
ax.set_xticks(non_zero_bar_positions)
ax.set_xticklabels(categories, rotation=45, fontsize=30, fontweight='bold')
ax.set_title('Clashscores (Top 5 Models)', fontsize=40, fontweight='bold', pad=20)
plt.yticks(fontsize=33)
plt.tight_layout()

# Save the figure using the directory specified in the config file
output_directory = config["Supp_Figures_directory"]
output_filename = 'Figure_S14_C_3.jpeg'
output_path = os.path.join(output_directory, output_filename)
plt.savefig(output_path, format='jpeg', dpi=1000)
print(f"Figure saved to: {output_path}")

# Display the plot
plt.show()


# In[ ]:


import matplotlib.pyplot as plt
import numpy as np
import os

# Data for TB and TF of AlphaFold and ColabFold
categories = ['High', 'Medium', 'Acceptable', 'Incorrect']
# AlphaFold values
af3_values  = [58.0, 28.0, 10.0, 4.0]
af_tb_values = [40.0, 30.0, 20.0, 10.0]
af_tf_values = [50.0, 20.0, 20.0, 10.0]
# ColabFold values (Relaxed)
re_tb_values = [48.0, 28.0, 22.0, 2.0]
re_tf_values = [46.0, 28.0, 24.0, 2.0]

# Positions of the bar-groups
barWidth = 0.15  # Narrower bar width to fit all bars
r1 = np.arange(len(categories))               # Positions for the first set of bars
r2 = [x + barWidth for x in r1]                # Positions for the second set
r3 = [x + barWidth for x in r2]                # Positions for the third set
r4 = [x + barWidth for x in r3]                # Positions for the fourth set
r5 = [x + barWidth for x in r4]                # Positions for the fifth set

plt.figure(figsize=(10, 8), facecolor='white')  # Set figure background to white
plt.bar(r1, af3_values, color='#FF8E44', width=barWidth, edgecolor='grey', label='AlphaFold3')
plt.bar(r2, af_tb_values, color='#006400', width=barWidth, edgecolor='grey', label='AlphaFold Template_based')
plt.bar(r3, af_tf_values, color='#b2d8b2', width=barWidth, edgecolor='grey', label='AlphaFold Template_Free')
plt.bar(r4, re_tb_values, color='#800080', width=barWidth, edgecolor='grey', label='Relaxed Template_based')
plt.bar(r5, re_tf_values, color='#D8BFD8', width=barWidth, edgecolor='grey', label='Relaxed Template_Free')

# Add xticks on the middle of the group bars and adjust font sizes
plt.xticks([r + 1.5 * barWidth for r in range(len(categories))], categories, rotation=0, fontsize=22)
plt.ylabel('Percentage of DockQ scores(%)', fontsize=23)
plt.xlabel('DockQ score categories', fontsize=23)
plt.tick_params(axis='y', labelsize=18)  # Adjust y-axis font size

# Create legend and add title with larger font
plt.title('AFM, Relaxed-AFM, and AF3 DockQ scores (top 5 models)', fontsize=19, fontweight='bold', pad=20)

# Save the figure to the directory specified in the config file
output_directory = config["Supp_Figures_directory"]
output_filename = 'Figure_S14_A_1.jpeg'
output_path = os.path.join(output_directory, output_filename)
plt.savefig(output_path, format='jpeg', dpi=1000, facecolor='white')

plt.show()


# In[ ]:


import matplotlib.pyplot as plt
import numpy as np
import os

# Data for TB and TF of AlphaFold and ColabFold
categories = ['High', 'Medium', 'Acceptable', 'Incorrect']
# AlphaFold values
af3_values  = [60.0, 30.0, 0.0, 10.0]
af_tb_values = [40.0, 30.0, 20.0, 10.0]
af_tf_values = [50.0, 30.0, 20.0, 0.0]
# ColabFold values
re_tb_values = [50.0, 30.0, 20.0, 0.0]
re_tf_values = [40.0, 40.0, 10.0, 10.0]

# Positions of the bar-groups
barWidth = 0.15  # Narrower bar width to fit all bars
r1 = np.arange(len(categories))               # Positions for the first set of bars
r2 = [x + barWidth for x in r1]                # Positions for the second set
r3 = [x + barWidth for x in r2]                # Positions for the third set
r4 = [x + barWidth for x in r3]                # Positions for the fourth set
r5 = [x + barWidth for x in r4]                # Positions for the fifth set

plt.figure(figsize=(10, 8), facecolor='white')  # Set figure background to white
plt.bar(r1, af3_values, color='#FF8E44', width=barWidth, edgecolor='grey', label='AlphaFold3')
plt.bar(r2, af_tb_values, color='#006400', width=barWidth, edgecolor='grey', label='AlphaFold Template_based')
plt.bar(r3, af_tf_values, color='#b2d8b2', width=barWidth, edgecolor='grey', label='AlphaFold Template_Free')
plt.bar(r4, re_tb_values, color='#800080', width=barWidth, edgecolor='grey', label='Relaxed Template_based')
plt.bar(r5, re_tf_values, color='#D8BFD8', width=barWidth, edgecolor='grey', label='Relaxed Template_Free')

# Add xticks on the middle of the group bars and adjust font sizes
plt.xticks([r + 1.5 * barWidth for r in range(len(categories))], categories, rotation=0, fontsize=22)
plt.ylabel('Percentage of DockQ scores(%)', fontsize=23)
plt.xlabel('DockQ score categories', fontsize=23)
plt.tick_params(axis='y', labelsize=18)  # Adjust y-axis font size

# Create legend and add title with larger font
plt.title('AFM, Relaxed-AFM, and AF3 DockQ scores (first-ranked models)', fontsize=19, fontweight='bold', pad=20)

# Save the figure using the directory specified in the config file
output_directory = config["Supp_Figures_directory"]
output_filename = 'Figure_S14_A_2.jpeg'
output_path = os.path.join(output_directory, output_filename)
plt.savefig(output_path, format='jpeg', dpi=1000, facecolor='white')

plt.show()


# In[ ]:




