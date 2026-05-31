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

# Define the directories for each category using paths from config
first_ranked_dockq_dir = config["First_ranked_DockQ"]
af3_data_dir = config["AF3_Data"]

# Just update this part to use the config paths
directories = {
    'AFM-TB': os.path.join(first_ranked_dockq_dir, 'AFM-TB.xlsx'),
    'AFM-TF': os.path.join(first_ranked_dockq_dir, 'AFM-TF.xlsx'),
    'CF-TB': os.path.join(first_ranked_dockq_dir, 'CF-TB.xlsx'),
    'CF-TF': os.path.join(first_ranked_dockq_dir, 'CF-TF.xlsx'),
    'AF3': os.path.join(first_ranked_dockq_dir, 'AF3.xlsx'),
}

# Print paths to verify
print("Using these file paths:")
for category, path in directories.items():
    print(f"{category}: {path}")
    print(f"  File exists: {os.path.exists(path)}")

# List to store individual DataFrames for each category
category_dfs = {cat: [] for cat in directories}

# Loop through all files in each directory
for category, file_path in directories.items():
    if os.path.exists(file_path):
        try:
            df = pd.read_excel(file_path)
            df['Category'] = category
            category_dfs[category] = df
            print(f"Successfully loaded data for {category}")
        except Exception as e:
            print(f"Error loading {file_path}: {e}")
    else:
        print(f"Warning: File not found: {file_path}")

# Combine DataFrames for each category
if any(len(df) > 0 for df in category_dfs.values()):
    # Only combine non-empty DataFrames
    dfs_to_combine = [df for df in category_dfs.values() if len(df) > 0]
    combined_df = pd.concat(dfs_to_combine, ignore_index=True)

    # Print a summary of the combined data
    print(f"\nCombined data shape: {combined_df.shape}")
    print(f"Categories present: {combined_df['Category'].unique()}")

    # Ensure the relevant columns are numeric
    parameters = ['Fnat', 'iRMS', 'LRMS', 'DockQ', 'IOU']
    for param in parameters:
        if param in combined_df.columns:
            combined_df[param] = pd.to_numeric(combined_df[param], errors='coerce')
        else:
            print(f"Warning: Column '{param}' not found in the data")

    # Define the new names for categories
    category_order_renamed = ['AFM-TB', 'AFM-TF', 'CF-TB', 'CF-TF', 'AF3']
    category_mapping = {
        'AFM-TB': 'AFM-TB',
        'AFM-TF': 'AFM-TF',
        'CF-TB': 'CF-TB',
        'CF-TF': 'CF-TF',
        'AF3': 'AF3'
    }

    # Updating the category names in the DataFrame
    combined_df['Category'] = combined_df['Category'].map(category_mapping)

    # Set font sizes for each part of the plot
    def plot_box_plots(output_directory, title_fontsize=22, xlabel_fontsize=20, ylabel_fontsize=20, xticks_fontsize=18, yticks_fontsize=20):
        plot_colors = ['lightgreen', 'lightblue', 'lightcoral', 'lightpink', 'lightgrey']

        # Create the output directory if it doesn't exist
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)
            print(f"Created output directory: {output_directory}")

        available_parameters = [p for p in parameters if p in combined_df.columns]
        if not available_parameters:
            print("No parameter columns found in the data. Cannot create plots.")
            return

        # Plotting boxplots for each parameter
        for idx, param in enumerate(available_parameters, start=1):
            color = plot_colors[idx % len(plot_colors)]  # Get color safely

            plt.figure(figsize=(12, 8))

            # Drop rows with NaN values in the parameter column
            df_to_plot = combined_df.dropna(subset=[param])

            if df_to_plot.empty:
                print(f"No valid data for parameter '{param}'. Skipping this plot.")
                plt.close()
                continue

            # Get available categories in the filtered data
            available_categories = df_to_plot['Category'].unique()
            filtered_order = [cat for cat in category_order_renamed if cat in available_categories]

            if not filtered_order:
                print(f"No categories with valid data for parameter '{param}'. Skipping this plot.")
                plt.close()
                continue

            # Creating the boxplot with specified order and new names
            box = sns.boxplot(x='Category', y=param, data=df_to_plot, order=filtered_order, palette=[color])

            # Change the color of the median line to red
            for line in box.lines:
                # Lines 4n+2 (n=0,1,2,...) are the median lines
                if box.lines.index(line) % 6 == 4:
                    line.set_color('red')
                    line.set_linewidth(2)

            # Set appropriate labels for iRMSD and LRMSD
            ylabel = 'iRMSD' if param == 'iRMS' else 'LRMSD' if param == 'LRMS' else param

            plt.title(f'Distribution of {ylabel} (first-ranked models)', fontsize=title_fontsize)
            plt.xlabel('', fontsize=0)  # Hide the x-axis label
            plt.ylabel(ylabel, fontsize=ylabel_fontsize)
            plt.xticks(fontsize=xticks_fontsize, rotation=45)
            plt.yticks(fontsize=yticks_fontsize)
            plt.grid(True)

            # Save the plot as a JPEG image with name ending in Fig_S1 to Fig_S5
            plot_path = os.path.join(output_directory, f'{param}_Fig_S{idx}.jpeg')
            plt.savefig(plot_path, dpi=1000, format='jpeg', bbox_inches='tight')
            print(f"Saved plot to {plot_path}")

            # Display the plot
            plt.show()

            # Print statistics for each category and parameter
            print(f"Statistics for {param}:")
            for category in filtered_order:
                data = df_to_plot[df_to_plot['Category'] == category][param].dropna()
                if len(data) > 0:
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
                else:
                    print(f"{category}: No data available")
            print("\n")

    # Usage example
    output_directory = config["Supp_Figures_directory"]
    print(f"Output directory: {output_directory}")
    plot_box_plots(output_directory, title_fontsize=22, xlabel_fontsize=20, ylabel_fontsize=20, xticks_fontsize=18, yticks_fontsize=20)
else:
    print("No data could be loaded from any file. Please check the file paths.")


# In[ ]:


import os
import pandas as pd
import plotly.graph_objects as go
import json

# Load configuration and define a global variable
def load_config():
    global config
    with open('config_Distributions.json') as config_file:
        config = json.load(config_file)

# Call the function to load the configuration
load_config()

# Define the directory paths from the configuration file
directory_path = config["First_ranked_DockQ"]
output_directory = config["Main_Figures_directory"]

# Define the file names, including the new "GT" dataset
file_names = ['AFM-TB.xlsx', 'AFM-TF.xlsx', 'AF3.xlsx', 'CF-TB.xlsx', 'CF-TF.xlsx', 'GT']

# Load the data from the first sheet of each file (excluding the GT dataset)
dfs = {file_name: pd.read_excel(os.path.join(directory_path, file_name), sheet_name='Sheet1') for file_name in file_names[:-1]}  # Exclude GT

# Select the columns for the radar plot
columns_for_radar = ['Fnat', 'iRMS', 'LRMS', 'IOU', 'DockQ']

# Replace 'Not Found' with NaN and convert to numeric for all datasets
for df in dfs.values():
    for column in columns_for_radar:
        df[column] = pd.to_numeric(df[column], errors='coerce')

# Calculate the median values for all datasets
median_values_actual = {file_name: df[columns_for_radar].median() for file_name, df in dfs.items()}

# Add the GT dataset with ideal values directly
median_values_actual['GT'] = pd.Series({'Fnat': 1, 'iRMS': 0, 'LRMS': 0, 'IOU': 1, 'DockQ': 1})

# Create a DataFrame to hold the median values
median_df = pd.DataFrame(median_values_actual).T

# Initialize the scaled DataFrame
scaled_median_df = pd.DataFrame(index=median_df.index, columns=columns_for_radar)

# Scale the median values for each parameter
scale_min = 0.25
scale_max = 1.0

for column in columns_for_radar:
    min_val = median_df[column].min()
    max_val = median_df[column].max()

    if column in ['iRMS', 'LRMS']:
        # For iRMSD and LRMSD, lower is better
        scaled_median_df[column] = scale_max - (median_df[column] - min_val) / (max_val - min_val) * (scale_max - scale_min)
    else:
        # For DockQ, IOU, and Fnat, higher is better
        scaled_median_df[column] = scale_min + (median_df[column] - min_val) / (max_val - min_val) * (scale_max - scale_min)

# The GT dataset should be 1.0 across all parameters after scaling
scaled_median_df.loc['GT'] = [1.0] * len(columns_for_radar)

# Reorder the datasets by their scaled values
reordered_data = {name: [0]*len(columns_for_radar) for name in scaled_median_df.index}
for i, category in enumerate(columns_for_radar):
    sorted_datasets = scaled_median_df[category].sort_values(ascending=False)
    for dataset_name in sorted_datasets.index:
        reordered_data[dataset_name][i] = sorted_datasets[dataset_name]

# Prepare the categories and values for the radar plot
categories = ['Fnat', 'iRMSD', 'LRMSD', 'IOU', 'DockQ']

# Create the radar plot without the legend
fig = go.Figure()

# Define colors for each dataset, ensuring GT has maximum contrast (black)
colors = {
    'AFM-TB.xlsx': 'rgba(31, 119, 180, 1.0)',   # Blue (Tableau blue)
    'AFM-TF.xlsx': 'rgba(255, 127, 14, 1.0)',     # Orange (Tableau orange)
    'AF3.xlsx': 'rgba(44, 160, 44, 1.0)',         # Green (Tableau green)
    'CF-TB.xlsx': 'rgba(214, 39, 40, 1.0)',       # Red (Tableau red)
    'CF-TF.xlsx': 'rgba(148, 103, 189, 1.0)',      # Purple (Tableau purple)
    'GT': 'rgba(0, 0, 0, 1.0)'                    # Black for GT dataset
}

# Add the datasets to the radar plot without filled lines
for dataset_name, values in reordered_data.items():
    fig.add_trace(go.Scatterpolar(
        r=values + [values[0]],  # Complete the loop
        theta=categories + [categories[0]],  # Complete the loop
        fill=None,  # No filled lines
        opacity=1.0,  # Fully visible lines
        name=dataset_name,
        line=dict(color=colors[dataset_name], width=3)  # Thicker lines
    ))

# Update layout with larger angular axis label font size and set label color to black
fig.update_layout(
    polar=dict(
        radialaxis=dict(
            visible=True,
            range=[0, 1],  # Set the range for scaled values
            showticklabels=False  # Hide the tick labels on the radial axis
        ),
        angularaxis=dict(
            tickfont=dict(size=20, color='black')  # Larger font size and black color for parameter labels
        )
    ),
    showlegend=False,  # Hide the legend in this plot
    title=dict(
        text='Median Values of DockQ Parameters',
        font=dict(size=22, color='black'),  # Title font settings
        xanchor='center',  # Center the title horizontally
        yanchor='top',
        x=0.5
    )
)

# Show the radar plot without the legend
fig.show()

# Save the radar plot without the legend as an image
output_image_path = os.path.join(output_directory, 'Fig_3a.jpeg')
fig.write_image(output_image_path, format='png', scale=10, engine='kaleido')

print(f'Radar plot without legend saved at {output_image_path}')

# Create a separate plot just for the legend
'''legend_fig = go.Figure()

# Add empty traces with the correct colors and names to build the legend
for dataset_name, color in colors.items():
    legend_fig.add_trace(go.Scatterpolar(
        r=[None],  # Empty trace
        theta=[None],  # Empty trace
        name=dataset_name,
        line=dict(color=color, width=3)
    ))

legend_fig.update_layout(
    showlegend=True,
    legend=dict(
        font=dict(size=18, color='black'),  # Legend font settings (black color)
        orientation="h",  # Horizontal legend
        xanchor='center',
        x=0.5,
        yanchor='bottom',
        y=1.02
    ),
    margin=dict(l=0, r=0, t=0, b=0),  # Remove margins for the legend plot
    polar=dict(
        radialaxis=dict(visible=False),  # Hide the polar grid in legend
        angularaxis=dict(visible=False)  # Hide the angular labels in legend
    )
)'''

# (Optional) Save the legend as a separate image if needed
# legend_image_path = os.path.join(output_directory, 'Legend.png')
# legend_fig.write_image(legend_image_path, format='png', scale=10, engine='kaleido')
# print(f'Legend saved separately at {legend_image_path}')


# In[ ]:


import os
import pandas as pd
import plotly.graph_objects as go
import json

# Load configuration and define a global variable
def load_config():
    global config
    with open('config_Distributions.json') as config_file:
        config = json.load(config_file)

# Call the function to load the configuration
load_config()

# Define the directory paths from the configuration file
directory_path = config["First_ranked_DockQ"]
output_directory = config["Supp_Figures_directory"]

# Define the file names, including the new "GT" dataset
file_names = ['AFM-TB.xlsx', 'AFM-TF.xlsx', 'AF3.xlsx', 'CF-TB.xlsx', 'CF-TF.xlsx', 'GT']

# Load the data from the first sheet of each file (excluding the GT dataset)
dfs = {file_name: pd.read_excel(os.path.join(directory_path, file_name), sheet_name='Sheet1') 
       for file_name in file_names[:-1]}  # Exclude GT

# Select the columns for the radar plot
columns_for_radar = ['Fnat', 'iRMS', 'LRMS', 'IOU', 'DockQ']

# Replace 'Not Found' with NaN and convert to numeric for all datasets
for df in dfs.values():
    for column in columns_for_radar:
        df[column] = pd.to_numeric(df[column], errors='coerce')

# Calculate the median values for all datasets
median_values_actual = {file_name: df[columns_for_radar].median() for file_name, df in dfs.items()}

# Add the GT dataset with ideal values directly
median_values_actual['GT'] = pd.Series({'Fnat': 1, 'iRMS': 0, 'LRMS': 0, 'IOU': 1, 'DockQ': 1})

# Create a DataFrame to hold the median values
median_df = pd.DataFrame(median_values_actual).T

# Initialize the scaled DataFrame
scaled_median_df = pd.DataFrame(index=median_df.index, columns=columns_for_radar)

# Scale the median values for each parameter
scale_min = 0.25
scale_max = 1.0

for column in columns_for_radar:
    min_val = median_df[column].min()
    max_val = median_df[column].max()

    if column in ['iRMS', 'LRMS']:
        # For iRMS and LRMS, lower is better
        scaled_median_df[column] = scale_max - (median_df[column] - min_val) / (max_val - min_val) * (scale_max - scale_min)
    else:
        # For DockQ, IOU, and Fnat, higher is better
        scaled_median_df[column] = scale_min + (median_df[column] - min_val) / (max_val - min_val) * (scale_max - scale_min)

# The GT dataset should be 1.0 across all parameters after scaling
scaled_median_df.loc['GT'] = [1.0] * len(columns_for_radar)

# Reorder the datasets by their scaled values (this step uses the original order)
reordered_data = {name: [0]*len(columns_for_radar) for name in scaled_median_df.index}
for i, category in enumerate(columns_for_radar):
    sorted_datasets = scaled_median_df[category].sort_values(ascending=False)
    for dataset_name in sorted_datasets.index:
        reordered_data[dataset_name][i] = sorted_datasets[dataset_name]

# Prepare the categories and values for the radar plot
categories = ['Fnat', 'iRMSD', 'LRMSD', 'IOU', 'DockQ']

# --- New mapping for ordering and colors ---
# Define a mapping from file keys to simplified names (removing '.xlsx')
name_mapping = {
    'AFM-TB.xlsx': 'AFM-TB',
    'AFM-TF.xlsx': 'AFM-TF',
    'AF3.xlsx': 'AF3',
    'CF-TB.xlsx': 'CF-TB',
    'CF-TF.xlsx': 'CF-TF',
    'GT': 'GT'
}
# Create a reverse mapping (simplified name -> file key)
reverse_mapping = {v: k for k, v in name_mapping.items()}

# Define the new order and new color codes with the updated color codes:
new_order = ['GT', 'AFM-TB', 'AFM-TF', 'CF-TB', 'CF-TF', 'AF3']
colors_zero = ['#08306b', '#08519c', '#2171b5', '#4292c6', '#6baed6', '#9ecae1']
new_colors = dict(zip(new_order, colors_zero))

# Create the radar plot
fig = go.Figure()

# Loop over the new_order to add traces in the desired order with the new colors.
for dataset in new_order:
    file_key = reverse_mapping[dataset]  # Get the original file key
    values = reordered_data[file_key]
    fig.add_trace(go.Scatterpolar(
         r=values + [values[0]],  # Complete the loop
         theta=categories + [categories[0]],  # Complete the loop
         fill=None,  # No filled lines
         opacity=1.0,  # Fully visible lines
         name=dataset,
         line=dict(color=new_colors[dataset], width=3)
    ))

fig.update_layout(
    paper_bgcolor='white',   # Set white background for the paper
    plot_bgcolor='white',    # Set white background for the plot area
    polar=dict(
        radialaxis=dict(
            visible=True,
            range=[0, 1],  # Set the range for scaled values
            showticklabels=False  # Hide the tick labels
        ),
        angularaxis=dict(
            tickfont=dict(size=14, color='black')  # Font settings for the angular labels
        )
    ),
    showlegend=False,  # Hide the legend in this plot
    title=dict(
        text='Median Values of DockQ Parameters',
        font=dict(size=22, color='black'),
        xanchor='center',
        yanchor='top',
        x=0.5
    )
)

# Save the radar plot without the legend as an image
output_image_path = os.path.join(output_directory, 'Fig_S6_no_legend.jpeg')
fig.write_image(output_image_path, format='png', scale=10, engine='kaleido')

print(f'Radar plot without legend saved at {output_image_path}')

# Create a separate figure for the horizontal legend with enough space
legend_fig = go.Figure()

# Add traces with the correct colors and labels to build the legend
for dataset in new_order:
    legend_fig.add_trace(go.Scatter(
        x=[None],
        y=[None],
        mode='lines',
        line=dict(color=new_colors[dataset], width=3),
        name=dataset
    ))

legend_fig.update_layout(
    paper_bgcolor='white',
    plot_bgcolor='white',
    showlegend=True,
    width=1600,
    height=250,
    legend=dict(
        font=dict(size=22, color='black'),
        orientation="h",
        xanchor='center',
        x=0.5,
        yanchor='bottom',
        y=0.5,
        traceorder='normal'
    ),
    margin=dict(l=0, r=0, t=0, b=0),
    xaxis=dict(visible=False),
    yaxis=dict(visible=False)
)

legend_image_path = os.path.join(output_directory, 'Fig_S6_horizontal_legend.png')
legend_fig.write_image(legend_image_path, format='png', scale=10, engine='kaleido')

print(f'Legend saved separately at {legend_image_path}')


# In[ ]:




