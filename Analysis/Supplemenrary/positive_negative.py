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


import matplotlib.pyplot as plt
import seaborn as sns
import json
import os

# Load configuration and define a global variable
def load_config():
    global config
    with open('config_Distributions.json') as config_file:
        config = json.load(config_file)

# Call the function to load the configuration
load_config()

# Data for plotting
software_names = ['AFM-Score', 'PyRosetta', 'FoldX-Stability', 'FoldX-Interaction', 
                  'HADDOCK-emscore', 'HADDOCK-mdscore', 'GNN_Dove', 'Vina', 
                  'Vinardo', 'DeepRank-GNN-esm']
percentages = [55.00, 51.67, 49.15, 38.98, 60.00, 68.33, 55.93, 39.29, 44.64, 70.0]

# Creating the bar chart with percentage values above each bar, all bars have the same color, and larger font size for labels
plt.figure(figsize=(12, 8), facecolor='white')  # Set the background to white
bars = sns.barplot(x=software_names, y=percentages, color="blue")  # All bars with the same color

# Adding percentage values above each bar with larger font size
for bar in bars.patches:
    label = f'{bar.get_height():.2f}%'    
    plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), label, 
             ha='center', va='bottom', fontsize=15)

plt.title('Positive Spearman correlations in TB models', fontsize=23, fontweight='bold', pad=20)
plt.ylabel('Percentage of positive correlation (%)', fontsize=18, fontweight='bold')
plt.xlabel('Methods', fontsize=18, fontweight='bold')
plt.yticks(fontsize=16)
plt.xticks(rotation=90, fontsize=14)  # Rotate labels for better readability
plt.ylim(0, 100)  # Space above bars for labels
plt.grid(axis='y')

# Define the desired save path using the config file directory
desired_save_path = os.path.join(config["Supp_Figures_directory"], "Figure_S16_a.jpeg")

# Ensure the directory exists
os.makedirs(os.path.dirname(desired_save_path), exist_ok=True)

# Save the figure in JPEG format with 1000 DPI
plt.savefig(desired_save_path, format='jpeg', dpi=1000, bbox_inches='tight', facecolor='white')

plt.show()


# In[ ]:


import matplotlib.pyplot as plt
import seaborn as sns
import json
import os

# Load configuration and define a global variable
def load_config():
    global config
    with open('config_Distributions.json') as config_file:
        config = json.load(config_file)

# Call the function to load the configuration
load_config()

# Data for plotting
software_names = ['AFM-Score', 'PyRosetta', 'FoldX-Stability', 'FoldX-Interaction', 
                  'HADDOCK-emscore', 'HADDOCK-mdscore', 'GNN_Dove', 'Vina', 
                  'Vinardo', 'DeepRank-GNN-esm']
percentages = [63.33, 56.67, 57.63, 52.54, 56.67, 56.67, 67.8, 50, 41.07, 60.00]

# Creating the bar chart with percentage values above each bar,
# all bars have the same color, and larger font size for labels
plt.figure(figsize=(12, 8), facecolor='white')  # Set the background to white
bars = sns.barplot(x=software_names, y=percentages, color="brown")  # All bars with the same color

# Adding percentage values above each bar with larger font size
for bar in bars.patches:
    label = f'{bar.get_height():.2f}%'
    plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), label, 
             ha='center', va='bottom', fontsize=15)

plt.title('Positive Spearman correlations in TF models', fontsize=23, fontweight='bold', pad=20)
plt.ylabel('Percentage (%)', fontsize=18, fontweight='bold')
plt.yticks(fontsize=16)
plt.xticks(rotation=90, fontsize=14)  # Rotate labels for better readability
plt.ylim(0, 100)  # Ensure space above bars for the labels
plt.grid(axis='y')

# Define the desired path for saving the plot in JPEG format using the config file directory
desired_save_path = os.path.join(config["Supp_Figures_directory"], "Figure_S16_b.jpeg")

# Ensure the directory exists
os.makedirs(os.path.dirname(desired_save_path), exist_ok=True)

# Save the figure in JPEG format with 1000 DPI
plt.savefig(desired_save_path, format='jpeg', dpi=1000, bbox_inches='tight', facecolor='white')

plt.show()


# In[ ]:




