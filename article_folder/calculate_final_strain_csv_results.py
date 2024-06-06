import pandas as pd

# Load the datasets
confindr_data = pd.read_csv("confindr_final_results.csv")
nanomgt_data = pd.read_csv("nanomgt_final_results.csv")

# Combine the datasets
data = pd.concat([confindr_data, nanomgt_data])

# Define the group mapping
species_to_group = {
    'ecoli': 'group 1',
    'salmonella': 'group 1',
    'klebsiella': 'group 2',
    'staph': 'group 2',
    'campylobacter': 'group 2'
}

# Map the species to the group
data['group'] = data['specie'].map(species_to_group)

# Group the data by the new groups along with tool, maf, batch, and depth
grouped_data = data.groupby(['group', 'tool', 'maf', 'depth', 'batch']).agg({
    'precision': 'mean',
    'recall': 'mean',
    'f1score': 'mean'
}).reset_index()

# Output the result to a new CSV file
grouped_data.to_csv("grouped_final_results.csv", index=False)

print("Grouped data has been saved to 'grouped_final_results.csv'.")
