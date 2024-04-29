import pandas as pd

# Load the data from both CSV files
confindr_data = pd.read_csv("csv_files/confindr_final_individual_results.csv")
nanomgt_data = pd.read_csv("csv_files/nanomgt_final_all_individual_results.csv")

# Combine the data from both tools into one DataFrame
combined_data = pd.concat([confindr_data, nanomgt_data], ignore_index=True)

# Define the relevant sequencing IDs for ecoli, salmonella, and staph aureus
ecoli_ids = ['SRR25689478', 'SRR26036455', 'ERR12533301']
salmonella_ids = ['SRR28399428', 'SRR27136088', 'SRR27755678']
staph_ids = ['SRR28370694', 'ERR8958843']

# Filter data for the specified depth, specie, and sequencing ID
filtered_data = combined_data[
    (combined_data['depth'] == 220) &
    (combined_data['sequencing_id'].isin(ecoli_ids + salmonella_ids + staph_ids)) &
    (combined_data['specie'].isin(['ecoli', 'salmonella_enterica', 'staph_aureus']))
]

# Rename species for uniformity
species_map = {
    'ecoli': 'Escherichia coli',
    'salmonella_enterica': 'Salmonella enterica',
    'staph_aureus': 'Staphylococcus aureus'
}
filtered_data['specie'] = filtered_data['specie'].map(species_map)

# Save the filtered data to a new CSV file
filtered_data.to_csv("csv_files/filtered_individual_co_occurence_noise.csv", index=False)

# Print the head of the filtered DataFrame to verify
print(filtered_data.head())
