import pandas as pd

# Load the CSV file into a DataFrame
df = pd.read_csv('csv_files/average_performance_results.csv')

# Define the MAF thresholds and their corresponding minimum batch values
maf_thresholds = {
    0.01: 1,
    0.02: 2,
    0.03: 3,
    0.04: 4,
    0.05: 5,
}

# Mapping of species short names to full names
species_mapping = {
    'ecoli': 'Escherichia coli',
    'staph': 'Staphylococcus aureus',
    'campylobacter': 'Campylobacter jejuni',
    'salmonella': 'Salmonella enterica',
    'klebsiella': 'Klebsiella pneumoniae'
}

# Initialize a list to store the results
results = []

# Iterate through each species
for short_name, full_name in species_mapping.items():
    # Filter the DataFrame based on the species
    specie_df = df[df['specie'] == short_name]

    # Iterate through each MAF threshold and calculate the average scores
    for maf, min_batch in maf_thresholds.items():
        # Filter the DataFrame based on MAF and batch values
        filtered_df = specie_df[(specie_df['maf'] == maf) & (specie_df['batch'] >= min_batch)]

        for tool in ['Confindr', 'NanoMGT']:
            # Filter the DataFrame based on the tool
            tool_df = filtered_df[filtered_df['tool'] == tool]

            if not tool_df.empty:
                # Calculate the average precision, recall, and f1score
                avg_precision = tool_df['precision'].mean()
                avg_recall = tool_df['recall'].mean()
                avg_f1score = tool_df['f1score'].mean()
            else:
                avg_precision = avg_recall = avg_f1score = None

            # Store the results in the list
            results.append([full_name, maf, tool, avg_precision, avg_recall, avg_f1score])

# Create a DataFrame from the results list
results_df = pd.DataFrame(results, columns=['Species', 'MAF', 'Tool', 'Precision', 'Recall', 'F1 Score'])

# Print the DataFrame as CSV
print(results_df.to_csv(index=False))

results_df.to_csv('csv_files/specie_maf_performance_for_plots.csv', index=False)

