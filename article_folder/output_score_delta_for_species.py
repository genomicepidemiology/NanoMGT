import pandas as pd

# Load the CSV file
data = pd.read_csv("csv_files/average_performance_by_specie.csv")

# Convert MAF to float if it's read as string
data['maf'] = data['maf'].astype(float)

# Filter rows where the batch number is greater than 100 times the MAF value
data['batch'] = data['batch'].astype(int)
data = data[data['batch'] > data['maf'] * 100]

# Define MAF values of interest
maf_values = [0.01, 0.02, 0.03, 0.04, 0.05]

# Calculate differences for each metric for each specie and each MAF
results = []
for specie in data['specie'].unique():
    for maf in maf_values:
        # Filter the data for the current specie and MAF
        filtered_data = data[(data['specie'] == specie) & (data['maf'] == maf)]

        if not filtered_data.empty:
            # Separate the data by tool
            confindr_data = filtered_data[filtered_data['tool'] == 'Confindr']
            nanomgt_data = filtered_data[filtered_data['tool'] == 'NanoMGT']

            # Calculate differences if both tools have data
            if not confindr_data.empty and not nanomgt_data.empty:
                precision_diff = confindr_data['precision'].mean() - nanomgt_data['precision'].mean()
                recall_diff = confindr_data['recall'].mean() - nanomgt_data['recall'].mean()
                f1_diff = confindr_data['f1score'].mean() - nanomgt_data['f1score'].mean()
                results.append((specie, maf, precision_diff, recall_diff, f1_diff))
            else:
                results.append((specie, maf, 'Data not available for comparison', 'Data not available for comparison',
                                'Data not available for comparison'))
        else:
            results.append((specie, maf, 'No data for this MAF', 'No data for this MAF', 'No data for this MAF'))

# Print results
for result in results:
    print(
        f"Specie: {result[0]}, MAF: {result[1]}, Precision Difference: {result[2]}, Recall Difference: {result[3]}, F1-Score Difference: {result[4]}")
