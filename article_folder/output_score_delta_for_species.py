import pandas as pd

# Load the CSV file
data = pd.read_csv("csv_files/average_performance_by_specie.csv")

# Convert MRD to float if it's read as string
data['mrd'] = data['mrd'].astype(float)

# Filter rows where the batch number is greater than 100 times the MRD value
data['batch'] = data['batch'].astype(int)
data = data[data['batch'] > data['mrd'] * 100]

# Define MRD values of interest
mrd_values = [0.01, 0.02, 0.03, 0.04, 0.05]

# Calculate differences for each metric for each specie and each MRD
results = []
for specie in data['specie'].unique():
    for mrd in mrd_values:
        # Filter the data for the current specie and MRD
        filtered_data = data[(data['specie'] == specie) & (data['mrd'] == mrd)]

        if not filtered_data.empty:
            # Separate the data by tool
            confindr_data = filtered_data[filtered_data['tool'] == 'Confindr']
            nanomgt_data = filtered_data[filtered_data['tool'] == 'NanoMGT']

            # Calculate differences if both tools have data
            if not confindr_data.empty and not nanomgt_data.empty:
                precision_diff = confindr_data['precision'].mean() - nanomgt_data['precision'].mean()
                recall_diff = confindr_data['recall'].mean() - nanomgt_data['recall'].mean()
                f1_diff = confindr_data['f1score'].mean() - nanomgt_data['f1score'].mean()
                results.append((specie, mrd, precision_diff, recall_diff, f1_diff))
            else:
                results.append((specie, mrd, 'Data not available for comparison', 'Data not available for comparison',
                                'Data not available for comparison'))
        else:
            results.append((specie, mrd, 'No data for this MRD', 'No data for this MRD', 'No data for this MRD'))

# Print results
for result in results:
    print(
        f"Specie: {result[0]}, MRD: {result[1]}, Precision Difference: {result[2]}, Recall Difference: {result[3]}, F1-Score Difference: {result[4]}")
