import pandas as pd

# Load the data from both CSV files
confindr_data = pd.read_csv("csv_files/confindr_final_results.csv")
nanomgt_data = pd.read_csv("csv_files/nanomgt_final_results.csv")

# Add a column to each dataframe to identify the tool before combining
confindr_data['tool'] = 'Confindr'
nanomgt_data['tool'] = 'NanoMGT'

# Combine the data from both tools into one DataFrame
combined_data = pd.concat([confindr_data, nanomgt_data], ignore_index=True)

# Calculate the averages of precision, recall, and f1score based on maf, batch, depth, and tool
average_performance = combined_data.groupby(['tool', 'maf', 'batch', 'depth']).agg({
    'precision': 'mean',
    'recall': 'mean',
    'f1score': 'mean'
}).reset_index()

# Write the results to a new CSV file with the specified headers
average_performance.to_csv("csv_files/average_performance_by_depth.csv", index=False)

# Print out the head of the average_performance DataFrame to verify
print(average_performance.head())
