import pandas as pd

# Load the data from both CSV files
confindr_data = pd.read_csv("csv_files/confindr_final_results.csv")
nanomgt_data = pd.read_csv("csv_files/nanomgt_final_results.csv")

# Combine the data from both tools into one DataFrame
combined_data = pd.concat([confindr_data, nanomgt_data], ignore_index=True)

# Calculate the averages of precision, recall, and f1score
average_performance = combined_data.groupby(['tool', 'specie', 'maf', 'batch']).agg({
    'precision': 'mean',
    'recall': 'mean',
    'f1score': 'mean'
}).reset_index()

# Write the results to a new CSV file
average_performance.to_csv("csv_files/average_performance_results.csv", index=False)

# Print out the head of the average_performance DataFrame to verify
print(average_performance.head())
