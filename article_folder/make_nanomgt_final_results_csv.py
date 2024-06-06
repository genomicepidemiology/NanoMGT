import os
import csv
from collections import defaultdict

# Directory where the main folders are located
base_dir = '/home/people/malhal/test/new_nanomgt_results'

# Header for the new CSV file
csv_header = ['tool', 'specie', 'maf', 'batch', 'precision', 'recall', 'f1score']

# File to save the consolidated results
output_file = '/home/people/malhal/test/new_nanomgt_results/nanomgt_final_results.csv'

# Dictionary to hold the data for averaging
data = defaultdict(list)

# Start processing each folder in the base directory
for main_folder in range(1, 6):
    folder_path = f"{base_dir}/final_nanomgt_performance_output_{main_folder}"
    maf = f"0.0{main_folder}"  # Construct the MAF based on folder index

    # Check if the directory exists
    if os.path.exists(folder_path):
        # Iterate over each species folder inside the main folder
        for species_folder in os.listdir(folder_path):
            species_path = os.path.join(folder_path, species_folder)

            # Extract batch number
            batch = species_folder.split('_')[-2].replace('batch', '')

            # Build the path to the target CSV file
            csv_file = os.path.join(species_path, 'top_result.csv')

            # Process the CSV if it exists
            if os.path.isfile(csv_file):
                with open(csv_file, mode='r') as infile:
                    reader = csv.DictReader(infile)
                    for row in reader:
                        # Check if parameters are empty or missing necessary values
                        if not row['Parameters'] or row['F1 Score'] == '':
                            precision = recall = f1score = 0
                        else:
                            precision = float(row['Precision'])
                            recall = float(row['Recall'])
                            f1score = float(row['F1 Score'])

                        # Fix the indexing issue for 'ecoli'
                        specie = species_folder.split('_')[1] if 'ecoli' not in species_folder else 'ecoli'

                        # Append data for averaging
                        data[(specie, maf, batch)].append((precision, recall, f1score))

# Write averaged results to file
with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(csv_header)

    for (specie, maf, batch), values in data.items():
        avg_precision = sum(v[0] for v in values) / len(values)
        avg_recall = sum(v[1] for v in values) / len(values)
        avg_f1score = sum(v[2] for v in values) / len(values)

        # Write the processed row to the output CSV
        writer.writerow(['NanoMGT', specie, maf, batch, avg_precision, avg_recall, avg_f1score])

print(f"Data has been consolidated and averaged into {output_file}.")
