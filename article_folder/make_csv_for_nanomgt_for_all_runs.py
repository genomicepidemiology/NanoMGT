import os
import csv

# Directory where the main folders are located
base_dir = '/home/people/malhal/test/new_nanomgt_results'

# Header for the new CSV file
csv_header = ['tool', 'depth', 'specie', 'sequencing_id', 'maf', 'batch', 'precision', 'recall', 'f1score']

# File to save the consolidated results
output_file = '/home/people/malhal/test/new_nanomgt_results/nanomgt_final_all_individual_results.csv'

# Start processing each folder in the base directory
results = []

for main_folder in range(1, 6):
    folder_path = f"{base_dir}/final_nanomgt_performance_output_{main_folder}"
    maf = f"0.0{main_folder}"  # Construct the MAF based on folder index

    # Check if the directory exists
    if os.path.exists(folder_path):
        # Iterate over each species folder inside the main folder
        for species_folder in os.listdir(folder_path):
            species_path = os.path.join(folder_path, species_folder)

            # Extract the depth, batch from the folder name
            parts = species_folder.split('_')
            depth = parts[0].replace('depth', '')
            batch = parts[-2].replace('batch', '')

            # Determine the sequencing ID based on species
            if 'ecoli' in species_folder:
                sequencing_id = parts[2]
                specie = 'ecoli'
            else:
                sequencing_id = parts[3]
                specie = '_'.join(parts[1:3])  # Joining parts to get complete species name if not 'ecoli'

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

                        # Store each record directly
                        results.append([
                            'NanoMGT', depth, specie, sequencing_id, maf, batch,
                            precision, recall, f1score
                        ])

# Write results to file
with open(output_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(csv_header)
    writer.writerows(results)

print(f"Data has been consolidated into {output_file}.")
