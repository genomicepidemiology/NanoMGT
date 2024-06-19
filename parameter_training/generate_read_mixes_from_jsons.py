import os
import json
import pandas as pd
import random

# Load genome sizes from JSON
with open('genome_sizes.json', 'r') as f:
    genome_sizes = json.load(f)

# List of JSON files to load
json_files = [
    'campylobacter_jejuni_validation.json',
    'escherichia_coli_validation.json',
    'klebsiella_pneumoniae_validation.json',
    'listeria_monocytogenes_training.json',
    'salmonella_enterica_training.json',
    'staphylococcus_aureus_training.json'
]

# Initialize dictionaries to store sequencing IDs and batch information
species_sequencing_ids = {}
batch_info = {}

# Load data from JSON files
for file in json_files:
    species_name = file.replace('_validation.json', '').replace('_training.json', '').replace('_', ' ')
    with open(file, 'r') as f:
        data = json.load(f)
        species_sequencing_ids[species_name] = []
        for item in data:
            species_sequencing_ids[species_name].extend(item.values())
        batch_info[species_name] = data

# Path to the folder containing FASTQ files
fastq_folder_path = '/home/people/malhal/data/new_nanomgt/sup_data'  # Replace this with the actual path

# Function to calculate the average read length for a given FASTQ file
def calculate_average_read_length(file_path):
    total_length = 0
    total_reads = 0
    with open(file_path, 'r') as f:
        for i, line in enumerate(f):
            if i % 4 == 1:  # Sequence lines are every 4th line starting from the second line
                total_length += len(line.strip())
                total_reads += 1
    if total_reads == 0:
        return 0
    return total_length / total_reads

# Merge function definition
def merge_files(batch_files, output_directory, species_name, batch_id):
    output_file_path = os.path.join(output_directory, f'{species_name}_{batch_id}.fastq')
    with open(output_file_path, 'w') as outfile:
        for fname in batch_files:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
    print(f"Merged file created: {output_file_path}")

# Function to count nucleotides in a FASTQ file
def count_nucleotides_in_fastq(file_path):
    nucleotide_count = 0
    with open(file_path, 'r') as fastq_file:
        for i, line in enumerate(fastq_file):
            if i % 4 == 1:  # Sequence lines are every 4th line starting from the second line
                nucleotide_count += len(line.strip())
    return nucleotide_count

# Function to simulate batches
def simulate_batches(species, ids, genome_size, batch_info, output_directory_base, num_batches=10, sample_complexity=3):
    random.seed(42)  # For reproducibility

    # Adjust sample complexity for species with fewer IDs
    if len(ids) <= 2:
        sample_complexity = 2

    # Convert species name to lowercase with underscores for file naming
    species_name_for_files = species.lower().replace(' ', '_')

    # Iterate over different depths for each organism
    depths = [220]
    for depth in depths:
        output_directory = os.path.join(output_directory_base, f'depth{depth}_{species_name_for_files}')
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

        output_prefix = f'depth{depth}_{species_name_for_files}'
        desired_total_nucleotides = genome_size * depth

        nucleotide_distribution_summary = []  # Collect nucleotide distribution information

        for batch_info_item in batch_info:
            batch_id = batch_info_item['batch_id']
            minority_abundance = batch_info_item['minority_abundance']

            selected_ids = []
            percentages = []

            # Prepare the IDs and percentages for the current batch
            for group, ids_dict in minority_abundance.items():
                for id, percentage in ids_dict.items():
                    selected_ids.append(id)
                    percentages.append(percentage)

            subsampled_files = []
            for id, percentage in zip(selected_ids, percentages):
                num_nucleotides = int(desired_total_nucleotides * percentage / 100)
                input_file = os.path.join(fastq_folder_path, f'{id.strip()}.fastq')
                output_file = os.path.join(output_directory, f'{output_prefix}_batch{batch_id}_{id.strip()}.fastq')

                # Check if input file exists before processing
                if not os.path.exists(input_file):
                    print(f"Error: Input file {input_file} does not exist.")
                    continue

                # Calculate the average read length for this file
                avg_read_length = calculate_average_read_length(input_file)
                num_reads = int(num_nucleotides // avg_read_length)
                command = f'seqtk sample -s100 {input_file} {num_reads} > {output_file}'
                print(f"Running command: {command}")
                os.system(command)
                subsampled_files.append(output_file)

            # Merge the subsampled files for this batch
            merge_files(subsampled_files, output_directory, species_name_for_files, batch_id)

            # Count nucleotides in the merged file
            merged_file_path = os.path.join(output_directory, f'{species_name_for_files}_{batch_id}.fastq')
            merged_nucleotide_count = count_nucleotides_in_fastq(merged_file_path)

            # Avoid division by zero
            if merged_nucleotide_count == 0:
                print(f"Warning: Merged nucleotide count is zero for {merged_file_path}. Skipping.")
                continue

            # Collect information for nucleotide distribution summary
            for subsampled_file in subsampled_files:
                file_nucleotide_count = count_nucleotides_in_fastq(subsampled_file)
                percentage_of_total = (file_nucleotide_count / merged_nucleotide_count) * 100
                nucleotide_distribution_summary.append(f'{subsampled_file}: {percentage_of_total:.2f}% of nucleotides')

        print(f"Subsampling, merging, and nucleotide counting complete for {output_prefix}.")

        # Print nucleotide distribution summary
        with open(f'{output_directory}/{output_prefix}_nucleotide_distribution_summary.txt', 'w') as f:
            f.write("Nucleotide Distribution Summary:\n")
            for summary in nucleotide_distribution_summary:
                f.write(f"{summary}\n")

print(genome_sizes)
for species, ids in species_sequencing_ids.items():
    formatted_species = species[0].upper() + species[1:]
    genome_size = genome_sizes.get(formatted_species, None)
    if genome_size:
        # Create output directory based on JSON file name
        for item in json_files:
            if item.startswith(species.split(' ')[0]):
                json_file_name = item
                break
        output_directory_base = 'simulated_' + json_file_name.replace('.json', '')
        simulate_batches(species, ids, genome_size, batch_info[species], output_directory_base)

print("Simulation complete.")
