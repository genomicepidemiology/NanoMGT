import os
import sys
import json
import random
import subprocess

# Define the directory containing the JSON files
json_directory = '/home/people/malhal/data/new_nanomgt/mixed_data'
read_dir = '/home/people/malhal/data/new_nanomgt/mixed_data'

# Automatically find all training and validation JSON files
json_files = [f for f in os.listdir(json_directory) if f.endswith('training.json') or f.endswith('validation.json')]

# Load genome sizes from JSON
with open('genome_sizes.json', 'r') as f:
    genome_sizes = json.load(f)

# Initialize dictionaries to store sequencing IDs and batch information
species_sequencing_ids = {}
batch_info = {"training": {}, "validation": {}}

# Function to calculate the average read length for a given FASTQ file
def calculate_average_read_length(file_path):
    total_length = 0
    total_reads = 0
    with open(file_path, 'r') as f:
        for i, line in enumerate(f):
            if i % 4 == 1:  # Sequence lines are every 4th line starting from the second line
                total_length += len(line.strip())
                total_reads += 1
    return total_length / total_reads if total_reads > 0 else 0

# Calculate average read lengths for all FASTQ files in the directory

average_read_lengths = {}
fastq_files = [f for f in os.listdir(read_dir) if f.endswith('.fastq')]
for fastq_file in fastq_files:
    file_path = os.path.join(read_dir, fastq_file)
    average_read_lengths[fastq_file] = calculate_average_read_length(file_path)
    print('average_read_lengths', average_read_lengths[fastq_file])

# Load data from JSON files and categorize into training or validation
for file in json_files:
    type_key = "validation" if "validation" in file else "training"
    species_name = file.replace('_validation.json', '').replace('_training.json', '').replace('_', ' ')
    with open(os.path.join(json_directory, file), 'r') as f:
        data = json.load(f)
        if not data:
            print(f"{file} is empty, skipping...")
            batch_info[type_key][species_name] = []
        else:
            species_sequencing_ids[species_name] = [item for sublist in data for item in sublist.values()]
            batch_info[type_key][species_name] = data

# Define the function to simulate batches, merge files, and clean up intermediate files
def simulate_batches(type_key, species, ids, genome_size, batch_info, output_directory_base):
    random.seed(42)  # For reproducibility
    sample_complexity = 2 if len(ids) <= 2 else 3
    species_name_for_files = species.lower().replace(' ', '_')
    depth = 220

    print (output_directory_base)
    output_directory = os.path.join(output_directory_base)
    os.makedirs(output_directory, exist_ok=True)

    for batch_info_item in batch_info[type_key][species]:
        batch_id = batch_info_item['batch_id']
        minority_abundance = batch_info_item['minority_abundance']
        selected_ids, percentages = zip(*( (id, pct) for group in minority_abundance.values() for id, pct in group.items() ))
        subsampled_files = []

        for id, percentage in zip(selected_ids, percentages):
            num_nucleotides = int(genome_size * depth * percentage / 100)
            file_name = f'{id}.fastq'
            input_file = os.path.join(read_dir, file_name)
            output_file = os.path.join(output_directory, f'{species_name_for_files}_{batch_id}_{id}.fastq')

            if not os.path.exists(input_file) or file_name not in average_read_lengths:
                print(f"Error: Input file {input_file} does not exist or average read length is unknown.")
                continue

            avg_read_length = average_read_lengths[file_name]
            num_reads = int(num_nucleotides / avg_read_length)
            command = f'seqtk sample -s100 {input_file} {num_reads} > {output_file}'
            subprocess.run(command, shell=True)
            subsampled_files.append(output_file)

        merged_file_path = os.path.join(output_directory, f'{species_name_for_files}_{batch_id}.fastq')
        with open(merged_file_path, 'w') as outfile:
            for fname in subsampled_files:
                with open(fname) as infile:
                    outfile.write(infile.read())
                os.remove(fname)
        print(f"Merged file created: {merged_file_path}")

# Process each species
# Process each species
for type_key in ['training', 'validation']:
    for species, ids in species_sequencing_ids.items():
        species_underscored = species.split()[0] + '_' + species.split()[1]
        formatted_species = species[0].upper() + species[1:]
        genome_size = genome_sizes.get(formatted_species, None)
        if genome_size:
            output_directory_base = f'simulated_{species_underscored}_{type_key}'
            simulate_batches(type_key, species, ids, genome_size, batch_info, output_directory_base)


print("Simulation complete.")
