import pandas as pd
import os
import sys

# Merge function definition
def merge_files(batch_files, batch, output_directory, output_prefix):
    output_file_path = os.path.join(output_directory, f'{output_prefix}_batch{batch}_merged.fastq')
    with open(output_file_path, 'w') as outfile:
        for fname in batch_files:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

# Function to count nucleotides in a FASTQ file
def count_nucleotides_in_fastq(file_path):
    nucleotide_count = 0
    with open(file_path, 'r') as fastq_file:
        for i, line in enumerate(fastq_file):
            if i % 4 == 1:  # Sequence lines are every 4th line starting from the second line
                nucleotide_count += len(line.strip())
    return nucleotide_count

# Parameters
genome_sizes = {
    'ecoli': 5000000,
    'staph_aureus': 3000000,
    'campylobacter_jejuni': 1700000,
    'salmonella_enterica': 5000000,
    'klebsiella_pneumoniae': 5500000
}

# Mapping of isolates to organisms
isolate_to_organism = {
    'SRR25689478': 'ecoli',
    'SRR26036455': 'ecoli',
    'ERR12533301': 'ecoli',
    'ERR8958737': 'klebsiella_pneumoniae',
    'SRR27348733': 'klebsiella_pneumoniae',
    'SRR27638397': 'campylobacter_jejuni',
    'SRR27710526': 'campylobacter_jejuni',
    'SRR28399428': 'salmonella_enterica',
    'SRR27136088': 'salmonella_enterica',
    'SRR27755678': 'salmonella_enterica',
    'SRR28370694': 'staph_aureus',
    'ERR8958843': 'staph_aureus'
}

# Average read lengths for each file
average_read_lengths = {
    'ERR8958843.fastq': 7798.34,
    'SRR26036455.fastq': 2445.88,
    'SRR27136088.fastq': 7780.01,
    'SRR27638397.fastq': 3516.32,
    'SRR27710526.fastq': 3907.54,
    'SRR27755678.fastq': 2891.15,
    'SRR28370694.fastq': 4247.29,
    'SRR28399428.fastq': 7268.81,
    'SRR25689478.fastq': 5753.44,
    'SRR27348733.fastq': 8129.34,
    'ERR8958737.fastq': 8661.8,
    'ERR12533301.fastq': 7730.87
}

# File paths
file_paths = [
    'SRR27638397_majority_batches.csv',
    'SRR27710526_majority_batches.csv',
    'SRR25689478_majority_batches.csv',
    'SRR26036455_majority_batches.csv',
    'ERR12533301_majority_batches.csv',
    'ERR8958737_majority_batches.csv',
    'SRR27348733_majority_batches.csv',
    'SRR27136088_majority_batches.csv',
    'SRR27755678_majority_batches.csv',
    'SRR28399428_majority_batches.csv',
    'ERR8958843_majority_batches.csv',
    'SRR28370694_majority_batches.csv'
]

for file_path in file_paths:
    print(file_path)
    description_file = f'simulated_batches/{file_path}'

    # Load the samples description
    df = pd.read_csv(description_file)
    df.columns = [col.strip() for col in df.columns]  # Strip leading/trailing spaces from column names

    # Iterate over different depths for each organism
    depths = [120, 170, 220]
    for depth in depths:
        file_name = file_path.split('.')[0]
        output_directory = f'depth{depth}_{file_name}'
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

        output_prefix = f'depth{depth}_{file_name}'

        # Extract the organism from the first ID in the dataframe
        first_id = df.iloc[0]['ID1'].strip()
        organism = isolate_to_organism[first_id]
        genome_size = genome_sizes[organism]
        desired_total_nucleotides = genome_size * depth

        nucleotide_distribution_summary = []  # Collect nucleotide distribution information

        for index, row in df.iterrows():
            batch = row['Batch']
            ids = [row[col] for col in df.columns if 'ID' in col]
            percentages = [float(row[col].rstrip('%')) / 100 for col in df.columns if 'Percentage' in col]

            subsampled_files = []
            for id, percentage in zip(ids, percentages):
                num_nucleotides = int(desired_total_nucleotides * percentage)
                # Adjust to where the fastq files are.
                input_file = f'/home/people/malhal/data/new_nanomgt/{id.strip()}.fastq'
                output_file = os.path.join(output_directory, f'{output_prefix}_batch{batch}_{id.strip()}.fastq')

                # Check if input file exists before processing
                if not os.path.exists(input_file):
                    print(f"Error: Input file {input_file} does not exist.")
                    continue

                # Use the correct average read length for this file
                avg_read_length = average_read_lengths.get(f'{id.strip()}.fastq', 10000)  # Fallback to 10000 if not found
                num_reads = num_nucleotides // avg_read_length
                command = f'seqtk sample -s100 {input_file} {num_reads} > {output_file}'
                os.system(command)
                subsampled_files.append(output_file)

            # Merge the subsampled files for this batch
            merge_files(subsampled_files, batch, output_directory, output_prefix)

            # Count nucleotides in the merged file
            merged_file_path = os.path.join(output_directory, f'{output_prefix}_batch{batch}_merged.fastq')
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

# Note: This script assumes all necessary input files and directories are correctly set up and accessible.
