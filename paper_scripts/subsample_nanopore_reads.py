import pandas as pd
import os


# Function to count nucleotides in a FASTQ file
def count_nucleotides_in_fastq(file_path):
    nucleotide_count = 0
    with open(file_path, 'r') as fastq_file:
        for i, line in enumerate(fastq_file):
            if i % 4 == 1:  # Sequence lines are every 4th line starting from the second line
                nucleotide_count += len(line.strip())
    return nucleotide_count


# Parameters
description_file = 'simulated_batches/ecoli_2s.csv'  # Path to your file
genome_size = 5000000  # Genome size
total_depth = 300  # Total depth
desired_total_nucleotides = genome_size * total_depth
output_prefix = 'ecoli_2s'  # Prefix for output files
output_directory = 'simulated_ecoli_2s'  # Directory to store output files

# Average read lengths for each file
average_read_lengths = {
    'ERR8958843.fastq': 7798.34,
    'SRR23387317.fastq': 17704,
    'SRR26643493.fastq': 12211.8,
    'SRR26899121.fastq': 2339.21,
    'SRR26899125.fastq': 1329.76,
    'SRR26899129.fastq': 2704.13,
    'SRR27136088.fastq': 7780.01,
    'SRR27167517.fastq': 1191.5,
    'SRR27638397.fastq': 3516.32,
    'SRR27710526.fastq': 3907.54,
    'SRR27710531.fastq': 4850.71,
    'SRR27755678.fastq': 2891.15,
    'SRR28370694.fastq': 4247.29,
    'SRR28399428.fastq': 7268.81,
    'SRR28399430.fastq': 4996.18,
    'SRR25689478.fastq': 5753.44
}

# Ensure output directory exists
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Load the samples description
df = pd.read_csv(description_file)


# Function to merge subsampled files
def merge_files(batch_files, batch):
    output_file_path = os.path.join(output_directory, f'{output_prefix}_batch{batch}_merged.fastq')
    with open(output_file_path, 'w') as outfile:
        for fname in batch_files:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
            # Optionally, remove the subsampled file after merging
            # os.remove(fname)


nucleotide_distribution_summary = []  # Collect nucleotide distribution information

# Iterate through each batch to subsample based on nucleotide count
total_nucleotide_counts = {}  # Dictionary to hold total nucleotide counts for each batch

for index, row in df.iterrows():
    batch = row['Batch']

    # Collect IDs and percentages
    ids = []
    percentages = []
    for col in df.columns:
        if 'ID' in col:
            ids.append(row[col])
        elif 'Percentage' in col:
            percentages.append(row[col].rstrip('%'))

    subsampled_files = []  # Keep track of subsampled files for merging
    for id, percentage in zip(ids, percentages):
        id = id.strip()
        # Calculate the number of nucleotides for this sample
        num_nucleotides = int(desired_total_nucleotides * (float(percentage) / 100))
        # Retrieve the average read length for this sample ID
        file_name = f'{id}.fastq'
        avg_read_length = average_read_lengths.get(file_name, 10000)  # Default to 10000 if not found
        # Estimate the number of reads to achieve the desired nucleotide count
        num_reads = num_nucleotides // avg_read_length
        input_file = file_name
        output_file = os.path.join(output_directory, f'{output_prefix}_batch{batch}_{id}.fastq')

        # Construct and execute seqtk command
        command = f'seqtk sample -s100 {input_file} {num_reads} > {output_file}'
        os.system(command)

        subsampled_files.append(output_file)

    # Merge the subsampled files for this batch
    merged_file_path = os.path.join(output_directory, f'{output_prefix}_batch{batch}_merged.fastq')
    merge_files(subsampled_files, batch)

    # Count nucleotides in the merged file
    merged_nucleotide_count = count_nucleotides_in_fastq(merged_file_path)
    total_nucleotide_counts[batch] = merged_nucleotide_count

    # Collect information for nucleotide distribution summary
    for subsampled_file in subsampled_files:
        file_nucleotide_count = count_nucleotides_in_fastq(subsampled_file)
        percentage_of_total = (file_nucleotide_count / merged_nucleotide_count) * 100
        nucleotide_distribution_summary.append(f'{subsampled_file}: {percentage_of_total:.2f}% of nucleotides')

print("Subsampling, merging, and nucleotide counting complete.")

# Print nucleotide distribution summary
print("\nNucleotide Distribution Summary:")
for summary in nucleotide_distribution_summary:
    print(summary)
