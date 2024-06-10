import os
import subprocess

# List of folders to process
path = "/home/people/malhal/data/new_nanomgt/genomes/genomes/simulated_reads/renamed_reads/"
folders = os.listdir(path)

# Loop through each folder
for folder in folders:
    if folder.startswith("depth"):
        # Get all 'merged.fastq' files in the folder
        fastq_files = [f for f in os.listdir(path + folder) if f.endswith('merged.fastq')]

        # Process each file
        for file in fastq_files:
            output_name = file[:-6]  # Removes the '.fastq' part from the file name for the output directory
            input_file_path = os.path.join(path, folder, file)

            # Construct the command
            command = f"~/NanoMGT/bin/nanomgt --nanopore {input_file_path} --o {output_name} --threads 8 --db_dir ~/nanomgt_db/"

            os.system(command)

import os
import subprocess

# List of folders to process
path = "/home/people/malhal/data/new_nanomgt/genomes/genomes/simulated_reads/renamed_reads/"
folders = os.listdir(path)

# Loop through each folder
for folder in folders:
    if folder.startswith("depth"):
        # Get all 'merged.fastq' files in the folder
        fastq_files = [f for f in os.listdir(path + folder) if f.endswith('merged.fastq')]

        # Process each file
        for file in fastq_files:
            output_name = file[:-6]  # Removes the '.fastq' part from the file name for the output directory
            input_file_path = os.path.join(path, folder, file)

            # Construct the command
            command = f"~/NanoMGT/bin/nanomgt --nanopore {input_file_path} --o {output_name} --threads 8 --db_dir ~/nanomgt_db/"

            os.system(command)


"""
# Use this for aligning individual isolates to make the variant position map used to benchmark
path = "/home/people/malhal/data/new_nanomgt/genomes/genomes/simulated_reads/renamed_reads/"

# Align individual isolates to make maps of variance
fastq_files = [f for f in os.listdir(path) if f.endswith('.fastq')]

# Process each file
for file in fastq_files:
    output_name = file[:-6]  # Removes the '.fastq' part from the file name for the output directory
    input_file_path = os.path.join(path, file)

    # Construct the command
    command = f"~/NanoMGT/bin/nanomgt --nanopore {input_file_path} --o {output_name} --threads 8 --db_dir ~/nanomgt_db/"

    os.system(command)

"""