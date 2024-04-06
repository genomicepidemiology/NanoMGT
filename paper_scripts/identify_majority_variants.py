import os
import subprocess
import sys

# Path to the directory
directory_path = '/home/people/malhal/data/new_nanomgt'

# Placeholder command to run for each .fastq file
placeholder_command = 'echo Processing file'

# List all files in the directory
files_in_directory = os.listdir(directory_path)

# Filter for .fastq files
fastq_files = [file for file in files_in_directory if file.endswith('.fastq')]

# Iterate over each .fastq file and run the placeholder command
for fastq_file in fastq_files:
    # Construct the full path to the file
    full_path = os.path.join(directory_path, fastq_file)
    name = fastq_file.split('.')[0]

    # Run the placeholder command with the full path of the file
    # For example, 'echo Processing file /path/to/file.fastq'
    cmd = '~/NanoMGT/bin/nanomgt --nanopore {} --db_dir ~/nanomgt_db/ --o {} --threads 4'.format(full_path, name)

    # Execute the command
    process = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Optionally, print the output of the command
    print(process.stdout.decode())

print("Processing complete.")