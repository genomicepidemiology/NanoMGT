import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

# Define the path where the FASTQ files are located
path = "/home/people/malhal/data/new_nanomgt/clean_data/120/120_validation_reads"

# List all fastq files in the directory
fastq_files = [f for f in os.listdir(path) if f.endswith('.fastq')]

# Function to process each file
def process_file(file):
    output_name = file[:-6]  # Removes the '.fastq' part from the file name for the output directory
    input_file_path = os.path.join(path, file)

    # Construct the command to run NanoMGT with specified parameters
    command = f"~/NanoMGT/bin/nanomgt --nanopore {input_file_path} --o {output_name} " \
              f"--threads 2 --db_dir ~/nanomgt_db/ --cor 0 --np 0 --dp 0 --pp 0 --ii 0 --maf 0.05"
    subprocess.run(command, shell=True)

# Create a ThreadPoolExecutor to run processes in parallel
with ThreadPoolExecutor(max_workers=10) as executor:
    # Map the process_file function to each FASTQ file
    executor.map(process_file, fastq_files)

print("All processes are complete.")
