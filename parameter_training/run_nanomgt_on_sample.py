import os
import subprocess

# Use this for aligning individual isolates to make the variant position map used to benchmark
path = "/home/people/malhal/data/new_nanomgt/sup_data"

# Align individual isolates to make maps of variance
fastq_files = [f for f in os.listdir(path) if f.endswith('.fastq')]

# Process each file
for file in fastq_files:
    output_name = file[:-6]  # Removes the '.fastq' part from the file name for the output directory
    input_file_path = os.path.join(path, file)

    # Normal run
    #command = f"~/NanoMGT/bin/nanomgt --nanopore {input_file_path} --o {output_name} --threads 8 --db_dir ~/nanomgt_db/"
    # With all params set to 0 for characterization of isolate noise.
    command = f"~/NanoMGT/bin/nanomgt --nanopore {input_file_path} --o {output_name}" \
              f" --threads 8 --db_dir ~/nanomgt_db/ --cor 0 --np 0 --dp 0 --pp 0 --ii 0"
    os.system(command)
