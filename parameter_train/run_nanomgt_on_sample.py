import os
import subprocess

# List of folders to process
# Modify this script with the correct path to the data
#path = "/some/path/to/data/"
path = '/home/people/malhal/test/training_data_set/'
files = os.listdir(path)
fastq_files = [f for f in os.listdir(path) if f.endswith('.fastq')]


# Loop through each folder
for file in fastq_files:
    # Get all 'merged.fastq' files in the folder

    # Process each filex
    output_name = file[:-6]  # Removes the '.fastq' part from the file name for the output directory
    input_file_path = os.path.join(path, file)

    # Construct the command
    # This runs NanoMGT once with the default parameters to produce the required alignment.
    # The alignment can then be used for new training of parameters.
    # Modify this script with the correct paths and the database to be used.
    command = f"~/NanoMGT/bin/nanomgt --nanopore {input_file_path} --o {output_name} --threads 8 --db_dir ~/nanomgt_db/"
    print (command)
    os.system(command)
