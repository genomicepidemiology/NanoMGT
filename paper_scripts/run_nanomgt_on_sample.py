import os
import subprocess

# List of folders to process
path = "/home/people/malhal/data/new_nanomgt/"
folders = os.listdir(path)

# Loop through each folder
for folder in folders:
    if folder.startswith("seed"):

        # Get all 'merged.fastq' files in the folder
        fastq_files = [f for f in os.listdir(path + folder) if f.endswith('merged.fastq')]

        # Process each file
        for file in fastq_files:
            input_file_path = os.path.join(path, folder, file)
            output_name = file[:-6]  # Removes the '.fastq' part from the file name for the output directory

            # Construct the command
            command = f"~/NanoMGT/bin/nanomgt --nanopore {input_file_path} --o {output_name} --threads 8 --db_dir ~/nanomgt_db/"

            # Run the command
            try:
                subprocess.run(command, shell=True, check=True)
                print(f"Processed {file} successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error processing {file}: {e}")

    print("Processing completed.")
