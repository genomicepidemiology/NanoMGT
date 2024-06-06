import os

def construct_command(file_path, output_dir):
    # Remove the .fastq extension to get the output name
    base_name = os.path.basename(file_path)
    output_name = base_name[:-6]  # Strip '.fastq' from the end

    #Adjust these parameters
    # Construct the command with all parameters set to 0
    #command = f"~/NanoMGT/bin/nanomgt --nanopore {file_path} --o {output_dir}/{output_name} --threads 8 --db_dir ~/nanomgt_db/" \
    #          f" --maf 0.03 --cor 0 --np 0 --dp 0 --ii 0 --pp 0"

    # Type co-occurencing noise at default setting for 0.03 maf
    command = f"~/NanoMGT/bin/nanomgt --nanopore {file_path} --o {output_dir}/{output_name} --threads 8 --db_dir ~/nanomgt_db/" \
             f" --maf 0.05"
    return command

def main():
    directory = "/home/people/malhal/data/new_nanomgt"
    #output_directory = "majority_variants"  # Define or adapt your output directory path
    output_directory = "co_noise"  # Define or adapt your output directory path

    # List all files in the directory
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)

        # Check if the file ends with .fastq
        if file_path.endswith(".fastq"):
            # Generate and print the command
            command = construct_command(file_path, output_directory)
            os.system(command)

if __name__ == "__main__":
    main()
