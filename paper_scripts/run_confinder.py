import os
import subprocess

# Define the path to the directory containing the folders
data_dir_path = "/home/people/malhal/data/new_nanomgt/confindr_data_dir"

# List all entries in the directory and filter out the folders
folders = [name for name in os.listdir(data_dir_path)
           if os.path.isdir(os.path.join(data_dir_path, name))]

bf = [0.01, 0.02, 0.03, 0.04, 0.05]
# Iterate over the list of folders and run the echo command for each
for rate in bf:
    os.system(f"mkdir /home/people/malhal/data/new_nanomgt/confindr_results/bf_{rate}")
    for folder in folders:
        cmd = f"confindr -i {data_dir_path}/{folder}/ -o /home/people/malhal/data/new_nanomgt/confindr_results/bf_{rate}/{folder} -d /home/people/malhal/contamErase/benchmarking/confindr/rmlst_db/ --rmlst -b 3 -bf {rate} -dt Nanopore -q 14"
        os.system(cmd)
        #subprocess.run(command, shell=True)