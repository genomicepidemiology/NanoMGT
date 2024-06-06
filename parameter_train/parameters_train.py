import os
import sys
import gzip
import time
import pandas as pd
import numpy as np
import csv
import multiprocessing

from Bio import SeqIO
from itertools import combinations
import concurrent.futures
from itertools import product


sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')] + sys.path



# Perhaps this path won't work, if not look into adjusting your path
from nanomgt.nanopore_variantcaller import train_parameters

# List of folders to process
# Modify this script with the correct path to the data
#path = "/some/path/to/data/"
path = '/home/people/malhal/test/training_data_set/'
files = os.listdir(path)
fastq_files = [f for f in os.listdir(path) if f.endswith('.fastq')]

# Training values
# Use INT here, they will get divided by 100 later
maf = 1

#This represents a gridsearch of the parameters
#dp = [0.1, 0.2, 0.3]
#np = [0.1. 0.2, 0.3]
#pp = [0.1, 0.2, 0.3]
#ii = [0.1, 0.2, 0.3]
#cor = [0.1, 0.2, 0.3]
#For fine-tuning consider a wide range for 1-2 parameters, and then a narrow range for the other parameters to limit the search space
#
cor = [0.1, 0.2, 0.3]
dp = [0.1]
np = [0.1]
pp = [0.1]
ii = [0.1]

parameters = {
    'cor_interval': cor,
    'iteration_increase_interval': ii,
    'pp_interval': pp,
    'np_interval': np,
    'dp_interval': dp
}

cpus = cpu_count_mp = multiprocessing.cpu_count()
cpus = int(cpus/2) #Use half capacity.



def run_jobs_in_parallel(max_workers, new_output_folder, alignment_folder, maf, parameters):
    # Fixed parameters
    min_n = 3
    proxi = 5
    dp_window = 15

    #First grid search
    cor_interval = parameters['cor_interval']
    iteration_increase_interval = parameters['iteration_increase_interval']
    pp_interval = parameters['pp_interval']
    np_interval = parameters['np_interval']
    dp_interval = parameters['dp_interval']

    # Best score initialization
    best_score = 0
    best_params = None
    top_precision = None
    top_recall = None
    top_tp = None
    top_fp = None
    top_fn = None

    # Create all combinations of parameters
    all_params = list(product(cor_interval, iteration_increase_interval, pp_interval, np_interval, dp_interval))
    total_combinations = len(all_params)
    print(f"Total number of parameter combinations: {total_combinations}")

    results_filename = new_output_folder + "/all_results.csv"
    top_result_filename = new_output_folder + "/top_result.csv"

    all_results = []

    processed_combinations = 0

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Create a future for each parameter combination
        futures_to_params = {
            executor.submit(
                train_parameters, maf, alignment_folder, min_n,
                combo[0], new_output_folder, combo[1], proxi, dp_window,
                combo[2], combo[3], combo[4]
            ): combo for combo in all_params
        }

        # Process completed futures
        for future in concurrent.futures.as_completed(futures_to_params):
            params = futures_to_params[future]
            processed_combinations += 1
            #if processed_combinations % 1 == 0:
            print(f"Processed {processed_combinations}/{total_combinations} combinations.")
            try:
                f1, parameter_string, precision, recall, tp, fp, fn = future.result()

                # Write each result to the CSV
                all_results.append([f1, parameter_string, precision, recall, tp, fp, fn])

                if f1 > best_score:
                    best_score = f1
                    best_params = parameter_string
                    top_precision = precision
                    top_recall = recall
                    top_tp = tp
                    top_fp = fp
                    top_fn = fn
            except Exception as exc:
                print(f"Generated an exception: {exc}")

    # Write all results to a CSV

    print (all_results)

    with open(results_filename, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['F1 Score', 'Parameters', 'Precision', 'Recall', 'TP', 'FP', 'FN'])
        writer.writerows(all_results)

    # Write the top result to its own CSV
    with open(top_result_filename, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['F1 Score', 'Parameters', 'Precision', 'Recall', 'TP', 'FP', 'FN'])
        writer.writerow([best_score, best_params, top_precision, top_recall, top_tp, top_fp, top_fn])

output_training_folder = 'training_output_{}'.format(maf)
os.makedirs(output_training_folder, exist_ok=True)

# Loop through each folder
for file in fastq_files:
    # Get all 'merged.fastq' files in the folder

    # Process each filex
    output_name = file[:-6]  # Removes the '.fastq' part from the file name for the output directory
    input_file_path = os.path.join(path, file)

    # This is folder in which the run_nanomgt_on_sample.py script produced folders with alignments.
    alignment_folder = '/home/people/malhal/test/training_test/{}/'.format(output_name)
    new_output_folder = output_training_folder + '/' + output_name + '/'
    os.makedirs(new_output_folder, exist_ok=True)

    run_jobs_in_parallel(cpus, new_output_folder, alignment_folder, maf / 100, parameters)

