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
import argparse

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')] + sys.path

from nanomgt import nanopore_variantcaller as nvc

# List of folders to process
# Modify this script with the correct path to the data
path = '/home/people/malhal/test/training_data_set/'
files = os.listdir(path)
fastq_files = [f for f in os.listdir(path) if f.endswith('.fastq')]

# Training values
# Use INT here, they will get divided by 100 later
maf = 1

# This represents a gridsearch of the parameters
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
cpus = int(cpus / 2)  # Use half capacity.

def train_parameters(maf, results_folder, min_n, cor, new_output_folder,
                    iteration_increase, proxi, dp_window, pp, np, dp):
    arguments = argparse.Namespace()
    arguments.maf = maf
    arguments.output = results_folder
    arguments.min_n = min_n
    arguments.cor = cor
    arguments.new_output = new_output_folder
    arguments.iteration_increase = iteration_increase
    arguments.proxi = proxi
    arguments.dp = dp
    arguments.pp = pp
    arguments.np = np
    arguments.dp_window = dp_window

    # Build a consensus dictionary from alignment results
    consensus_dict = nvc.build_consensus_dict(os.path.join(arguments.output, 'rmlst_alignment.res'),
                                          os.path.join(arguments.output, 'rmlst_alignment.mat'))

    confirmed_mutation_dict = nvc.derive_mutation_positions(consensus_dict, arguments)

    # Perform biological validation of mutations
    bio_validation_dict = nvc.bio_validation_mutations(consensus_dict, os.path.join(results_folder, 'specie.fsa'))
    # Co-occurrence analysis until convergence
    confirmed_mutation_dict, co_occurrence_tmp_dict, iteration_count =\
        nvc.co_occurrence_until_convergence(arguments, confirmed_mutation_dict,
                                        consensus_dict, {}, bio_validation_dict)


    # Format and output the results
    format_output(new_output_folder, confirmed_mutation_dict, consensus_dict, bio_validation_dict,
                  co_occurrence_tmp_dict)

    print (arguments.output)
    print (results_folder)

    sample = arguments.output.split('/')[-1]

    print (sample)
    print (sample)

    minor_mutation_expected = benchmark_analysis_result(sample)

    minor_mutation_results = convert_mutation_dict_to_object(confirmed_mutation_dict)

    print (len(minor_mutation_expected), len(minor_mutation_results))

    precision, recall, f1, tp, fp, fn = calculate_metrics(minor_mutation_expected, minor_mutation_results)


    precision = 1
    recall = 1
    f1 = 1
    tp = 1
    fp = 1
    fn = 1

    parameter_string = f"maf_{maf}_cor_{cor}_pp_{pp}_np_{np}_dp_{dp}_iteration_increase_{iteration_increase}"



    return f1, parameter_string, precision, recall, tp, fp, fn

def benchmark_analysis_result(sample):
    #batch = sample.split('_')[-2]
    #print (sample)
    print (sample)
    batch_id = int(sample.split('_')[-2][5:])
    #print (batch_id)
    #print(f"Batch ID: {batch_id}")
    specie = sample.split('_')[1:-5]
    batch_csv = "_".join(sample.split('_')[1:-2]) + ".csv"
    data = load_data('/home/people/malhal/data/new_nanomgt/simulated_batches/' + batch_csv)
    # Change this batch_id to test different batches
    # batch_id = 10
    #print(f"Highest percentage ID for batch {batch_id}: {find_highest_percentage_id(batch_id, data)}")

    top_id, minor = find_highest_percentage_id(batch_id, data)
    # Use names and batch ID to get the correct mutation map
    if 'salmonella_enterica' in sample:
        map_file_1 = '/home/people/malhal/data/new_nanomgt/majority_variants/minor_variant_maps/major_{}_minor_{}.txt'.format(
            top_id, minor[0])
        map_file_2 = '/home/people/malhal/data/new_nanomgt/majority_variants/minor_variant_maps/major_{}_minor_{}.txt'.format(
            top_id, minor[1])
        mutation_map = load_mutations_from_files([map_file_1, map_file_2])
    if 'ecoli' in sample:
        map_file_1 = '/home/people/malhal/data/new_nanomgt/majority_variants/minor_variant_maps/major_{}_minor_{}.txt'.format(
            top_id, minor[0])
        map_file_2 = '/home/people/malhal/data/new_nanomgt/majority_variants/minor_variant_maps/major_{}_minor_{}.txt'.format(
            top_id, minor[1])
        mutation_map = load_mutations_from_files([map_file_1, map_file_2])
    if 'staph_aureus' in sample:
        map_file_1 = '/home/people/malhal/data/new_nanomgt/majority_variants/minor_variant_maps/major_{}_minor_{}.txt'.format(
            top_id, minor[0])
        mutation_map = load_mutations_from_files([map_file_1])
    if 'campylobacter_jejuni' in sample:
        map_file_1 = '/home/people/malhal/data/new_nanomgt/majority_variants/minor_variant_maps/major_{}_minor_{}.txt'.format(
            top_id, minor[0])
        mutation_map = load_mutations_from_files([map_file_1])
    #print ('My expected output')
    #for item in mutation_map:
    #    print(item, mutation_map[item])

    return mutation_map


def format_output(new_output_folder, confirmed_mutation_dict, consensus_dict, bio_validation_dict, co_occurrence_tmp_dict):
    """
    Format and print the output of confirmed mutations with additional information.

    Args:
        confirmed_mutation_dict (dict): A dictionary containing confirmed mutations for alleles.
        consensus_dict (dict): A dictionary containing consensus information for alleles.
        bio_validation_dict (dict): A dictionary containing biological validation data for genes.

    Returns:
        None
    """
    with open(new_output_folder + '/minor_mutations.csv', 'w') as outfile:
        header = 'Gene,Position,MajorityBase,MutationBase,MutationDepth,TotalDepth,GeneLength,MutationComment,CoOccurrence'
        print(header, file=outfile)

        for allele in confirmed_mutation_dict:
            for mutation in zip(confirmed_mutation_dict[allele][0], confirmed_mutation_dict[allele][1]):
                position = mutation[0].split('_')[0]
                mutation_base = mutation[0].split('_')[1]
                mutation_depth = mutation[1]
                majority_base = consensus_dict[allele][1][int(position) - 1]
                total_depth = sum(consensus_dict[allele][0][int(position) - 1])
                biological_existence = nvc.check_single_mutation_existence(bio_validation_dict, allele, mutation[0])
                gene_length = len(consensus_dict[allele][1])
                if mutation[0] in co_occurrence_tmp_dict[allele]:
                    co_occurrence = 'Yes'
                else:
                    co_occurrence = 'No'

                if biological_existence:
                    print('{},{},{},{},{},{},{},{},{}'.format(allele, position, majority_base, mutation_base,
                                                              mutation_depth, total_depth, gene_length,
                                                              'Mutation seen in database', co_occurrence), file=outfile)
                else:
                    print('{},{},{},{},{},{},{},{},{}'.format(allele, position, majority_base, mutation_base,
                                                              mutation_depth, total_depth, gene_length,
                                                              'Novel mutation', co_occurrence), file=outfile)


def convert_mutation_dict_to_object(mutation_dict):
    mutation_object = {}
    for allele in mutation_dict:
        gene = allele.split('_')[0]
        mutation_object[gene] = set()
        for mutation in mutation_dict[allele][0]:
            mutation_object[gene].add(mutation)
    return mutation_object

def calculate_metrics(expected_mutations, actual_mutations):
    """ Calculate precision, recall, and F1 score for the predicted mutations. """
    # Initialize variables to calculate sum of metrics across all genes
    sum_precision, sum_recall, sum_f1 = 0, 0, 0
    genes_counted = 0

    total_tp = 0
    total_fp = 0
    total_fn = 0

    # Iterate over expected mutations by gene to calculate metrics
    for gene, expected_set in expected_mutations.items():
        if gene in actual_mutations:
            actual_set = actual_mutations[gene]
            # True Positives (TP): Mutations correctly predicted
            tp = len(expected_set & actual_set)
            # False Positives (FP): Mutations incorrectly predicted (not in expected but in actual)
            fp = len(actual_set - expected_set)
            # False Negatives (FN): Mutations missed (in expected but not in actual)
            fn = len(expected_set - actual_set)

            total_tp += tp
            total_fp += fp
            total_fn += fn

    precision = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else 0
    recall = total_tp / (total_tp + total_fn) if (total_tp + total_fn) > 0 else 0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

    return precision, recall, f1, total_tp, total_fp, total_fn


def run_jobs_in_parallel(max_workers, new_output_folder, alignment_folder, maf, parameters):
    # Fixed parameters
    min_n = 3
    proxi = 5
    dp_window = 15

    # First grid search
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
            print(f"Processed {processed_combinations}/{total_combinations} combinations.")
            try:
                result = future.result()
                f1, parameter_string, precision, recall, tp, fp, fn = result
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

    # Process each file
    output_name = file[:-6]  # Removes the '.fastq' part from the file name for the output directory
    input_file_path = os.path.join(path, file)

    # This is folder in which the run_nanomgt_on_sample.py script produced folders with alignments.
    alignment_folder = '/home/people/malhal/test/training_test/{}/'.format(output_name)
    new_output_folder = output_training_folder + '/' + output_name
    os.makedirs(new_output_folder, exist_ok=True)

    #run_jobs_in_parallel(cpus, new_output_folder, alignment_folder, maf / 100, parameters)
    train_parameters(maf / 100, alignment_folder, 3, 0.1, new_output_folder, 1, 5, 0.1, 0.1, 0.1, 0.1)
