import os
import sys
import gzip
import time
import pandas as pd
import numpy as np
import csv
import re
import multiprocessing
import json

from Bio import SeqIO
from itertools import combinations
import concurrent.futures
from itertools import product
import argparse
from collections import defaultdict

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')] + sys.path

from nanomgt import nanopore_variantcaller as nvc

alignment_results_path = '/home/people/malhal/test/training_test/'
maps_path = '/home/people/malhal/test/training_test/maps/'
simulated_batches_csv_path = '/home/people/malhal/test/training_test/data/simulated_batches/'
files = os.listdir(alignment_results_path)
folders = [f for f in os.listdir(alignment_results_path)]

maf_interval = [1, 2, 3, 4, 5]

cor_interval_search = [0.1, 0.3, 0.5, 0.7]
dp_interval_search = [0.1, 0.3, 0.5]
np_interval_search = [1.5, 2.5, 3.5]
pp_interval_search = [0.01, 0.3, 0.5, 0.7]
ii_interval_search = [0.01, 0.2, 0.40, 0.60]

parameters_interval_search = {
    'cor_interval': cor_interval_search,
    'iteration_increase_interval': ii_interval_search,
    'pp_interval': pp_interval_search,
    'np_interval': np_interval_search,
    'dp_interval': dp_interval_search
}

cpus = 40

def train_parameters(maf, results_folder, min_n, cor, new_output_folder, maps_path, simulated_batches_csv_path,
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

    consensus_dict = nvc.build_consensus_dict(os.path.join(arguments.output, 'rmlst_alignment.res'),
                                          os.path.join(arguments.output, 'rmlst_alignment.mat'))

    confirmed_mutation_dict = nvc.derive_mutation_positions(consensus_dict, arguments)

    bio_validation_dict = nvc.bio_validation_mutations(consensus_dict, os.path.join(results_folder, 'specie.fsa'))

    confirmed_mutation_dict, co_occurrence_tmp_dict, iteration_count =\
        nvc.co_occurrence_until_convergence(arguments, confirmed_mutation_dict,
                                        consensus_dict, {}, bio_validation_dict)

    format_output(new_output_folder, confirmed_mutation_dict, consensus_dict, bio_validation_dict,
                  co_occurrence_tmp_dict)

    sample = arguments.output.split('/')[-1]

    minor_mutation_expected = benchmark_analysis_result(sample, simulated_batches_csv_path, maps_path)

    minor_mutation_results = convert_mutation_dict_to_object(confirmed_mutation_dict)

    print (len(minor_mutation_expected), len(minor_mutation_results))

    precision, recall, f1, tp, fp, fn = calculate_metrics(minor_mutation_expected, minor_mutation_results)

    parameter_string = f"maf_{maf}_cor_{cor}_pp_{pp}_np_{np}_dp_{dp}_iteration_increase_{iteration_increase}"

    return f1, parameter_string, precision, recall, tp, fp, fn

def load_data(filepath):
    return pd.read_csv(filepath, skipinitialspace=True)

def find_highest_percentage_id(batch_id, df):
    batch_data = df[df['Batch'] == batch_id]
    if batch_data.empty:
        return "Batch ID not found."
    highest_percentage = 0
    highest_percentage_id = ''
    all = []
    for index, row in batch_data.iterrows():
        if 'Percentage3' in df.columns:
            for i in range(1, 4):
                percentage = int(row[f'Percentage{i}'][:-1])
                if percentage > highest_percentage:
                    highest_percentage = percentage
                    highest_percentage_id = row[f'ID{i}']
                all.append((row[f'ID{i}']))
        else:
            for i in range(1, 3):
                percentage = int(row[f'Percentage{i}'][:-1])
                if percentage > highest_percentage:
                    highest_percentage = percentage
                    highest_percentage_id = row[f'ID{i}']
                all.append((row[f'ID{i}']))

    minor = []

    for item in all:
        if item != highest_percentage_id:
            minor.append(item)

    return highest_percentage_id, minor

def load_mutations_from_files(file_paths):
    mutations_dict = {}
    for file_path in file_paths:
        with open(file_path, 'r') as file:
            lines = file.readlines()
            for i in range(0, len(lines), 2):
                gene_id = lines[i].strip()
                mutations = set(lines[i + 1].strip().split(','))

                if gene_id in mutations_dict:
                    mutations_dict[gene_id].update(mutations)
                else:
                    mutations_dict[gene_id] = mutations

    return mutations_dict

def get_number_of_columns(dataframe):
    return dataframe.shape[1]

def benchmark_analysis_result(sample, batch_csv_path, maps_path):
    batch_id = int(sample.split('_')[-2][5:])
    batch_csv = batch_csv_path + "_".join(sample.split('_')[1:-2]) + ".csv"
    data = load_data(batch_csv)
    sample_number = int((get_number_of_columns(data) - 1) / 2)

    top_id, minor = find_highest_percentage_id(batch_id, data)

    map_files = []
    for i in range(sample_number-1):
        map_file = f'{maps_path}/major_{top_id}_minor_{minor[i]}.txt'
        map_files.append(map_file)

    mutation_map = load_mutations_from_files(map_files)

    return mutation_map

def format_output(new_output_folder, confirmed_mutation_dict, consensus_dict, bio_validation_dict, co_occurrence_tmp_dict):
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
    sum_precision, sum_recall, sum_f1 = 0, 0, 0
    genes_counted = 0

    total_tp = 0
    total_fp = 0
    total_fn = 0

    for gene, expected_set in expected_mutations.items():
        if gene in actual_mutations:
            actual_set = actual_mutations[gene]
            tp = len(expected_set & actual_set)
            fp = len(actual_set - expected_set)
            fn = len(expected_set - actual_set)

            total_tp += tp
            total_fp += fp
            total_fn += fn

    precision = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else 0
    recall = total_tp / (total_tp + total_fn) if (total_tp + total_fn) > 0 else 0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

    return precision, recall, f1, total_tp, total_fp, total_fn

def run_jobs_in_parallel(max_workers, new_output_folder, alignment_folder, maf, parameters, maps_path, simulated_batches_csv_path):
    min_n = 3
    proxi = 5
    dp_window = 15

    cor_interval = parameters['cor_interval']
    iteration_increase_interval = parameters['iteration_increase_interval']
    pp_interval = parameters['pp_interval']
    np_interval = parameters['np_interval']
    dp_interval = parameters['dp_interval']

    best_score = 0
    best_params = None
    top_precision = None
    top_recall = None
    top_tp = None
    top_fp = None
    top_fn = None

    all_params = list(product(cor_interval, iteration_increase_interval, pp_interval, np_interval, dp_interval))
    total_combinations = len(all_params)
    print(f"Total number of parameter combinations: {total_combinations}")

    results_filename = new_output_folder + "/all_results.csv"
    top_result_filename = new_output_folder + "/top_result.csv"

    all_results = []

    processed_combinations = 0

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures_to_params = {
            executor.submit(
                train_parameters, maf, alignment_folder, min_n,
                combo[0], new_output_folder, maps_path, simulated_batches_csv_path, combo[1], proxi, dp_window,
                combo[2], combo[3], combo[4]
            ): combo for combo in all_params
        }

        for future in concurrent.futures.as_completed(futures_to_params):
            params = futures_to_params[future]
            processed_combinations += 1
            print(f"Processed {processed_combinations}/{total_combinations} combinations.")
            try:
                result = future.result()
                f1, parameter_string, precision, recall, tp, fp, fn = result
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

    with open(results_filename, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['F1 Score', 'Parameters', 'Precision', 'Recall', 'TP', 'FP', 'FN'])
        writer.writerows(all_results)

    with open(top_result_filename, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['F1 Score', 'Parameters', 'Precision', 'Recall', 'TP', 'FP', 'FN'])
        writer.writerow([best_score, best_params, top_precision, top_recall, top_tp, top_fp, top_fn])

def extract_parameters(param_string):
    param_pattern = r'([a-z_]+)_([0-9.]+)'
    return dict(re.findall(param_pattern, param_string))

def calculate_best_parameters(file_name):
    df = pd.read_csv(file_name)

    parameters_data = defaultdict(lambda: {'f1_scores': [], 'best_f1_score': float('-inf'), 'best_params': {}})

    for index, row in df.iterrows():
        f1_score = row['F1 Score']
        parameters = extract_parameters(row['Parameters'])

        for param, value in parameters.items():
            key = f"{param}_{value}"
            parameters_data[key]['f1_scores'].append(f1_score)

            if f1_score > parameters_data[key]['best_f1_score']:
                parameters_data[key]['best_f1_score'] = f1_score
                parameters_data[key]['best_params'] = {param: float(value)}

    best_params_per_file = {}
    for key, data in parameters_data.items():
        param_name, param_value = key.rsplit('_', 1)
        param_value = float(param_value)

        if param_name not in best_params_per_file:
            best_params_per_file[param_name] = (param_value, data['best_f1_score'])
        elif data['best_f1_score'] > best_params_per_file[param_name][1]:
            best_params_per_file[param_name] = (param_value, data['best_f1_score'])

    return {param: value for param, (value, _) in best_params_per_file.items()}

output_training_folder = 'nanomgt_training_output'
os.makedirs(output_training_folder, exist_ok=True)

for maf in maf_interval:
    os.makedirs(output_training_folder + '/maf_' + str(maf), exist_ok=True)
    for folder in folders:
        if folder.startswith('depth'):
            batch_id = int(folder.split('_')[-2][5:])
            if maf >= batch_id:
                input_file_path = os.path.join(alignment_results_path, folder)
                alignment_folder = '/home/people/malhal/test/training_test/{}'.format(folder)
                new_output_folder = output_training_folder + '/' + 'maf_' + str(maf) + '/' + folder
                os.makedirs(new_output_folder, exist_ok=True)
                run_jobs_in_parallel(cpus, new_output_folder, alignment_folder, maf / 100, parameters_interval_search, maps_path, simulated_batches_csv_path)

all_best_params = defaultdict(list)

for maf in maf_interval:
    print(f"maf_{maf}")
    average_best_params = {}
    for folder in folders:
        if folder.startswith('depth'):
            batch_id = int(folder.split('_')[-2][5:])
            if maf >= batch_id:
                new_output_folder = output_training_folder + '/' + 'maf_' + str(maf) + '/' + folder
                results_filename = new_output_folder + "/all_results.csv"
                best_params = calculate_best_parameters(results_filename)
                for param, value in best_params.items():
                    all_best_params[param].append(value)
    for param, values in all_best_params.items():
        average_value = sum(values) / len(values)
        average_best_params[param] = average_value
        print(f"Average of best {param}: {average_value:.4f}")
    output_file_path = os.path.join(output_training_folder, "{}_average_best_params.json".format('maf_' + str(maf)))
    with open(output_file_path, 'w') as json_file:
        json.dump(average_best_params, json_file, indent=4)

    print(f"Averages saved to {output_file_path}")
