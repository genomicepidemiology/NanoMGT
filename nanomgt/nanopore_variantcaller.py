import os
import sys
import gzip
import time
import json
import numpy as np
import argparse
from Bio import SeqIO
from nanomgt import kma
from nanomgt import util
from itertools import combinations
from scipy.interpolate import UnivariateSpline

from nanomgt.nanopore_mutations import parse_sam_and_find_mutations
from nanomgt.nanopore_mutations import extract_alignment
from nanomgt.nanopore_mutations import create_mutation_vector
import nanomgt.curve_functions as cf

def nanopore_metagenomics_variantcaller(arguments):
    """
    Conducts metagenomics variant calling for Nanopore sequencing data. This function orchestrates
    the whole process from setting up the output directory, running KMA alignments, extracting reads,
    and identifying significant variants.

    Args:
        arguments: Parsed command-line arguments containing parameters and file paths.
    """
    # Set up output directory and verify input file

    auto_cor, auto_iteration_increase, auto_pp, auto_np, auto_dp = load_parameters(arguments.maf)

    arguments = initialize_parameters(arguments, auto_cor, auto_iteration_increase, auto_pp, auto_np, auto_dp)

    set_up_output_and_check_input(arguments)

    # Run KMA alignment for bacteria mapping
    kma.KMARunner(arguments.nanopore,
                  os.path.join(arguments.output, "bacteria_mapping"),
                  os.path.join(arguments.db_dir, "bac_db"),
                  f"-mem_mode -Sparse -ss c -t {arguments.threads}").run()

    # Identify the highest scoring bacterial template
    highest_scoring_template = highest_scoring_hit(os.path.join(arguments.output, "bacteria_mapping.spa"))
    primary_specie = ' '.join(highest_scoring_template.split()[1:3])

    # Produce a species-specific KMA database
    produce_specie_specific_kma_db(primary_specie,
                                   os.path.join(arguments.db_dir, 'rmlst.fsa'),
                                   os.path.join(arguments.db_dir, 'rmlst_scheme.txt'),
                                   arguments.output)

    # Run KMA alignment for initial rMLST
    kma.KMARunner(arguments.nanopore,
                  os.path.join(arguments.output, "initial_rmlst_alignment"),
                  os.path.join(arguments.output, 'specie_db'),
                  f"-t {arguments.threads} -ID 10 -ont -md 1.5 -matrix -eq {arguments.q_score} -mct 0.75 -sam 2096 > {os.path.join(arguments.output, 'initial_rmlst_alignment.sam')}").run()

    # Index top hits from the initial RMLST alignment
    index_top_hits_db(arguments.output)


    # Run KMA alignment for rMLST
    kma.KMARunner(arguments.nanopore,
                  os.path.join(arguments.output, "rmlst_alignment"),
                  os.path.join(arguments.output, 'top_hits_db'),
                  f"-t {arguments.threads} -ID 10 -ont -md 1.5 -eq {arguments.q_score} -matrix -mct 0.75 -sam 2096 > {os.path.join(arguments.output, 'rmlst_alignment.sam')}").run()

    
    os.system(f'gunzip {os.path.join(arguments.output, "rmlst_alignment.mat.gz")}')

    # Build a consensus dictionary from alignment results
    consensus_dict = build_consensus_dict(os.path.join(arguments.output, 'rmlst_alignment.res'),
                                          os.path.join(arguments.output, 'rmlst_alignment.mat'))

    print_majority_alelles(consensus_dict, arguments)

    if arguments.majority_alleles_only: #End the program if only majority alleles are requested
        sys.exit()

    # Adjust the consensus dictionary based on individual quality scores
    # Currently not doing anything.
    consensus_dict, read_positions_blacklisted_dict = adjust_consensus_dict_for_individual_qscores(consensus_dict, os.path.join(arguments.output, 'rmlst_alignment.sam'), arguments.nanopore, arguments.q_score)

    # Derive mutation positions from consensus data
    confirmed_mutation_dict = derive_mutation_positions(consensus_dict, arguments)

    # Perform biological validation of mutations
    bio_validation_dict = bio_validation_mutations(consensus_dict, os.path.join(arguments.output, 'specie.fsa'))

    # Co-occurrence analysis until convergence
    confirmed_mutation_dict, co_occurrence_tmp_dict, iteration_count = co_occurrence_until_convergence(arguments, confirmed_mutation_dict, consensus_dict, read_positions_blacklisted_dict, bio_validation_dict)

    for item in confirmed_mutation_dict:
        print(item, confirmed_mutation_dict[item])

    print_minor_variants(confirmed_mutation_dict, consensus_dict, arguments.output)
    ## Format and output the results
    format_output(arguments, confirmed_mutation_dict, consensus_dict, bio_validation_dict, co_occurrence_tmp_dict)

    # Write majority sequences to file
    with open(os.path.join(arguments.output, 'majority_seqs.fasta'), 'w') as f:
        for allele in consensus_dict:
            print(f'>{allele}', file=f)
            print(consensus_dict[allele][1], file=f)

    sys.exit()


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
    consensus_dict = build_consensus_dict(os.path.join(arguments.output, 'rmlst_alignment.res'),
                                          os.path.join(arguments.output, 'rmlst_alignment.mat'))

    confirmed_mutation_dict = derive_mutation_positions(consensus_dict, arguments)

    # Perform biological validation of mutations
    #bio_validation_dict = bio_validation_mutations(consensus_dict, os.path.join(results_folder, 'specie.fsa'))
    # Co-occurrence analysis until convergence
    #confirmed_mutation_dict, co_occurrence_tmp_dict, iteration_count =\
    #    co_occurrence_until_convergence(arguments, confirmed_mutation_dict,
    #                                    consensus_dict, {}, bio_validation_dict)

    # Format and output the results
    #format_output(new_output_folder, confirmed_mutation_dict, consensus_dict, bio_validation_dict,
    #              co_occurrence_tmp_dict)

    sample = results_folder.split('/')[-1]

    minor_mutation_expected = benchmark_analysis_result(sample)

    minor_mutation_results = convert_mutation_dict_to_object(confirmed_mutation_dict)

    precision, recall, f1, tp, fp, fn = calculate_metrics(minor_mutation_expected, minor_mutation_results)

    parameter_string = f"maf_{maf}_cor_{cor}_pp_{pp}_np_{np}_dp_{dp}_iteration_increase_{iteration_increase}"

    precision = 1
    recall = 1
    f1 = 1
    tp = 1
    fp = 1
    fn = 1
    parameter_string = 'test'

    return f1, parameter_string, precision, recall, tp, fp, fn
def initialize_parameters(arguments, auto_cor, auto_iteration_increase, auto_pp, auto_np, auto_dp):
    if arguments.cor == 'auto':
        arguments.cor = auto_cor
    else:
        arguments.cor = float(arguments.cor)

    if arguments.iteration_increase == 'auto':
        arguments.iteration_increase = auto_iteration_increase
    else:
        arguments.iteration_increase = float(arguments.iteration_increase)

    if arguments.pp == 'auto':
        arguments.pp = auto_pp
    else:
        arguments.pp = float(arguments.pp)

    if arguments.np == 'auto':
        arguments.np = auto_np
    else:
        arguments.np = float(arguments.np)

    if arguments.dp == 'auto':
        arguments.dp = auto_dp
    else:
        arguments.dp = float(arguments.dp)

    return arguments


def load_parameters(maf_value):

    # Names of the parameters and their corresponding JSON files
    parameters = {
        'cor': cf.load_cor(),
        'iteration': cf.load_ii(),
        'pp': cf.load_pp(),
        'np': cf.load_np(),
        'dp': cf.load_dp()
    }

    splines = {}

    # Load splines from JSON files
    for param, filename in parameters.items():
        splines[param] = load_spline_from_json(parameters[param])

    # Calculate and print the parameters values for the given MAF
    results = {}
    for param, spline in splines.items():
        value = calculate_parameter_value(spline, maf_value)
        results[param] = value


    return results['cor'], results['iteration'], results['pp'], results['np'], results['dp']

def load_spline_from_json(data):
    # Convert keys back to floats and sort by keys
    x_values = np.array(sorted(map(float, data.keys())))
    y_values = np.array([data[str(x)] for x in x_values])
    spline = UnivariateSpline(x_values, y_values, s=None)
    return spline

def calculate_parameter_value(spline, maf):
    """ Calculate the parameter value for a given maf using the spline. """
    value = spline(maf)
    # Ensure that the value is non-negative
    return max(value, 0)


def print_majority_alelles(consensus_dict, arguments):
    nucleotides = ['A', 'C', 'G', 'T']

    #TBD Change this maybe, not really a vcf output. Maybe simply to a fasta file?
    with open(arguments.output + '/majority_variants.vcf', 'w') as file:
        # Print header lines for VCF format and describe the INFO field contents
        print("##fileformat=VCFv4.2", file=file)
        print("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">", file=file)
        print("##INFO=<ID=MD,Number=1,Type=Integer,Description=\"Depth of Majority Variant\">", file=file)
        print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", file=file)

        for item in consensus_dict:
            allele = item
            depths = consensus_dict[item][0]
            for i, sequence_position in enumerate(depths):
                nucleotide_frequencies = dict(zip(nucleotides, sequence_position))
                majority_nucleotide = max(nucleotide_frequencies, key=nucleotide_frequencies.get)
                total_depth = sum(nucleotide_frequencies.values())
                majority_depth = nucleotide_frequencies[majority_nucleotide]

                # Assuming the allele identifier as CHROM and additional fields as described
                chrom = allele
                pos = i + 1
                ref = consensus_dict[item][1][i]  # Assuming this is the correct reference base
                alt = majority_nucleotide
                qual = 'NA'
                filter_status = "PASS"

                # Constructing the INFO field content
                info = f"DP={total_depth};MD={majority_depth}"

                # Printing the line with the INFO field included
                print(f"{chrom}\t{pos}\t{allele}\t{ref}\t{alt}\t{qual}\t{filter_status}\t{info}", file=file)


def print_minor_variants(confirmed_mutation_dict, consensus_dict, output_path):
    with open(f'{output_path}/minor_variants.vcf', 'w') as file:
        # Print header lines for VCF format
        print("##fileformat=VCFv4.2", file=file)
        print("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">", file=file)
        print("##INFO=<ID=MD,Number=1,Type=Integer,Description=\"Depth of Majority Variant\">", file=file)
        print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", file=file)

        for item in confirmed_mutation_dict:
            allele = item
            mutations = confirmed_mutation_dict[item][0]
            depths = confirmed_mutation_dict[item][1]

            # Ensure allele exists in consensus_dict before proceeding
            if allele not in consensus_dict:
                continue

            sequence = consensus_dict[allele][1]  # The nucleotide sequence for the allele

            for mutation, depth in zip(mutations, depths):
                position_str, alt = mutation.split('_')
                pos = int(position_str)
                chrom = allele
                total_depth = sum(consensus_dict[allele][0][pos - 1])

                # Use the position from the mutation to get the reference nucleotide from the sequence
                # Note: Assuming that the positions in mutations are 1-based indexing
                try:
                    ref = sequence[pos - 1]  # Adjusting for zero-based indexing
                except IndexError:
                    # If the position is out of range, set REF as unknown
                    ref = "N"

                qual = 'NA'  # Placeholder for quality
                filter_status = "PASS"
                info = f"DP={total_depth};MD={depth}"


                # Printing the line with the INFO field included
                print(f"{chrom}\t{pos}\t{allele}\t{ref}\t{alt}\t{qual}\t{filter_status}\t{info}", file=file)

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


def highest_scoring_hit(file_path):
    """
    Identifies and returns the highest scoring template from a tab-separated file.

    This function reads through a specified file, assuming it to be tab-separated,
    and identifies the row with the highest score in a designated score column.
    It returns the template corresponding to this highest score.

    Args:
        file_path (str): The path to the file to be read. Assumes a specific format where
                         the score is in the third column and the template in the first.

    Returns:
        str: The identifier of the highest scoring template.
    """

    highest_score = 0
    highest_scoring_template = ""

    with open(file_path, 'r') as file:
        next(file)  # Skip the header line
        for line in file:
            columns = line.split('\t')
            try:
                # Extract score and compare to find the highest
                score = int(columns[2])  # Score is expected in the 3rd column
                if score > highest_score:
                    highest_score = score
                    highest_scoring_template = columns[0]  # Template is expected in the 1st column
            except ValueError:
                # Skip line if score is not an integer or line is malformed
                continue

    return highest_scoring_template

def co_occurrence_until_convergence(arguments, confirmed_mutation_dict, consensus_dict, read_positions_blacklisted_dict, bio_validation_dict):
    """
    Executes an iterative process to identify co-occurring mutations in reads until convergence is achieved.
    The process adjusts the parameters for correlation and density penalty in each iteration and checks for
    the stabilization of the number of mutations, indicating convergence.

    Args:
        arguments (object): An object containing various parameters for the function, including output path and penalty adjustments.
        confirmed_mutation_dict (dict): A dictionary of confirmed mutations.
        consensus_dict (dict): A dictionary containing consensus data.
        read_positions_blacklisted_dict (dict): A dictionary with blacklisted read positions.
        bio_validation_dict (dict): A dictionary used for biological validation.

    Returns:
        tuple: A tuple containing the updated dictionary of confirmed mutations, temporary co-occurrence data, and iteration count.
    """

    print('loading reads_mutation_dict...')

    reads_mutation_dict = parse_sam_and_find_mutations(arguments.output + '/rmlst_alignment.sam',
                                                       confirmed_mutation_dict,
                                                       consensus_dict,
                                                       read_positions_blacklisted_dict)
    #print (reads_mutation_dict['SRR27755678.258255 BACT000001_1153'])
    print('reads_mutation_dict loaded')

    current_count = count_mutations_in_mutations_dict(confirmed_mutation_dict)
    iteration_count = 0
    original_cor = arguments.cor
    original_dp = arguments.dp
    start_time = time.time()


    with open(arguments.output + '/convergence_results.txt', 'w') as convergence_file:
        print('Interations,Mutations', file=convergence_file)

        while True:
            # Adjust the correlation and density penalty parameters for each iteration
            arguments.cor += original_cor * arguments.iteration_increase # Increase of 20% per iteration
            arguments.dp += original_dp * arguments.iteration_increase  # Increase of 20% per iteration

            # Perform upper co-occurring mutations analysis
            confirmed_mutation_dict, co_occurrence_tmp_dict = upper_co_occuring_mutations_in_reads(
                arguments,
                confirmed_mutation_dict,
                consensus_dict,
                read_positions_blacklisted_dict,
                bio_validation_dict,
                reads_mutation_dict
            )

            new_count = count_mutations_in_mutations_dict(confirmed_mutation_dict)
            iteration_count += 1
            print(f'Iteration: {iteration_count}', file=sys.stderr)
            print(f'Mutations: {new_count}', file=sys.stderr)
            print(f'{iteration_count},{new_count}', file=convergence_file)
            #for item in co_occurrence_tmp_dict:
            #    print(f'{item},{co_occurrence_tmp_dict[item]}', file=sys.stderr)

            # Check for convergence: no change in mutation count
            if new_count == current_count:
                break
            current_count = new_count

    end_time = time.time()
    print(f'Time taken for all iterations: {end_time - start_time} seconds', file=sys.stderr)

    return confirmed_mutation_dict, co_occurrence_tmp_dict, iteration_count



def count_mutations_in_mutations_dict(mutation_dict):
    """
    Counts the total number of mutations present in a dictionary of mutations.

    Args:
        mutation_dict (dict): A dictionary where each key maps to a list containing mutation data,
                              with the mutation information in the first element of the list.

    Returns:
        int: The total count of mutations in the dictionary.
    """
    count = 0
    for key in mutation_dict:
        count += len(mutation_dict[key][0])
    return count


def bio_validation_mutations(consensus_dict, fsa_file):
    """
    Generates a dictionary of mutations for biological validation.

    Args:
        consensus_dict (dict): A dictionary containing consensus information for alleles.
        fsa_file (str): The path to a FASTA file containing sequence data.

    Returns:
        dict: A dictionary mapping genes to sets of mutation strings.
    """
    # Derive correct length headers for sequences
    correct_length_dict = derive_correct_length_headers(consensus_dict, fsa_file)
    mutation_dict = {}

    # Iterate through each gene and sequence to construct mutation strings
    for gene, (_, sequences) in correct_length_dict.items():
        mutation_set = set()
        for sequence in sequences:
            for position, base in enumerate(sequence, start=1):
                mutation_set.add(f"{position}_{base}")
        mutation_dict[gene] = mutation_set

    return mutation_dict



def index_top_hits_db(output_directory):
    """
    Processes initial RMLST alignment results to identify and index top hits. This function reads two files:
    'initial_rmlst_alignment.res' to identify top hits and 'specie_db.name' to match these hits with sequence numbers.
    It then creates a FASTA file of these top hits and indexes it for use with the KMA (K-mer alignment) tool.

    Args:
        output_directory (str): The directory where the input files are located and where the output files will be saved.

    Note:
        This function assumes specific formats for the 'initial_rmlst_alignment.res' and 'specie_db.name' files.
    """
    # Read and process initial alignment results to identify top hits
    initial_alignment_file = os.path.join(output_directory, 'initial_rmlst_alignment.res')
    top_hits = {}
    with open(initial_alignment_file, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                elements = line.strip().split('\t')
                allele, score = elements[:2]
                gene = allele.split('_')[0]

                # Update the top hits for each gene
                if gene not in top_hits or score > top_hits[gene][1]:
                    top_hits[gene] = [allele, score, 0]

    # Match the top hits with their respective sequence number
    species_name_file = os.path.join(output_directory, 'specie_db.name')
    with open(species_name_file, 'r') as file:
        for sequence_number, line in enumerate(file, start=1):
            allele = line.strip()
            gene = allele.split('_')[0]

            if gene in top_hits and top_hits[gene][0] == allele:
                top_hits[gene][-1] = sequence_number

    # Prepare the sequences for FASTA conversion and indexing
    sequence_numbers = ','.join(str(top_hits[gene][-1]) for gene in top_hits)

    # Generate FASTA file from the sequences
    fasta_file = os.path.join(output_directory, 'top_hits.fsa')
    specie_db = os.path.join(output_directory, 'specie_db')
    os.system(f'kma seq2fasta -seqs {sequence_numbers} -t_db {specie_db} > {fasta_file}')

    # Index the generated FASTA file
    top_hits_db = os.path.join(output_directory, 'top_hits_db')
    os.system(f'kma index -i {fasta_file} -o {top_hits_db} 2>/dev/null')

def adjust_consensus_dict_for_individual_qscores(consensus_dict, sam_file, fastq_file, q_score):
    """
    Adjusts the consensus dictionary based on individual quality scores derived from SAM and FASTQ files.
    This adjustment involves blacklisting positions with low quality scores and updating the consensus
    dictionary with mutations from high-quality positions.

    Args:
        consensus_dict (dict): Dictionary containing initial consensus data.
        sam_file (str): Path to the SAM file containing alignment information.
        fastq_file (str): Path to the FASTQ file containing sequence data and quality scores.

    Returns:
        tuple: A tuple containing the adjusted consensus dictionary and a dictionary of blacklisted positions.
    """

    # Blacklist low-quality positions based on quality scores from the FASTQ file
    #black_listed_positions = blacklist_positions(fastq_file, quality_threshold=q_score)
    black_listed_positions = blacklist_positions(fastq_file, quality_threshold=1) #Decided not to use quality threshold. Revisit another time. Invidivial positions are not being blacklisted due to their respective quality scores.

    # Initialize the adjusted consensus dictionary
    adjusted_consensus_dict = {}

    # Setting up the structure of the adjusted consensus dictionary
    for allele, (positions, allele_seq) in consensus_dict.items():
        # Initialize nucleotide counts for each position
        adjusted_consensus_dict[allele] = [[[0, 0, 0, 0, 0, 0] for _ in positions], allele_seq]

    # Process alignments from the SAM file
    with open(sam_file, 'r') as file:
        for alignment in file:
            if alignment.startswith('@'):  # Skip header lines
                continue

            # Extract relevant data from the alignment
            qname, _, rname, pos_str, _, cigar_str, _, _, _, seq = alignment.strip().split('\t')[:10]
            read_id = qname.split(' ')[0]
            position = int(pos_str)
            template_seq = adjusted_consensus_dict[rname][1]

            # Process only full-length reads
            if position == 1 and len(seq) >= len(template_seq):
                aligned_ref, aligned_query = extract_alignment(template_seq, seq, cigar_str)
                mutation_vector = create_mutation_vector(aligned_ref, aligned_query)

                # Update the consensus dictionary with mutation data
                for i, mutation in enumerate(mutation_vector):
                    if i not in black_listed_positions.get(read_id, []):
                        adjusted_consensus_dict[rname][0][i][mutation_to_index(mutation)] += 1

    return adjusted_consensus_dict, black_listed_positions

def mutation_to_index(mutation):
    """
    Converts a mutation character to its corresponding index.

    Args:
        mutation (str): A single character representing a mutation.

    Returns:
        int: The index corresponding to the mutation.
    """
    nucleotide_list = ['A', 'C', 'G', 'T', 'N', '-']
    return nucleotide_list.index(mutation)

#Here
def blacklist_positions(fastq_file, quality_threshold):
    """
    Generates a dictionary of low-quality positions for each read in a FASTQ file.

    Args:
        fastq_file (str): The path to the FASTQ file containing sequencing data and quality scores.
        quality_threshold (int): A threshold below which the quality score is considered too low.

    Returns:
        dict: A dictionary mapping each read ID to a list of blacklisted (low-quality) positions.
    """
    blacklist_dict = {}

    # Parsing each read in the FASTQ file
    for record in SeqIO.parse(fastq_file, "fastq"):
        # Identifying positions with quality scores below the threshold
        blacklist = [pos for pos, score in enumerate(record.letter_annotations["phred_quality"]) if score < quality_threshold]

        # Adding identified positions to the blacklist dictionary
        blacklist_dict[record.id] = blacklist

    return blacklist_dict

def format_output_for_plots(arguments, confirmed_mutation_dict, consensus_dict, bio_validation_dict, co_occurrence_tmp_dict):
    """
    Format and print the output of confirmed mutations with additional information.

    Args:
        confirmed_mutation_dict (dict): A dictionary containing confirmed mutations for alleles.
        consensus_dict (dict): A dictionary containing consensus information for alleles.
        bio_validation_dict (dict): A dictionary containing biological validation data for genes.

    Returns:
        None
    """
    with open(arguments.output + '/minor_mutations.csv', 'w') as outfile:
        header = 'Gene,Position,MajorityBase,MutationBase,MutationDepth,TotalDepth,GeneLength,MutationComment,CoOccurrence'
        print(header, file = outfile)

        for allele in confirmed_mutation_dict:
            for mutation in zip(confirmed_mutation_dict[allele][0], confirmed_mutation_dict[allele][1]):
                position = mutation[0].split('_')[0]
                mutation_base = mutation[0].split('_')[1]
                mutation_depth = mutation[1]
                majority_base = consensus_dict[allele][1][int(position) - 1]
                total_depth = sum(consensus_dict[allele][0][int(position) - 1])
                biological_existence = check_single_mutation_existence(bio_validation_dict, allele, mutation[0])
                gene_length = len(consensus_dict[allele][1])
                if mutation[0] in co_occurrence_tmp_dict[allele]:
                    co_occurrence = 'Yes'
                else:
                    co_occurrence = 'No'

                if biological_existence:
                    print('{},{},{},{},{},{},{},{},{}'.format(allele, position, majority_base, mutation_base, mutation_depth, total_depth, gene_length, 'Mutation seen in database', co_occurrence), file = outfile)
                else:
                    print('{},{},{},{},{},{},{},{},{}'.format(allele, position, majority_base, mutation_base, mutation_depth, total_depth, gene_length, 'Novel mutation', co_occurrence), file = outfile)

def format_output(arguments, confirmed_mutation_dict, consensus_dict, bio_validation_dict, co_occurrence_tmp_dict):
    """
    Format and print the output of confirmed mutations with additional information.

    Args:
        confirmed_mutation_dict (dict): A dictionary containing confirmed mutations for alleles.
        consensus_dict (dict): A dictionary containing consensus information for alleles.
        bio_validation_dict (dict): A dictionary containing biological validation data for genes.

    Returns:
        None
    """
    with open(arguments.output + '/minor_mutations.csv', 'w') as outfile:
        header = 'Gene,Position,MajorityBase,MutationBase,MutationDepth,TotalDepth,GeneLength,MutationComment,CoOccurrence'
        print(header, file = outfile)

        for allele in confirmed_mutation_dict:
            for mutation in zip(confirmed_mutation_dict[allele][0], confirmed_mutation_dict[allele][1]):
                position = mutation[0].split('_')[0]
                mutation_base = mutation[0].split('_')[1]
                mutation_depth = mutation[1]
                majority_base = consensus_dict[allele][1][int(position) - 1]
                total_depth = sum(consensus_dict[allele][0][int(position) - 1])
                biological_existence = check_single_mutation_existence(bio_validation_dict, allele, mutation[0])
                gene_length = len(consensus_dict[allele][1])
                if mutation[0] in co_occurrence_tmp_dict[allele]:
                    co_occurrence = 'Yes'
                else:
                    co_occurrence = 'No'

                if biological_existence:
                    print('{},{},{},{},{},{},{},{},{}'.format(allele, position, majority_base, mutation_base, mutation_depth, total_depth, gene_length, 'Mutation seen in database', co_occurrence), file = outfile)
                else:
                    print('{},{},{},{},{},{},{},{},{}'.format(allele, position, majority_base, mutation_base, mutation_depth, total_depth, gene_length, 'Novel mutation', co_occurrence), file = outfile)

def extract_mapped_rmlst_read(output_directory, nanopore_fastq):
    """
    Extract and trim mapped rMLST reads from an initial alignment file.

    Args:
        output_directory (str): The directory where output files will be saved.
        nanopore_fastq (str): The path to the nanopore FASTQ file.

    Returns:
        None
    """
    read_set = set()
    total = 0

    # Extract read IDs from the initial rMLST alignment file
    with open(output_directory + '/initial_rmlst_alignment.sam', 'r') as sam_file:
        for line in sam_file:
            if not line.startswith('@'):
                read_set.add(line.split('\t')[0])
                total += 1

    # Write the extracted read IDs to a text file
    with open(output_directory + '/rmlst_reads.txt', 'w') as f:
        for item in read_set:
            f.write(item + '\n')

    print ('Number of reads extracted: ', len(read_set))
    print ('Total number of alignments in sam file: ', total)

    # Use seqtk to extract the corresponding reads from the nanopore FASTQ file
    os.system('seqtk subseq {} {} > {}'.format(nanopore_fastq, output_directory + '/rmlst_reads.txt',
                                               output_directory + '/trimmed_rmlst_reads.fastq'))


def derive_mutation_positions(consensus_dict, arguments):
    """
    Derive mutation positions and their depths from a consensus dictionary.

    Args:
        consensus_dict (dict): A dictionary containing consensus information for alleles.
        arguments: Arguments containing parameters for mutation derivation.

    Returns:
        dict: A dictionary containing derived mutation positions and depths for each allele.
    """
    all_confirmed_mutation_dict = {}

    for allele, allele_data in consensus_dict.items():
        all_confirmed_mutation_dict[allele] = [[], []]

        for i in range(len(allele_data[0])):
            positions = allele_data[0][i][:4]
            max_number = max(positions)
            index_of_max = positions.index(max_number)
            nucleotide_index = ['A', 'C', 'G', 'T']

            for t in range(len(positions)):
                if t != index_of_max:
                    if positions[t] >= arguments.min_n:
                        total_depth = sum(positions)
                        relative_depth = positions[t] / total_depth

                        if relative_depth >= arguments.maf - (arguments.maf * arguments.cor):
                            # Only consider mutations with minimum depth >= 2
                            if nucleotide_index[t] != 'N' and nucleotide_index[t] != '-': #We don't consider there SNVs
                                all_confirmed_mutation_dict[allele][0].append(
                                    '{}_{}'.format(i + 1, nucleotide_index[t]))
                                all_confirmed_mutation_dict[allele][1].append(positions[t])

    return all_confirmed_mutation_dict



def upper_co_occuring_mutations_in_reads(arguments, confirmed_mutation_dict, consensus_dict,
                                         read_positions_blacklisted_dict, bio_validation_dict, reads_mutation_dict):
    """
    Filter and adjust confirmed mutations based on co-occurrence, depth, and biological validation.

    Args:
        arguments: Arguments containing parameters for filtering.
        confirmed_mutation_dict (dict): A dictionary containing confirmed mutations for alleles.
        consensus_dict (dict): A dictionary containing consensus information for alleles.
        read_positions_blacklisted_dict (dict): A dictionary of blacklisted positions.
        bio_validation_dict (dict): A dictionary containing biological validation data for genes.

    Returns:
        dict: A filtered and adjusted mutation dictionary for alleles.
    """

    co_occurrence_tmp_dict = {}
    co_occurrence_matrix_dict = {}
    for allele in confirmed_mutation_dict:
        mutation_list = confirmed_mutation_dict[allele][0]
        num_mutations = len(mutation_list)
        if num_mutations > 1:
            co_occurrence_matrix = [[0] * num_mutations for _ in range(num_mutations)]
            for read in reads_mutation_dict:
                read_allele = read.split(' ')[1]
                if read_allele == allele:
                    read_mutations = reads_mutation_dict[read]
                    valid_mutations = [mutation for mutation in read_mutations if mutation in mutation_list]
                    #if allele == 'BACT000001_1153':
                    #    if 'SRR27755678.258255' in read:
                    #        print (read)
                    #        print (valid_mutations)
                    #        print (read_mutations)
                    #        print (mutation_list)
                    if len(valid_mutations) > 1:
                        for i in range(len(valid_mutations)):
                            for j in range(i + 1, len(valid_mutations)):
                                mutation1 = mutation_list.index(valid_mutations[i])
                                mutation2 = mutation_list.index(valid_mutations[j])
                                co_occurrence_matrix[mutation1][mutation2] += 1
                                co_occurrence_matrix[mutation2][mutation1] += 1
            co_occurrence_matrix_dict[allele] = [co_occurrence_matrix, mutation_list]
            #if allele == 'BACT000001_1153':
            #    print (allele)
            #    for i in range(len(co_occurrence_matrix)):
            #        print (mutation_list[i], co_occurrence_matrix[i])
    adjusted_mutation_dict = {}
    for allele in confirmed_mutation_dict:
        co_occurrence_tmp_dict[allele] = []
        if allele in co_occurrence_matrix_dict:
            adjusted_mutation_dict[allele] = [[], []]
            matrix = co_occurrence_matrix_dict[allele][0]
            mutation_list = co_occurrence_matrix_dict[allele][1]
            depth_list = confirmed_mutation_dict[allele][1]
            for i in range(len(matrix)):
                row = matrix[i]
                mutation = mutation_list[i]
                position = int(mutation.split('_')[0])
                position_depth = sum(consensus_dict[allele][0][position - 1])
                mutation_depth = depth_list[i]
                proxi_mutations = find_mutations_proximity_specific_mutation(mutation_list, mutation, arguments.proxi)
                density_mutations = find_mutations_proximity_specific_mutation(mutation_list, mutation, arguments.dp_window)
                biological_existence = check_single_mutation_existence(bio_validation_dict, allele, mutation)

                #TBD LOOK INTO THIS
                mutation_threshold = position_depth * arguments.maf
                co_occurrence_list = check_mutation_co_occurrence(row, mutation_list, mutation,
                                                                 position_depth, arguments.cor, arguments.pp, arguments.maf, proxi_mutations, mutation_depth)
                #if allele == 'BACT000001_1153':
                #    print (mutation, co_occurrence_list)
                #    print (mutation_threshold, mutation_depth)
                if co_occurrence_list != []:
                    if mutation not in co_occurrence_tmp_dict[allele]:
                        co_occurrence_tmp_dict[allele].append(mutation)
                    for item in co_occurrence_list:
                        if item not in co_occurrence_tmp_dict[allele]:
                            co_occurrence_tmp_dict[allele].append(item)
                    mutation_threshold = mutation_threshold - position_depth * arguments.maf * arguments.cor

                if not biological_existence:
                    mutation_threshold = mutation_threshold + arguments.np * position_depth * arguments.maf
                if proxi_mutations != []:
                    mutation_threshold = mutation_threshold + arguments.pp * position_depth * arguments.maf
                if density_mutations != []:
                    mutation_threshold = mutation_threshold + arguments.dp * position_depth * arguments.maf * len(density_mutations)
                if mutation_depth >= mutation_threshold:
                    adjusted_mutation_dict[allele][0].append(confirmed_mutation_dict[allele][0][i])
                    adjusted_mutation_dict[allele][1].append(confirmed_mutation_dict[allele][1][i])


        else:
            adjusted_mutation_dict[allele] = [[], []]
            if confirmed_mutation_dict[allele][0] != []:
                mutation = confirmed_mutation_dict[allele][0][0]
                position = int(mutation.split('_')[0])
                position_depth = sum(consensus_dict[allele][0][position - 1])
                mutation_threshold = position_depth * arguments.maf
                depth = confirmed_mutation_dict[allele][1][0]
                biological_existence = check_single_mutation_existence(bio_validation_dict, allele, mutation)
                if not biological_existence:
                    mutation_threshold = mutation_threshold + (arguments.np-1) * position_depth * arguments.maf

                if depth >= mutation_threshold:
                    adjusted_mutation_dict[allele][0].append(confirmed_mutation_dict[allele][0][0])
                    adjusted_mutation_dict[allele][1].append(confirmed_mutation_dict[allele][1][0])
    return adjusted_mutation_dict, co_occurrence_tmp_dict



def check_mutation_co_occurrence(list_of_mutation_co_occurrence, mutation_list, mutation,
                                 position_depth, correlation_coefficient, proximity_penalty, maf, proximity_mutations, mutation_depth):
    """
    Check for co-occurrence of a mutation with other mutations in a list.

    Args:
        list_of_mutation_co_occurrence (list): A list of co-occurrence counts for each mutation.
        mutation_list (list): A list of mutations.
        mutation (str): The mutation for which co-occurrence is being checked.
        position_depth (int): The depth at which the mutation occurs.
        correlation_coefficient (float): The correlation coefficient used for threshold calculation.
        proximity_penalty (float): The penalty factor for mutations within proximity.
        maf (float): The mutation rate difference.
        proximity_mutations (list): A list of mutations within proximity.

    Returns:
        list: A list of mutations that co-occur with the given mutation.
    """
    if mutation not in mutation_list:
        # Should never happen
        return []  # No co-occurrence and not in proximity

    co_threshold = mutation_depth * (1/2)
    if co_threshold < 3:
        co_threshold = 3

    # Find the index of the mutation in the mutation list
    mutation_index = mutation_list.index(mutation)

    co_occurrence_list = []
    # Check if the co-occurrence count of the mutation with any other mutation is above the threshold
    for i, count in enumerate(list_of_mutation_co_occurrence):
        if mutation_list[i] in proximity_mutations:
            # Add penalty for proximity to make it harder to get the co-occurrence reward
            # for mutations within the proximity
            co_threshold = co_threshold * proximity_penalty
        if i != mutation_index and count >= co_threshold:
            co_occurrence_list.append(mutation_list[i])

    return co_occurrence_list

def check_single_mutation_existence(bio_validation_dict, allele, specific_mutation):
    """
    Check if a specific mutation exists for a given allele and gene in a biological validation dictionary.

    Args:
        bio_validation_dict (dict): A dictionary containing biological validation data for genes.
        allele (str): The allele for which the existence of a specific mutation is checked.
        specific_mutation (str): The specific mutation to check for.

    Returns:
        bool: True if the specific mutation exists for the allele, False otherwise.
    """
    gene = allele.split('_')[0]

    if specific_mutation in bio_validation_dict.get(gene, []):
        return True

    return False

def find_mutations_proximity_specific_mutation(mutations, specific_mutation, proximity):
    """
    Find mutations that are in proximity (within a certain number of positions) to a specific mutation.

    Args:
        mutations (list): A list of mutation strings in the format "position_base".
        specific_mutation (str): The specific mutation to which proximity is determined.
        proximity (int): The maximum number of positions for mutations to be considered in proximity.

    Returns:
        list: A list of mutations that are in proximity to the specific mutation.
    """
    specific_mutation_pos = int(specific_mutation.split('_')[0])
    proximity_mutations = []

    # Split mutations into position and base, and convert positions to integers
    split_mutations = [(int(mutation.split('_')[0]), mutation) for mutation in mutations]

    for pos, mutation in split_mutations:
        # Check if the mutation is within 'proximity' positions of the specific mutation
        if abs(pos - specific_mutation_pos) <= proximity:
            if mutation != specific_mutation:
                proximity_mutations.append(mutation)

    return proximity_mutations

def derive_correct_length_headers(consensus_dict, fsa_file):
    """
    Derive correct length headers and sequences from a FASTA file based on a consensus dictionary.

    Args:
        consensus_dict (dict): A dictionary containing consensus information for alleles.
        fsa_file (str): The path to the input FASTA file.

    Returns:
        dict: A dictionary mapping gene names to correct length sequences.
    """
    correct_length_dict = {}

    for allele in consensus_dict:
        gene = allele.split('_')[0]
        correct_length_dict[gene] = [len(consensus_dict[allele][0]), []]

    sequence = ''
    gene = None

    with open(fsa_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if gene is not None:
                    if sequence != '':
                        if len(sequence) == correct_length_dict[gene][0]:
                            correct_length_dict[gene][1].append(sequence)
                header = line.strip()[1:]
                allele = header
                gene = allele.split('_')[0]
                sequence = ''
            else:
                sequence += line.strip()

    if gene is not None:
        if sequence != '':
            if len(sequence) == correct_length_dict[gene][0]:
                correct_length_dict[gene][1].append(sequence)

    return correct_length_dict

def produce_specie_specific_kma_db(specie, fsa_file, scheme_file, output_directory):
    """
    Produce a specie-specific KMA database using the provided FASTA file and scheme file.

    Args:
        specie (str): The target specie for which the KMA database is being created.
        fsa_file (str): The path to the input FASTA file.
        scheme_file (str): The path to the scheme file containing gene information.
        output_directory (str): The directory where the output KMA database will be created.

    Returns:
        None
    """
    gene_set = set()
    t = 0

    with open(scheme_file, 'r') as f:
        for line in f:
            if line.startswith('rST'):
                headers = line.strip().split('\t')[1:54]
            else:
                if line.strip().split('\t') != ['']:
                    if line.strip().split('\t')[55] == specie:
                        t += 1
                        for i in range(len(headers)):
                            allele = headers[i] + '_' + line.strip().split('\t')[i + 1]
                            gene_set.add(allele)

    # Create a specie-specific FASTA file with the selected genes
    produce_specie_fsa_file(fsa_file, gene_set, output_directory)

    # Create a KMA database from the specie-specific FASTA file
    os.system('kma index -i {}/specie.fsa -o {}/specie_db 2>/dev/null'.format(output_directory, output_directory))


def produce_specie_fsa_file(fsa_file, gene_set, output_directory):
    """
    Produce a specie-specific FASTA file containing sequences for genes in the given gene set.

    Args:
        fsa_file (str): The path to the input FASTA file.
        gene_set (set): A set containing gene IDs to include in the output.
        output_directory (str): The directory where the output file will be saved.

    Returns:
        None
    """
    output_file = output_directory + '/specie.fsa'

    with open(output_file, 'w') as outfile:
        with open(fsa_file, 'r') as f:
            write_sequence = False  # A flag to indicate whether the sequence should be written to the output
            for line in f:
                if line.startswith('>'):
                    # Check if the gene_id (without '>') is in the gene_set
                    gene_id = line.strip().split()[0][1:]
                    write_sequence = gene_id in gene_set
                # Write the line (header or sequence) if write_sequence is True
                if write_sequence:
                    outfile.write(line)

#Here
def build_consensus_dict(res_file, mat_file):
    """
    Build a consensus dictionary from result and matrix files.

    Args:
        res_file (str): The name of the result file.
        mat_file (str): The name of the matrix file.

    Returns:
        dict: A dictionary containing consensus information for alleles.
    """
    consensus_dict = {}

    # Read and process the result file
    with open(res_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                line = line.strip().split('\t')
                allele = line[0]
                length = int(line[3])
                consensus_dict[allele] = [[], '']
                for i in range(length):
                    consensus_dict[allele][0].append([0, 0, 0, 0, 0, 0])  # [A, C, G, T, N, -]

    # Read and process the matrix file
    with open(mat_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#'):
                allele = line[1:].strip()
                index = 0
            elif line != '':
                line = line.split('\t')
                if line[0] != '-':  # Excludes read gaps. Reconsider?
                    line = line[1:]
                    for i in range(len(line)):
                        consensus_dict[allele][0][index][i] += int(line[i])
                    index += 1

    # Generate consensus sequences for alleles
    for allele in consensus_dict:
        for position in consensus_dict[allele][0]:
            consensus_dict[allele][1] += 'ACGT'[position[:4].index(max(position[:4]))]

    return consensus_dict



def number_of_bases_in_file(filename):
    """
    Calculate the total number of bases in a FASTA or FASTQ file.

    Args:
        filename (str): The name of the input file.

    Returns:
        int: The total number of bases in the file.
    """
    gzipped, file_type = determine_file_type(filename)

    if file_type == 'fasta':
        total_bases = 0
        with open(filename, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    total_bases += len(line.strip())
        return total_bases

    elif file_type == 'fastq':
        total_bases = 0
        line_count = 1
        with (gzip.open(filename, 'r') if gzipped else open(filename, 'r')) as f:
            for line in f:
                if line_count == 2:
                    total_bases += len(line.strip())
                line_count += 1
                if line_count == 5:
                    line_count = 1
        return total_bases

def eval_bacteria_results(results, total_bacteria_aligning_bases):
    """
    Evaluate bacterial alignment results to determine primary and candidate results.

    Args:
        results (list of dict): List of dictionaries containing bacterial alignment results.
        total_bacteria_aligning_bases (int): Total number of bases aligning to bacteria.

    Returns:
        tuple: A tuple containing the primary result and a dictionary of candidate results.
    """
    primary = results[0]['#Template']
    candidate_dict = dict()

    for item in results:
        # Calculate the number of bases hit based on depth
        bases_hit = int(item['Template_length']) * float(item['Depth'].strip())

        # Calculate the template identity
        template_id = float(item['Template_Identity'])

        # Calculate relative template depth
        relative_template_depth = bases_hit / total_bacteria_aligning_bases

        # Check if the result qualifies as a candidate
        if relative_template_depth > 0.01 or template_id > 20.00:
            candidate_dict[item['#Template']] = [relative_template_depth, int(item['Template_length'])]

    return primary, candidate_dict




def determine_file_type(file):
    """
    Determine the file type and whether it is gzipped.

    Args:
        file (str): The name of the input file.

    Returns:
        tuple: A tuple containing a boolean indicating if the file is gzipped and a string indicating the file type.
    """
    gzipped = False
    file_type = None

    # Check if the file has a '.gz' extension, indicating it is gzipped
    if file.endswith('.gz'):
        gzipped = True
        file = file[:-3]  # Remove the '.gz' extension

    # Check the file type based on its extension
    if file.endswith('.fastq') or file.endswith('.fq'):
        file_type = 'fastq'
    elif file.endswith('.fasta') or file.endswith('.fa') or file.endswith('.fna') or file.endswith('.fsa'):
        file_type = 'fasta'

    return gzipped, file_type


def sort_lines_by_score(filename):
    """
    Sort lines in a tab-delimited file by the 'Score' column in descending order.

    Args:
        filename (str): The name of the input file.

    Returns:
        list: A list of dictionaries representing the sorted data.
    """
    data = []

    # Read data from the input file
    with open(filename, 'r') as file:
        lines = file.readlines()
        header = lines[0].strip().split('\t')

        # Parse and store data as dictionaries
        for line in lines[1:]:
            values = line.strip().split('\t')
            data.append(dict(zip(header, values)))

    # Sort the data by the 'Score' column in descending order
    sorted_data = sorted(data, key=lambda x: int(x['Score']), reverse=True)

    return sorted_data

def benchmark_analysis_result(sample):
    #batch = sample.split('_')[-2]
    #print (sample)
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





def set_up_output_and_check_input(arguments):
    """
    Set up the output directory and check the existence of input files.

    Args:
        arguments: Parsed command-line arguments.
    """
    # Create the output directory if it doesn't exist
    if not os.path.exists(arguments.output):
        os.makedirs(arguments.output)

    # Check if a nanopore input file is provided and if it exists
    if arguments.nanopore is not None:
        if not os.path.isfile(arguments.nanopore):
            print('Input file does not exist')
            sys.exit()