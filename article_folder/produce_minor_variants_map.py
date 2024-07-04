import re
import sys
import os
from Bio import pairwise2
from Bio.Seq import Seq


def calculate_proximity_density(mutations):
    # Extract base positions and mutation values from the set
    mutation_positions = {int(mut.split('_')[0]): mut for mut in mutations}

    # Calculate density for each mutation
    def density(mutation):
        position = int(mutation.split('_')[0])
        count = 0
        for other_position in mutation_positions.keys():
            if abs(other_position - position) <= 15 and other_position != position:
                count += 1
        return count

    # Calculate total density score
    total_density = sum(density(mut) for mut in mutations)

    # Return the total density score and the number of mutations
    return total_density, len(mutations)

def check_single_mutation_existence(bio_validation_dict, gene, specific_mutation):
    """
    Check if a specific mutation exists for a given allele and gene in a biological validation dictionary.

    Args:
        bio_validation_dict (dict): A dictionary containing biological validation data for genes.
        allele (str): The allele for which the existence of a specific mutation is checked.
        specific_mutation (str): The specific mutation to check for.

    Returns:
        bool: True if the specific mutation exists for the allele, False otherwise.
    """
    if specific_mutation in bio_validation_dict.get(gene, []):
        return True

    return False


def derive_correct_length_headers(ref_dict, fsa_file):
    """
    Derive correct length headers and sequences from a FASTA file based on a consensus dictionary.

    Args:
        ref_dict (dict): A dictionary containing consensus information for alleles.
        fsa_file (str): The path to the input FASTA file.

    Returns:
        dict: A dictionary mapping gene names to correct length sequences.
    """
    correct_length_dict = {}

    for gene in ref_dict:
        correct_length_dict[gene] = [len(ref_dict[gene]), []]

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


def bio_validation_mutations(query_dict, fsa_file):
    """
    Generates a dictionary of mutations for biological validation.

    Args:
        query_dict (dict): A dictionary containing consensus information for alleles.
        fsa_file (str): The path to a FASTA file containing sequence data.

    Returns:
        dict: A dictionary mapping genes to sets of mutation strings.
    """
    # Derive correct length headers for sequences
    correct_length_dict = derive_correct_length_headers(query_dict, fsa_file)
    mutation_dict = {}

    # Iterate through each gene and sequence to construct mutation strings
    for gene, (_, sequences) in correct_length_dict.items():
        mutation_set = set()
        for sequence in sequences:
            for position, base in enumerate(sequence, start=1):
                mutation_set.add(f"{position}_{base}")
        mutation_dict[gene] = mutation_set

    return mutation_dict

def create_mutation_vector(aligned_ref, aligned_query):
    """
    Creates a mutation vector based on the aligned reference and query sequences.

    Parameters:
    - aligned_ref (str): Aligned portion of the reference sequence.
    - aligned_query (str): Aligned portion of the query sequence.

    Returns:
    - list[str]: A mutation vector representing the alignment with respect to the original reference sequence.
    """
    mutation_vector = []

    # Loop through the aligned sequences
    for i in range(len(aligned_ref)):
        ref_nt = aligned_ref[i]
        query_nt = aligned_query[i]

        if ref_nt == "-":
            # Skip insertions in the reference
            continue
        elif query_nt == "-":
            # Represent deletions in the query with "-"
            mutation_vector.append("-")
        else:
            # Otherwise, represent the nucleotide from the query
            mutation_vector.append(query_nt)

    return mutation_vector


def identify_mutations(mutation_vector, reference_sequence):
    """
    Identify all mutation positions from a mutation vector compared to the reference.

    Parameters:
    - mutation_vector (list[str]): The mutation vector.
    - reference_sequence (str): The original reference sequence.

    Returns:
    - list[str]: A list where each mutation is described as a string in the format "POSITION_NUCLEOTIDE".
    """
    mutations = []

    # Ensure that mutation vector and reference have equal lengths
    if len(mutation_vector) != len(reference_sequence):
        print (len(mutation_vector), len(reference_sequence))
        raise ValueError("The mutation vector and reference sequence must have the same length.")

    for i in range(len(mutation_vector)):
        if mutation_vector[i] != reference_sequence[i] and mutation_vector[i] != "-":
            mutations.append('{}_{}'.format(i+1, mutation_vector[i]))

    return mutations


def align_and_identify_mutations(seq1, seq2):
    # seq1 = qeury, seq2 = ref
    match_score = 2  # Positive score for matches, encourages alignment to match characters
    mismatch_score = -3  # Negative score for mismatches, penalizes alignment for mismatching characters
    gap_open = -5  # Negative score for opening a gap, makes the algorithm less willing to introduce gaps
    gap_extend = -2  # Negative score for extending a gap, penalizes extending gaps over single mismatches
    # Perform local alignment with custom scoring
    alignments = pairwise2.align.localms(seq1, seq2, match_score, mismatch_score, gap_open, gap_extend)

    # Assuming we take the first alignment (typically the highest scoring)
    alignment_a, alignment_b, score, begin, end = alignments[0]

    # Identify point mutations, excluding gaps
    #point_mutations = []
    #for i in range(len(alignment_a)):
    #    if alignment_a[i] != alignment_b[i]:
    #        point_mutations.append((i + 1, alignment_a[i], alignment_b[i]))  # i+1 for 1-based indexing

    return alignment_a, alignment_b

def load_majority_seq_from_vcf(vcf_file):
    ref_dict = dict()
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            chrom = line.split('\t')[0].split('_')[0]
            if chrom not in ref_dict:
                ref_dict[chrom] = ''
            if line.split('\t')[3] != '-':
                ref_dict[chrom] += line.split('\t')[3]
    return ref_dict

def load_majority_seqs_from_fasta_file(fasta_file):
    ref_dict = dict()
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                chrom = line[1:].strip().split('_')[0]
                ref_dict[chrom] = ''
            else:
                ref_dict[chrom] += line.strip()
    return ref_dict


def find_close_mutations(mutations):
    # Initialize a set to store mutations that are close to each other
    close_mutations = set()

    # Get the total number of mutations
    num_mutations = len(mutations)

    # Iterate through each mutation
    for i in range(num_mutations):
        # Extract the position of the first mutation
        pos1 = int(mutations[i].split('_')[0])

        # Compare it with each subsequent mutation in the list
        for j in range(i + 1, num_mutations):
            # Extract the position of the second mutation
            pos2 = int(mutations[j].split('_')[0])

            # Check if the positions are within 5 of each other
            if abs(pos1 - pos2) <= 5:
                # If so, add both mutations to the set
                close_mutations.add(mutations[i])
                close_mutations.add(mutations[j])

    # Return the set of close mutations
    return close_mutations

def load_ids_from_directory(directory):
    ids = []
    for filename in os.listdir(directory):
        if filename.startswith("ERR") or filename.startswith("SRR"):
            ids.append(filename.split('_')[1])
    return ids

#ids = load_ids_from_directory(directory)


txt_files = os.listdir('/home/people/malhal/test/mixed_isolates_nanomgt/variant_maps')

for item in txt_files:
    if item.endswith('txt'):
        basename = item[:-4]
        k = basename.split('_')[1]
        j = basename.split('_')[3]
        print (j, k)
        mutations_in_proximity = 0
        novel_mutation_count = 0
        total_mutations = 0
        proximity_mutations = []
        total_proximity_density = 0
        with open('major_{}_minor_{}.txt'.format(k, j), 'w') as write_file:
            query_dict = load_majority_seqs_from_fasta_file(f'/home/people/malhal/test/mixed_isolates_nanomgt/{j}/majority_seqs.fasta')
            ref_dict = load_majority_seqs_from_fasta_file(f'/home/people/malhal/test/mixed_isolates_nanomgt//{k}/majority_seqs.fasta')
            bio_validation_dict = bio_validation_mutations(query_dict, '/home/people/malhal/test/mixed_isolates_nanomgt/' + k + '/specie.fsa')
            for key in query_dict:
                mutations = []
                if key not in ref_dict:
                    print (f"Reference sequence not found for {key}")
                    sys.exit()
                query = query_dict[key]
                ref = ref_dict[key]
                alignment_query, alignment_ref = align_and_identify_mutations(query, ref)
                mutation_vector = create_mutation_vector(alignment_ref, alignment_query)
                mutations = identify_mutations(mutation_vector, ref)


                if mutations != []:
                    total_mutations += len(mutations)
                    print (key, file = write_file)
                    print (",".join(mutations), file = write_file)
                    proxi_mutations = find_close_mutations(mutations)
                    mutations_in_proximity += len(proxi_mutations)
                    if len(proxi_mutations) > 0:
                        proximity_density, num_mutations = calculate_proximity_density(proxi_mutations)
                        total_proximity_density += proximity_density
                    for mutation in mutations:
                        biological_existence = check_single_mutation_existence(bio_validation_dict, key,
                                                                               mutation)
                        if not biological_existence:
                            novel_mutation_count += 1
                            print (key, mutation)
            print (f"Total mutations: {total_mutations}")
            print (f"Novel mutations: {novel_mutation_count}")
            print (f"Mutations in proximity: {mutations_in_proximity}")
            print (f"Total density: {total_proximity_density}")
            if mutations_in_proximity > 0:
                print (f"Average proximity density: {total_proximity_density/mutations_in_proximity}")
