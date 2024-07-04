import re
import os
from Bio import pairwise2
from Bio.Seq import Seq

def calculate_proximity_density(mutations):
    mutation_positions = {int(mut.split('_')[0]): mut for mut in mutations}

    def density(mutation):
        position = int(mutation.split('_')[0])
        count = 0
        for other_position in mutation_positions.keys():
            if abs(other_position - position) <= 15 and other_position != position:
                count += 1
        return count

    total_density = sum(density(mut) for mut in mutations)
    return total_density, len(mutations)

def check_single_mutation_existence(bio_validation_dict, gene, specific_mutation):
    return specific_mutation in bio_validation_dict.get(gene, [])

def derive_correct_length_headers(ref_dict, fsa_file):
    correct_length_dict = {}
    for gene in ref_dict:
        correct_length_dict[gene] = [len(ref_dict[gene]), []]

    sequence = ''
    gene = None
    with open(fsa_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if gene is not None and sequence != '':
                    if len(sequence) == correct_length_dict[gene][0]:
                        correct_length_dict[gene][1].append(sequence)
                header = line.strip()[1:]
                gene = header.split('_')[0]
                sequence = ''
            else:
                sequence += line.strip()

    if gene is not None and sequence != '':
        if len(sequence) == correct_length_dict[gene][0]:
            correct_length_dict[gene][1].append(sequence)

    return correct_length_dict

def bio_validation_mutations(query_dict, fsa_file):
    correct_length_dict = derive_correct_length_headers(query_dict, fsa_file)
    mutation_dict = {}
    for gene, (_, sequences) in correct_length_dict.items():
        mutation_set = set()
        for sequence in sequences:
            for position, base in enumerate(sequence, start=1):
                mutation_set.add(f"{position}_{base}")
        mutation_dict[gene] = mutation_set
    return mutation_dict

def create_mutation_vector(aligned_ref, aligned_query):
    mutation_vector = []
    for i in range(len(aligned_ref)):
        ref_nt = aligned_ref[i]
        query_nt = aligned_query[i]
        if ref_nt == "-":
            continue
        elif query_nt == "-":
            mutation_vector.append("-")
        else:
            mutation_vector.append(query_nt)
    return mutation_vector

def identify_mutations(mutation_vector, reference_sequence):
    mutations = []
    if len(mutation_vector) != len(reference_sequence):
        print(len(mutation_vector), len(reference_sequence))
        raise ValueError("The mutation vector and reference sequence must have the same length.")
    for i in range(len(mutation_vector)):
        if mutation_vector[i] != reference_sequence[i] and mutation_vector[i] != "-":
            mutations.append('{}_{}'.format(i + 1, mutation_vector[i]))
    return mutations

def align_and_identify_mutations(seq1, seq2):
    match_score = 2
    mismatch_score = -3
    gap_open = -5
    gap_extend = -2
    alignments = pairwise2.align.localms(seq1, seq2, match_score, mismatch_score, gap_open, gap_extend)
    alignment_a, alignment_b, _, _, _ = alignments[0]
    return alignment_a, alignment_b

def load_majority_seqs_from_fasta_file(fasta_file):
    ref_dict = {}
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                chrom = line[1:].strip().split('_')[0]
                ref_dict[chrom] = ''
            else:
                ref_dict[chrom] += line.strip()
    return ref_dict

def find_close_mutations(mutations):
    close_mutations = set()
    num_mutations = len(mutations)
    for i in range(num_mutations):
        pos1 = int(mutations[i].split('_')[0])
        for j in range(i + 1, num_mutations):
            pos2 = int(mutations[j].split('_')[0])
            if abs(pos1 - pos2) <= 5:
                close_mutations.add(mutations[i])
                close_mutations.add(mutations[j])
    return close_mutations

def load_ids_from_directory(directory):
    ids = []
    for filename in os.listdir(directory):
        if filename.startswith("ERR") or filename.startswith("SRR"):
            ids.append(filename.split('_')[1])
    return ids

def calculate_stats(txt_files, path_prefix, specific_ids=None):
    total_mutations = 0
    total_proximity_mutations = 0
    total_novel_mutations = 0
    total_cooccurring_mutations = 0
    total_proximity_density = 0
    total_pairs = 0
    file_count = 0

    for item in txt_files:
        if item.endswith('txt'):
            basename = item[:-4]
            k = basename.split('_')[1]
            j = basename.split('_')[3]
            if specific_ids and (k not in specific_ids or j not in specific_ids):
                continue
            file_count += 1
            print(j, k)
            with open(f'{path_prefix}/variant_maps/{item}', 'r') as read_file:
                query_dict = load_majority_seqs_from_fasta_file(f'{path_prefix}/{j}/majority_seqs.fasta')
                ref_dict = load_majority_seqs_from_fasta_file(f'{path_prefix}/{k}/majority_seqs.fasta')
                bio_validation_dict = bio_validation_mutations(query_dict, f'{path_prefix}/{k}/specie.fsa')
                for key in query_dict:
                    mutations = []
                    if key not in ref_dict:
                        print(f"Reference sequence not found for {key}")
                        sys.exit()
                    query = query_dict[key]
                    ref = ref_dict[key]
                    alignment_query, alignment_ref = align_and_identify_mutations(query, ref)
                    mutation_vector = create_mutation_vector(alignment_ref, alignment_query)
                    mutations = identify_mutations(mutation_vector, ref)

                    if mutations:
                        total_mutations += len(mutations)
                        total_pairs += 1
                        proxi_mutations = find_close_mutations(mutations)
                        total_proximity_mutations += len(proxi_mutations)
                        cooccurring_mutations = sum(1 for mut in mutations if len(mut.split(',')) > 1)
                        total_cooccurring_mutations += cooccurring_mutations
                        if proxi_mutations:
                            proximity_density, _ = calculate_proximity_density(proxi_mutations)
                            total_proximity_density += proximity_density
                        for mutation in mutations:
                            if not check_single_mutation_existence(bio_validation_dict, key, mutation):
                                total_novel_mutations += 1

    total_percentage_proximity = (total_proximity_mutations / total_mutations) * 100 if total_mutations > 0 else 0
    total_percentage_novel = (total_novel_mutations / total_mutations) * 100 if total_mutations > 0 else 0
    total_percentage_cooccurring = (total_cooccurring_mutations / total_mutations) * 100 if total_mutations > 0 else 0
    average_proximity_density = (total_proximity_density / total_proximity_mutations) if total_proximity_mutations > 0 else 0
    average_mutations_per_pair = total_mutations / total_pairs if total_pairs > 0 else 0

    return {
        'total_mutations': total_mutations,
        'total_novel_mutations': total_novel_mutations,
        'total_proximity_mutations': total_proximity_mutations,
        'total_cooccurring_mutations': total_cooccurring_mutations,
        'total_proximity_density': total_proximity_density,
        'total_percentage_proximity': total_percentage_proximity,
        'total_percentage_novel': total_percentage_novel,
        'total_percentage_cooccurring': total_percentage_cooccurring,
        'average_proximity_density': average_proximity_density,
        'average_mutations_per_pair': average_mutations_per_pair,
        'file_count': file_count
    }

def print_stats(stats, label):
    print(f"Stats for {label}:")
    print(f"Total mutations: {stats['total_mutations']}")
    print(f"Total novel mutations: {stats['total_novel_mutations']} ({stats['total_percentage_novel']:.2f}%)")
    print(f"Total mutations in proximity: {stats['total_proximity_mutations']} ({stats['total_percentage_proximity']:.2f}%)")
    print(f"Total co-occurring mutations: {stats['total_cooccurring_mutations']} ({stats['total_percentage_cooccurring']:.2f}%)")
    print(f"Total proximity density: {stats['total_proximity_density']}")
    print(f"Average proximity density: {stats['average_proximity_density']:.2f}")
    print(f"Average mutations per pair: {stats['average_mutations_per_pair']:.2f}")
    print(f"Number of .txt files: {stats['file_count']}")
    print()

txt_files = os.listdir('/home/people/malhal/test/mixed_isolates_nanomgt/variant_maps')

# Define the data sets
clean_data_set = {
    "Campylobacter jejuni": [
        "SRR27638397",
        "SRR27710526"
    ],
    "Escherichia coli": [
        "SRR28370668",
        "SRR25689478",
        "SRR26036455",
        "SRR28789754",
        "SRR26036458"
    ],
    "Klebsiella pneumoniae": [
        "ERR8958810",
        "ERR8958737",
        "SRR27348733",
        "SRR29213739",
        "SRR24833081"
    ],
    "Staphylococcus aureus": [
        "SRR28370694",
        "SRR28370638",
        "SRR25865495"
    ],
    "Salmonella enterica": [
        "SRR27755684",
        "SRR27755678",
        "SRR26899146"
    ],
    "Listeria monocytogenes": [
        "SRR27755667",
        "SRR27755674",
        "SRR26899103",
        "SRR25999202"
    ],
    "Pseudomonas aeruginosa": [
        "SRR24833086",
        "ERR8958866"
    ]
}

contaminated_data_set = {
    "Campylobacter jejuni": [
        "SRR26899118",
        "SRR26353490",
        "SRR27710532"
    ],
    "Escherichia coli": [
        "SRR24837710",
        "SRR24837712",
        "SRR24834173",
        "SRR28789463",
        "SRR28789469",
        "SRR28800580"
    ],
    "Klebsiella pneumoniae": [],
    "Staphylococcus aureus": [
        "ERR8958848",
        "SRR25890190",
        "ERR8958843"
    ],
    "Salmonella enterica": [
        "SRR27136090",
        "SRR28399428",
        "SRR27136088"
    ],
    "Listeria monocytogenes": [],
    "Pseudomonas aeruginosa": []
}

# Flatten the ids
clean_ids = [item for sublist in clean_data_set.values() for item in sublist]
contaminated_ids = [item for sublist in contaminated_data_set.values() for item in sublist]

# Calculate stats for each set
original_stats = calculate_stats(txt_files, '/home/people/malhal/test/mixed_isolates_nanomgt')
clean_stats = calculate_stats(txt_files, '/home/people/malhal/test/mixed_isolates_nanomgt', clean_ids)
contaminated_stats = calculate_stats(txt_files, '/home/people/malhal/test/mixed_isolates_nanomgt', contaminated_ids)

# Print the results
print_stats(original_stats, 'original')
print_stats(clean_stats, 'clean')
print_stats(contaminated_stats, 'contaminated')
