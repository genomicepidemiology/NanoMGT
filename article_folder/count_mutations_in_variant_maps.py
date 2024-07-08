import os
import sys

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

def check_single_mutation_existence(bio_validation_dict, gene, specific_mutation):
    return specific_mutation in bio_validation_dict.get(gene, [])

clean_ids = [
    "SRR27638397", "SRR27710526", "SRR28370668", "SRR25689478", "SRR26036455", "SRR28789754", "SRR26036458",
    "ERR8958810", "ERR8958737", "SRR27348733", "SRR29213739", "SRR24833081", "SRR28370694", "SRR28370638",
    "SRR25865495", "SRR27755684", "SRR27755678", "SRR26899146", "SRR27755667", "SRR27755674", "SRR26899103",
    "SRR25999202", "SRR24833086", "ERR8958866"
]

contaminated_ids = [
    "SRR26899118", "SRR26353490", "SRR27710532", "SRR24837710", "SRR24837712", "SRR24834173", "SRR28789463",
    "SRR28789469", "SRR28800580", "ERR8958848", "SRR25890190", "ERR8958843", "SRR27136090", "SRR28399428", "SRR27136088"
]

sets = [clean_ids, contaminated_ids]
path = '/home/people/malhal/test/mixed_isolates_nanomgt/variant_maps'

def count_mutations(file_path, bio_validation_dict):
    total_mutations = 0
    proximity_mutations = 0
    co_occuring_mutations = 0
    novel_mutations = 0

    with open(file_path, 'r') as file:
        print(file_path)
        lines = file.readlines()

        for line in lines:
            if not line.startswith('BACT'):
                mutations = line.strip().split(',')
                gene_mutations = []
                gene_mutations.append([int(mut.split('_')[0]) for mut in mutations])
                total_mutations += len(mutations)
                if len(mutations) > 1:
                    co_occuring_mutations += 1

                for mut in mutations:
                    if not check_single_mutation_existence(bio_validation_dict, gene, mut):
                        novel_mutations += 1

                for positions in gene_mutations:
                    positions.sort()
                    for i in range(len(positions) - 1):
                        if positions[i+1] - positions[i] <= 5:
                            proximity_mutations += 1
            else:
                gene = line.strip()

    return total_mutations, proximity_mutations, co_occuring_mutations, novel_mutations

total_counts = {"clean": {"files": 0, "mutations": 0, "proximity": 0, "co_occuring": 0, "novel": 0},
                "contaminated": {"files": 0, "mutations": 0, "proximity": 0, "co_occuring": 0, "novel": 0}}

for data_set, label in zip(sets, ["clean", "contaminated"]):
    for j in data_set:
        for k in data_set:
            file_path = os.path.join(path, f'major_{k}_minor_{j}.txt')
            if os.path.exists(file_path):
                path_prefix = '/home/people/malhal/test/mixed_isolates_nanomgt'
                query_dict = load_majority_seqs_from_fasta_file(f'{path_prefix}/{j}/majority_seqs.fasta')
                bio_validation_dict = bio_validation_mutations(query_dict, f'{path_prefix}/{k}/specie.fsa')
                total_counts[label]["files"] += 1
                mutations, proximity, co_occuring, novel = count_mutations(file_path, bio_validation_dict)
                total_counts[label]["mutations"] += mutations
                total_counts[label]["proximity"] += proximity
                total_counts[label]["co_occuring"] += co_occuring
                total_counts[label]["novel"] += novel

# Print the results
for label in ["clean", "contaminated"]:
    print(f"{label.capitalize()} Data Set:")
    print(f"  Number of .txt files: {total_counts[label]['files']}")
    print(f"  Total mutations: {total_counts[label]['mutations']}")
    print(f"  Proximity mutations: {total_counts[label]['proximity']}")
    print(f"  Co-occurring mutations: {total_counts[label]['co_occuring']}")
    print(f"  Novel mutations: {total_counts[label]['novel']}")
    print()

# Output for the entire data set (clean + contaminated)
print("Combined Data Set:")
total_files = total_counts["clean"]["files"] + total_counts["contaminated"]["files"]
total_mutations = total_counts["clean"]["mutations"] + total_counts["contaminated"]["mutations"]
total_proximity = total_counts["clean"]["proximity"] + total_counts["contaminated"]["proximity"]
total_co_occuring = total_counts["clean"]["co_occuring"] + total_counts["contaminated"]["co_occuring"]
total_novel = total_counts["clean"]["novel"] + total_counts["contaminated"]["novel"]
print(f"  Number of .txt files: {total_files}")
print(f"  Total mutations: {total_mutations}")
print(f"  Proximity mutations: {total_proximity}")
print(f"  Co-occurring mutations: {total_co_occuring}")
print(f"  Novel mutations: {total_novel}")
