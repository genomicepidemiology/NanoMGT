import re
import sys
import json
from Bio import pairwise2
from Bio.Seq import Seq

alignment_result_path = "/home/people/malhal/test/sup_nanomgt_data"

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
    if specific_mutation in bio_validation_dict.get(gene, []):
        return True
    return False


def derive_correct_length_headers(ref_dict, fsa_file):
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
        raise ValueError("The mutation vector and reference sequence must have the same length.")

    for i in range(len(mutation_vector)):
        if mutation_vector[i] != reference_sequence[i] and mutation_vector[i] != "-":
            mutations.append('{}_{}'.format(i+1, mutation_vector[i]))

    return mutations


def align_and_identify_mutations(seq1, seq2):
    match_score = 2
    mismatch_score = -3
    gap_open = -5
    gap_extend = -2
    alignments = pairwise2.align.localms(seq1, seq2, match_score, mismatch_score, gap_open, gap_extend)
    alignment_a, alignment_b, score, begin, end = alignments[0]
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


# Load species sequencing IDs from file
with open('species_sequencing_ids.json', 'r') as f:
    species_sequencing_ids = json.load(f)

for species, sequencing_ids in species_sequencing_ids.items():
    for i in range(len(sequencing_ids)):
        for j in range(len(sequencing_ids)):
            if i != j:
                major_id = sequencing_ids[i]
                minor_id = sequencing_ids[j]
                print(f'{major_id}, {minor_id}')
                print(f'major_{minor_id}_minor_{major_id}.txt')

                mutations_in_proximity = 0
                novel_mutation_count = 0
                total_mutations = 0
                proximity_mutations = []
                total_proximity_density = 0

                with open(f'major_{minor_id}_minor_{major_id}.txt', 'w') as write_file:
                    query_dict = load_majority_seqs_from_fasta_file(f'{alignment_result_path}/{major_id}/majority_seqs.fasta')
                    ref_dict = load_majority_seqs_from_fasta_file(f'{alignment_result_path}/{minor_id}/majority_seqs.fasta')
                    bio_validation_dict = bio_validation_mutations(query_dict, f'{alignment_result_path}/{minor_id}/specie.fsa')

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
                            print(key, file=write_file)
                            print(",".join(mutations), file=write_file)
                            proxi_mutations = find_close_mutations(mutations)
                            mutations_in_proximity += len(proxi_mutations)
                            if proxi_mutations:
                                proximity_density, num_mutations = calculate_proximity_density(proxi_mutations)
                                total_proximity_density += proximity_density
                            for mutation in mutations:
                                biological_existence = check_single_mutation_existence(bio_validation_dict, key, mutation)
                                if not biological_existence:
                                    novel_mutation_count += 1
                                    print(key, mutation)

                    print(f"Total mutations: {total_mutations}")
                    print(f"Novel mutations: {novel_mutation_count}")
                    print(f"Mutations in proximity: {mutations_in_proximity}")
                    print(f"Total density: {total_proximity_density}")
                    if mutations_in_proximity > 0:
                        print(f"Average proximity density: {total_proximity_density / mutations_in_proximity}")

