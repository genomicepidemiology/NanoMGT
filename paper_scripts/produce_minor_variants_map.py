import re
import sys
from Bio import pairwise2
from Bio.Seq import Seq

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


ecoli_sequencing_ids = ['SRR25689478', 'SRR26899125'] # BACT000018 {225, 228}

# Staphylococcus aureus
staph_aureus_sequencing_ids = ['SRR28370694', 'ERR8958843', 'SRR27167517'] # BACT000030 {690, 693}

# Campylobacter jejuni
campylobacter_jejuni_sequencing_ids = ['SRR27638397', 'SRR27710526'] # BACT000051 {426: SRR27710526, 444: SRR26899121, SRR27638397}

# Salmonella enterica
salmonella_enterica_sequencing_ids = ['SRR28399428', 'SRR27136088', 'SRR27755678'] #All share rqual lengths

combined_list_of_ids = [ecoli_sequencing_ids, staph_aureus_sequencing_ids, campylobacter_jejuni_sequencing_ids, salmonella_enterica_sequencing_ids]

print (combined_list_of_ids)


#for item in combined_list_of_ids:
for item in combined_list_of_ids:
    for j in item:
        for k in item:
            if j != k:
                print (j, k)
                with open('major_{}_minor_{}.txt'.format(k, j), 'w') as write_file:
                    print ('major_{}_minor_{}'.format(k, j))
                    query_dict = load_majority_seqs_from_fasta_file(f'/home/people/malhal/data/new_nanomgt/majority_variants/{j}/majority_seqs.fasta')
                    ref_dict = load_majority_seqs_from_fasta_file(f'/home/people/malhal/data/new_nanomgt/majority_variants/{k}/majority_seqs.fasta')
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
                            print (key, file = write_file)
                            print (",".join(mutations), file = write_file)
