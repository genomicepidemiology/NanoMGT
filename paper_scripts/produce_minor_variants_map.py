import re
import sys
from Bio import pairwise2
from Bio.Seq import Seq

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


ecoli_sequencing_ids = ['SRR28399430', 'SRR26899125'] # BACT000018 {225, 228}

# Staphylococcus aureus
staph_aureus_sequencing_ids = ['SRR28370694', 'ERR8958843', 'SRR27167517'] # BACT000030 {690, 693}

# Acinetobacter baumannii
acinetobacter_baumannii_sequencing_ids = ['SRR26643493', 'SRR27710531', 'SRR23387317'] #All share rqual lengths

# Campylobacter jejuni
campylobacter_jejuni_sequencing_ids = ['SRR27638397', 'SRR27710526'] # BACT000051 {426: SRR27710526, 444: SRR26899121, SRR27638397}

# Salmonella enterica
salmonella_enterica_sequencing_ids = ['SRR28399428', 'SRR27136088', 'SRR27755678'] #All share rqual lengths

for j in salmonella_enterica_sequencing_ids:
    for k in salmonella_enterica_sequencing_ids:
        if j != k:
            with open('major_{}_minor_{}.txt'.format(k, j), 'w') as write_file:
                print ('major_{}_minor_{}'.format(k, j))
                query_dict = load_majority_seq_from_vcf(f'/home/people/malhal/data/new_nanomgt/majority_variants/{j}/majority_variants.vcf')
                ref_dict = load_majority_seq_from_vcf(f'/home/people/malhal/data/new_nanomgt/majority_variants/{k}/majority_variants.vcf')

                for key in query_dict:
                    mutations = []
                    if key not in ref_dict:
                        print (f"Reference sequence not found for {key}")
                        sys.exit()
                    query = query_dict[key]
                    ref = ref_dict[key]
                    alignment_query, alignment_ref = align_and_identify_mutations(query, ref)


                    for i in range(len(alignment_query)):
                        # Check if there is a gap in the reference or a mismatch
                        # Check if there is a gap in the reference or a mismatch
                        if alignment_ref[i] != alignment_query[i]:
                            # For gaps in the reference, we do not increment index_ref
                            if alignment_query[i] != '-':
                                mutations.append(f"{i+1}_{alignment_query[i]}")
                    if mutations != []:
                        print (key, file = write_file)
                        print (",".join(mutations), file = write_file)