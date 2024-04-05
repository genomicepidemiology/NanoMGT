from Bio import pairwise2
from Bio.Seq import Seq


# Function to perform alignment and identify point mutations with custom scoring
def align_and_identify_mutations(seq1, seq2):
    match_score = 2  # Positive score for matches, encourages alignment to match characters
    mismatch_score = -3  # Negative score for mismatches, penalizes alignment for mismatching characters
    gap_open = -5  # Negative score for opening a gap, makes the algorithm less willing to introduce gaps
    gap_extend = -2  # Negative score for extending a gap, penalizes extending gaps over single mismatches
    # Perform local alignment with custom scoring
    alignments = pairwise2.align.localms(seq1, seq2, match_score, mismatch_score, gap_open, gap_extend)

    # Assuming we take the first alignment (typically the highest scoring)
    alignment_a, alignment_b, score, begin, end = alignments[0]
    print(alignment_a)
    print(alignment_b)

    # Identify point mutations, excluding gaps
    point_mutations = []
    for i in range(len(alignment_a)):
        if alignment_a[i] != alignment_b[i]:
            point_mutations.append((i + 1, alignment_a[i], alignment_b[i]))  # i+1 for 1-based indexing

    return point_mutations


# Example sequences
query = "AGTAGCTAGCTAGCTTTAGCTAGCAAG"
ref = "AGTAGCTAGCTACTTTAGCTAGCAAG"


# Align the sequences and identify point mutations with custom scoring
mutations = align_and_identify_mutations(query, ref)

# Print the results
print("Point mutations (excluding gaps):")
for position, nucleotide_a, nucleotide_b in mutations:
    print(f"Position: {position}, {nucleotide_b} -> {nucleotide_a}")  #from ref to query
