import re
import sys
from Bio import pairwise2
from Bio.Seq import Seq

def align_sequences(seq1, seq2):
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


def extract_alignment(reference, query, cigar_string):
    """
    Extracts the alignment between reference and query based on the provided CIGAR string.

    Parameters:
    - reference (str): The reference sequence.
    - query (str): The query sequence.
    - cigar_string (str): The CIGAR string.

    Returns:
    - str, str: The portions of the reference and query sequences that align, respectively.
    """
    # Parse the CIGAR string
    cigar_tuples = re.findall(r'(\d+)([MIDNSHP=X])', cigar_string)

    aligned_ref = ""
    aligned_query = ""

    ref_index = 0
    query_index = 0

    # Loop through CIGAR tuples and build the aligned sequences
    for length, operation in cigar_tuples:
        length = int(length)

        if operation == "M":  # Match or mismatch
            aligned_ref += reference[ref_index:ref_index + length]
            aligned_query += query[query_index:query_index + length]
            ref_index += length
            query_index += length
        elif operation == "I":  # Insertion
            aligned_ref += "-" * length
            aligned_query += query[query_index:query_index + length]
            query_index += length
        elif operation == "D":  # Deletion
            aligned_ref += reference[ref_index:ref_index + length]
            aligned_query += "-" * length
            ref_index += length
        elif operation == "S":  # Soft clipping
            query_index += length
        elif operation == "=":  # Match
            aligned_ref += reference[ref_index:ref_index + length]
            aligned_query += query[query_index:query_index + length]
            ref_index += length
            query_index += length
        elif operation == "X":  # Mismatch
            aligned_ref += reference[ref_index:ref_index + length]
            aligned_query += query[query_index:query_index + length]
            ref_index += length
            query_index += length
        # Optionally handle other CIGAR operations (N, H, P) as per requirements

    return aligned_ref, aligned_query

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


def identify_mutations(mutation_vector, reference_sequence, gene_mutations, read_id, read_positions_blacklisted_dict):
    """
    Identify all mutation positions from a mutation vector compared to the reference.

    Parameters:
    - mutation_vector (list[str]): The mutation vector.
    - reference_sequence (str): The original reference sequence.

    Returns:
    - list[str]: A list where each mutation is described as a string in the format "POSITION_NUCLEOTIDE".
    """
    mutations = []
    if len(mutation_vector) == len(reference_sequence):
        #TBD implement check for high amount of mutations is super bad coding, write better as a function when you have time.
        alignment_query = ''.join(mutation_vector)
        alignment_ref = reference_sequence
        index = 0

        for i in range(len(alignment_query)):
            # Check if there is a gap in the reference or a mismatch
            # Check if there is a gap in the reference or a mismatch
            if alignment_ref[i] != alignment_query[i]:
                # For gaps in the reference, we do not increment index_ref
                if alignment_query[i] != '-':
                    mutations.append(f"{index + 1}_{alignment_query[i]}")
            if alignment_query[i] != '-':
                index += 1

        if len(mutations)/len(reference_sequence) > 0.05: #5% mutations, likely a need for alignment
            mutations = []
            alignment_query, alignment_ref = align_sequences(''.join(mutation_vector), reference_sequence)
            index = 0

            for i in range(len(alignment_query)):
                # Check if there is a gap in the reference or a mismatch
                # Check if there is a gap in the reference or a mismatch
                if alignment_ref[i] != alignment_query[i]:
                    # For gaps in the reference, we do not increment index_ref
                    if alignment_query[i] != '-':
                        mutations.append(f"{index + 1}_{alignment_query[i]}")
                if alignment_query[i] != '-':
                    index += 1
    else:
        alignment_query, alignment_ref = align_sequences(''.join(mutation_vector), reference_sequence)

        index = 0

        for i in range(len(alignment_query)):
            # Check if there is a gap in the reference or a mismatch
            # Check if there is a gap in the reference or a mismatch
            if alignment_ref[i] != alignment_query[i]:
                # For gaps in the reference, we do not increment index_ref
                if alignment_query[i] != '-':
                    mutations.append(f"{index + 1}_{alignment_query[i]}")
            if alignment_query[i] != '-':
                index += 1
    return mutations

def parse_sam_and_find_mutations(sam_file_path, confirmed_mutation_dict, consensus_dict, read_positions_blacklisted_dict):
    """
    Parses a SAM file, extracts necessary information and finds mutations in each read.

    Parameters:
    - sam_file_path (str): The path to the SAM file.

    Returns:
    - dict: A dictionary where keys are read names and values are lists of mutation strings.

    """

    mutations_dict = {}
    with open(sam_file_path, 'r') as sam_file:
        for line in sam_file:
            # Skip header lines
            if line.startswith('@'):
                continue
            # Extract relevant columns: [QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ]
            cols = line.strip().split('\t')
            qname, flag, rname, pos, mapq, cigar_str, rnext, pnext, tlen, seq = cols[:10]
            read_id = qname.split(' ')[0]

            # Convert string columns to appropriate types
            pos = int(pos)
            tlen = int(tlen)

            #Should be start pos of the alignment and not of the read
            #TBD could we optimize this is still work with partial alignments? I guess we dont need to have full coverage.
            if pos == 1 and len(seq) >= tlen:
                majority_seq = consensus_dict[rname][1]
                #TBD consider if this majority_seq is a problem if it contains gaps
                # Obtaining the alignment using your function
                aligned_ref, aligned_query = extract_alignment(majority_seq[pos-1:pos-1+tlen], seq, cigar_str)
                print (aligned_ref)
                print (aligned_query)
                sys.exit()
                mutation_vector = create_mutation_vector(aligned_ref, aligned_query)
                if 'BACT000001_' in rname:
                    print (mutation_vector)
                mutations = identify_mutations(mutation_vector, majority_seq[pos-1:pos-1+tlen], confirmed_mutation_dict[rname][0], read_id, read_positions_blacklisted_dict)
                name = read_id + ' ' + rname
                mutations_dict[name] = mutations

    return mutations_dict

def parse_fsa_get_references(fsa_file_path):
    """
    Parse the FSA file and extract all unique reference names.

    :param fsa_file_path: Path to the FSA file.
    :return: Set of reference names.
    """
    references = dict()
    with open(fsa_file_path, 'r') as fsa_file:
        for line in fsa_file:
            # Skip header lines
            if line.startswith(">"):
                gene = line.strip()[1:].split('_')[0]
                references[gene] = ""
            else:
                references[gene] += line.strip()
    return references


def parse_sam_get_references(sam_file_path):
    """
    Parse the SAM file and extract all unique reference names.

    :param sam_file_path: Path to the SAM file.
    :return: Set of reference names.
    """
    references = set()
    with open(sam_file_path, 'r') as sam_file:
        for line in sam_file:
            # Skip header lines
            if line.startswith("@"):
                continue

            fields = line.split("\t")
            if len(fields) > 2:
                references.add(fields[2])  # The reference name is in the 3rd field
    return references


def load_references_from_fasta(fasta_file_path, references):
    """
    Load reference sequences from a FASTA file that are in the provided reference names set.

    :param fasta_file_path: Path to the FASTA file.
    :param references: Set of reference names to load.
    :return: Dictionary with reference names as keys and sequences as values.
    """
    ref_sequences = {}
    with open(fasta_file_path, 'r') as fasta_file:
        current_ref_name = None
        current_sequence = []
        for line in fasta_file:
            line = line.strip()
            # Check if the line is a header line
            if line.startswith(">"):
                # If we were building a sequence, add it to the dictionary
                if current_ref_name and (current_ref_name in references):
                    ref_sequences[current_ref_name] = "".join(current_sequence)

                # Get the new reference name
                current_ref_name = re.sub(r"^>\s*", "", line).split()[0]  # Remove ">" and extract the first word

                # Reset the sequence
                current_sequence = []
            else:
                # Add the line to the sequence
                current_sequence.append(line.upper())

        # Add the last sequence
        if current_ref_name and (current_ref_name in references):
            ref_sequences[current_ref_name] = "".join(current_sequence)
    return ref_sequences