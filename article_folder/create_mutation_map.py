import pandas as pd


# Function to load VCF and return the desired dictionary
def load_vcf_to_dict(vcf_file):
    # Load the VCF file, skipping the metadata lines
    vcf_df = pd.read_csv(vcf_file, comment='#', sep='\t', header=0,
                         names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'])

    # Initialize an empty dictionary to store the data
    chrom_ref_dict = {}

    # Iterate over the dataframe to populate the dictionary
    for index, row in vcf_df.iterrows():
        # If the CHROM is not in the dict, initialize it with the REF value
        if row['CHROM'] not in chrom_ref_dict:
            chrom_ref_dict[row['CHROM']] = row['REF']
        else:
            # Append the REF value to the existing string for this CHROM
            chrom_ref_dict[row['CHROM']] += row['REF']

    return chrom_ref_dict


ecoli_sequencing_ids = ['SRR28399430', 'SRR26899125', 'SRR26899129'] # BACT000018 {225, 228}

# Staphylococcus aureus
staph_aureus_sequencing_ids = ['SRR28370694', 'ERR8958843', 'SRR27167517'] # BACT000030 {690, 693}

# Campylobacter jejuni
campylobacter_jejuni_sequencing_ids = ['SRR27638397', 'SRR26899121', 'SRR27710526'] # BACT000051 {426: SRR27710526, 444: SRR26899121, SRR27638397}

# Salmonella enterica
salmonella_enterica_sequencing_ids = ['SRR28399428', 'SRR27136088', 'SRR27755678'] #All share rqual lengths

list_of_sequencing_ids = [ecoli_sequencing_ids, staph_aureus_sequencing_ids, acinetobacter_baumannii_sequencing_ids,
                            campylobacter_jejuni_sequencing_ids, salmonella_enterica_sequencing_ids]

for sequencing_ids in list_of_sequencing_ids:
    print (sequencing_ids)
    specie_genes = {}

    for id in sequencing_ids:

        # Path to your VCF file
        vcf_file_path = id + '/majority_variants.vcf'

        # Load the VCF file and get the dictionary
        chrom_ref_sequence = load_vcf_to_dict(vcf_file_path)

        # Print the dictionary to verify
        for chrom, ref_seq in chrom_ref_sequence.items():
            gene = chrom.split('_')[0]
            if gene not in specie_genes:
                specie_genes[gene] = set()
                specie_genes[gene].add(len(ref_seq))
            else:
                specie_genes[gene].add(len(ref_seq))
    for gene, lengths in specie_genes.items():
        print(gene, lengths)
