import os
import subprocess
import math

def ensure_dir(directory):
    """Ensure directory exists."""
    if not os.path.exists(directory):
        os.makedirs(directory)

def find_folders_and_files(data_dir):
    target_folders = []
    for root, dirs, files in os.walk(data_dir):
        if 'merged' in root:
            for file in files:
                if file.endswith('.fastq'):
                    target_folders.append((root, file))
    return target_folders

def map_reads(folders, reference_dir):
    maf_values = [0.01, 0.02, 0.03, 0.04, 0.05]

    for folder, fastq_file in folders:
        name = fastq_file.split('.fastq')[0]
        depth = int(name.split('_')[0][-3:])  # Assuming depth is part of the file name like 'depth220'
        # Extract the sequencing ID
        parts = folder.split('_')
        seq_id = [part for part in parts if part.startswith('SRR') or part.startswith('ERR')][0]
        output_dir = os.path.join(os.getcwd(), f'{name}_output')
        ensure_dir(output_dir)
        print(f"Output directory created at: {output_dir}")

        reference_path = os.path.join(reference_dir, f'{seq_id}.fasta')
        fastq_path = os.path.join(folder, fastq_file)
        bam_output = os.path.join(output_dir, f'{name}_output.bam')
        minimap_cmd = f"minimap2 -ax map-ont {reference_path} {fastq_path} | samtools sort -o {bam_output}"
        subprocess.run(minimap_cmd, shell=True, check=True)
        print(f"Minimap2 and Samtools processing completed. Output: {bam_output}")

        # Indexing the BAM file
        subprocess.run(['samtools', 'index', bam_output], check=True)

        # Generate consensus HDF file
        consensus_hdf = os.path.join(output_dir, f"{name}_consensus.hdf")
        medaka_consensus_cmd = f"medaka consensus {bam_output} {consensus_hdf}"
        print(f"Running Medaka consensus...")
        subprocess.run(medaka_consensus_cmd, shell=True, check=True)
        print(f"Medaka consensus completed. Output: {consensus_hdf}")

        vcf_output = os.path.join(output_dir, f'variants.vcf')
        medaka_variant_cmd = f"medaka variant {reference_path} {consensus_hdf} {vcf_output}"
        print(f"Executing Medaka variant for HDF: {consensus_hdf}...")
        subprocess.run(medaka_variant_cmd, shell=True, check=True)
        print(f"Medaka variant completed. VCF Output: {vcf_output}")

data_directory = '/home/people/malhal/data/new_nanomgt/confindr_data_dir'
reference_directory = '/home/people/malhal/test/new_nanomgt_results/references'

folders_files = find_folders_and_files(data_directory)
if not folders_files:
    print("No FASTQ files found matching the criteria.")
else:
    map_reads(folders_files, reference_directory)
