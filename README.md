# NanoMGT

NanoMGT is a Python-based tool designed for metagenomic variant calling using Nanopore sequencing data. It provides efficient methods for aligning sequences and identifying genetic variations in metagenomic samples.

## Features

- Metagenomic variant calling tailored for metagenomic Nanopore sequencing data of the same species.
- Identifies minority position variant calls which originates from strains which are not the most abundant.

## Installation

To install NanoMGT, clone the repository and install the required dependencies:

```bash
conda install -c genomicepidemiology nanomgt
```

## Downloading Ribosomal MLST Database from PubMLST

Follow these steps to download the Ribosomal MLST genome and Ribosomal MLST locus/sequence databases from PubMLST:

1. **Create a PubMLST Account**
   - Sign up for an account at [PubMLST](https://pubmlst.org/).

2. **Request Database Access**
   - Once registered, request access to the Ribosomal MLST genome and the Ribosomal MLST locus/sequence databases.

3. **Download Ribosomal MLST Profiles and Allele Sequences**
   - Go to the Ribosomal MLST database at [PubMLST rMLST Seqdef](https://pubmlst.org/bigsdb?db=pubmlst_rmlst_seqdef).
   - In the section on the right, download the Ribosomal MLST profiles and Allele sequences.
   - To download Allele sequences:
     ```bash
     # Navigate to Downloads > Allele sequences > All loci folder.
     # Click on the 53 separate download buttons for each locus.
     ```

4. **Concatenate Loci Files and Save Profiles**
   - After downloading all 53 loci files, concatenate them into one file named `rmlst.fsa` and save the profiles in a file called `rmlst_scheme.txt`.
     ```bash
     cat BACT000*.fsa > rmlst.fsa
     # Ensure 'rmlst_scheme.txt' is saved in the same directory.
     ```

5. **Download Homology Reduced Bacteria Database**
   - Download the homology reduced bacteria database from the CGE server:
     ```bash
     wget https://cge.food.dtu.dk/services/nanomgt/nanomgt_db.tar.gz
     tar -xf nanomgt_db.tar.gz
     ```

6. **Organize Files**
   - Place `rmlst.fsa` and `rmlst_scheme.txt` in the same folder as the homology reduced bacteria database:
     ```bash
     mv rmlst.fsa rmlst_scheme.txt /path/to/homology_reduced_bacteria_db
     ```

Replace `/path/to/homology_reduced_bacteria_db` with the actual path to the homology reduced bacteria database folder.

In the final rmlst_db folder there should be the following files:

```bash
rmlst.fsa
rmlst_scheme.txt
bac_db.name
bac_db.seq.b
bac_db.comp.b
bac_db.len.b
```

## Usage

```bash
nanomgt --nanopore /complete/path/to/input.fastq --o any_output_name --threads <int, default:4> --maf <float, default:0.05> --db_dir /path/to/database_as_described_above/
```

## License

Apache License 2.0

## Authors
Malte B. Hallgren
