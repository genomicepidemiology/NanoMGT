# Most important scripts used for the nanoMGT article experiments:

1. simulate_batches.py to simulate the batch files
2. subsample_nanopore.py to subsample the batch files with the fastq files download with ENA
3. produce_minor_variants_map.py is used to determine the expected true positive set between the isolates' concensus sequences.
4. run_nanomgt_to_id_noise_and_majority_seqs.py to identify noise with each input isolate.
4. run_nanomgt_on_sample to produce nanoMGT results for each multistrain sample.
5. run_cf.py to run confindr (requires some manual work of setting of the Confidnr database and put each input fastq into its own directory.)
6. run_longshot.py to run longshot (requires some manuel setup of the rMLST gene references with samtools and bwa)
7. run_medaka.py to run medaka (requires some manuel setup of the rMLST gene references with samtools and bwa)
8. The rest of the scripts are used to generate the figures, tables and calculate F1, precision and recall scores in the article.