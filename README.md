# PlasEval experiments

## Data

This branch contains the input and output of the experiments with 54 *E. coli* bacterial samples. 

Short read assemblies for these samples were used as input for the plasmid binning tools. 
Contigs from hybrid (short reads and long reads) assemblies for these samples that were at most 500,000 bp long were used as the ground truth plasmids.
Short read contigs were mapped onto the hybrid "plasmid" contigs using BLAST to obtain the ground truth plasmid bins. 
The output for the BLAST mapping between the short read and hybrid assemblies for each sample has also been provided. 

Short reads and hybrid assemblies are available at: https://zenodo.org/records/10785151

## Scripts
Two scripts have been provided for formatting the results of the binning tools used in the experiments as well the results of the BLAST mapping to generate the ground truth bins.

#### Usage

1. The following commands format the binning results of the specified tool in the PlasEval input format.
```
cd scripts
python format_binning_results --tool TOOL --assembly PATH_TO_ASSEMBLY_GFA --results PATH_TO_BINNING_RESULTS_FILE --outdir PATH_TO_OUTPUT_FOLDER --outfile OUTPUT_FILENAME
```
Here, `TOOL` is the code used for the methods (`pbf` for PlasBin-flow, `hy` for HyAsP, `mob` for MOB-recon and `gp` for gplas2). `PATH_TO_ASSEMBLY_GFA` is the location of the short read assembly file in GFA format. `PATH_TO_BINNING_RESULTS_FILE` is the location of the binning results file for the chosen method. 

2. The following commands format the BLAST output to generate the ground truth bins in the PlasEval input format.
```
cd scripts
python ground_truth_mapping.py --references REFERENCE_FASTA --assembly ASSEMBLY_FASTA --mapping BLAST_OUTPUT_TSV --outdir OUT_DIR --outfile OUT_FILENAME --min_pid PID_THR --min_cov COV_THR --max_len LEN_THR
```
Here, `REFERENCE_FASTA` and `ASSEMBLY_FASTA` are the FASTA files for the hybrid and short read assemblies respectively. `BLAST_OUTPUT_TSV` is the location of the BLAST mapping output file in (outfmt 6) TSV format. `PID_THR` is the threshold used for percent identity of a particular match between a short read contig and a hybrid contig (reference sequence). `COV_THR` is the threshold used for proportion of short read contig cumulatively covered by reference sequences and `LEN_THR` is the maximum length of a hybrid contig for it to be considered a plasmid. 

