# PlasEval experiments

This branch contains data and scripts required to reproduce the experiments with 54 *E. coli* bacterial samples. 

The purpose of the experiments is to illustrate the use of PlasEval in assessing the performance of 4 plasmid binning tools:
<a href="https://github.com/cchauve/PlasBin-flow">PlasBin-flow</a>,
<a href="https://github.com/cchauve/HyAsP">HyAsP</a>,
<a href="https://github.com/phac-nml/mob-suite">MOB-recon</a>,
<a href="https://gitlab.com/mmb-umcu/gplas2">gplas2</a>.

## Data

Short read assemblies for these samples were used as input for the plasmid binning tools. 
Contigs from hybrid (short reads and long reads) assemblies were used as the ground truth plasmids.
Short read contigs were mapped onto the hybrid "plasmid" contigs using BLAST to obtain the ground truth plasmid bins. 
The output for the BLAST mapping between the short read and hybrid assemblies for each sample has also been provided. 

Short reads and hybrid assemblies, as well as the raw results of the plasmid binning tools PlasBin-flow, HyAsP, MOB-recon and gplas2 applied to the short reads assemblies, are available at: https://zenodo.org/records/10785151

## Experiments

## Ground truth

To generate ground truth plasmid bins from the hybrid assemblies, contigs that were longer than 500,000bp were labelled as chromosomes,
irrespective of circularity and circular contigs shorter than 500,000bp were labelled as plasmids. In the chosen samples, all hybrid contigs were either labelled as chromosome or plasmids.

The the short reads assembly contigs were mapped, using BLAST, against the hybrid assemblies. The BLAST mapping results are provided in this repo.
BLAST hits with identity below 95% or covering less than 80% of the short reads contig were discarded, and he remaining hits were used to define ground truth plasmid bins.
To generate these ground truth plasmid bins, the following command was used:

```
cd scripts
python ground_truth_mapping.py --references REFERENCE_FASTA --assembly ASSEMBLY_FASTA --mapping BLAST_OUTPUT_TSV --outdir OUT_DIR --outfile OUT_FILENAME --min_pid PID_THR --min_cov COV_THR --max_len LEN_THR
```
Here, `REFERENCE_FASTA` and `ASSEMBLY_FASTA` are the FASTA files for the hybrid and short read assemblies respectively. `BLAST_OUTPUT_TSV` is the location of the BLAST mapping output file in (outfmt 6) TSV format. `PID_THR` is the threshold used for percent identity of a particular match between a short read contig and a hybrid contig (reference sequence). `COV_THR` is the threshold used for proportion of short read contig cumulatively covered by reference sequences and `LEN_THR` is the maximum length of a hybrid contig for it to be considered a plasmid. 

## Reformatting plasmid bins

The following commands format the plasmid bins results of the four considered tools in the PlasEval input format.

```
cd scripts
python format_binning_results --tool TOOL --assembly PATH_TO_ASSEMBLY_GFA --results PATH_TO_BINNING_RESULTS_FILE --outdir PATH_TO_OUTPUT_FOLDER --outfile OUTPUT_FILENAME
```
Here, `TOOL` is the code used for the methods (`pbf` for PlasBin-flow, `hy` for HyAsP, `mob` for MOB-recon and `gp` for gplas2). `PATH_TO_ASSEMBLY_GFA` is the location of the short read assembly file in GFA format. `PATH_TO_BINNING_RESULTS_FILE` is the location of the binning results file for the chosen method. 

## Running PlasEval

The instructions to run PlasEval on the reformatted plasmid bins are provided in the README of the master branch of this repo.
