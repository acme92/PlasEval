# PlasEval

## Evaluation script
### Input
1. A file with predicted contig chains with or without orientations. The file should contain one chain per line.<br/>
Format with orientation:<br/>
plasmid_1;27+,28-,12+<br/>
plasmid_2;13-,1+,4+<br/>
Format without orientation:<br/>
plasmid_1;27,28,12<br/>
plasmid_2;13,1,4<br/>

2. A file that maps a contig to a plasmid, chromosome or ambiguous sequence with one contig per line. <br/>
Format:<br/>
short_read_contig_id,short_read_contig_length,mapped_against_hybrid_contig_id,hybrid_contig_length,number_of_residue_matches,is_circular,label

3. Path to output directory and name of output file.

#### Usage
```
python evaluate_sample.py --pred prediction_file --map mapping_file --out out_dir --res out_file --amb 0/1 --ori 0/1
```

Additional arguments
```
--amb		Argument to indicate if contigs mapped ambiguous sequences should be considered plasmidic. If yes, then pass value 1, else 0.
--ori		Argument to indicate if predicted contig chains have oriented contigs. If yes, then pass value 1, else 0.                           
```

### Output
A text file that contains the following information about the evaluation of a particular prediction.<br/>
1. No. of references chromosomes and plasmids and their respective lengths,<br/>
2. No. of predicted plasmids and their respective lengths,<br/>
3. Matrix with proportion of reference plasmid covered by predicted plasmid,<br/>
4. Matrix with proportion of predicted plasmid covered by reference plasmid,<br/>
5. Overall covered proportion of predicted plasmids,<br/>
6. Overall covered proportion of reference plasmids,<br/>
7. Pairs of predicted and reference plasmids with covered proprtion over a specific threshold in both directions,<br/>
8. Precision, recall and F1 score<br/>

## Comparison script
### Input
1. Two .tsv files with predicted / reference plasmids as sets of contigs. The file should contain one chain per line.<br/>
Format:<br/>
Plasmid_ID	Contig_ID 	Contig_Length<br/>

To generate the input in the desired format from the HyAsP output, we use the 'putative_plasmid_contigs.fasta' file. Format:<br/>
'>'Contig_ID|x_plasmid_ID<br/>
Contig_sequence

#### Usage
```
python format_input_v2.py --i input_dir/putative_plasmid_contigs.fasta --o output_dir --f formatted_file.tsv
```

2. Path to output directory and name of output file.

#### Usage
```
python plasmid_comparison_main.py --l left_plasmids.tsv --r right_plasmids.tsv --out out_dir --res out_file
```

### Output
A log file that contains the following information for each possible labelling combination. The combinations are listed in increasing order of the dissimilarity score.:<br\>
1. Combination number,<br/>
2. Dissimilarity score,<br/>
3. Total length of contigs that are unique to one set of plasmids and half the total length of contigs that are common but grouped differently in both sets of plasmids,<br/>
4. Labelling of contigs in both sets of plasmids,<br/>
5. Splits in both sets of plasmids
