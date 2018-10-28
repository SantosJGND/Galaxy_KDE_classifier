
## Summary Analysis.

The application of the following tools follows a first layer of analysis of local genomic windows.

Reference KDEs were extracted across data sets of contiguous bi-allelic markers. Sample-level statistics
were extracted at each window and stored in matrices *Blocks_Request* and *Blocks_profiles* (optional).

The generation of ideograms allows an initial, visual exploration of this output (1). This repository 
provides tools for a quantitative exploration of this data.

- *1* We performed an analysis of Core Asian rice variation (3K RG, Alexandrov *et al*. 2014). A classification
of that data was performed. Core Ideograms were stored on another repository: [Core Ideograms](https://imgur.com/a/lpD0r31).

### gff_Blocksmerge.py - combining information 

The script gff_Blocksmerge.py produces this merged data set of local classification and genic location data. it takes a gff file, 
and the user can chose to focus on particular elements (miRNA, genes or both), as well as the columns of information on each 
to extract, by name.

**Context.** The genome crawl divides the genome into a set of overlapping windows. Our downstream analysis queries these windows for associations. 
For various reasons it can be usefull to know which genes are present at each window along the genome. 

*Note:* this information can be easily included in the an product dash application.


### Summary_stats.py - Physical classification summary.

Calculate the extent of physical regions assigned to each class by selected individual for a set of pre-determined
parameters.

### Intermediate_labels_Exp.py - Classification distribution

Calculate the density of target classifications across genomic regions for a group of selected individuals.

*app ready*: this script outputs a dash application that combines an analysis of the distribution of classifications
with the coordinates of genic regions in the merged gff-to-blocks file.

### pVal_gff_overlap.py - KDE overlap by genes

Calculate a measure of distribution overlap across genes in merged gff file. Distribution overlap is calculated
as the sum of proportional differences in assignment to reference distributions. Overlap is measured using 
reference samples only. 

### Complement - Galaxy_summary_tools.py

## Output examples.

Summary output, accompanied of jupyter notebooks or, in the case of `Intermediate_labels_Exp.py` analysis, dash application.

- **/Summary_classification**: Group and individual classification summaries.

- **/Intermediate**: Genome-wide density of selected classifications across group of interest.

- **/Overlap_gff**: Distribution overlap measure at regions of interest.


## References

- Alexandrov N. (2017) Rice SNP-seek database update: new SNPs, indels, and queries. Nucl. Acids Res. 45(D1):D1075-D1081.