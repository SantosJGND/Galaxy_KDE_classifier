## Kernel classification.

Local haplotype classification and data exploration. This software was built with the purpose of comparing local 
genetic distributions of diverging biological populations. The comparison of distributions is accompanied by a set of 
visualisation tools we deemed useful for inference of evolutionary relationships.

_Note_: Application of this pipeline of analysis should take into considerations the particular qualities of the data set for which 
it was developped. The requirements for an informative output are only nearly met given the stable and high density of genetic markers,
 a state of near complete homozygosity, to add to the quality of the data.

### Kernel Density Estimation -- script: Kernel_mPLib3D_FM36_Galaxy.py 

The initial analysis is conducted following a machine learning approach 
of dimensionality reduction and kernel density estimation. 

This analysis is conducted on the basis of a crawling window of fixed size of genetic markers. Window size (-w), 
dimensionality reduction (-Dr), the number of axis to retain in feature space (-ncomp) and the level of heterozygosity allowed 
for (-het) are customizable. If PCA is chosen for dimensionality reduction, the user can chose to print the variance 
of each axis along the crawl (--VARprint, boolean).

Parallel to the estimation of local supervised distributions, an unsupervised clustering using the MeanShift algorithm 
(Comaniciu & Meer, 2002) is performed. KD estimates of each cluster identified are stored. These are necessary for downstream analysis 
but are not automatically printed.

While the software could in theory handle more than three reference populations, we advise this number to be limited to three. 
This is because, for plotting as well as targeted cluster analysis, we allowed for intermediate classes among the references 
provided. At K=3, this results in 8 classes in total (3 pure classes, 3 two-way intermediates and 1 three way intermediate). 
At K == 4, this would imply 6 two way intermediate classes. Unless reference populations are strongly structured, the resulting plots 
can become messy.

**Input**

- *[.geno](https://github.com/SantosJGND/Galaxy_KDE_classifier/blob/master/Simulation_related/test_.05_.2/test.geno) file*: genotype file, bi-allelic markers. Columns = Inds; Lines = SNPs; Codes: {Reference hom: 0; Het: 1; Alternative Hom: 2}.
- *[.bim](https://github.com/SantosJGND/Galaxy_KDE_classifier/blob/master/Simulation_related/test_.05_.2/test.bim) file*: SNP information file; plink format.
- *[.fam](https://github.com/SantosJGND/Galaxy_KDE_classifier/blob/master/Simulation_related/test_.05_.2/test.fam) file*: Individual information file; plink format.
- *[ref](https://github.com/SantosJGND/Galaxy_KDE_classifier/blob/master/Simulation_related/refs_sim.txt) file*: Individual to group file, tab delimited. lines: ID (present in .fam file); group code.
- *[admix](https://github.com/SantosJGND/Galaxy_KDE_classifier/blob/master/Simulation_related/admx_sim.txt) file*: Individual to group file, tab delimited. lines: ID (present in .fam file); group code.

**Output**

- *Blocks_Request*: reference-specific KDE estimated *p*-values across windows surveyed. 
- (optional) *Blocks_profiles*: unsupervised KDE estimated *p*-values across windows surveyed.

**Example files**

- Simulations: see [/Simulation_related/test_.05_.2](Galaxy_KDE_classifier/Simulation_related/test_.05_.2/Ideo_sample32_CHR01_Z1.6_bin9.png)
- Command line: see [Command_examples.txt](Command_examples.txt)

### Ideogram plotting -- script: Ideogram_plots_Galaxy.py

The first layer of analysis following the genome crawl is the visualization of local statistics pertaining 
to the distributions of the reference populations chosen. This step constitutes an exploration of variation 
in our data set. The purpose of this exploration can be to visualize possible exchanges of genetic material among 
reference material, the local composition of non-reference accessions included in the analysis, or to investigate 
the presence of differentiated (outlier) material.

For this purpose, the construction of _chromosome paintings_ allows for user defined thresholds of comparison and an 
outlier threshold to control the assignment into pure, intermediate and outlier classes.

A savitsky golay filter is automatically applied if the flag '--coarse' is not passed. We recommend the use of this flag.

**Input**

- [ID](https://github.com/SantosJGND/Galaxy_KDE_classifier/blob/master/Simulation_related/Focus_IDs.txt) file
- Blocks_request file

**Output**

- Ideogram plots (.png)

Example figure: [Ideogram cut](Fig_4.png), see [Summary Analysis](https://github.com/SantosJGND/Galaxy_KDE_classifier/tree/master/Summary_functions) section


## Galaxy integration

The two scripts described above have been integrated to the Galaxy platform ([Giardine *et al*. 2005](https://genome.cshlp.org/content/15/10/1451.short)). 
The portal provides an accessible interface to deploy various bioinformatics tools while ensuring reproducibility. Under the header *KDE_classifier*, an instance 
was developed for the deployment of a reference based classification of consecutive genomic window using KDE and its subsequent analysis in the form of ideograms.

The scripts for classification and plotting are available as separate tools:

- [KDE](http://galaxy.southgreen.fr/galaxy/root?tool_id=KDE1) 
- [Ideograms](http://galaxy.southgreen.fr/galaxy/root?tool_id=Ideogram)

In each case the parameters of the script are accessible and defaults are provided. 
Through Galaxy, these tools can be combined into a [workflow of analysis](http://galaxy.southgreen.fr/galaxy/u/acomte/p/reconstruction-of-mosaic-genomes) to facilitate reproduction.

## Summary Analysis - /Summary_functions

Ideograms provide a visual summary of reference KDE assignments across local data sets. This repository provides tools for 
a quantitative analysis of local KDE estimates.

## Downstream Analysis - /Downstream_functions

The study of reference distributions at local genomic regions consists of the first layer in a process of inference. This repository 
provides tools to complement supervised information with targeted measures of local genetic correlation and distance.

## References

- Comaniciu D and Meer P. 2002. Mean shift: A robust approach toward feature space analysis. IEEE Trans Patt Mach Int 24: 603-619.

- Fuller DQ, Sato YI, Castillo C, Qin L, Weisskopf AR, Kingwell-Banham EJ, Song J, Ahn SM and Van Etten, J. 2010. Consilience of genetics and archaeobotany in the entangled history of rice. Archae and Anthro Sci 2: 115-131.

- Giardine B, Riemer C, Hardison RC, Burhans R, Elnitski L, Shah P, Zhang Y, Blankenberg D, Albert I, Taylor J and Miller W. 2005. Galaxy: a platform for interactive large-scale genome analysis. Genome research, 15: 1451-1455.

## Supplementary Material

A series of [jupyter notebooks](https://github.com/SantosJGND/Genetic-data-analysis) was created in accompaniment to the development of this
pipeline.