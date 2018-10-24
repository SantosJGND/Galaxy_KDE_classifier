# Kernel classification.

Local haplotype classification and data exploration. This software was built with the purpose of comparing local 
genetic distributions of diverging biological populations. The comparison of distributions is accompanied by a set of 
visualisation tools we deemed useful for inference of evolutionary relationships.

_Note_: Application of this pipeline of analysis should take into considerations the particular qualities of the data set for which 
it was developped. The requirements for an informative output are only nearly met given the stable and high density of genetic markers,
 a state of near complete homozygosity, to add to the quality of the data.

 - Fuller DQ, Sato YI, Castillo C, Qin L, Weisskopf AR, Kingwell-Banham EJ, Song J, Ahn SM and Van Etten, J. 2010. Consilience of genetics and archaeobotany in the entangled history of rice. Archae and Anthro Sci 2: 115-131.
 

## Kernel Density Estimation -- script: Kernel_mPLib3D_FM36_Galaxy.py

The initial analysis is conducted following a machine learning approach 
of dimensionality reduction and kernel density estimation. The latter was selected for the lack of local 
priors in the convoluted history of modern domesticated populations (Fuller et al. 2010). 
The Kernel Density Estimation is used for the identification both of local mislabelling (genetic exchanges between 
reference populations) and the presence of outlier material among hybrids (possible introgressions from  
differentiated sources). Follow this series of series of [jupyter notebooks](https://github.com/Joaos3092/Genetic-data-analysis) 
for a step-by-step explanation on how KDE can be used for these purposes.

This analysis is conducted on the basis of a crawling window of fixed size of genetic markers. Window size, 
dimensionality reduction, the number of axis to retain in feature space and the level of heterozygosity allowed 
for are customizable. If PCA is chosen for dimensionality reduction, the user can chose to print the variance 
of each axis along the crawl (--VARprint, boolean).

Parallel to the estimation of local supervised distributions, an unsupervised clustering using the MeanShift algorithm 
(Comaniciu & Meer, 2002) is performed. KD estimates of each cluster identified are stored. These are necessary for downstream analysis 
but are not automatically printed.

While the software could in theory handle more than three reference populations, we advise this number to be limited to three. 
This is because, for plotting as well as targeted cluster analysis, we allowed for intermediate classes among the references 
provided. At K=3, this results in 8 classes in total (3 pure classes, 3 two-way intermediates and 1 three way intermediate). 
At K == 4, this would imply 6 two way intermediate classes. Unless reference populations are strongly structured, the resulting plots 
can become messy.


## Ideogram plotting -- script: Ideogram_plots_Galaxy.py

The first layer of analysis following the genome crawl should be the visualization of local statistics pertaining 
to the distributions of the reference populations provided. This step also constitutes an exploration of variation 
in our data set. The purpose of this exploration is to visualize possible exchanges of genetic material among 
reference material, the local composition of non-reference accessions included in the analysis, and to investigate 
the presence of differentiated (outlier) material.

For this purpose, the construction of the _chromosome painting_ along pure and intermediate classes allows for user 
defined thresholds of comparison and an outlier threshold.

A savitsky golay filter is automatically applied if the flag '--coarse' is not passed. We recommend the use of this flag.

# Downstream Analysis

## -- script: code_track.py

Because we allow for intermediate classes for plotting and to guide other steps of the analysis, population labels are recoded. 
In the following section, we shall discuss how to perform a targeted analysis of cluster associations on the data extracted 
during the genome crawl. In order to direct this code the user must know the code-label map produced. code_track.py takes the 
reference file used for the supervised analysis to replicate the codes produced and the colors associated with each code in the 
ideogram plots.

## -- script: gff_Blocksmerge.py combining information 

The genome crawl divides the genome into a set of overlapping windows. Our downstream analysis queries these windows for associations. 
For various reasons it can be usefull to know which genes are present at each window along the genome. The script gff_Blocksmerge.py 
produces this merged data set. it takes a gff file, and the user can chose to focus on particular elements (miRNA, genes or both), as 
well as the columns of information on each to extract, by name.

This information can be easily included in the end product dash application.

## -- script: targeted_analysis.py

The endpoint of association analysis in population genetics. Estimates have been drawn genome wide of the proximity of each accession to each 
reference population in play. However, reference populations can be large and structured, in particular older populations that have amassed 
considerable genetic variation through the years, be it through mutation or the incorporation of novel genetic material from cryptic sources.

The idea is this: from a purely informational point of view, our supervised classification has served the purpose of simplifying our genomic 
landscape into a handful of signals we can read. We can now pinpoint contradictory signals along the genomes of particular accessions which would, 
in a genome-wide analysis, have resulted in intermediate or even outlier positions. Separate analyses of cleaner groups of signals can inform us 
on the structure of a given group in more depth.

The script *targeted_analysis.py* takes one or multiple target labels (read code_track.py output, column _code_) and a list of individuals the 
user wishes to focus on. The script proceeds to identify the windows of the genome crawl at which the focal accessions are classed to the code 
provided. Profiles of clusters these individuals are associated with at these windows are then queried (_--MSprint_ must have been toggled 
when first running the crawl). Individual-to-profile association is done using a minimum threshold defined by the user. 

Total profiles extracted are reduced through principal component analysis and Kmeans clustering performed in feature space.

Through the *--app* option, a python dash application is produced for the exploration of associations revealed in the profiles extracted and 
their location along the genome of focal accessions. Procfile, requirements and gitignore files are automatically produced so that the 
appplication can be readily set up on any server as well as locally.

- Comaniciu D and Meer P. 2002. Mean shift: A robust approach toward feature space analysis. IEEE Trans Patt Mach Int 24: 603-619.
