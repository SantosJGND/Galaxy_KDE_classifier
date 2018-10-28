## Downstream Analysis

### -- script: code_track.py

Because we allow for intermediate classes for plotting and to guide other steps of the analysis, population labels are recoded. 
In the following section, we shall discuss how to perform a targeted analysis of cluster associations on the data extracted 
during the genome crawl. In order to direct this code the user must know the label-to-code map produced. code_track.py takes the 
reference file used for the supervised analysis to replicate the codes produced and the colors associated with each code in the 
ideogram plots.

### -- script: targeted_analysis.py

The script *targeted_analysis.py* takes one or multiple target labels (read code_track.py output, column _code_) and a list of individuals the 
user wishes to focus on. The script proceeds to identify the windows of the genome crawl at which the focal accessions are classed to the code 
provided. Profiles of clusters these individuals are associated with at these windows are then queried (_--MSprint_ must have been toggled 
when first running the crawl). Individual-to-profile association is done using a minimum threshold defined by the user. 

Total profiles extracted are reduced through principal component analysis and Kmeans clustering performed in feature space.

Through the *--app* option, a python dash application is produced for the exploration of associations revealed in the profiles extracted and 
their location along the genome of focal accessions. Procfile, requirements and gitignore files are automatically produced so that the 
appplication can be readily set up on any server as well as locally.

**Context.** The endpoint of association analysis in population genetics. Following local classificaiton, estimates are available of the genome-wide 
proximity of each accession to each reference population in play. However, reference populations can be large and structured, in particular 
older populations that might have amassed considerable genetic variation through the years, be it through mutation or the incorporation 
of novel genetic material from cryptic sources.

Our supervised classification has served the purpose of simplifying our genomic landscape into a handful of signals we can read. We can now 
pinpoint contradictory signals along the genomes of particular accessions which would, in a genome-wide analysis, have resulted in intermediate 
and even outlier positions. Separate analyses of cleaner groups of signals can inform us on the structure of a given group in more depth.
