## Downstream Analysis

**Warning** The development of this branch has been moved to an independent repository.

Link to be provided soon.

### -- script: code_track.py

**Context.** The endpoint of association analysis in population genetics. Following local classificaiton, estimates are available of the genome-wide 
proximity of each accession to each reference population in play. However, reference populations can be large and structured, in particular 
older populations that might have amassed considerable genetic variation through the years, be it through mutation or the incorporation 
of novel genetic material from cryptic sources.

Our supervised classification has served the purpose of simplifying our genomic landscape into a handful of signals we can read. We can now 
pinpoint contradictory signals along the genomes of particular accessions which would, in a genome-wide analysis, have resulted in intermediate 
and even outlier positions. Separate analyses of cleaner groups of signals can inform us on the structure of a given group in more depth.

Because we allow for intermediate classes for plotting and to guide other steps of the analysis, population labels are recoded. 
In the following section, we shall discuss how to perform a targeted analysis of cluster associations on the data extracted 
during the genome crawl. In order to direct this code the user must know the label-to-code map produced. code_track.py takes the 
reference file used for the supervised analysis to replicate the codes produced and the colors associated with each code in the 
ideogram plots.

### script: targeted_analysis.py

The script *targeted_analysis.py* takes one or multiple target labels (read code_track.py output, column _code_) and a list of individuals the 
user wishes to focus on. The script proceeds to identify the windows of the genome crawl at which the focal accessions are classed to the code 
provided. Profiles of clusters these individuals are associated with at these windows are then queried (_--MSprint_ must have been toggled 
when first running the crawl). Individual-to-profile association is done using a minimum threshold defined by the user. 

Total profiles extracted are reduced through principal component analysis and Kmeans clustering performed in feature space.

**Input**

- `.geno` file, position. file name must include string `chr` followed by chromosome integer.

- `--fam`  accompanying .fam file to .geno. Used here because the function `read_reds` was not adapted to this layer.

- `--ref`  reference accession file. Must be same as used for supervised analysis. 

- `--focus` one column file. Which accessions to quer for selected classification.

**Parameters**: `--CHR`: direct which chromosome to run analyses on. Will use all .geno 
files available if not provided; `--start`: integer. Position to begin query for windows, 
in base pairs; `--end`  integer. position to end query for windows, in base pairs. `--target`: integer. 
Target classification to query. see `code_track.py`; `--id`: string. Name your analysis, 
incorporated in filenames; `--height` and `--width`: Output ideogram size, in inches; 
`--threshold`: *p*-value comparison threshold for intermediate classification; `--outlier`: lower 
*p*-value threshold for outlier classification; `--coarse`: avoid smoothing. savgol filter filter
used other wise. `--plot`: controls ideogram output. `--reduc`: reduces PCA projections for jupyter 
plots.

**Output**

- *_comp_*: queried cluster PCA projections.
- *_request_*: accession projections.
- *_coordinates_*: genome coordinates of clusters queried.
- *_Profile_*: accession level average likelihood across cluster groups (Kmeans).

Through the *--app* option, a python dash application is produced for the exploration of associations revealed in the profiles extracted and 
their location along the genome of focal accessions. Procfile, requirements and gitignore files are automatically produced so that the 
appplication can be readily set up on any server as well as locally.

- [command example](Cluster_target_command.txt)

#### Output analysis:

- targeted KDE summary: Vector cluster analysis against accession structure, average likelihood, geographic origin and genomic coordinates. 

> [notebook](https://nbviewer.jupyter.org/github/SantosJGND/Galaxy_KDE_classifier/blob/master/Downstream_functions/Analyses_Jsubtrop_self_KDE/Targeted_analysis_plot.ipynb)


### script: Gen_distances_Layers.py

The script *Gen_distances_Layers.py* adds another layer of information to the analysis conducted through the 
targeted analysis of cluster membership provided by *targeted_analysis.py*. The user inputs the coordinates file output
of the targeted analysis and indicates which cluster group to focus on. Local genomic windows are queried by 
cluster queried. Local haplotypes assigned to queried cluster are identified by maximum *p*-value. Centroid
of target haplotypes is calculated in PCA feature space. Centroid of reference accessions in `ref` file, excluding
target haplotypes. Euclidean distances of all accessions to target and reference centroids are calculated.

**Input**
- *Blocks_profiles*: output of *Galaxy_KDE.py*
- `--genoSuf`: suffix of geno files included in coordinates file. must be followed by int(chr)+ '.txt'
- `--ref`  reference accession file. Must be same as used for supervised analysis. 
- `--focus` one column file. Includes accessions if absent from `ref`.
- `--coords` cluster coordinates file.
- `--code` Kmeans cluster code in coordinates cluster.
- `--fam`  accompanying .fam file to .geno. Used here because the function `read_reds` was not yet adapted to this layer.
- `--bim` plink .bim file.

`--id`: string. Name your analysis. Incorporated in filenames; `--height` and `--width`: Output ideogram size, in inches; 
`--plot`: controls ideogram output. `--reduc`: reduces PCA projections for jupyter plots.


**Output**

- *_comp_*: queried distances PCA projections.
- *_request_*: accession projections.
- *_coordinates_*: genome coordinates of clusters queried.
- *_Profile_*: accession level average likelihood across cluster groups (Kmeans).
- *layer_analysis.p*: Cluster analyses using different cluster methods. pickle serialization.


- [command example](Distance_command_command.txt)

#### Output analysis:

- KDE targeted distance analysis.

> [notebook](https://nbviewer.jupyter.org/github/SantosJGND/Galaxy_KDE_classifier/blob/master/Downstream_functions/JapanKorea_tropical_Rdist/model_Rdist.ipynb)


## Example test studies. 
*python applications running on [heroku](https://www.heroku.com/platform) servers.*

Apps were modified to read more than one output. 

Targeted studies combining supervised haplotype classification with mean shift cluster analysis. 

**warning**: Heroku servers are 512 Mb RAM and take one user at a time only.

- Genome-wide association to a target reference: [App_global](https://cbasmati-japonica-signals.herokuapp.com/)

- Chromosome specific association to target references with regional focus: [App_regional](https://cbasmati-chr08-examples.herokuapp.com/)

- Association to target references at gene of interest: [App_gene](https://sh4-gene.herokuapp.com/)

- Individual specific association to target references: [App_ind](https://iris-313-12074.herokuapp.com/)

