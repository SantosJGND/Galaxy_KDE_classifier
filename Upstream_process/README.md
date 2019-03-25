## Upstream process

The KDE pipeline of analysis uses a combination of input formats. The 
genotype data is read in .geno format, but accession and SNP information is read 
from plink files `.fam` and `.bim` respectively. 

This repository proposes a script, `Process_plink.py`, to automate the initial 
step of formatting.

This script re-purposes filtering commands from the plink toolkit to pre-process
data, generating a new plink file set. Then, if the argument `--vcf` is given, 
the filtered data is reformated to vcf. Finally, if the `--geno` argument is also
provided, a .geno file is generated. This final step uses the sNMF function 
*vcf2geno*. The `--vcf` option must have been toggled.

For the script to work, the software sNMF (>= v.1.2) and plink (>= v.1.9) must be
installed.

[command example](process_command.txt)

