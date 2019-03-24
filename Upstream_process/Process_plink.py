import subprocess
import numpy as np
import itertools as it
import time
import os

from Kernel_tools import IDread

import argparse
parser = argparse.ArgumentParser()

parser.add_argument("--plinkfile",type=str,help = "plink file to use. Full path")

parser.add_argument("--Nsample",type=int,default= 0,help = "Number of SNPs to sample")

parser.add_argument("--chr",type=str,default= '',help = "chromosomes to sample, comma sep for multiple")

parser.add_argument("--thin",type= float,default= 0,help = "percentage of SNPs to keep")

parser.add_argument("--region",type= str,default= '',help = "physical region to extract")

parser.add_argument("--maf",type= float,default= 0,help = "min allele freq filter")

parser.add_argument("--miss",type= float,default= 0,help = "miss data freq filter. locus.")

parser.add_argument("--mind",type= float,default= 0,help = "miss data freq filter. ind. ")

parser.add_argument("--hwe",type= float,default= 0,help = "hwe significance filter")

parser.add_argument("--vcf",action= "store_true", help= "if given reformats to .vcf")
### 
parser.add_argument("--geno",action= "store_true", help= "if given and vcf == True reformats vcf to .geno")
### 


args = parser.parse_args()


###
plink_file= args.plinkfile
Nsample= args.Nsample
thin= args.thin
region= args.region

if region:
    region= region.split(',')
    region= [int(x) for x in region]

###
Home= time.strftime("%d-%m-%Y")

print('Processing date: ' + Home)
print('source bfile: ' + plink_file)

Bimfile= plink_file + '.bim'

Miss, Gindex = IDread(Bimfile)

chromosomes= Miss.keys()

if args.chr:
    chromosomes= args.chr.split()
    chromosomes= [int(x) for x in chromosomes]

Absent= [x for x in chromosomes if x not in Gindex.keys()]

if Absent:
    chromosomes= [x for x in chromosomes if x in Gindex.keys()]
    print('chromosomes {} absent from .bim file. Will proceed without them.')
    print('Will return error if no chr index is found in .bim file.')


if Nsample:
    Sample= {
        z: list(sorted(np.random.choice(list(Gindex[z].keys()),Nsample,replace= False))) for z in chromosomes
    }

    extract_files= {}

    for z in Sample.keys():
        
        snp_list= 'extract_chr{}.txt'.format(str(z).zfill(2))
        
        output= open(snp_list,'w')
        output.write('\n'.join([str(int(x)) for x in Sample[z]]))
        output.close()
        extract_files[z] = snp_list


for Chr in chromosomes:
    out_plink= {}

    plink_args= {
        '--bfile': plink_file,
        '--chr': str(Chr)
    }

    if Nsample:
        outfile= 'Extract_Chr{}_{}'.format(Chr,Nsample)
        out_plink[Chr]= outfile
        
        plink_args['--extract']= extract_files[Chr]
        plink_args['--make-bed --out']= outfile

    if thin:
        outfile= 'Extract_Chr{}_t{}'.format(Chr,thin)
        out_plink[Chr]= outfile
        plink_args['--thin']= str(thin)
        plink_args['--make-bed --out']= outfile

    if region:
        
        outfile= 'Extract_Chr{}_R{}-{}'.format(Chr,region[0],region[1])
        #region= [int(x + min(Gindex[Chr].keys()) - Miss[Chr][0][0]) for x in region]
        
        out_plink[Chr]= outfile
        plink_args['--from-bp']= str(region[0])
        plink_args['--to-bp']= str(region[1])
        plink_args['--make-bed --out']= outfile

    ### form command
    print('SNPs in chromosome: {}'.format(len(Miss[Chr])))
    
    if len(plink_args) > 1:
        plink_command= ['plink']
        plink_command.extend(list(it.chain(*[[v,g] for v,g in plink_args.items()])))

        ## load plink
        ## 
        plink_house= 'bioinfo/plink/1.90b3k'
        call_plink= ['module','load',plink_house]
        subprocess.run(' '.join(call_plink), check=True, shell=True)

        ### run extraction
        subprocess.run(' '.join(plink_command), check=True, shell=True)

    else:
        outfile= plink_file
        out_plink[Chr]= outfile
        
    ###
    ### Filter
    filt_args= {
        '--bfile': out_plink[Chr],
    }

    filfile= 'Filtered_Chr{}'.format(Chr)

    if args.maf:
        filt_args['--maf']= str(args.maf)
        filfile += '_maf{}'.format(args.maf)

    if args.mind:
        filt_args['--mind']= str(args.mind)
        filfile += '_mind{}'.format(args.mind)

    if args.miss:
        filt_args['--geno']= str(args.miss)
        filfile += '_miss{}'.format(args.miss)

    if args.hwe:
        filt_args['--hwe']= str(args.hwe)
        filfile += '_hwe{}'.format(args.hwe)


    ### Filter if necessary
    if len(filt_args) > 1:
        ### form command
        print('hi')
        filt_args['--make-bed --out']= filfile
        plink_command= ['plink']
        plink_command.extend(list(it.chain(*[[v,g] for v,g in filt_args.items()])))

        ##run filter
        subprocess.run(' '.join(plink_command), check=True, shell=True)
        
        out_plink[Chr]= filfile
        

    ###
    ### convert to vcf
    if args.vcf:
        for z in out_plink.keys():
            
            vcf_convert= ['plink','--bfile',out_plink[z],'--recode','vcf','--out',out_plink[z]]
            subprocess.run(' '.join(vcf_convert), check=True, shell=True)


    ### convert vcf to geno if args.geno is given
    sNMF_house= 'bioinfo/sNMF_CL/1.2'

    if args.vcf and args.geno:
        call_snmf= ['module','load',sNMF_house]
        subprocess.run(' '.join(call_snmf), check=True, shell=True)
        
        for z in out_plink.keys():
            run_snmf= ['vcf2geno',out_plink[z] + '.vcf',out_plink[z] + '.geno']
            subprocess.run(' '.join(run_snmf), check=True, shell=True)
        
