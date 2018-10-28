# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 15:41:52 2018

@author: jgarcia
"""


from Kernel_tools import *
########## START HERE #############

import os
import argparse

import collections
import time
import itertools as it
import numpy as np
import re

import pandas as pd

import Galaxy_summary_tools

parser = argparse.ArgumentParser()

### optional arguments
parser.add_argument("books",type=str,metavar= 'N',nargs= '+',
                    help = "Reference files to read. Any number can be given.")

parser.add_argument("--bim",help = "snp information bim format.")

parser.add_argument("--loci",action = "store_true",help = "write allele assignment by individual")

parser.add_argument("--summary",action = "store_true",help = "write classification summary.")

parser.add_argument("--focus",type= str,help = "IDs of accessions to plot.")

parser.add_argument("--id",type= str,help = "ID for this particular plot. if not given, name of last accession in focus.")

parser.add_argument("--CHR",type= int,help = "chromosome to draw ideogram of.")

parser.add_argument("--start",type= int,help = "where to begin, in markers. Only makes sense if --CHR is also used.")

parser.add_argument("--end",type= int,help = "where to end, in markers. Only makes sense if --CHR is also used.")

parser.add_argument("--out",type= str,default= '',help = "output directory")

parser.add_argument("--coarse",action='store_false',help= 'to smooth or not to smooth.')

parser.add_argument("--bin",default = 5,type= int,help = "smoothing parameter, must be uneven [savgol filter]")

parser.add_argument("--sg_order",default = 3,type= int,help = "staviksy golay filter order")

parser.add_argument("--outlier",type=float,default = 1e-3,help = "Outlier threshold")

parser.add_argument("--threshold",type = float,default = 1.6,help = "Intermediate classification threshold")

parser.add_argument("--chrom_height",type= float, default= 1, help= "height of ideograms")

parser.add_argument("--chrom_gap",type= float,default= 0,help= "gap between ideograms.")

parser.add_argument("--height",type= float, default= 10,help= "figure height, in inches.")

parser.add_argument("--width",type= float,default= 20,help= "figure width, in inches.")

parser.add_argument('--xticks',type= int,default= 100000,help= 'xticks on final ideogram')

args = parser.parse_args()


Home= args.out



####
#### read block and profile files
####

print('To begin reading from: ')
print(args.books)

Ref_profiles, Names, Out = read_3D_profiles_list(args.books)

Home= args.out

######

if args.focus:
    Focus = read_focus(args.focus)
else:
    Focus = Names

focus_indexes= [x for x in range(len(Names)) if Names[x] in Focus]


Blocks = Merge_class(Ref_profiles,focus_indexes,Out,args.threshold,args.bin,args.outlier)

print("Number of chromosomes selected for analysis: {0}".format(len(Blocks)))

####
####

if args.CHR:
    chromosomes= [args.CHR]
else:
    chromosomes= [x for x in Blocks.keys()]


if args.start:
    Blocks= {Chr:{x:Blocks[Chr][x] for x in Blocks[Chr].keys() if x >= args.start} for Chr in Blocks.keys()}
    if not args.CHR:
        print("start was selected with no CHR specification.")

if args.end:
    Blocks= {Chr:{x:Blocks[Chr][x] for x in Blocks[Chr].keys() if x <= args.end} for Chr in Blocks.keys()}
    if not args.CHR:
        print("end was selected with no CHR specification.")


chromosome_list = []

Ideo = []

for here in range(len(Focus)):
    Subject = Focus[here]
    
    chromosome_list.extend(['chr'+str(Chr)+ '_' + Subject for Chr in chromosomes])
    
    color_ref= ['red','yellow','blue','black','orange','purple','green','silver','red3','deepskyeblue','navy','chartreuse','darkorchid3','goldenrod2']
    Stock = [[['chr'+str(Chr)+ '_' + Subject,bl,Out[Chr][bl],color_ref[Blocks[Chr][bl][here] - 1]] for bl in sorted(Blocks[Chr].keys())] for Chr in chromosomes]
    Stock = [y for y in it.chain(*[z for z in it.chain(*[Stock])])]
    
    Ideo.extend(Stock)

###########################################
###########################################

# Read in ideogram.txt, downloaded from UCSC Table Browser
ideo = pd.DataFrame(Ideo,columns = ['chrom', 'start', 'end', 'gieStain'])

# Filter out chromosomes not in our list
ideo = ideo[ideo.chrom.apply(lambda x: x in chromosome_list)]

# Colors for different chromosome stains
color_lookup = {
    'red': [255, 0, 0],
    'yellow': [255, 255, 0],
    'blue': [0, 0, 255],
    'orange': [255, 165, 0],
    'green': [50, 205, 50],
    'black': [0, 0, 0],
    'purple': [128, 0, 128],
    'silver': [211, 211, 211],
}




ideo = compress_ideo(ideo,chromosome_list)

#ideo = ideo[(ideo.start > 5.06e6) & (ideo.start < 7.06e6)]

ideo['colors'] = ideo['gieStain'].apply(lambda x: tuple([round(y / float(255),1) for y in color_lookup[x]]))
# Add a new column for width
ideo['width'] = ideo.end - ideo.start

### Sum compressed window sizes by class and individual.
###
if args.summary:
    
    Compass= []
    
    for ind in ideo.ID.unique():
        for Chr in ideo.chrom.unique():
            for color in ideo.gieStain.unique():
                sub= ideo[(ideo.ID== ind) & (ideo.gieStain == color) & (ideo.chrom == Chr)]
                Compass.append([ind,Chr,color,sum(sub.width),np.mean(sub.width),len(sub)])
    
    Compass= pd.DataFrame(Compass,columns= ['ID','chrom','color','t_length','mean_size','N'])
    
    Output= open('Summary_assignments_CHR' + str(chromosomes[-1]).zfill(2) + '_Z' + str(args.threshold) + '_bin' + str(args.bin) + '.txt','w')
    
    Output.write('\t'.join(Compass.columns))
    Output.write('\n')
    
    for gene in range(len(Compass)):
        Output.write('\t'.join([str(c) for c in Compass.iloc[gene,:]]))
        Output.write('\n')
    
    Output.close()

### Output classification of individual loci based on compressed windows.
### Not advised. Rough compression not suited for locus-specific analysis.
### Better solution to be developped.

if args.bim and args.loci:
    print('loci classification requested.')
    
    Ind_home= Home + 'loci/'
    
    print('Merging..'
    )
    
    MissG, Gindex = BIMread(args.bim)
    
    guys= ideo.ID.unique()
    
    for guy in guys:
        
        for Chr in ideo.chrom.unique():
            snps= sorted([x for x in Gindex[Chr].keys()])
            
            sub= ideo[(ideo.ID==guy) & ideo.chrom==Chr]
            sub= sub.sort_values('start')
            
            ## this is cheating, but next Analysis will include every snp.
            sub.loc[(sub.start==max(sub.start)),'end']= max(snps)
            
            library= []
            included= []
            for step in range(len(sub)):
                    box= [x for x in snps if x >= sub.start.iloc[step] and x <= sub.end.iloc[step]]
                    included.extend(box)
                    truth= ' '.join([str(int(color_ref[x] == sub.gieStain.iloc[step])) for x in range(4)])
                    
                    truth= np.tile(truth,(len(box),1))
                            
                    library.extend(truth)
            
            library= np.array(library)
            
            filename= Ind_home + guy + 'CHR'+str(Chr).zfill(2)+'_loci_class.txt'
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            Output= open(filename,'w')
            
            for ellen in library:
                Output.write(ellen[0] + '\n')
            
            Output.close()

print('Done.')