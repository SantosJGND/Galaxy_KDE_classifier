
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 15:25:38 2017

@author: jgarcia
"""

from Kernel_tools import *
from Galaxy_Ideogram_tools import *

import os
import argparse
parser = argparse.ArgumentParser()

### optional arguments
parser.add_argument("books",type=str,metavar= 'N',nargs= '+',
                    help = "Reference files to read. Any number can be given.")

parser.add_argument("--focus",type= str,help = "IDs of accessions to plot.")

parser.add_argument("--id",type= str,help = "ID for this particular plot. if not given, name of last accession in focus.")

parser.add_argument("--CHR",type= int,help = "chromosome to draw ideogram of.")

parser.add_argument("--start",type= int,help = "where to begin, in markers. Only makes sense if --CHR is also used.")

parser.add_argument("--end",type= int,help = "where to end, in markers. Only makes sense if --CHR is also used.")

parser.add_argument("--out",type= str,default= '',help = "output directory")

parser.add_argument("--clust",type= str,default= 'None',help = "Do you wish to cluster ideograms based on similarity? provide method for scipy linkage function")

parser.add_argument("--coarse",action='store_false',help= 'to smooth or not to smooth.')

parser.add_argument("--bornes",type= int,default= 0,help = "bornes")

parser.add_argument("--square",action= 'store_true',help= 'Draw rectangle around region of interest. only makes sense when "bornes" argument is used too.')

parser.add_argument("--bin",default = 5,type= int,help = "smoothing parameter, must be uneven [savgol filter]")

parser.add_argument("--sg_order",default = 3,type= int,help = "staviksy golay filter order")

parser.add_argument("--outlier",type=float,default = 1e-3,help = "Outlier threshold")

parser.add_argument("--threshold",type = float,default = 1.6,help = "Intermediate classification threshold")

parser.add_argument("--chrom_height",type= float, default= 1, help= "height of ideograms")

parser.add_argument("--chrom_gap",type= float,default= 0,help= "gap between ideograms.")

parser.add_argument("--fontsize",type= float,default= 5,help= "label font size.")

parser.add_argument("--height",type= float, default= 10,help= "figure height, in inches.")

parser.add_argument("--width",type= float,default= 20,help= "figure width, in inches.")

parser.add_argument('--xticks',type= int,default= 100000,help= 'xticks on final ideogram')

args = parser.parse_args()


Home= args.out

if len(Home) > 0:
    Home= Home + '/'

########## Complementary files.

def read_focus(index_file):
    indxs = []
    
    Input = open(index_file,'r')
    for line in Input:
        line = line.split()
        indxs.append(line[0])
    
    Input.close()
    
    return indxs



import collections
import time
import itertools as it
import numpy as np
import re


def recursively_default_dict():
        return collections.defaultdict(recursively_default_dict)

####
#### read block and profile files
####

print('To begin reading from: ')
print(args.books)

Ref_profiles, Names, Out = read_3D_profiles(args.books)


######
######

if args.focus:
    Focus = read_focus(args.focus)
else:
    Focus = Names


Absent= [x for x in Focus if x not in Names]

if len(Absent) > 0:
    print('The following individuals were not found in the files provided: {}'.format(Absent))
    print('Analysis will proceed without them.')
    Focus= [x for x in Focus if x in Names]

focus_indexes = [Names.index(x) for x in Focus]


if args.coarse:
    print('statistics will be smoothed, using a savgol filter of order {0} and a bin size of {1}'.format(args.sg_order,args.bin))

Blocks, N_pops = Merge_class(Ref_profiles,focus_indexes,Out,args.threshold,args.bin,args.outlier,args.coarse,args.sg_order)

print("Number chromosomes selected: {0}".format(len(Blocks)))
####
####

if args.CHR:
    chromosomes= [args.CHR]
else:
    chromosomes= Blocks.keys()


if args.start:
    Blocks= {Chr:{x:Blocks[Chr][x] for x in Blocks[Chr].keys() if Out[Chr][x] >= args.start - args.bornes} for Chr in Blocks.keys()}
    if not args.CHR:
        print("start was selected with no CHR specification.")

if args.end:
    Blocks= {Chr:{x:Blocks[Chr][x] for x in Blocks[Chr].keys() if x <= args.end + args.bornes} for Chr in Blocks.keys()}
    if not args.CHR:
        print("end was selected with no CHR specification.")


#################### define color refs:
color_ref= ['red','yellow','blue','black','orange','purple','green','silver','red3','deepskyeblue','navy','chartreuse','darkorchid3','goldenrod2']

prim_colors= ['red','yellow','blue']

if N_pops <= 3:
    step= prim_colors[:N_pops]
    step.extend(color_ref[3:])
    color_ref= step


#################### PROCEED
####################

chromosome_list = []

Ideo = []

for here in range(len(Focus)):
    Subject = Focus[here]
    
    chromosome_list.extend(['chr'+str(Chr)+ '_' + Subject for Chr in chromosomes])
    Stock = [[['chr'+str(Chr)+ '_' + Subject,bl,Out[Chr][bl],color_ref[Blocks[Chr][bl][here] - 1]] for bl in sorted(Blocks[Chr].keys())] for Chr in chromosomes]
    Stock = [y for y in it.chain(*[z for z in it.chain(*[Stock])])]
    
    Ideo.extend(Stock)


###########################################
###########################################
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
import pandas as pd


# Height of each ideogram
chrom_height = args.chrom_height

# Spacing between consecutive ideograms
chrom_spacing = args.chrom_gap

# Height of the gene track. Should be smaller than `chrom_spacing` in order to
# fit correctly
gene_height = 0.0

# Padding between the top of a gene track and its corresponding ideogram
gene_padding = 0.0


# Keep track of the y positions for ideograms and genes for each chromosome,
# and the center of each ideogram (which is where we'll put the ytick labels)
ybase = 0
chrom_ybase = {}
gene_ybase = {}
chrom_centers = {}

# Iterate in reverse so that items in the beginning of `chromosome_list` will
# appear at the top of the plot
for chrom in chromosome_list[::-1]:
    chrom_ybase[chrom] = ybase
    chrom_centers[chrom] = ybase + chrom_height / 2.
    gene_ybase[chrom] = ybase - gene_height - gene_padding
    ybase += chrom_height + chrom_spacing

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

# Add a new column for colors
ideo = compress_ideo(ideo,chromosome_list, Out)

ideo['colors'] = ideo['gieStain'].apply(lambda x: tuple([round(y / float(255),1) for y in color_lookup[x]]))
# Add a new column for width
ideo['width'] = ideo.end - ideo.start

# Width, height (in inches)
figsize = (args.width, args.height)

fig = plt.figure(figsize=figsize)
ax = fig.add_subplot(111)

if args.square:
    from matplotlib.patches import Rectangle
    someX= args.start
    someY= -10
    
    currentAxis= plt.gca()
    currentAxis.add_patch(Rectangle((someX,someY),
    args.end-args.start,len(chromosome_list) + 20, 
    fill= None, alpha= 1, linewidth= 5))

# Now all we have to do is call our function for the ideogram data...
print("adding ideograms...")
for collection in chromosome_collections(ideo, chrom_ybase, chrom_height, edgecolors=None, linewidths= 0):
    ax.add_collection(collection)

# Axes tweaking
ax.set_xticks([x for x in range(min(ideo.start),max(ideo.end),int(args.xticks))])
plt.xticks(fontsize = 5,rotation = 90)
ax.tick_params(axis = 'x',pad = 10)

ax.tick_params(axis='y', which='major', pad=30)
ax.set_yticks([chrom_centers[i] for i in chromosome_list])
ax.set_yticklabels(chromosome_list, fontsize = args.fontsize)
ax.axis('tight')

if args.id:
    Subject= args.id

plt.savefig(Home + 'Ideo_' + Subject +'_CHR' + str(chromosomes[-1]).zfill(2)+'_st' + str(min(ideo.start)) + '_Z' +str(args.threshold)+ '_bin'+ str(args.bin)+'.png',bbox_inches = 'tight')

print('Done.')
