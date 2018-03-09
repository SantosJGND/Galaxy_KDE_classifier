# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 14:12:22 2018

@author: jgarcia
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 15:17:06 2018

@author: jgarcia
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 15:25:38 2017

@author: jgarcia
"""

from Kernel_tools import *
########## START HERE #############

import os
import argparse
parser = argparse.ArgumentParser()

### optional arguments
parser.add_argument("books",type=str,metavar= 'N',nargs= '+',
                    help = "Reference files to read. Any number can be given.")

parser.add_argument("--focus",type= str,help = "IDs of accessions to plot.")

parser.add_argument("--CHR",type= int,help = "chromosome to draw ideogram of.")

parser.add_argument("--start",type= int,help = "where to begin, in markers. Only makes sense if --CHR is also used.")

parser.add_argument("--end",type= int,help = "where to end, in markers. Only makes sense if --CHR is also used.")

parser.add_argument("--out",type= str,default= '',help = "output directory")

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



def Merge_class(Ref_profiles,focus_indicies,Out,Diff_threshold,BIN,X_threshold):
    Blocks_genome = recursively_default_dict()
    
    for CHR in Ref_profiles.keys():
        print(CHR)
        Points = sorted(Ref_profiles[CHR].keys())
        Likes = Ref_profiles[CHR]
        N_pops= len(Likes[[x for x in Likes.keys()][0]])
        
        print("number of reference populations: {0}".format(N_pops))
        Likes = {x:[Likes[bl][x] for bl in sorted(Likes.keys())] for x in range(N_pops)}
        Likes = {x:np.array(Likes[x]) for x in Likes.keys()}
        
        
        Topo = []
        
        #range_Parents = [x + Aro.shape[0] for x in range(Daddy.shape[0])]
        #range_Crossed = [x for x in range(Aro.shape[0])]
        
        for acc in focus_indicies:
            Guys = np.array([Likes[x][:,acc] for x in range(N_pops)])
            Guys = np.nan_to_num(Guys)
            
            
            Test = [int(x < X_threshold) for x in np.amax(np.array(Guys),axis = 0)]
            Test = savgol_filter(Test,BIN,3,mode = "nearest")
            Test = [round(x) for x in Test]
            
            Guys = [[[y,0][int(y<=X_threshold)] for y in x] for x in Guys]
            Guys = [savgol_filter(x,BIN,args.sg_order,mode = "nearest") for x in Guys]
            #
            Guys = np.array(Guys).T
            
            maxim = np.argmax(Guys,axis = 1)
            where_X = [x for x in range(Guys.shape[0]) if Test[x] == 1]
            #where_X = [x for x in range(Guys.shape[0]) if len([c for c in Guys[x,:3] if c <= .0001]) == 3]
            Consex = [x for x in it.combinations(range(N_pops),2)]
            if Consex:
                for h in range(len(maxim)):
                    CL = []
                    for j in Consex:
                        Diff = Guys[h,j]
                        if maxim[h] not in j or len([x for x in Diff if x < X_threshold]) > 0:
                            continue
                        if max(Diff) <= X_threshold:
                            Diff = 0
                        else:
                            #Diff = int(len([x for x in Diff if x > 1e-5]) == 1)
                            Diff = abs(max(Diff)) / abs(min(Diff))
                            Diff = int(Diff > Diff_threshold)
                        
                        #print(Diff)
                        if Diff == 0:
                            CL.append(j)
                    
                    if len(CL) == 2:
                        maxim[h] = 7
                    if len(CL) == 1:
                        maxim[h] = sum(CL[0]) + N_pops
            
            maxim[where_X] = N_pops
            
            if not Consex:
                for h in range(len(maxim)):
                    maxim[h] = int(10*Guys[h,0])    
            
            #Similar = [sum(Guys[g,:] == np.amax(Guys,axis=1)[g]) for g in range(Guys.shape[0])]
            #Similar = [sum(Guys[g,:] == 0) for g in range(Guys.shape[0])]
            #Similar = [int(maxim[x] > 3 or sum(Guys[x,:] == 0) > 2 and maxim[x] != 3) + 1 for x in range(len(maxim))]
            Similar = [int(sum(Guys[x,:] <= 0) == Guys.shape[1] and maxim[x] != 3) + 1 for x in range(len(maxim))]
            
            peaks = [Points[x] for x in range(len(Points)) if Similar[x] == 1]
            where_peaks = [x for x in range(len(Points)) if Similar[x] == 1]
            
            d = 'none'
            if peaks:
                for l in range(len(Similar)):
                    if Similar[l] == 1:
                        if d == 'none' and l > 0:
                            maxim[:l] = [maxim[l] for x in range(l)]
                        d = l
                    if Similar[l] > 1:
                        if d != 'none':
                            #if max(Guys[l,:]) > 0:
                            #    Close = [x for x in range(Guys.shape[1]) if Guys[l,x] == max(Guys[l,:])]
                            #    maxim[l] = sum(Close) + 3
                            #else:
                            Distances = [abs(peaks[x] - Points[l]) for x in range(len(peaks))]
                            #print(maxim[l],maxim[where_peaks[Distances.index(min(Distances))]],Similar[l])
                            maxim[l] = maxim[where_peaks[Distances.index(min(Distances))]]                 
            
            ###
            ### PARK nber: 1
            ###    
            
            Topo.append(maxim + 1)
        
        
        Topo = np.array(Topo).T
        
        Clove = {CHR:{Points[x]:Topo[x,] for x in range(len(Points))}}
        
        Blocks_genome.update(Clove)
    
    return Blocks_genome


def read_3D_profiles(File_list):
    
    s0 = time.time()
    
    Blocks= recursively_default_dict()
#    Profiles= recursively_default_dict()
    Out = recursively_default_dict()
    Names= []        
    
    for File in File_list:
        
        ##### read blocks file
        Input= open(File,'r')
        d= 0
        
        for line in Input:
            line= line.split()
            
            if d== 0:
                line=line[4:]
                Names= line
            else:
                Blocks[int(line[0])][int(float(line[1]))][int(line[3])]= [float(x) for x in line[4:]]
                if line[3] == '0':
                    Out[int(line[0])][int(float(line[1]))] = int(float(line[2]))
            d += 1
        
        Input.close()
        #    
    s1= time.time()
    
    print(str(s1-s0) + ' seconds elapsed.')
    print(str((s1-s0) / float(60)) + ' minutes elapsed.')
    return Blocks, Names, Out



####
#### read block and profile files
####

print('To begin reading from: ')
print(args.books)

Ref_profiles, Names, Out = read_3D_profiles(args.books)


######
###### deciding who you're going to obe looking at;
###### in the future, replace with input file only, it should be easier.

if args.focus:
    Focus = read_focus(args.focus)
else:
    Focus = Names

focus_indexes= [x for x in range(len(Names)) if Names[x] in Focus]


Blocks = Merge_class(Ref_profiles,focus_indexes,Out,args.threshold,args.bin,args.outlier)

print("Number chromosomes selected: {0}".format(len(Blocks)))
####
####

if args.CHR:
    chromosomes= [args.CHR]
else:
    chromosomes= Blocks.keys()


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
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
import pandas as pd


# Here's the function that we'll call for each dataframe (once for chromosome
# ideograms, once for genes).  The rest of this script will be prepping data
# for input to this function
#
def chromosome_collections(df, y_positions, height,  **kwargs):
    """
    Yields BrokenBarHCollection of features that can be added to an Axes
    object.
    Parameters
    ----------
    df : pandas.DataFrame
        Must at least have columns ['chrom', 'start', 'end', 'color']. If no
        column 'width', it will be calculated from start/end.
    y_positions : dict
        Keys are chromosomes, values are y-value at which to anchor the
        BrokenBarHCollection
    height : float
        Height of each BrokenBarHCollection
    Additional kwargs are passed to BrokenBarHCollection
    """
    del_width = False
    if 'width' not in df.columns:
        del_width = True
        df['width'] = df['end'] - df['start']
    for chrom, group in df.groupby('chrom'):
        print(chrom)
        yrange = (y_positions[chrom], height)
        xranges = group[['start', 'width']].values
        yield BrokenBarHCollection(
            xranges, yrange, facecolors=group['colors'], **kwargs)
    if del_width:
        del df['width']


def compress_ideo(df,chromosome_list):
    
    new_set = []
    
    for CHR in range(len(chromosome_list)):
        
        Chr = int(re.search('chr(.+?)_',chromosome_list[CHR]).group(1))
        sub = df[df.chrom == chromosome_list[CHR]]
        Coordinates = sorted(sub.start)
        Size = sub.shape[0]
        start = min(df.start)
        First = sub.gieStain.iloc[0]
        for index in range(len(Coordinates)):
            row = sub[sub.start == Coordinates[index]]
            if index == 0:
                continue
            if index == (Size - 1):
                if row.gieStain.iloc[0] == First:
                    new_set.append([chromosome_list[CHR],start,Out[Chr][max(sub.start)],First])
                else:
                    new_set.append([chromosome_list[CHR],start,Out[Chr][max(sub.start)],First])
                    First = row.gieStain.iloc[0]
                    start = row.start.iloc[0]
                    new_set.append([chromosome_list[CHR],start,Out[Chr][max(sub.start)],First])
            else:
                if row.gieStain.iloc[0] == First:
                    continue
                else:
                    new_set.append([chromosome_list[CHR],start,row.start.iloc[0]-1,First])
                    First = row.gieStain.iloc[0]
                    start = row.start.iloc[0]
        
    new_set = pd.DataFrame(new_set,columns = ['chrom', 'start', 'end', 'gieStain'])
    return new_set




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

ideo = compress_ideo(ideo,chromosome_list)

#ideo = ideo[(ideo.start > 5.06e6) & (ideo.start < 7.06e6)]

ideo['colors'] = ideo['gieStain'].apply(lambda x: tuple([round(y / float(255),1) for y in color_lookup[x]]))
# Add a new column for width
ideo['width'] = ideo.end - ideo.start

# Width, height (in inches)
figsize = (args.width, args.height)

fig = plt.figure(figsize=figsize)
ax = fig.add_subplot(111)

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
ax.set_yticklabels(chromosome_list, fontsize = 5)
ax.axis('tight')

plt.savefig(Home + 'Ideo_' + Subject +'_CHR' + str(chromosomes[-1]).zfill(2)+'_st' + str(min(ideo.start)) + '_Z' +str(args.threshold)+ '_bin'+ str(args.bin)+'.png',bbox_inches = 'tight')

print('Done.')
