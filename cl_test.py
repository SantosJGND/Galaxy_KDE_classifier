# -*- coding: utf-8 -*-
"""
Created on Wed Feb 07 17:45:43 2018

@author: jgarcia
"""

from Kernel_tools import *
import collections
import time
import itertools as it
import numpy as np
import re

import sys
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("--fam",type=str,help = "Fam file, ind id's, same as for kernel analysis")

#args = parser.parse_args()
print(sys.argv)

flag = -1
X_flag = -1
D_flag = -1
bin_flag= -1

if '--fam' not in sys.argv:
    flag = -1
    print('Missing fam file')
    sys.exit()
    
else:
    flag = sys.argv.index('--fam')
    FAM = sys.argv[flag + 1]




def recursively_default_dict():
        return collections.defaultdict(recursively_default_dict)


Home= "/work/jgarcia/PO_lunch/COMP/"
FAM = "/gs7k1/home/jgarcia/olde_home/FILTER_GENO/F3G/NG_001.fam"
Ident = '/gs7k1/home/jgarcia/PO_Core/Rice_Identifiers.txt'
Geno_Q = "/gs7k1/home/jgarcia/olde_home/FILTER_GENO/SplitGENO/GENO_total.4.Q"

Parents = 3

def Merge_class(Ref_profiles,focus_indicies,Out,Diff_threshold,BIN,X_threshold):
    Blocks_genome = recursively_default_dict()
    
    for CHR in Ref_profiles.keys():
        print(CHR)
        Points = sorted(Out[CHR].keys())
        Likes = Ref_profiles[CHR]
        Likes = {x:np.array(Likes[x]) for x in Likes.keys()}
        
        Topo = []
        
        #range_Parents = [x + Aro.shape[0] for x in range(Daddy.shape[0])]
        #range_Crossed = [x for x in range(Aro.shape[0])]
        
        for acc in focus_indicies:
            Guys = np.array([Likes[x][:,acc] for x in range(Parents)])
            Guys = np.nan_to_num(Guys)
            
            
            Test = [int(x < X_threshold) for x in np.amax(np.array(Guys),axis = 0)]
            Test = savgol_filter(Test,BIN,3,mode = "nearest")
            Test = [round(x) for x in Test]
            
            Guys = [[[y,0][int(y<=X_threshold)] for y in x] for x in Guys]
            Guys = [savgol_filter(x,BIN,3,mode = "nearest") for x in Guys]
            #    
            Guys = np.array(Guys).T
            
            maxim = np.argmax(Guys,axis = 1)
            where_X = [x for x in range(Guys.shape[0]) if Test[x] == 1]
            #where_X = [x for x in range(Guys.shape[0]) if len([c for c in Guys[x,:3] if c <= .0001]) == 3]
            Consex = [x for x in it.combinations(range(Parents),2)]
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
                        maxim[h] = sum(CL[0]) + Parents
                maxim[where_X] = Parents
            
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


def read_3D_profiles(Home):
    
    s0 = time.time()
    
    Blocks= recursively_default_dict()
#    Profiles= recursively_default_dict()
    Out = recursively_default_dict()
    Names= []
    
    request= 'Request'
        
    
    for Chr in range(1,13):
        Blocks[Chr]= {x:[] for x in range(Parents)}
#        Profiles[Chr]= recursively_default_dict()
        
        start= '1'
        
        blocks_file= Home + 'Blocks_'+request+'_st' + start + '_CHR' + str(Chr).zfill(2) + '.txt'
#        profiles_file = Home + 'Blocks_profiles_st'+start+'_CHR' + str(Chr).zfill(2) + '.txt'
            
        ##### read blocks file
        Input= open(blocks_file,'r')
        d= 0
        
        for line in Input:
            line= line.split()
            
            if d== 0:
                line=line[4:]
                Names= line
            else:
                Blocks[int(line[0])][int(line[3])].append([float(x) for x in line[4:]])
                if line[3] == '0':
                    Out[int(line[0])][int(line[1])] = int(line[2])
            d += 1
        
        Input.close()
        
#        ##### read profiles file
#        
#        Input= open(profiles_file,'r')
#        d= 0
#        
#        for line in Input:
#            line= line.split()
#            
#            if d== 0:
#                line=line[3:]
#                Names= line
#            else:
#                Profiles[int(line[0])][int(line[1])][int(line[2])]= [float(x) for x in line[3:]]
#            
#            d += 1
#        
#        Input.close()
#    
    s1= time.time()
    
    print(str(s1-s0) + ' seconds elapsed.')
    print(str((s1-s0) / float(60)) + ' minutes elapsed.')
    return Blocks, Names, Out



####
#### read block and profile files
####


Ref_profiles, Names, Out = read_3D_profiles(Home)


######
###### deciding who you're going to obe looking at;
###### in the future, replace with input file only it should be easier.
Fam = FAMread(FAM)


Geneal = Geneo_CORE(Fam,Ident)
Geneo = OriginbySNMF(Geno_Q,.7)
Excluded_varieties = ["IRIS_313-10032","IRIS_313-11765","IRIS_313-9978"]
Aro = [x for x in Geneo[2] if Fam[x] not in Excluded_varieties and x not in Geneal]
Geneo = OriginbySNMF(Geno_Q,.8)
Geneo = {x:[y for y in Geneo[x] if y in Geneal] for x in Geneo.keys()}
Geneo[2].extend(Aro)


####
#### some parameters
####
Diff_threshold = 2
X_threshold= 5e-3
BIN = 15

#Focus = ['CX59', 'CX65', 'CX67', 'CX104', 'CX143', 'CX149', 'IRIS_313-8268', 'IRIS_313-8326', 'IRIS_313-8385', 'IRIS_313-8656', 'IRIS_313-8712', 'IRIS_313-8747', 'IRIS_313-8813', 'IRIS_313-9083', 'IRIS_313-9172', 'IRIS_313-9601', 'IRIS_313-9629', 'IRIS_313-10670', 'IRIS_313-10851', 'IRIS_313-10868', 'IRIS_313-10926', 'IRIS_313-10933', 'IRIS_313-11021', 'IRIS_313-11022', 'IRIS_313-11023', 'IRIS_313-11026', 'IRIS_313-11218', 'IRIS_313-11258', 'IRIS_313-11268', 'IRIS_313-11289', 'IRIS_313-11293', 'IRIS_313-11451', 'IRIS_313-11564', 'IRIS_313-11567', 'IRIS_313-11625', 'IRIS_313-11627', 'IRIS_313-11630', 'IRIS_313-11632', 'IRIS_313-11743', 'IRIS_313-11825']
Focus = [Fam[x] for x in Geneo[2]]
#Focus = ['CX110']
Focus = Names
focus_indexes= [x for x in range(len(Names)) if Names[x] in Focus]

#Focus = ['IRIS_313-9800','IRIS_313-10073','IRIS_313-10074','IRIS_313-10078','IRIS_313-10080','IRIS_313-10094','IRIS_313-10067']

Blocks = Merge_class(Ref_profiles,focus_indexes,Out,Diff_threshold,BIN,X_threshold)


####
#

chromosomes = [7]

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
        Coordinates = sorted(Out[Chr].keys())
        sub = df[df.chrom == chromosome_list[CHR]]
        Size = sub.shape[0]
        start = 0
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
chrom_height = 1

# Spacing between consecutive ideograms
chrom_spacing = .05

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

ideo = ideo[(ideo.start > 5.06e6) & (ideo.start < 7.06e6)]

ideo['colors'] = ideo['gieStain'].apply(lambda x: tuple([round(y / float(255),1) for y in color_lookup[x]]))
# Add a new column for width
ideo['width'] = ideo.end - ideo.start

# Width, height (in inches)
figsize = (5, 30)

fig = plt.figure(figsize=figsize)
ax = fig.add_subplot(111)

# Now all we have to do is call our function for the ideogram data...
print("adding ideograms...")
for collection in chromosome_collections(ideo, chrom_ybase, chrom_height, edgecolors=None, linewidths= 0):
    ax.add_collection(collection)


# Axes tweaking
ax.set_xticks([x for x in range(0,max(ideo.end),2000000)])
plt.xticks(fontsize = 10,rotation = 90)
ax.tick_params(axis = 'x',pad = 10)

ax.tick_params(axis='y', which='major', pad=30)
ax.set_yticks([chrom_centers[i] for i in chromosome_list])
ax.set_yticklabels(chromosome_list, fontsize = 5)
ax.axis('tight')

plt.savefig('Ideo_' + Subject + '.png',bbox_inches = 'tight')


