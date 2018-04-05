# -*- coding: utf-8 -*-
"""
Created on Wed Apr 04 21:10:47 2018

@author: jgarcia
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 15:30:52 2018

@author: jgarcia
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 14:36:11 2018

@author: jgarcia
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Dec 06 18:53:36 2017

@author: jgarcia
"""

from Kernel_tools import *
import collections
import time
import re

import os, sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("books",type=str,metavar= 'N',nargs= '+',
                    help = "Reference files to read. Any number can be given.")

parser.add_argument("--focus",type= str,help = "reference accessions indexes in genofile.")

parser.add_argument("--id",type= str,help = "Name your analysis")

parser.add_argument("--target",type= str,help = "target population")

parser.add_argument("--ref",type= str,help = "reference accessions indexes in genofile.")

parser.add_argument("--CHR",type= int,help = "chromosome to draw ideogram of.")

parser.add_argument("--start",type= int,help = "where to begin, in markers. Only makes sense if --CHR is also used.")

parser.add_argument("--end",type= int,help = "where to end, in markers. Only makes sense if --CHR is also used.")

parser.add_argument("--fam",help = "accession name file. same order as in geno.")

parser.add_argument("--plot",action= "store_true",help = "if given reduces cluster points using KDE extraction per cluster identified.")

parser.add_argument("--info",type= str,help = "optional information file on accessions analysed. requires columns: ID and label, with IDs and reference populations.")

parser.add_argument("--app",action= "store_true",help = "if given writes dash application and accompanying files.")

parser.add_argument("--reduc",action= "store_true",help = "if given prints cluster stats.")

parser.add_argument("--coarse",action='store_false',help= 'to smooth or not to smooth.')

parser.add_argument("--bin",default = 5,type= int,help = "smoothing parameter, must be uneven [savgol filter]")

parser.add_argument("--sg_order",default = 3,type= int,help = "staviksy golay filter order")

parser.add_argument("--outlier",type=float,default = 1e-4,help = "Outlier threshold")

parser.add_argument("--ms",type=float,default = .1,help = "cluster profile selection threshold")

parser.add_argument("--shared",type=float,default = .1,help = "Proportion of shared ancestry at given locus.")

parser.add_argument("--threshold",type = float,default = 2,help = "Intermediate classification threshold")

parser.add_argument("--chrom_height",type= float, default= 1, help= "height of ideograms")

parser.add_argument("--chrom_gap",type= float,default= 0,help= "gap between ideograms.")

parser.add_argument("--height",type= float, default= 10,help= "figure height, in inches.")

parser.add_argument("--width",type= float,default= 20,help= "figure width, in inches.")

parser.add_argument('--xticks',type= int,default= 1000000,help= 'xticks on final ideogram')


args = parser.parse_args()


def recursively_default_dict():
        return collections.defaultdict(recursively_default_dict)


Home= 'Analyses_' + args.id


FAM = "/gs7k1/home/jgarcia/olde_home/FILTER_GENO/F3G/NG_001.fam"
Ident = '/gs7k1/home/jgarcia/PO_Core/Rice_Identifiers.txt'
Geno_Q = "/gs7k1/home/jgarcia/olde_home/FILTER_GENO/SplitGENO/GENO_total.4.Q"
File_home= '/work/jgarcia/PO_lunch/COMP/'

########## Complementary files.

def read_focus(index_file):
    indxs = []
    
    Input = open(index_file,'r')
    for line in Input:
        line = line.split()
        indxs.append(line[0])
    
    Input.close()
    
    return indxs


def read_refs(index_file):
    indxs = recursively_default_dict()
    
    Input = open(index_file,'r')
    for line in Input:
        line = line.split()
        indxs[int(line[0])][int(line[1])] = []
    
    Input.close()
    
    indxs = {gop:[x for x in indxs[gop].keys()] for gop in indxs.keys()}
    
    return indxs, [x for x in indxs.keys()]


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
        
        
        for acc in focus_indicies:
            Guys = np.array([Likes[x][:,acc] for x in range(N_pops)])
            Guys = np.nan_to_num(Guys)
            Guys = [[[y,0][int(y<=X_threshold)] for y in x] for x in Guys]
            
            Test = [int(x < X_threshold) for x in np.amax(np.array(Guys),axis = 0)]
            
            if args.coarse:
                
                Guys = [savgol_filter(x,BIN,args.sg_order,mode = "nearest") for x in Guys]
                Test = savgol_filter(Test,BIN,args.sg_order,mode = "nearest")
                Test = [round(x) for x in Test]
            
            #
            Guys = np.array(Guys).T
            
            maxim = np.argmax(Guys,axis = 1)
            where_X = [x for x in range(Guys.shape[0]) if Test[x] == 1]
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
#                            Diff = int(len([x for x in Diff if x <= Diff_threshold]) == 1)
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
            
            ###
            ### PARK nber: 1
            ###
            
            Topo.append(maxim + 1)
        
        
        Topo = np.array(Topo).T
        
        Clove = {CHR:{Points[x]:Topo[x,] for x in range(len(Points))}}
        
        Blocks_genome.update(Clove)
    
    return Blocks_genome



def read_books(books):
    library= []
    
    for shelf in books:
        card= shelf.split('/')
        
        cover= card[-1].split('_')
        Chr= int([re.findall(r'\d+',i)[0] for i in cover if re.search('CHR',i)][0])
        tag= cover[1]
        start= cover[2]
        library.append([shelf,tag,start,Chr])
    
    library= pd.DataFrame(library,columns= ['file','tag','start','Chr'])
    return library
        
    


def read_3D_profiles(Books):
    
    s0 = time.time()
    
    Blocks= recursively_default_dict()
    Profiles= recursively_default_dict()
    Out = recursively_default_dict()
    Names= []
    
    
    for Chr in Books.Chr.unique():
        
        Shelf= Books[(Books.Chr== Chr)]
        
        Blocks[Chr]= recursively_default_dict()
        Profiles[Chr]= recursively_default_dict()
        
        start= '2'
        
        for start in Shelf.start.unique():
            
            Series= Shelf[(Shelf.start== start)]
            
            if len(Series) != 2:
                print('{0} files were provided for analysis Chr: {1} - start: {2}. 2 files are expected.'.format(len(Series),Chr,start))
                break
            
            blocks_file= Series.loc[Series['tag'] == 'Request','file']
            blocks_file= blocks_file.iloc[0]
            
            profiles_file = Series.loc[Series['tag'] == 'profiles','file']
            profiles_file= profiles_file.iloc[0]
            ##### read blocks file
            Input= open(blocks_file,'r')
            d= 0
            
            for line in Input:
                line= line.split()
                
                if d== 0:
                    line=line[4:]
                    Names= line
                else:
                    Blocks[int(line[0])][int(line[1])][int(line[3])]= [float(x) for x in line[4:]]
                    if line[3] == '0':
                        Out[int(line[0])][int(line[1])] = int(line[2])
                d += 1
            
            Input.close()
            
            ##### read profiles file
            
            Input= open(profiles_file,'r')
            d= 0
            
            for line in Input:
                line= line.split()
                
                if d== 0:
                    line=line[3:]
                    Names= line
                else:
                    Profiles[int(line[0])][int(line[1])][int(line[2])]= [float(x) for x in line[3:]]
                
                d += 1
            
            Input.close()
    
    s1= time.time()
    
    print(str(s1-s0) + ' seconds elapsed.')
    print(str((s1-s0) / float(60)) + ' minutes elapsed.')
    return Blocks, Profiles, Names, Out




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
        
        yrange = (y_positions[chrom], height)
        xranges = group[['start', 'width']].values
        yield BrokenBarHCollection(
            xranges, yrange, facecolors=group['colors'], **kwargs)
    if del_width:
        del df['width']


def compress_ideo(df,chromosome_list):
    
    new_set = []
    
    for CHR in range(len(chromosome_list)):
        print(chromosome_list[CHR])
        Chr = int(re.search('chr(.+?)_',chromosome_list[CHR]).group(1))
        sub = df[df.chrom == chromosome_list[CHR]]
        Coordinates = sorted(sub.start)
        Size = sub.shape[0]
        start = min(sub.start)
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




def plot_ideo(Blocks,CHR,Focus,label,ideo_height,ideo_spacing,height,width,Home,ID):
    chromosomes = CHR
    
    chromosome_list = []
    
    Ideo = []
    
    for here in range(len(Focus)):
        Subject = Focus[here]
        
        chromosome_list.extend(['chr'+str(Chr)+ '_' + Subject for Chr in chromosomes])
        
        color_ref= ['red','yellow','blue','black','orange','purple','green','silver','red3','deepskyeblue','navy','chartreuse','darkorchid3','goldenrod2']
        color_ref= ['white','red']
        
        Stock = [[['chr'+str(Chr)+ '_' + Subject,bl,Out[Chr][bl],color_ref[Blocks[Chr][bl][here]]] for bl in sorted(Blocks[Chr].keys())] for Chr in chromosomes]
        Stock = [y for y in it.chain(*[z for z in it.chain(*[Stock])])]
        
        Ideo.extend(Stock)
    
    
    
    # Height of each ideogram
    chrom_height = ideo_height
    
    # Spacing between consecutive ideograms
    chrom_spacing = ideo_spacing
    
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
        'white': [255,255,255]
    }
    
    # Add a new column for colors
    
    ideo = compress_ideo(ideo,chromosome_list)
    
    #ideo = ideo[(ideo.start > 5.06e6) & (ideo.start < 7.06e6)]
    
    ideo['colors'] = ideo['gieStain'].apply(lambda x: tuple([round(y / float(255),1) for y in color_lookup[x]]))
    # Add a new column for width
    ideo['width'] = ideo.end - ideo.start
    
    # Width, height (in inches)
    figsize = (width, height)
    
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    
    # Now all we have to do is call our function for the ideogram data...
    print("adding ideograms...")
    for collection in chromosome_collections(ideo, chrom_ybase, chrom_height, edgecolors=None, linewidths= 0):
        ax.add_collection(collection)
    
    
    # Axes tweaking
    ax.set_xticks([x for x in range(min(ideo.start),max(ideo.end),int(args.xticks))])
    plt.xticks(fontsize = 10,rotation = 90)
    ax.tick_params(axis = 'x',pad = 10)
    
    ax.tick_params(axis='y', which='major', pad=30)
    ax.set_yticks([chrom_centers[i] for i in chromosome_list])
    ax.set_yticklabels(chromosome_list, fontsize = 5)
    ax.axis('tight')
    filename= Home+ '/Ideo_id.'+ ID +'_label.' + str(label) +'_CHR' + str(chromosomes[-1]).zfill(2)+'.png'
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    plt.savefig(filename,bbox_inches = 'tight')



def plot_ideo2(Blocks,CHR,Focus,label,ideo_height,ideo_spacing,height,width,Home,ID):
    chromosomes = CHR
    
    chromosome_list = []
    all_blocks= [x for x in it.chain(*[Blocks[y].keys() for y in chromosomes])]
    
    Ideo = []
    
    for here in range(len(Focus)):
        Subject = Focus[here]
        
        chromosome_list.extend(['chr'+str(Chr)+ '_' + Subject for Chr in chromosomes])
        
        color_ref= ['red','yellow','blue','black','orange','purple','green','silver','red3','deepskyeblue','navy','chartreuse','darkorchid3','goldenrod2']
        color_ref= ['white','red']
        
        Stock = [[['chr'+str(Chr)+ '_' + Subject,bl,Out[Chr][bl],color_ref[Blocks[Chr][bl][here]]] for bl in sorted(Blocks[Chr].keys()) \
        if Blocks[Chr][bl][here] != 0] for Chr in chromosomes]
        Stock = [y for y in it.chain(*[z for z in it.chain(*[Stock])])]
        ## very lazy way of making sure matplotlib plots everyting, including whites.It hasn't been doing so despite the ax.set_xlim().. dunno why.
        
        Ideo.extend(Stock)
    
    
    
    # Height of each ideogram
    chrom_height = ideo_height
    
    # Spacing between consecutive ideograms
    chrom_spacing = ideo_spacing
    
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
        'white': [255,255,255]
    }
    
    # Add a new column for colors
    
    #ideo = compress_ideo(ideo,chromosome_list)
    
    #ideo = ideo[(ideo.start > 5.06e6) & (ideo.start < 7.06e6)]
    
    ideo['colors'] = ideo['gieStain'].apply(lambda x: tuple([round(y / float(255),1) for y in color_lookup[x]]))
    # Add a new column for width
    ideo['width'] = ideo.end - ideo.start
    
    # Width, height (in inches)
    figsize = (width, height)
    
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    
    # Now all we have to do is call our function for the ideogram data...
    print("adding ideograms...")
    for collection in chromosome_collections(ideo, chrom_ybase, chrom_height, edgecolors=None, linewidths= 0):
        ax.add_collection(collection)
    
    
    # Axes tweaking
    
    ax.set_xlim(min(all_blocks),max(all_blocks))
    ax.set_xticks([x for x in range(min(all_blocks),max(all_blocks),int(args.xticks))])
    plt.xticks(fontsize = 10,rotation = 90)
    ax.tick_params(axis = 'x',pad = 10)
    
    ax.tick_params(axis='y', which='major', pad=30)
    ax.set_yticks([chrom_centers[i] for i in chromosome_list])
    ax.set_yticklabels(chromosome_list, fontsize = 5)
    
    filename= Home+ '/Ideo_id.'+ ID +'_label.' + str(label) +'_CHR' + str(chromosomes[-1]).zfill(2)+'.png'
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    plt.savefig(filename,bbox_inches = 'tight')
    plt.close(fig)




print('to begin reading from: ')
print(args.books)

Library= read_books(args.books)

print('library:')
print(Library.sort_values(by= ['Chr','start']))


Ref_profiles, Profiles, Names, Out = read_3D_profiles(Library)


######
###### deciding who you're going to obe looking at;
###### in the future, replace with input file only it should be easier.
Fam = FAMread(args.fam)


####
#### some parameters
####

Diff_threshold = args.threshold
X_threshold= args.outlier

refs_lib, Parents = read_refs(args.ref)


if args.focus:
    Focus = read_focus(args.focus)
else:
    Focus = Names

focus_indexes= [x for x in range(len(Focus))]

Absent= [x for x in Focus if x not in Names]

if len(Absent) > 0:
    print('The following individuals were not found in the files provided: {}'.format(Absent))
    print('Analysis will proceed without them.')
    Focus= [x for x in Focus if x in Names]

compound_reference = [Names.index(x) for x in Focus]

Blocks = Merge_class(Ref_profiles,compound_reference,Out,args.threshold,args.bin,args.outlier)

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


####

threshold= args.shared
MS_threshold= args.ms


def codes(Parents):
    Store= recursively_default_dict()
    
    # Parent labels
    for x in range(len(Parents)):
        Store[x + 1]= [Parents[x]]
    
    ## outlier label
    Store[len(Parents) + 1]= [x for x in Parents]
    
    ## Combination labels:
    combine= [x for x in it.combinations(Parents,2)]
    for i in range(len(combine)):
        Store[len(Parents) + 2 + i]= [x for x in combine[i]]
    
    if len(Parents) >=3:
        ## three way label    
        Store[len(Store) + 1]= [x for x in Parents]
    
    return Store


code_back= codes(Parents)


if args.target:
    target= [int(x) for x in args.target.split(',')]
else:
    target= [x for x in code_back.keys()]

print('target: {0}'.format(target))

##### Selecting windows to focus on:

#### Chose_profiles: automatically chose clusters with at least one included 
#### focus accession of 'target' color at a threshold >= target_threshold
Chose_profiles = {CHR:{bl:[y for y in Profiles[CHR][bl].keys() if sum([int(Profiles[CHR][bl][y][z] >= MS_threshold) \
for z in [compound_reference[x] for x in focus_indexes if Blocks[CHR][bl][x] in target]]) >= 1] \
for bl in Blocks[CHR].keys() if \
len([x for x in focus_indexes if Blocks[CHR][bl][x] in target]) / float(len(Focus)) >= threshold} \
for CHR in Blocks.keys() if CHR in chromosomes}

Coordinates = [[[[CHR,bl,Out[CHR][bl],x] for x in Chose_profiles[CHR][bl]] for bl in sorted(Chose_profiles[CHR].keys())] for CHR in sorted(Chose_profiles.keys())]
Coordinates = [z for z in it.chain(*[y for y in it.chain([x for x in it.chain(*Coordinates)])])]
Coordinates= np.array(Coordinates)

Clover= [[[Profiles[CHR][bl][x] for x in Chose_profiles[CHR][bl]] for bl in sorted(Chose_profiles[CHR].keys())] for CHR in sorted(Chose_profiles.keys())]
Clover= [z for z in it.chain(*[y for y in it.chain(*Clover)])]
Clover= np.array(Clover)

####
#### Pre processing and dimensionality reduction of matrix
#### of selected clusters.
####


from sklearn import preprocessing

Clover = np.nan_to_num(Clover)
preProc_Clover = Clover

print('Clover shape: ', Clover.shape)

Clover = preprocessing.scale(Clover,axis = 1)
#

print("Clover shape: ", Clover.shape)


reefer= [g for g in it.chain(*[refs_lib[y] for y in sorted(refs_lib.keys())])]

if not args.focus:
    mary= [Names.index(x) for x in Names if Fam[x] not in reefer]
else:
    mary= [Names.index(x) for x in Names if Fam[x] not in reefer and x in Focus]

reefer= [Names.index(Fam[x]) for x in reefer]

Subset= [x for x in it.chain(*[mary,reefer])]
Trend= np.repeat([0,1],[len(mary),len(reefer)])


## apply pca to reference accessions, transform the rest.
variation_focus= [Names.index(Fam[x]) for x in it.chain(*[refs_lib[z] for z in list(set(it.chain(*[code_back[y] for y in target])))])]

### PCA
pca = PCA(n_components=5, whiten=False).fit(Clover[:,variation_focus].T)
X_se = pca.transform(Clover[:,Subset].T)
COMPS = pca.components_.T*np.sqrt(pca.explained_variance_)


###############################################################################
########################### PAINTING SHIT!! ###################################
###############################################################################

## 
## CLUSTER EIGENVALUES
##

### Clustering on decomposition
    
#    bandwidth = estimate_bandwidth(COMPS, quantile=0.1)
#    if bandwidth==0:
#        bandwidth = 0.1
#    
#    ms = MeanShift(bandwidth=bandwidth, bin_seeding=False, cluster_all=True, min_bin_freq=1)
#    ms.fit(COMPS)
#    labels1 = ms.labels_
#    label_select = {y:[x for x in range(len(labels1)) if labels1[x] == y] for y in sorted(list(set(labels1)))}
#    
### HDBSCAN
#    from sklearn.cluster import DBSCAN
#    db = DBSCAN(min_samples=35).fit(COMPS)
#    labels1= db.labels_
#    label_select = {y:[x for x in range(len(labels1)) if labels1[x] == y] for y in sorted(list(set(labels1))) if y != -1}

#    ### Ward

#    from sklearn.cluster import AgglomerativeClustering
#    clustering = AgglomerativeClustering(linkage='ward', n_clusters=10)
#    clustering.fit(COMPS)
#    labels1= clustering.labels_
#    label_select = {y:[x for x in range(len(labels1)) if labels1[x] == y] for y in sorted(list(set(labels1)))}
#    
#    #
#    ### K-means
from sklearn.cluster import KMeans
kmeans = KMeans(n_clusters=10, random_state=0).fit(COMPS)
labels1 = kmeans.labels_
label_select = {y:[x for x in range(len(labels1)) if labels1[x] == y] for y in sorted(list(set(labels1)))}



############################################################################
############################################################################
############ Bring out ideograms ##########################################

if args.plot == True:
    Ideo_home= Home + '/Ideos'
    
    Blancs= {aim:{Chr:{bl:[0]*len(Focus) for bl in Blocks[Chr].keys()} for Chr in Blocks.keys()} for aim in list(set(labels1))}
    
    target_block= {Chr:{bl:[int(x in target) for x in Blocks[Chr][bl]] for bl in Blocks[Chr].keys()} for Chr in chromosomes}
    
    
    for Chr in chromosomes:
        plot_ideo2(target_block,[Chr],Focus,'all',args.chrom_height,args.chrom_gap,args.height,args.width,Ideo_home,args.id)
    
    for n in range(len(Coordinates)):
        site= Coordinates[n]
        trigger= [x for x in range(len(Focus)) if Profiles[site[0]][site[1]][site[3]][compound_reference[x]] >= MS_threshold and Blocks[site[0]][site[1]][x] in target]
        
        for v in trigger:
            Blancs[labels1[n]][site[0]][site[1]][v] = 1
    
    for aim in Blancs.keys():
        if len(Focus)== 1:
            plot_ideo(Blancs[aim],Blocks.keys(),Focus,aim,args.chrom_height,args.chrom_gap,args.height,args.width,Ideo_home,args.id)
        else:
            for Chr in chromosomes:
                plot_ideo2(Blancs[aim],[Chr],Focus,aim,args.chrom_height,args.chrom_gap,args.height,args.width,Ideo_home,args.id)


###
### DELETE Ref_profiles AND Profiles (trying to save space)
###
Ref_profiles= {}
Profiles= {}


###############################################################################
#### Average normalized likelihhod among clustered eigenvectors by haplotype #####
###############################################################################


Cameo = []

for cramp in sorted(label_select.keys()):
    Clamp = np.mean(preProc_Clover[label_select[cramp],:],axis = 0)
    Fry = [Clamp[x] for x in Subset]
    Cameo.append(Fry)

Cameo = np.array(Cameo).T


###########################################################################
### cosine of the clustered eigenvectors with haplotype coordinates ######## DEPRECATED
###########################################################################

#cos_threshold = .6
#
#
#from numpy import dot
#from numpy.linalg import norm
#
#SpaceY = recursively_default_dict()
#
#for g in label_select.keys():
#    Green = COMPS[label_select[g],:]
#    SpaceY[g] = [mean(dot(X_se[x],Green.T)/norm(Green,axis=1)/norm(X_se[x])) for x in range(X_se.shape[0])]
#
#Globe = np.array([SpaceY[x] for x in sorted(SpaceY.keys())]).T
#
#

######## Reducing the number of cluster profiles to print:


new_labs= labels1

if args.reduc:
    new_labs= []
    reduced_comp= []
    Size= 2000
    
    params = {'bandwidth': np.linspace(np.min(COMPS), np.max(COMPS),30)}
    grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0)
    
    for lab in label_select.keys():
        N_prop= round(len(label_select[lab]) * Size / float(sum([len(x) for x in label_select.values()])))
        
        if len(label_select[lab]) <= 3:
            reduced_comp.extend(COMPS[label_select[lab],:])
            new_labs.extend([lab]*N_prop)
            continue
        
        grid.fit(COMPS[label_select[lab],:])    
        
        kde = grid.best_estimator_
        
        new_data = kde.sample(N_prop, random_state=0)
        
        reduced_comp.extend(new_data)
        new_labs.extend([lab]*N_prop)
    
    
    COMPS= np.array(reduced_comp)



### print


CHR = [x for x in chromosomes][-1]
start= args.id

### Writting files for dash application;

if args.app:
    try:
        import cPickle as pickle
    except ImportError:  # python 3.x
        import pickle
    
    with open('app_files.p', 'rb') as fp:
        data = pickle.load(fp)
    
    elements= {
    'ID=': start,
    'Where=': CHR,
    'ref=': target[0],
    'pop_refs=': [x for x in it.chain(*[['unlabelled'],[str(x) for x in Parents]])]
    }
    
    for x in range(len(data['app.py'])):
        line= data['app.py'][x]
        for el in elements.keys():
            if re.search(el,line):
                line= line.split('=')
                if el == 'pop_refs=':
                    new_line= '{}= {}'.format(line[0],elements[el]) + '\n'
                else:
                    new_line= '{}= "{}"'.format(line[0],elements[el]) + '\n'
                data['app.py'][x]= new_line
                line= new_line
    
    for file_name in data:
        filename= ''.join([Home,'/',file_name])
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        Output= open(filename,'w')
        for line in data[file_name]:
            Output.write(line)
        Output.close()
        
    if args.info:
        ordercore= list(open(args.info,'r'))
        filename= ''.join([Home,'/Order_core.txt'])
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        Output= open(filename,'w')
        for line in ordercore:
            Output.write(line)
        Output.close()
    
    else:
        crafty= []
        ## set up an accession to group map.
        reference= recursively_default_dict()
        for gp in refs_lib.keys():
            for acc in refs_lib[gp]:
                reference[Fam[acc]]= Parents.index(gp) + 1
        for acc in Names:
            if acc not in reference.keys():
                reference[acc]= 0
                
        for race in Names:
            crafty.append([race,elements['pop_refs='][reference[race]],reference[race]])
        crafty= pd.DataFrame(crafty,columns= ['ID','label','code'])
        filename= ''.join([Home,'/Order_core.txt'])
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        crafty.to_csv(filename,sep= '\t')



filename= Home + "/Profile_" + str(target[0]) + "_CHR" + str(CHR) + "."+str(start) + ".txt"
os.makedirs(os.path.dirname(filename), exist_ok=True)

Output = open(filename,"w")

for title in sorted(label_select.keys()):
    Output.write("G" + str(title) + "\t")

Output.write("\n")

for Future in range(Cameo.shape[0]):
    for Machine in range(Cameo.shape[1]):
        Output.write(str(Cameo[Future,Machine]) + "\t")
    Output.write("\n")

Output.close()




filename= Home + "/Profile_coordinates_" + str(target[0]) + "_CHR" + str(CHR) + "."+str(start) + ".txt"
os.makedirs(os.path.dirname(filename), exist_ok=True)

Output = open(filename,"w")

Output.write('\t'.join(['chrom','start','end','cluster','label']))
Output.write('\n')

for axe in range(Coordinates.shape[0]):
    Output.write('\t'.join([str(x) for x in Coordinates[axe,:]]) + '\t')
    Output.write(str(labels1[axe]) + '\n')

Output.close()



filename= Home + "/DIM_private_"+str(target[0])+"_request_CHR" + str(CHR) + "."+str(start)+".txt"
os.makedirs(os.path.dirname(filename), exist_ok=True)

Output = open(filename,"w")

for entry in range(X_se.shape[0]):
    Output.write(str(Trend[entry]) + "\t")
    
    Output.write(Names[Subset[entry]] + "\t")
    for scythe in range(X_se.shape[1]):
        Output.write(str(X_se[entry,scythe]) + "\t")
    Output.write("\n")

Output.close()



filename= Home + "/DIM_private_"+str(target[0])+"_comp_CHR" + str(CHR) + "."+str(start)+".txt"
os.makedirs(os.path.dirname(filename), exist_ok=True)

Output = open(filename,"w")
Output.write('0\t' + '\t'.join([str(x) for x in pca.explained_variance_ratio_]) + '\t\n')

for entry in range(COMPS.shape[0]):
    Output.write(str(new_labs[entry]) + "\t")
    for scythe in range(COMPS.shape[1]):
        Output.write(str(COMPS[entry,scythe]) + "\t")
    Output.write("\n")

Output.close()

print('Done.')