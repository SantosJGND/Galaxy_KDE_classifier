import collections
import time
import itertools as it
import numpy as np
import re
import pandas as pd

import collections

def recursively_default_dict():
        return collections.defaultdict(recursively_default_dict)

from matplotlib.collections import BrokenBarHCollection

## created 28-10-2018
## SantosJGND

##########################
## Load data

def read_3D_profiles(File_list):
    '''
    Reads Blocks files of reference KDE derived p-values across genomic windows.
    Returns dictionary: {CHR: {Windows: {refs: [p-value list]}}}
    '''
    s0 = time.time()
    
    Blocks= recursively_default_dict()
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



#######################
### KDE-based classification across windows. 


def Merge_class(Ref_profiles,focus_indicies,Out,Diff_threshold,BIN,X_threshold,coarse,sg_order):
    '''
    Receives dictionary of reference p-values by individual by window and chromosome.
    Return dictionary of individual classifications according to the parameters provided.
    '''
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
            Guys = [[[y,0][int(y<=X_threshold)] for y in x] for x in Guys]
            
            Test = [int(x <= X_threshold) for x in np.amax(np.array(Guys),axis = 0)]     
            
            if coarse:
                
                Guys = [savgol_filter(x,BIN,sg_order,mode = "nearest") for x in Guys]
                Test = savgol_filter(Test,BIN,sg_order,mode = "nearest")
                Test = [round(x) for x in Test]
            
            #
            Guys = np.array(Guys).T
            
            maxim = np.argmax(Guys,axis = 1)
            where_X = [x for x in range(Guys.shape[0]) if Test[x] == 1]
            
            #
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
                            Diff = abs(max(Diff)) / abs(min(Diff))
                            Diff = int(Diff > Diff_threshold)
                        
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
            
            
            Topo.append(maxim + 1)
            
        Topo = np.array(Topo).T
        
        Clove = {CHR:{Points[x]:Topo[x,] for x in range(len(Points))}}
        
        Blocks_genome.update(Clove)
    
    return Blocks_genome


#############################
#### Ideogram construction
###

### 
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


def compress_ideo(df,chromosome_list, Out):
    '''
    Merge neighboring windows of the same class individual-wise. Returns pandas df.
    '''
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
                    new_set.append([chromosome_list[CHR],start,Out[Chr][max(df.start)],First])
                else:
                    new_set.append([chromosome_list[CHR],start,Out[Chr][max(df.start)],First])
                    First = row.gieStain.iloc[0]
                    start = row.start.iloc[0]
                    new_set.append([chromosome_list[CHR],start,Out[Chr][max(df.start)],First])
            else:
                if row.gieStain.iloc[0] == First:
                    continue
                else:
                    new_set.append([chromosome_list[CHR],start,row.start.iloc[0]-1,First])
                    First = row.gieStain.iloc[0]
                    start = row.start.iloc[0]
    
    new_set = pd.DataFrame(new_set,columns = ['chrom', 'start', 'end', 'gieStain'])
    return new_set



################################################################################
################################################################################
#### Deprecated.
################################################################################

#####
##### Clustering accessions based on classification
## to use direct similarity (proportion of differences).
from scipy.spatial.distance import pdist, squareform

##### functions retrieved in https://gmarti.gitlab.io/ml/2017/09/07/how-to-sort-distance-matrix.html
### used here for the order only.
def seriation(Z,N,cur_index):
    '''
        input:
            - Z is a hierarchical tree (dendrogram)
            - N is the number of points given to the clustering process
            - cur_index is the position in the tree for the recursive traversal
        output:
            - order implied by the hierarchical tree Z
            
        seriation computes the order implied by a hierarchical tree (dendrogram)
    '''
    
    if cur_index < N:
        return [cur_index]
    else:
        left = int(Z[cur_index-N,0])
        right = int(Z[cur_index-N,1])
        return (seriation(Z,N,left) + seriation(Z,N,right))

def compute_serial_matrix(dist_mat,method="ward"):
    '''
        input:
            - dist_mat is a distance matrix
            - method = ["ward","single","average","complete"]
        output:
            - seriated_dist is the input dist_mat,
              but with re-ordered rows and columns
              according to the seriation, i.e. the
              order implied by the hierarchical tree
            - res_order is the order implied by
              the hierarhical tree
            - res_linkage is the hierarhical tree (dendrogram)
        
        compute_serial_matrix transforms a distance matrix into 
        a sorted distance matrix according to the order implied 
        by the hierarchical tree (dendrogram)
    '''
    
    N = len(dist_mat)
    flat_dist_mat = squareform(dist_mat)
    res_linkage = linkage(flat_dist_mat, method=method,preserve_input=True)
    res_order = seriation(res_linkage, N, N + N-2)
    seriated_dist = np.zeros((N,N))
    a,b = np.triu_indices(N,k=1)
    seriated_dist[a,b] = dist_mat[ [res_order[i] for i in a], [res_order[j] for j in b]]
    seriated_dist[b,a] = seriated_dist[a,b]
    
    return seriated_dist, res_order, res_linkage

###
### Deprecated final order based on classification clustering.
### currently returning ideograms in order provided in Focus file.
###
'''
Blocks_matrix= [[Blocks[c][w] for w in Blocks[c].keys()] for c in Blocks.keys()]
Blocks_matrix= np.array(Blocks_matrix).T

def Calc_diff(Mat,indexes):
    
    Dist_dict= recursively_default_dict()
    
    for clev in indexes:
        dude= Mat[clev,:]
        
        for other in indexes:
            dist= [int(Mat[clev,x] == Mat[other,x]) for x in range(Mat.shape[1])]
            dist= sum(dist) / float(len(dist))
            Dist_dict[clev][other]= dist
    
    return Dist_dict

Blocks_distances= Calc_diff(Blocks_matrix,focus_indexes)

Blocks_diff= np.array([Blocks_distances[guy] for guy in Blocks_distances.keys()])
Blocks_diff= squareform(Blocks_diff)

ordered_dist_mat, res_order, res_linkage = compute_serial_matrix(Blocks_distances,args.clust)

'''
