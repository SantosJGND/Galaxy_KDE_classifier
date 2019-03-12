import numpy as np
import pandas as pd
import itertools as it

import scipy

from sklearn.neighbors import KernelDensity
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV
from sklearn.cluster import estimate_bandwidth
from sklearn.cluster import MeanShift, estimate_bandwidth

from sklearn.metrics.pairwise import pairwise_distances
from sklearn.metrics.pairwise import euclidean_distances

import re
import matplotlib.pyplot as plt

from matplotlib.collections import BrokenBarHCollection

import collections

def recursively_default_dict():
        return collections.defaultdict(recursively_default_dict)


######################################################################################
### Load data 


def FAMread(Famfile):
    '''
    reads plink .fam file for names of accesssions in geneo file.
    '''
    File = open(Famfile,"r")
    Inds = recursively_default_dict()
    d = 0
    for line in File:
        line = line.split()
        Inds[d] = line[0]
        Inds[line[0]] = d
        d += 1
    File.close()
    return Inds




def read_refs(index_file,Fam_lib):
    '''
    ref file indexes individuals to population code.
    '''
    indxs = recursively_default_dict()
    Absent= []
    
    Input = open(index_file,'r')
    for line in Input:
        line = line.split()
        if line[1] not in Fam_lib.keys():
            Absent.append(line[1])
            continue
        indxs[int(line[0])][Fam_lib[line[1]]] = []
    
    Input.close()
    
    indxs = {gop:[x for x in indxs[gop].keys()] for gop in indxs.keys()}
    
    return indxs, [x for x in sorted(indxs.keys())], Absent




def BIMread(bimFile):
    '''
    reads .bim file from plink genomic data.
    returns dictionary of {geno_snp_index: locus}
    '''
    File = open(bimFile,"r")
    #Nsnps = {x:recursively_default_dict() for x in range(1,13)}
    Nsnps = recursively_default_dict()
    Gindex = recursively_default_dict()
    for x in range(1,13):
        Nsnps[x] = recursively_default_dict()
    d = 0
    CHR = 0
    for line in File:
        line = line.split()
        if d == 0:
            CHR = int(line[0])
        if d != 0 and int(line[0]) != CHR:
            d = 0
            CHR = int(line[0])
        Nsnps[CHR][d] = [float(line[3]),line[4],line[5]]
        if float(line[3]) in Gindex[CHR].keys():
            Gindex[CHR][float(line[3]) + .5] = [d]
        else:
            Gindex[CHR][float(line[3])] = [d]
        d += 1
    File.close()
    return Nsnps,Gindex




def read_geno_books(books):
    '''
    reads files. required pattern: _chr(i)
    where i = chromosome number.
    Tag will be string preceding underscore.
    '''
    library= []
    
    for shelf in books:
        card= shelf.split('/')
        
        cover= card[-1].split('_')
        
        Chr= int([re.findall(r'\d+',i)[0] for i in cover if re.search('chr',i)][0])
        
        tag= cover[0]
        
        library.append([shelf,tag,Chr])
    
    library= pd.DataFrame(library,columns= ['file','tag','Chr'])
    return library



def Gen_rand(Snp_lib,chromosomes,n,L):
    from random import randint
    
    Seen = {x:recursively_default_dict() for x in range(1,13)}
    
    for i in range(n):
        CHR= np.random.choice(chromosomes,1)[0]
        snp1= randint(0,len(Snp_lib[CHR]) - L)
        snp2= snp1 + L - 1
        
        snp1= Snp_lib[CHR][snp1][0]
        snp2= Snp_lib[CHR][snp2][0]
        
        print('positions {} through {} CHR {} taken. Tlength= {}'.format(snp1,snp2,CHR, snp2 - snp1))
        
        
        Seen[CHR][snp1] = [snp2,'rand']
    return Seen



def Extract_to_dict(Genes,MissG,Whose,Books):
    
    SequenceStore = {CHR:{GEN:{fy:[] for fy in Whose} for GEN in Genes[CHR].keys()} for CHR in Genes.keys()}
    
    for i in range(Books.shape[0]):
        
        Chr= Books.loc[i,'Chr']
        Miss = MissG[Chr]
        CHRsize = len(Miss)
        
        Geno = Books.loc[i,'file']
        Geno = open(Geno,"r")
        
        index = 0
        
        for line in Geno:
            Codes = [0,1,2,0,0,0,0,0,0,np.nan]
            d = Miss[index][0]
            
            for start in Genes[Chr].keys():
                
                if d >= (start) and d <= Genes[Chr][start][0]:
                    for judas in SequenceStore[Chr][start].keys():
                        SequenceStore[Chr][start][judas].append(Codes[int(line[judas])])
            index += 1
        
        Geno.close()  
    
    return SequenceStore


############################ Stats ##############################################""

def JustTease(i,A1,A2,Well):
    Well[i] =sum([int(A1[x] != A2[x]) for x in range(len(A1))])

def CombDiffrecv3(L,Matrix,Plunder,Well):
    for x in range(len(Plunder)):
        JustTease(x,Matrix[Plunder[x][0]],Matrix[Plunder[x][1]],Well)

def Org_comb(L,Dict_diff):
    Set = {a:[0] for a in range(L-1)}
    for unit in Set.keys():
        start = sum(range(L-unit,L))
        end = start + (L-unit) -1
        for St in range(start,end):
            Set[unit].append(Dict_diff[St])
    Set[L] = [0]
    return [x for x in Set.values()]


def SecondHalf(L,differences_matrix,populations,total_populations):
    population_list = list(set(populations))
    SSTOT = 0
    SSTOT = sum([sum(x) for x in differences_matrix])
    
    SSTOT = float(SSTOT)/float(L)				
    
    SSWP_each = 0
    SSWP_divisor = 0
    SSWP = 0
    
    for population in population_list:
        SSWP_each = 0
        SSWP_divisor = 0
        for i in range(len(differences_matrix)):
            differences = differences_matrix[i]
            for x in range(len(differences)):
                
                if populations[i] == populations[x+i] and populations[i] == population:
                    SSWP_each = SSWP_each + differences[x]
                    SSWP_divisor = SSWP_divisor + 1
        SSWP_divisor = (2*SSWP_divisor+0.25)**0.5 - 0.5
        if SSWP_each != 0:
            SSWP += float(SSWP_each)/float(SSWP_divisor)
    
    SSAP = SSTOT - SSWP
    
    squared_count_sum = float(sum([populations.count(x)**2 for x in list(set(populations))]))
    
    total_samples = float(L)
    total_pops = float(total_populations)
    dfAP = total_populations - 1
    dfWP = L - total_populations
    MSAP = float(SSAP/dfAP)
    MSWP = float(SSWP/dfWP)
    N0 = float((total_samples - float(squared_count_sum/total_samples)) * float(1/(total_pops-1)))
    VAP = float((MSAP - MSWP)/N0)
    if VAP + MSWP == 0:
        PhiPT = 0
    else:
        PhiPT = float(VAP/(VAP + MSWP))
    return PhiPT



def findPhiPT(allele_profiles,populations,n_boot):
    '''
    allele_profiles: list of haplotype vectors (numeric, string doesnt matter).
    populations: list of population assignment of accessions in allele_profiles.
                -> same length as allele_profiles, same order.
    will treat NA's as alleles.. either pre-remove them or modify CombDiffrecv3
    '''
    different_populations = list(set(populations))
    population_list = different_populations
    total_populations = len(different_populations)
    
    #allele_profiles = [''.join([str(x) for x in y]) for y in allele_profiles]
    differences_matrix = recursively_default_dict()
    
    Its = [x for x in it.combinations(range(len(allele_profiles)),2)]
    Size = len(Its)
    CombDiffrecv3(Size,allele_profiles,Its,differences_matrix)
    differences_matrix = Org_comb(len(allele_profiles),differences_matrix)
    
    SSTOT = 0
    SSTOT = sum([sum(x) for x in differences_matrix])
    
    SSTOT = float(SSTOT)/float(len(allele_profiles))				
    
    PhiPT = SecondHalf(len(allele_profiles),differences_matrix,populations,total_populations)
    District = []
    for IT in range(n_boot):
        District.append(SecondHalf(len(allele_profiles),differences_matrix,random.sample(populations,len(populations))),total_populations)
    Sign = 0
    if District:
        if PhiPT >= 1.98 * (np.std(District) + np.mean(District)) or PhiPT <=  - 1.98 * (np.std(District) + np.mean(District)):
            Sign += 1
    return PhiPT,Sign




def Structure_profiles(feats,label_select,N,Bandwidth_split):
    
    Struct_dict= {}
    contract= [x for x in label_select.keys() if len(label_select[x]) > 1]
    #
    for sub in range(2,len(contract)+ 1):
        
        combs= [x for x in it.combinations(contract,sub)]
        
        for combi in combs:
            
            subsection= {x:label_select[x] for x in label_select.keys() if x in combi}
            subst_profiles= Distance_profiles(feats,subsection,N,Bandwidth_split)
            
            if subst_profiles:

                vectors= np.amax([x for x in subst_profiles.values()],axis= 0)

                Struct_dict[combi]= vectors
            else: print('empty')
    
    return Struct_dict

    

def Distance_profiles(feats,label_select,N,Bandwidth_split):
    
    Proxy_data= []
    label_select_labels= [z for z in it.chain(*[[x] * len(label_select[x]) for x in label_select.keys()])]
    Center_store= {}
    Proxy_indexes= {}
    distance_vecs= {}
    params = {'bandwidth': np.linspace(np.min(feats), np.max(feats),Bandwidth_split)}
    grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0)
    
    for lab in label_select.keys():
        if len(label_select[lab]) < 3:
            continue
        
        Quanted_set= feats[label_select[lab],:]
        
        grid.fit(Quanted_set)
        
        kde = grid.best_estimator_
        Extract= kde.sample(N)
        
        center= np.mean(Extract,axis= 0)
        Center_store[lab]= center
        Proxy_indexes[lab]= [x for x in range((len(Center_store) - 1) * N, len(Center_store) * N)]        
        #
        Return= Extract
        #
        Proxy_data.extend(Return)
    
    Proxy_data= np.array(Proxy_data)
    ### get distances to other centers:
    Distances_vectors= []
    if len(Center_store) == 1:
        print('empty')
        
        return {}
    
    else:
        for lab in Center_store.keys():

            Others= [z for z in it.chain(*[Proxy_indexes[x] for x in Proxy_indexes.keys() if x != lab])]

            distances= euclidean_distances(Center_store[lab].reshape(1,-1),Proxy_data[Others,:])
            #

            X_plot = np.linspace(0, 10, 1000)

            kde = KernelDensity(kernel='gaussian', bandwidth=0.15).fit(np.array(distances).reshape(-1,1))
            log_dens = kde.score_samples(np.array(X_plot).reshape(-1,1))
            log_dens= [np.exp(x) for x in log_dens]
            #
            distance_vecs[lab] = log_dens
    
    return distance_vecs




##########################################################################################
### Fst
### Calculate pairwise Fst based on frequency vectors selected.
### return total Fst
def return_fsts(vector_lib,pops):
    
    H= {pop: [1-(vector_lib[pop,x]**2 + (1 - vector_lib[pop,x])**2) for x in range(vector_lib.shape[1])] for pop in pops}
    Store= []
    for comb in it.combinations(pops,2):
        P= [sum([vector_lib[x,i] for x in comb]) / len(comb) for i in range(vector_lib.shape[1])]
        HT= [2 * P[x] * (1 - P[x]) for x in range(len(P))]
        Fst= np.mean([(HT[x] - np.mean([H[p][x] for p in comb])) / HT[x] for x in range(len(P))])
        
        Store.append([comb,Fst])
    
    ### total fst:
    P= [sum([vector_lib[x,i] for x in pops]) / len(pops) for i in range(vector_lib.shape[1])]
    HT= [2 * P[x] * (1 - P[x]) for x in range(len(P))]
    FST= np.mean([(HT[x] - np.mean([H[p][x] for p in pops])) / HT[x] for x in range(len(P))])
    
    return pd.DataFrame(Store,columns= ['pops','fst']),FST


### Dont return global Fst
def return_fsts2(freq_array):
    pops= range(freq_array.shape[0])
    H= {pop: [1-(freq_array[pop,x]**2 + (1 - freq_array[pop,x])**2) for x in range(freq_array.shape[1])] for pop in range(freq_array.shape[0])}
    Store= []

    for comb in it.combinations(H.keys(),2):
        P= [sum([freq_array[x,i] for x in comb]) / len(comb) for i in range(freq_array.shape[1])]
        HT= [2 * P[x] * (1 - P[x]) for x in range(len(P))]
        per_locus_fst= [[(HT[x] - np.mean([H[p][x] for p in comb])) / HT[x],0][int(HT[x] == 0)] for x in range(len(P))]
        per_locus_fst= np.nan_to_num(per_locus_fst)
        Fst= np.mean(per_locus_fst)

        Store.append([comb,Fst])
    
    return pd.DataFrame(Store,columns= ['pops','fst'])



###################################################################""""
### Local sampling correct (Notebook 8.)
################################################

### Data generation

#### Union - calculate 3D overlap from KDE estimation
#### Function to calculate overlap given coordinates matrix and dictionary of indicies.
def extract_profiles_union(global_data,target_ind_dict,threshold,P):
    ## estimate the bandwith
    params = {'bandwidth': np.linspace(np.min(global_data), np.max(global_data),20)}
    grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0)

    ## perform MeanShift clustering.
    combine= {}
    for bull in target_ind_dict.keys():
        grid.fit(global_data[target_ind_dict[bull],:])
        combine[bull]= grid.best_estimator_    

    Stats= recursively_default_dict()

    for combo in it.combinations(target_ind_dict.keys(),2):
        pop1= combo[0]
        pop2= combo[1]

        All_coords= [x for x in it.chain(*[target_ind_dict[z] for z in combo])]

        Quanted_set= global_data[All_coords,:]

        i_coords, j_coords, z_coords = np.meshgrid(np.linspace(min(Quanted_set[:,0]),max(Quanted_set[:,0]),P),
                              np.linspace(min(Quanted_set[:,1]),max(Quanted_set[:,1]),P),
                                np.linspace(min(Quanted_set[:,2]),max(Quanted_set[:,2]),P), indexing= 'ij')


        traces= [x for x in it.product(range(P),range(P),range(P))]

        background= np.array([i_coords,j_coords,z_coords])

        background= [background[:,c[0],c[1],c[2]] for c in traces]

        background=np.array(background)

        pop1_fist= combine[pop1].score_samples(background)
        #pop1_fist= np.exp(pop1_fist)
        P_dist_pop1= combine[pop1].score_samples(global_data[target_ind_dict[pop1],:])
        pop1_fist = scipy.stats.norm(np.mean(P_dist_pop1),np.std(P_dist_pop1)).cdf(pop1_fist)
        pop1_fist= [int(x >= threshold) for x in pop1_fist]
        
        pop2_fist= combine[pop2].score_samples(background)
        #pop2_fist= np.exp(pop2_fist)
        P_dist_pop2= combine[pop2].score_samples(global_data[target_ind_dict[pop2],:])
        pop2_fist = scipy.stats.norm(np.mean(P_dist_pop2),np.std(P_dist_pop2)).cdf(pop2_fist)
        pop2_fist= [int(x >= threshold) for x in pop2_fist]

        
        pop1_and_2= len([x for x in range(background.shape[0]) if pop1_fist[x] == 1 and pop2_fist[x] == 1])
        pop1_I_pop2= pop1_and_2 / float(sum(pop1_fist))
        pop2_I_pop1= pop1_and_2 / float(sum(pop2_fist))
        
        total_overlap= pop1_and_2 / float(sum(pop1_fist) + sum(pop2_fist) - pop1_and_2)
        
        empty_space= 1 - (sum(pop1_fist) + sum(pop2_fist) - pop1_and_2) / background.shape[0]
        
        Stats[combo][pop1]= pop1_I_pop2
        Stats[combo][pop2]= pop2_I_pop1
        Stats[combo]['empty']= empty_space
        Stats[combo]['PU']= total_overlap
        
    
    return Stats


#### reference KDE

def extract_profiles(global_data,target_ind_dict):
    ## estimate the bandwith
    params = {'bandwidth': np.linspace(np.min(global_data), np.max(global_data),20)}
    grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0)
    
    cluster_profiles= {x:[] for x in target_ind_dict.keys()}
    
    
    ## perform MeanShift clustering.
    combine= {}
    for bull in target_ind_dict.keys():
        Quanted_set= global_data[target_ind_dict[bull],:]
        grid.fit(Quanted_set)
        kde = grid.best_estimator_

        P_dist = kde.score_samples(Quanted_set)
        Fist = kde.score_samples(global_data)

        ## Normalizing log-likelihood estimates by those of the reference set and extracting their cdf.
        Fist = scipy.stats.norm(np.mean(P_dist),np.std(P_dist)).cdf(Fist)
        cluster_profiles[bull].append(Fist)

    
    return cluster_profiles


def extract_profiles_class(global_data,target_ind_dict):
    '''
    copy of the previous function. change of name to deal with local 
    function similarity.
    '''
    ## estimate the bandwith
    cluster_profiles= recursively_default_dict()
    params = {'bandwidth': np.linspace(np.min(global_data), np.max(global_data),20)}
    grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0)
    
    combine= {}
    for bull in target_ind_dict.keys():
        
        Quanted_set= global_data[target_ind_dict[bull],:]
        grid.fit(Quanted_set)
        kde = grid.best_estimator_
        
        P_dist = kde.score_samples(Quanted_set)
        Fist = kde.score_samples(global_data)

        ## Normalizing log-likelihood estimates by those of the reference set and extracting their cdf.
        Fist = scipy.stats.norm(np.mean(P_dist),np.std(P_dist)).cdf(Fist)
        cluster_profiles[bull]=Fist

    
    return cluster_profiles
    


#########################################################
### Ideogram Processing

def compress_ideo(df,Out,chromosome_list):
    
    new_set = []
    
    for CHR in range(len(chromosome_list)):
        
        Chr = int(re.search('Region_(.+?)_',chromosome_list[CHR]).group(1))
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



def compress_ideo_vII(df,Out,chromosome_list):
    
    new_set = []
    
    for CHR in range(len(chromosome_list)):
        
        Chr = int(re.search('Region_(.+?)_',chromosome_list[CHR]).group(1))
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


def return_ideogram(ideo, chromosome_list, Comparison_threshold, Outlier_threshold, out= True):
    # Height of each ideogram
    chrom_height = 1

    # Spacing between consecutive ideograms
    chrom_spacing = 0

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
    
    ideo['colors'] = ideo['gieStain'].apply(lambda x: tuple([round(y / float(255),1) for y in color_lookup[x]]))
    # Add a new column for width
    ideo['width'] = ideo.end - ideo.start

    # Width, height (in inches)
    figsize = (10, 30)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)

    # Now all we have to do is call our function for the ideogram data...
    print("adding ideograms...")
    for collection in chromosome_collections(ideo, chrom_ybase, chrom_height, edgecolors=None, linewidths= 0):
        ax.add_collection(collection)

    # Axes tweaking
    ax.set_xticks([x for x in range(min(ideo.start),max(ideo.end),int(10000))])
    ax.set_xticklabels([round(x / float(10000),3) for x in range(min(ideo.start),max(ideo.end),int(10000))])
    plt.xticks(fontsize = 5,rotation = 90)
    ax.tick_params(axis = 'x',pad = 10)

    ax.tick_params(axis='y', which='major', pad=30)
    ax.set_yticks([chrom_centers[i] for i in chromosome_list])
    ax.set_yticklabels(chromosome_list, fontsize = 5)
    ax.axis('tight')
    if out == True:
        plt.savefig('Ideo_step_' + '_OutlierTh' + str(Outlier_threshold) + '_Z' +str(Comparison_threshold)+ '.png',bbox_inches = 'tight')
    return fig


