# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 17:04:03 2018

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

from Kernel_tools import recursively_default_dict, read_focus, read_refs, FAMread, BIMread
from Galaxy_Ideogram_tools import chromosome_collections, plot_ideo, plot_ideo2, compress_ideo, Merge_class
from Downstream_tools import read_books, read_3D_profiles, read_coords, Distance_profiles

import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
import pandas as pd

import collections
import time
import re
import itertools as it
import numpy as np

from sklearn.neighbors import KernelDensity
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV
from sklearn.cluster import estimate_bandwidth
from sklearn.cluster import MeanShift, estimate_bandwidth

from sklearn.metrics.pairwise import pairwise_distances
from sklearn.metrics.pairwise import euclidean_distances

import os, sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("books",type=str,metavar= 'N',nargs= '+',
                    help = "Reference files to read. Any number can be given.")

parser.add_argument("--focus",type= str,help = "reference accessions indexes in genofile.")

parser.add_argument("--id",type= str,help = "Name your analysis")

parser.add_argument("--ref",type= str,help = "reference accessions indexes in genofile.")

parser.add_argument("--coords",type= str,help = "Blocks coordinates and cluster codes.")

parser.add_argument("--genoSuf",type= str,help = ".geno file suffix.")

parser.add_argument("--code",type= str,help = "Cluster code to focus on.")

parser.add_argument("--CHR",type= int,help = "chromosome to draw ideogram of.")

parser.add_argument("--fam",help = "accession name file. same order as in geno.")

parser.add_argument("--bim",help = "snp information bim format.")

parser.add_argument("--het",type = float,default = 5e-2, help = "Heterozygosity filter")

parser.add_argument("--start",type= int,help = "where to begin, in markers. Only makes sense if --CHR is also used.")

parser.add_argument("--end",type= int,help = "where to end, in markers. Only makes sense if --CHR is also used.")

parser.add_argument("--plot",action= "store_true",help = "if given reduces cluster points using KDE extraction per cluster identified.")

parser.add_argument("--reduc",action= "store_true",help = "if given prints cluster stats.")

parser.add_argument("--coarse",action='store_false',help= 'to smooth or not to smooth.')

parser.add_argument("--bin",default = 5,type= int,help = "smoothing parameter, must be uneven [savgol filter]")

parser.add_argument("--sg_order",default = 3,type= int,help = "staviksy golay filter order")

parser.add_argument("--outlier",type=float,default = 0.05,help = "Outlier threshold")

parser.add_argument("--ms",type=float,default = .1,help = "cluster profile selection threshold")

parser.add_argument('--ncomps',type= int,default= 5,help= 'PCA comps to retain.')

parser.add_argument("--chrom_height",type= float, default= 1, help= "height of ideograms")

parser.add_argument("--chrom_gap",type= float,default= 0,help= "gap between ideograms.")

parser.add_argument("--height",type= float, default= 10,help= "figure height, in inches.")

parser.add_argument("--width",type= float,default= 20,help= "figure width, in inches.")

parser.add_argument('--xticks',type= int,default= 1000000,help= 'xticks on final ideogram')


args = parser.parse_args()


Home= 'Analyses_' + args.id

########## Complementary files.

####
#### read coordinates
####

request_coords= read_coords(args.coords)

### sep code labels select:
cluster_code= [int(x) - 1 for x in args.code.split(',')]

coords_list= [z for z in it.chain(*[request_coords[z] for z in cluster_code])]
coords_list= np.array(coords_list)

####
####
if args.CHR:
    coords_list= coords_list[[x for x in range(coords_list.shape[0]) if coords_list[x,0] == args.CHR],:]

if args.start:
    coords_list= coords_list[[x for x in range(coords_list.shape[0]) if coords_list[x,1] >= args.start],:]
    if not args.CHR:
        print("start was selected with no CHR specification.")

if args.end:
    coords_list= coords_list[[x for x in range(coords_list.shape[0]) if coords_list[x,2] <= args.end],:]
    if not args.CHR:
        print("end was selected with no CHR specification.")

chromosomes= list(set(coords_list[:,0]))



### extract windows surveyed during genome crawl
print('to begin reading from: ')
print(args.books)

Library= read_books(args.books)

print('library:')
print(Library.sort_values(by= ['Chr','start']))

Ref_blocks, Ref_profiles, Names, Out = read_3D_profiles(Library)


######
###### deciding who you're going to obe looking at;
###### in the future, replace with input file only it should be easier.
Fam = FAMread(args.fam)

####
#### some parameters
####

X_threshold= args.outlier
n_comps= args.ncomps

refs_lib, Parents, absent_refs = read_refs(args.ref,Fam)

if args.focus:
    Focus = read_focus(args.focus)
else:
    Focus = Names

Absent= [x for x in Focus if x not in Names]

if len(Absent) > 0:
    print('The following individuals were not found in the files provided: {}'.format(Absent))
    print('Analysis will proceed without them.')
    Focus= [x for x in Focus if x in Names]

compound_reference = [Names.index(x) for x in Focus]


####
MissG, Gindex = BIMread(args.bim)


##### Selecting windows to focus on:

Whose= [Fam[y] for y in Names]

Basket= {Chr:{bl:{fy:[] for fy in Whose} for bl in list(set(coords_list[coords_list[:,0] == Chr,1]))} for Chr in chromosomes}

Biblio= {Chr:{bl:list(set([coords_list[x,3] for x in range(coords_list.shape[0]) if coords_list[x,0] == Chr and coords_list[x,1] == bl])) for bl in Basket[Chr].keys()} for Chr in chromosomes}

Parent_list= [Fam[x] for x in Focus if Fam[x] not in it.chain(*refs_lib.values())]
Parent_list.extend([y for y in it.chain(*refs_lib.values())])

Clover= []
Coordinates= []
Clusters_coords= []

Dist_vectors= []
Ind_labels= {CHR:{bl:recursively_default_dict() for bl in Biblio[CHR].keys()} for CHR in Biblio.keys()}

Distances= []
center_distances= []

Ref_stats= []
Ref_stats_lib= recursively_default_dict()

Empty= []

Whose_parents= [Whose.index(x) for x in Parent_list]
iu1 = np.triu_indices(len(Parent_list))

for CHR in chromosomes:
    
    Fam = FAMread(args.fam)
    GenoFile= args.genoSuf + str(CHR).zfill(2) + '.geno'
    
    #### Fam
    #### Miss
    #### GenoFile
    Miss= MissG[CHR]
    
    ######################
    ######################
    ## Heterozygosity filter
    Filter_Het = args.het
    
    Geno = open(GenoFile,"r")
    
    Win = 0
    Index = 0
    
    for line in Geno:
        Codes = [0,0,2,0,0,0,0,0,0,np.nan]
        d = Miss[Index][0]
        recycle= []
        
        if len([x for x in Whose if line[x] =='1']) / float(len(Whose)) > Filter_Het:
            Index += 1
            continue
        
        for bl in sorted(Basket[CHR].keys()):
            if d >= bl and d <= Out[CHR][bl]:
                for judas in Basket[CHR][bl].keys():
                    Basket[CHR][bl][judas].append(Codes[int(line[judas])])
                continue
            if d > Out[CHR][bl]:
                
                data= np.array([Basket[CHR][bl][x] for x in Whose])
                
                ###
                polyM= np.nansum(data,axis= 0)
                var_index= [x for x in range(len(polyM)) if polyM[x] > 0]

                print('{} polymorphic sites.'.format(str(len(var_index))))
                
                if len(var_index) <= n_comps:
                    Empty.append([CHR,bl])
                    recycle.append(bl)
                    continue
                
                ####
                data= np.nan_to_num(data[:,var_index])
                
                Nsnps= data.shape[1]
                
                pca = PCA(n_components=n_comps, whiten=False).fit(data)
                data = pca.transform(data)
                COMPS = pca.components_.T*np.sqrt(pca.explained_variance_)
                
                ### get cluster membership
                Bls= list(Ref_profiles[CHR][bl].keys())
                pVals= np.array([Ref_profiles[CHR][bl][y] for y in Bls])
                print(pVals.shape)
                
                max_vals= np.amax(pVals,axis= 0)
                max_indx= np.argmax(pVals,axis= 0)
                
                inlier= [x for x in range(pVals.shape[1]) if max_vals[x] >= args.outlier]
                
                BL_select= list(set([max_indx[x] for x in inlier]))
                
                if not BL_select:
                    Empty.append([CHR,bl])
                    recycle.append(bl)
                    continue
                
                BL_select= {
                    x: pVals[x] for x in BL_select
                    }
                
                BLextract= list(BL_select.keys())
                
                Assignment= {
                        Bls[b]: [x for x in inlier if max_indx[x] == b] for b in BLextract
                    }
                
                for cl in Biblio[CHR][bl]:
                    
                    #Who= [x for x in range(len(Whose)) if Ref_profiles[CHR][bl][cl][x] >= X_threshold]
                    Who= Assignment[cl]
                    in_focus= [x for x in Who if Fam[Whose[x]] in Focus]
                    accuracy= len(in_focus) / float(len(Focus))
                    
                    print('accuracy: {} %'.format(round(accuracy,2)))
                    #Refs_local= [x for x in Whose_parents if x not in Who]
                    Who_feat= data[Who,:]
                    Ref_feat= data[Whose_parents,:]
                    
                    #### Normalize by distance between local centroids (to compensate for bias in sampling number).
                    #### identify these clusters using MS.
                    #### use reference accessions NOT in the target cluster.
                    Dpool= data[[x for x in Whose_parents if x not in Who],:]
                    Pdistances= []
                    
                    bandwidth = estimate_bandwidth(Dpool, quantile=0.15)
                    if bandwidth <= 0:
                        bandwidth= .1
                    params = {'bandwidth': np.linspace(np.min(Dpool), np.max(Dpool),30)}
                    grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0)
                    
                    ## perform MeanShift clustering.
                    ms = MeanShift(bandwidth=bandwidth, bin_seeding=False, cluster_all=False, min_bin_freq=25)
                    ms.fit(Dpool)
                    labels1 = ms.labels_
                    label_select = {y:[x for x in range(len(labels1)) if labels1[x] == y] for y in sorted(list(set(labels1))) if y != -1}
                    
                    centers= [np.mean(Dpool[label_select[z],:],axis= 0) for z in label_select.keys()]
                    
                    #### Data set of evenly sampled data. ##
                    ## We'll generate 50 new observations from each cluster identified locally. ##
                    N= 50
                    Proxy_data= []
                    label_select_labels= [z for z in it.chain(*[[x] * len(label_select[x]) for x in label_select.keys()])]
                    Center_store= {}
                    Proxy_indexes= {}
                    distance_vecs= []
                    
                    for lab in label_select.keys():
                        if len(label_select[lab]) < 3:
                            continue
                            
                        Quanted_set= Dpool[label_select[lab],:]
                        
                        if np.max(pairwise_distances(Quanted_set,metric= 'euclidean')) <= 1e-3:
                            Extract= Quanted_set[np.random.choice(Quanted_set.shape[0],N),:]
                        else:
                            grid.fit(Quanted_set)
                            kde = grid.best_estimator_
                            Extract= kde.sample(N)
                            
                        center= np.mean(Extract,axis= 0)
                        Center_store[lab]= center
                        Proxy_indexes[lab]= [x for x in range((len(Center_store) - 1) * N, len(Center_store) * N)]        
                        #Return= pca.inverse_transform(Extract)
                        
                        #Return= data_now[np.random.choice(label_select[lab],N),:]
                        Proxy_data.extend(Extract)
                    
                    Proxy_data= np.array(Proxy_data)
                    ##### Get pairwise distances between centroids.
                    
                    for pair in it.combinations(label_select.keys(),2):
                        coordinates= [np.mean(Dpool[label_select[z],:],axis= 0) for z in pair]
                        coordinates= np.array(coordinates)
                        iu_control= np.triu_indices(2,1)
                        MS_pair_dist= pairwise_distances(coordinates,metric= 'euclidean')
                        MS_pair_dist= MS_pair_dist[iu_control][0]
                        Pdistances.append(MS_pair_dist)
                    ## 
                    
                    
                    reference_centroid= np.mean(centers,axis= 0)
                    
                    proxy_distances= pairwise_distances(reference_centroid.reshape(1,-1), Proxy_data,metric= 'euclidean')
                    distances_to_center= pairwise_distances(reference_centroid.reshape(1,-1), Ref_feat,metric= 'euclidean')[0]
                    self_distances= pairwise_distances(reference_centroid.reshape(1,-1), Who_feat, metric= 'euclidean')
                    
                    #distances_to_center= distances_to_center[0] / np.mean(Pdistances)
                    #distances_to_center= scipy.stats.norm(np.mean(proxy_distances),np.std(proxy_distances)).cdf(distances_to_center)[0]
                    #distances_to_center= [(x - np.mean(proxy_distances)) / np.std(proxy_distances) for x in distances_to_center[0]]
                    
                    centroid= np.mean(Who_feat,axis= 0)
                    distances_pairwise= pairwise_distances(centroid.reshape(1,-1), Ref_feat, metric= 'euclidean')[0]
                    
                    #distances_pairwise= [(x - np.mean(proxy_distances)) / np.std(proxy_distances) for x in distances_pairwise[0]]
                    #distances_pairwise= scipy.stats.norm(np.mean(proxy_distances),np.std(proxy_distances)).cdf(distances_pairwise)[0]
                    #distances_pairwise= distances_pairwise[0] / np.mean(Pdistances)
                    
                    ## still deciding why this should be done.
                    
                    #dist_vector= Distance_profiles(centroid,Dpool,5,label_select)
                    #Dist_vectors.append(dist_vector)
                    ##
                    
                    Distances.append(distances_pairwise)
                    distances_pairwise= [(x - np.mean(proxy_distances)) / np.std(proxy_distances) for x in distances_pairwise]
                    Clover.append(distances_pairwise)
                    print(np.array(Clover).shape)
                    
                    FC_stats= [np.mean(proxy_distances),np.std(proxy_distances), np.mean(self_distances), np.std(self_distances)]
                    
                    Coord= [[CHR,bl,Out[CHR][bl],x,cl] for x in Who]
                    
                    Ref_stats.append(FC_stats)
                    Ref_stats_lib[CHR][bl][cl]= FC_stats
                    
                    center_distances.append(distances_to_center)
                    
                    Ind_labels[CHR][bl][cl]= Who
                    Coordinates.extend(Coord)
                    Clusters_coords.append([CHR,bl,Out[CHR][bl],cl,Nsnps])
                
                recycle.append(bl)
                continue
        
        for bl in recycle:
            print('CHR: {}; Coords: {}; left: {}'.format(CHR,len(Clusters_coords),len(Basket[CHR])))
            del Basket[CHR][bl]
                
        Index += 1
    
    Geno.close()




#Profiles= np.array(Profiles)

Ref_stats= np.array(Ref_stats)
Clover= np.array([x for x in Clover])
Coordinates= np.array(Coordinates)
Clusters_coords= np.array(Clusters_coords)
Distances= np.array(Distances)
center_distances= np.array([x for x in center_distances])
#center_distances= np.array([x.reshape(1,-1)[0] for x in center_distances])

from sklearn import preprocessing


### initial step of  filtering.
### removing windows where variation between references resulted
### in normalized min and maximum distances above threshold. 
###

Trim_threshold= 20

trim_indexes= [x for x in range(Clover.shape[0]) if \
    min(Clover[x]) >= -Trim_threshold and max(Clover[x]) <= Trim_threshold and \
    min(center_distances[x]) >= -Trim_threshold and max(center_distances[x]) <= Trim_threshold]

Ref_stats= Ref_stats[trim_indexes]
Distances= Distances[trim_indexes]
Clover= Clover[trim_indexes]
center_distances= center_distances[trim_indexes]
Clusters_coords= Clusters_coords[trim_indexes]

####
####
####

####
#### Pre processing and dimensionality reduction of matrix
#### of selected clusters.
####

Clover = np.nan_to_num(Clover)
Distances= np.nan_to_num(Distances)
center_distances= np.nan_to_num(center_distances)

preProc_Clover = Clover



print('Clover shape: ', Clover.shape)

reefer= [g for g in it.chain(*[refs_lib[y] for y in sorted(refs_lib.keys())])]

if not args.focus:
    mary= [Names.index(x) for x in Names if Fam[x] not in reefer]
else:
    mary= [Names.index(x) for x in Names if Fam[x] not in reefer and x in Focus]

reefer= [Names.index(Fam[x]) for x in reefer]

Subset= [x for x in it.chain(*[mary,reefer])]
Trend= np.repeat([0,1],[len(mary),len(reefer)])


## apply pca to reference accessions, transform the rest.
variation_focus= [x for x in range(len(Parent_list)) if Fam[Parent_list[x]] not in Focus]

### PCA
pca = PCA(n_components=5, whiten=False).fit(Clover[:,variation_focus].T)
X_se = pca.transform(Clover.T)
COMPS81 = pca.components_.T*np.sqrt(pca.explained_variance_)

###############################################################################
###############################################################################
###############################################################################

###### Library on decompositions (deprecated):
from sklearn.decomposition import PCA, FactorAnalysis, FastICA
Library= recursively_default_dict()

#dr_names= ['PCA','FA','FastICA']
#Decoms= [PCA(n_components= 5),FactorAnalysis(n_components= 5)]

####
#### Library on clustering algorithms
Group_selection = "Kmeans"
DOOM = 0
cos_threshold = .6

from sklearn.cluster import DBSCAN
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import KMeans

dr_names= ['MeanShift','DBscan','Ward','Kmeans']

dr= PCA(n_components= 5).fit(Clover)
COMPS= dr.transform(Clover)

bandwidth = estimate_bandwidth(COMPS, quantile=0.15)
if bandwidth <= 0:
    bandwidth= .1

Decoms= [
    MeanShift(bandwidth=bandwidth, bin_seeding=False, cluster_all=True, min_bin_freq=350),
    DBSCAN(min_samples=35),
    AgglomerativeClustering(linkage='ward', n_clusters=10),
    KMeans(n_clusters=10, random_state=0)
]

###

label_store= []

for decom in range(len(Decoms)):
    
    kmeans = Decoms[decom].fit(COMPS)
    labels1 = kmeans.labels_
    
    label_select = {y:[x for x in range(len(labels1)) if labels1[x] == y] for y in sorted(list(set(labels1)))}
    
    Cameo = []
    
    for cramp in sorted(label_select.keys()):
        Clamp = np.mean(preProc_Clover[label_select[cramp],:],axis = 0)
        Fry = [x for x in Clamp]
        Cameo.append(Fry)
    
    Cameo = np.array(Cameo).T
    print(Cameo.shape)
    print(Cameo[:5,:5])
    COMPS= pd.DataFrame(COMPS,columns= ['pc{}'.format(x+1) for x in range(COMPS.shape[1])])
    COMPS['label']= labels1
    
    
    Library[dr_names[decom]]= {
        'features': COMPS,
        'KDE':pd.DataFrame(Cameo),
        'labels_l1': labels1
    }
    
    print(Cameo[:5,:5])
    labels_second_layer= [-1]*Distances.shape[0]
    
    for Cl in label_select.keys():
        if len(label_select[Cl]) <= len(label_select):
            continue
        ## retrieve distances to centroids selected
        New_comp= Distances[[x for x in label_select[Cl]]]
        ## identify problem rows.
        NAs_row= [sum(np.isnan(New_comp[x])) for x in range(New_comp.shape[0])]
        NAs_row= [x for x in range(len(NAs_row)) if NAs_row[x] == New_comp.shape[1]]
        ## remove problem rows.
        New_comp= New_comp[[x for x in range(New_comp.shape[0]) if x not in NAs_row]]
        ## retrieve distances to global center for centroids selected
        distance_to_center= center_distances[[x for x in label_select[Cl]]]
        distance_to_center= distance_to_center[[x for x in range(New_comp.shape[0]) if x not in NAs_row]]
        
        ### PCA on Distances
        pca= PCA(n_components= 3,whiten= False).fit(New_comp)
        new_feat= pca.transform(New_comp)
        clock = pca.components_.T*np.sqrt(pca.explained_variance_)
        
        ## Kmeans in feature space
        from sklearn.cluster import KMeans
        kmeans = KMeans(n_clusters=8, random_state=0).fit(new_feat)
        new_labels = kmeans.labels_
        another_label_select = {y:[x for x in range(len(new_labels)) if new_labels[x] == y] for y in sorted(list(set(new_labels)))}
        
        ### set up second layer of labels:
        for color in another_label_select.keys():
            for thot in range(len(another_label_select[color])):
                if another_label_select[color][thot] not in NAs_row:
                    indexed= label_select[Cl][another_label_select[color][thot]]
                    labels_second_layer[indexed]= color
        
        ### average distances across clustered profiles.
        Cameo = []
        center_means= []
        for cramp in sorted(another_label_select.keys()):
            Fry = np.mean(New_comp[another_label_select[cramp],:],axis = 0)
            Phillip= np.mean(distance_to_center[another_label_select[cramp],:],axis = 0)
            
            center_means.append(Phillip)
            Cameo.append(Fry)
        
        ### prepare data to save
        Cameo = np.array(Cameo).T
        center_means= np.array(center_means).T
        
        new_feat= pd.DataFrame(new_feat,columns= ['pc{}'.format(x+1) for x in range(new_feat.shape[1])])
        new_feat['label']= new_labels
        
        clock= pd.DataFrame(clock, columns= ['pc{}'.format(x+1) for x in range(clock.shape[1])])
        clock['id']= [Fam[x] for x in Parent_list]
                
        Library[dr_names[decom]]['stats'][Cl]= {
            'profiles': new_feat,
            'features': clock,
            'averages': pd.DataFrame(Cameo),
            'shapes': pd.DataFrame(center_means)
        }
    Library[dr_names[decom]]['labels_l2']= labels_second_layer
    


Library['Ref_stats']= Ref_stats
Library['Distances']= Distances
Library['centre_dists']= center_distances

Home= args.id

try:
    import cPickle as pickle
except ImportError:  # python 3.x
    import pickle

filename= Home+ '/layer_analysis.p'
os.makedirs(os.path.dirname(filename), exist_ok=True)

with open(filename, 'wb') as fp:
    pickle.dump(Library, fp, protocol=pickle.HIGHEST_PROTOCOL)



labels1= Library[Group_selection]['labels_l1']
label_select = {y:[x for x in range(len(labels1)) if labels1[x] == y] for y in sorted(list(set(labels1)))}


###############################################################################
#### Average normalized likelihhod among clustered eigenvectors by haplotype #####
###############################################################################

#### if Haplotypes are clustered instead of clusters:
#Clamp = np.array([np.mean(Clover[[x for x in range(Clover.shape[0]) if Globe[x,y] > cos_threshold],:],axis = 0) for y in range(Globe.shape[1])]).T
#Clamp = np.mean(Clover[[x for x in range(Clover.shape[0]) if Globe[x,1] > cos_threshold],:],axis = 0)
#Fry = [Trend[x] for x in range(len(Trend)) if Clamp[x] > .1]
####


Cameo = []

for cramp in sorted(label_select.keys()):
    Clamp = np.mean(preProc_Clover[label_select[cramp],:],axis = 0)
    #Fry = [Clamp[x] for x in range(len(Clamp))]
    Cameo.append(Clamp)

Cameo = np.array(Cameo).T

Cameo= Library[Group_selection]['KDE']


###########################################################################
### cosine of the clustered eigenvectors with haplotype coordinates ######## DEPRECATED
###########################################################################

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

######## Reducing the size of the cluster profiles:
new_labs= labels1

if args.reduc:
    new_labs= []
    reduced_comp= []
    Size= 2000
    
    params = {'bandwidth': np.linspace(np.min(COMPS81), np.max(COMPS81),20)}
    grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0)
    
    for lab in label_select.keys():
        N_prop= round(len(label_select[lab]) * Size / float(sum([len(x) for x in label_select.values()])))
        
        if len(label_select[lab]) <= 3:
            reduced_comp.extend(COMPS81[label_select[lab],:])
            new_labs.extend([lab]*(len(label_select[lab])))
            continue
        
        grid.fit(COMPS81[label_select[lab],:])    
        
        kde = grid.best_estimator_
        
        new_data = kde.sample(N_prop, random_state=0)
        
        reduced_comp.extend(new_data)
        new_labs.extend([lab]*N_prop)
    
    
    COMPS81= np.array(reduced_comp)


#############################################################################
#############################################################################
############################### PRINT #######################################
#############################################################################
## Profile - Cameo: Average haplotye normalized likelihood by cluster 
## DIM_private_request - X_se: Haplotype coordinates
## DIM_private_comp - COMPS: cluster coordinates.
##

print(COMPS81.shape)
print(len(new_labs))

CHR = [x for x in chromosomes][-1]
start= args.id
target= args.code

filename= Home + "/Profile_" + str(target[0]) + "_CHR" + str(CHR) + "."+str(start) + ".txt"
os.makedirs(os.path.dirname(filename), exist_ok=True)

Output = open(filename,"w")

for title in sorted(label_select.keys()):
    Output.write("G" + str(title) + "\t")

Output.write("\n")

for Future in range(Cameo.shape[0]):
    for Machine in range(Cameo.shape[1]):
        Output.write(str(Cameo.iloc[Future,Machine]) + "\t")
    Output.write("\n")

Output.close()




filename= Home + "/Profile_coordinates_" + str(target[0]) + "_CHR" + str(CHR) + "."+str(start) + ".txt"
os.makedirs(os.path.dirname(filename), exist_ok=True)

Output = open(filename,"w")

Output.write('\t'.join(['chrom','start','end','cluster','Nsnps','label']))
Output.write('\n')

for axe in range(Clusters_coords.shape[0]):
    Output.write('\t'.join([str(x) for x in Clusters_coords[axe,:]]) + '\t')
    Output.write(str(labels1[axe]) + '\n')

Output.close()



filename= Home + "/DIM_private_"+str(target[0])+"_request_CHR" + str(CHR) + "."+str(start)+".txt"
os.makedirs(os.path.dirname(filename), exist_ok=True)

Output = open(filename,"w")
#Output = open("DIM_private_JAP_intro.txt","w")
#Order_core = [y for y in it.chain(*[Geneo[x] for x in [2,5,1,3,4]])]

for entry in range(X_se.shape[0]):
##     Output.write(str(Trend[entry]) + "\t")
    #Output.write(str(labels1[entry]) + "\t")
    Output.write(Fam[Parent_list[entry]] + "\t")
    for scythe in range(X_se.shape[1]):
        Output.write(str(X_se[entry,scythe]) + "\t")
    Output.write("\n")

Output.close()



filename= Home + "/DIM_private_"+str(target[0])+"_comp_CHR" + str(CHR) + "."+str(start)+".txt"
os.makedirs(os.path.dirname(filename), exist_ok=True)

Output = open(filename,"w")
#Output = open("DIM_private_JAP_intro.txt","w")
#Order_core = [y for y in it.chain(*[Geneo[x] for x in [2,5,1,3,4]])]

#Comps = pca.components_.T
#Output.write('0\t' + '\t'.join([str(x) for x in pca.explained_variance_ratio_]) + '\t\n')

for entry in range(COMPS81.shape[0]):
    Output.write(str(new_labs[entry]) + "\t")
#    Output.write(str(Coordinates[entry,0]) + '\t')
#    Output.write(str(Coordinates[entry,1]) + '\t')
#    Output.write(str(Coordinates[entry,2]) + '\t')
    for scythe in range(COMPS81.shape[1]):
        Output.write(str(COMPS81[entry,scythe]) + "\t")
    Output.write("\n")

Output.close()


############################################################################
############################################################################
############ Bring out ideograms ##########################################
#Coordinates= Coordinates[[x for x in range(len(Coordinates)) if Coordinates[x,3]]]


if args.plot == True:
    Ideo_home= Home + '/Ideos'
    
    Blancs= {aim:{Chr:{bl:[-1]*len(Focus) for bl in Ref_profiles[Chr].keys()} for Chr in chromosomes} for aim in list(set(labels1))}
    
    target_block= {Chr:{bl:[-1]*len(Focus) for bl in Ref_profiles[Chr].keys()} for Chr in chromosomes}
    
    for window in Coordinates:
        #[CHR,bl,Out[CHR][bl],x,cl] for x in Who]
        #[np.mean(proxy_distances),np.std(proxy_distances), np.mean(self_distances), np.std(self_distances)]
        #Ref_stats_lib
        
        local_Refmean= Ref_stats_lib[window[0]][window[1]][window[4]][0]
        local_Refsd= Ref_stats_lib[window[0]][window[1]][window[4]][1]
        local_targetmean= Ref_stats_lib[window[0]][window[1]][window[4]][2]
        
        from scipy.stats import norm
        Dists_endpoints= norm.interval(0.95, loc=local_Refmean, scale=local_Refsd)
        
        local_col= 1
        
        if local_Refsd <= 1e-3:
            local_col= 4
        else:
            if local_targetmean < Dists_endpoints[0]:
                local_col= 2
            if local_targetmean > Dists_endpoints[1]:
                local_col= 3
        
        if Fam[Whose[window[3]]] in Focus:
            target_block[window[0]][window[1]][Focus.index(Fam[Whose[window[3]]])] = local_col
    
    for Chr in chromosomes:
        plot_ideo2(target_block,[Chr],Focus,Out,'all',args.chrom_height,args.chrom_gap,args.height,args.width,Ideo_home,args.id,args.xticks)
    
    for unit in range(len(labels1)):
        site= Clusters_coords[unit]
        for native in Ind_labels[site[0]][site[1]][site[3]]:
            if Fam[Whose[native]] in Focus:
                for labelo in list(Blancs.keys()):
                    Blancs[labelo][site[0]][site[1]][Focus.index(Fam[Whose[native]])] = int(labelo == labels1[unit])
    
    for aim in Blancs.keys():
        if len(Focus)== 1:
            plot_ideo2(Blancs[aim],Blocks.keys(),Focus,Out,aim,args.chrom_height,args.chrom_gap,args.height,args.width,Ideo_home,args.id,args.xticks)
        else:
            for Chr in chromosomes:
                plot_ideo2(Blancs[aim],[Chr],Focus,Out,aim,args.chrom_height,args.chrom_gap,args.height,args.width,Ideo_home,args.id,args.xticks)




